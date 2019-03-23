{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
module Taiji.Pipeline.RNASeq.DropSeq.Functions
    ( dropSeqMkIndex
    , dropSeqAlign
    , filterSortBam
    , dropSeqQuantification
    ) where

import           Bio.Data.Bam
import           Bio.Data.Bed
import           Bio.Data.Experiment
import           Bio.HTS
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc                       (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (CI, original)
import           Data.Either                          (either, fromLeft)
import qualified Data.HashMap.Strict                  as M
import qualified Data.HashSet                         as HS
import qualified Data.IntervalMap.Strict              as IM
import qualified Data.IntSet                          as S
import           Data.List                            (foldl')
import           Data.Maybe                           (fromJust)
import           Data.Monoid                          ((<>))
import           Data.Serialize                       (decode)
import qualified Data.Text                            as T
import           Scientific.Workflow
import           System.IO.Temp                       (withTempFile)
import           Text.Printf                          (printf)

import           Taiji.Pipeline.RNASeq.DropSeq.Config

dropSeqMkIndex :: DropSeqConfig config => [a] -> WorkflowConfig config [a]
dropSeqMkIndex input
    | null input = return input
    | otherwise = do
        genome <- asks _dropSeq_genome_fasta
        starIndex <- asks _dropSeq_star_index
        anno <- asks _dropSeq_annotation
        liftIO $ do
            _ <- starMkIndex "STAR" starIndex [genome] anno 100
            return input

dropSeqAlign :: DropSeqConfig config
         => RNASeq S (File '[Gzip] 'Fastq)
         -> WorkflowConfig config (RNASeq S (File '[] 'Bam))
dropSeqAlign input = do
    dir <- asks _dropSeq_output_dir >>= getPath
    idx <- asks _dropSeq_star_index
    let outputGenome = printf "%s/%s_rep%d_genome.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        f fl = starAlign outputGenome idx (Left fl) opt >>=
            return . fst . fromLeft undefined
    input & replicates.traverse.files %%~ liftIO . f
  where
    opt = defaultSTAROpts & starCores .~ 4 & starTranscriptome .~ Nothing

filterSortBam :: DropSeqConfig config
              => RNASeq S (File '[] 'Bam)
              -> WorkflowConfig config (RNASeq S (File '[NameSorted] 'Bam))
filterSortBam input = do
    dir <- asks _dropSeq_output_dir >>= getPath
    let output = printf "%s/%s_rep%d_filt_srt.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl ->
        withTempFile "./" "tmp_file." $ \tmp _ ->
            filterBam "./" tmp fl >>= sortBamByName "./" output )

dropSeqQuantification :: DropSeqConfig config
                      => RNASeq S (File '[NameSorted] 'Bam, File '[] 'Other)
                      -> WorkflowConfig config (RNASeq S [File '[] 'Tsv])
dropSeqQuantification input = do
    let dirname = asDir $ printf "%s_rep%d" (T.unpack $ input^.eid)
            (input^.replicates._1)
    dir <- asks ((<> "/Quantification/" <> dirname) . _dropSeq_output_dir) >>= getPath
    anno_f <- asks _dropSeq_annotation
    genes <- liftIO $ readGenes anno_f
    let annotation = bedToTree HS.union $ map (\Gene{..} ->
            (asBed geneChrom geneLeft geneRight :: BED3, HS.singleton geneName)) genes
    input & replicates.traverse.files %%~ liftIO . ( \(bam, bcStat) -> do
        bc <- S.fromList . fst . unzip . filter ((>(10000::Int)) . snd) .
            either error id . decode <$> B.readFile (bcStat^.location)
        hdr <- getBamHeader $ bam^.location
        runResourceT $ runConduit $ streamBam (bam^.location) .|
            quantify dir bc annotation hdr )

quantify :: MonadIO m
         => FilePath -> S.IntSet -> BEDTree (HS.HashSet (CI B.ByteString))
         -> BAMHeader -> ConduitT BAM o m [File '[] 'Tsv]
quantify dir bcSet annotation hdr = do
    (acc, result, cur) <- foldMC fun ([], M.empty, -1)
    r <- outputResult cur result
    return $ r : acc
  where
    fun (acc, result, cur) b
        | not (cellBc `S.member` bcSet) = return (acc, result, cur)
        | cellBc == cur = do
            result' <- update molBc b result
            return (acc, result', cur)
        | M.null result = do
            result' <- update molBc b M.empty
            return (acc, result', cellBc)
        | otherwise = do
            r <- outputResult cur result
            result' <- update molBc b M.empty
            return (r:acc, result', cellBc)
      where
        (cellBc, molBc) = getBC b
    getBC b = let [cellBc, molBc, _] = B.split ':' $ queryName b
              in (readInt cellBc, readInt molBc)
    update molBc bam result = do
        let bed = asBed (fromJust $ refName hdr bam) (fromIntegral $ startLoc bam)
                (fromIntegral $ endLoc bam) :: BED3
            f Nothing  = Just $ S.singleton molBc
            f (Just x) = Just $ S.insert molBc x
        return $ foldl' (\m x -> M.alter f x m) result $ concatMap HS.toList $
            IM.elems $ intersecting annotation bed
    outputResult cellBc result = liftIO $ do
        let output = dir ++ "/" ++ show cellBc ++ ".tsv"
            result' = fmap S.size result
            total = foldl' (+) 0 result'
        B.writeFile output $ B.unlines $ map
            (\(a,b) -> original a <> "\t" <> B.pack (show b)) $ M.toList result'
        return $ emptyFile & location .~ output
                           & info .~ [ ("cell id", T.pack $ show cellBc)
                                     , ("total reads", T.pack $ show total) ]
{-# INLINE quantify #-}
