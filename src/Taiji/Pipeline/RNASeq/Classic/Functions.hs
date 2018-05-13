{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE GADTs                 #-}
{-# LANGUAGE OverloadedLists       #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE TemplateHaskell       #-}
module Taiji.Pipeline.RNASeq.Classic.Functions
    ( rnaMkIndex
    , rnaDownloadData
    , rnaAlign
    , quantification
    , geneId2Name
    , mkTable
    ) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Download
import           Bio.Pipeline.NGS
import           Bio.Pipeline.NGS.RSEM
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Report
import           Bio.Pipeline.Utils
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc                       (readDouble)
import           Control.Lens
import           Control.Monad                        (forM)
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import           Data.Bifunctor                       (bimap, second)
import           Data.Bitraversable                   (bitraverse)
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (CI, mk, original)
import           Data.Coerce                          (coerce)
import           Data.Double.Conversion.ByteString    (toShortest)
import           Data.Either                          (lefts, rights)
import qualified Data.HashMap.Strict                  as M
import qualified Data.HashSet                         as S
import           Data.List
import           Data.List.Ordered                    (nubSort)
import           Data.Maybe                           (fromJust, fromMaybe,
                                                       mapMaybe)
import           Data.Monoid                          ((<>))
import           Data.Promotion.Prelude.List          (Elem)
import           Data.Singletons                      (SingI)
import qualified Data.Text                            as T
import           Scientific.Workflow
import           Text.Printf                          (printf)

import           Taiji.Pipeline.RNASeq.Classic.Config

type RNASeqWithSomeFile = RNASeq N [Either SomeFile (SomeFile, SomeFile)]

rnaMkIndex :: RNASeqConfig config => [a] -> WorkflowConfig config [a]
rnaMkIndex input
    | null input = return input
    | otherwise = do
        genome <- asks (fromJust . _rnaseq_genome_fasta)
        starIndex <- asks (fromJust . _rnaseq_star_index)
        anno <- asks (fromJust . _rnaseq_annotation)
        rsemIndex <- asks (fromJust . _rnaseq_rsem_index)
        liftIO $ do
            starMkIndex "STAR" starIndex [genome] anno 100
            rsemMkIndex rsemIndex anno [genome]
            return input

rnaAlign :: RNASeqConfig config
         => Bool  -- ^ Whether to output transcriptome
         -> RNASeq S ( Either (SomeTags 'Fastq)
                              (SomeTags 'Fastq, SomeTags 'Fastq) )
         -> WorkflowConfig config (
                RNASeq S ( Either (File '[] 'Bam, Maybe (File '[] 'Bam))
                    (File '[PairedEnd] 'Bam, Maybe (File '[PairedEnd] 'Bam))))
rnaAlign transcriptome input = do
    dir <- asks _rnaseq_output_dir >>= getPath
    idx <- asks (fromJust . _rnaseq_star_index)
    let outputGenome = printf "%s/%s_rep%d_genome.bam" dir (T.unpack $ input^.eid)
            (runIdentity (input^.replicates) ^. number)
        outputTranscriptome = printf "%s/%s_rep%d_transcriptome.bam" dir
            (T.unpack $ input^.eid) (runIdentity (input^.replicates) ^. number)
        opt = defaultSTAROpts & starCores .~ 4 & starTranscriptome .~
            if transcriptome then Just outputTranscriptome else Nothing
        f (Left fl) = if fl `hasTag` Gzip
            then let fl' = Left $ fromSomeTags fl ::
                        Either (File '[Gzip] 'Fastq) (File '[] 'Fastq, File '[] 'Fastq)
                 in starAlign outputGenome idx fl' opt
            else let fl' = Left $ fromSomeTags fl ::
                        Either (File '[] 'Fastq) (File '[] 'Fastq, File '[] 'Fastq)
                 in starAlign outputGenome idx fl' opt
        f (Right (f1, f2)) = if f1 `hasTag` Gzip && f2 `hasTag` Gzip
            then let fl' = Right (fromSomeTags f1, fromSomeTags f2) ::
                        Either (File '[] 'Fastq) (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)
                 in starAlign outputGenome idx fl' opt
            else let fl' = Right (fromSomeTags f1, fromSomeTags f2) ::
                        Either (File '[] 'Fastq) (File '[] 'Fastq, File '[] 'Fastq)
                 in starAlign outputGenome idx fl' opt
    input & replicates.traverse.files %%~ liftIO . f

rnaDownloadData :: RNASeqConfig config
                => [RNASeqWithSomeFile]
                -> WorkflowConfig config [RNASeqWithSomeFile]
rnaDownloadData dat = do
    dir <- asks _rnaseq_output_dir >>= getPath . (<> (asDir "/Download"))
    liftIO $ dat & traverse.replicates.traverse.files.traverse %%~ downloadFiles dir

quantification :: (RNASeqConfig config, SingI tags1, SingI tags2)
               => RNASeq S ( Either (File tags1 'Bam)
                                    (File tags2 'Bam) )
               -> WorkflowConfig config (RNASeq S
                    (File '[GeneQuant] 'Tsv, File '[TranscriptQuant] 'Tsv))
quantification input = do
    dir <- asks _rnaseq_output_dir >>= getPath
    idx <- asks (fromJust . _rnaseq_rsem_index)
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (runIdentity (input^.replicates) ^. number)
    input & replicates.traverse.files %%~ liftIO . either
        (\x -> rsemQuant output idx x opt)
        (\x -> rsemQuant output idx x opt)
  where
    opt = defaultRSEMOpts & rsemCores .~ 4

-- | Retrieve gene names
geneId2Name :: (RNASeqConfig config, Elem 'GeneQuant tags ~ 'True)
            => ([RNASeqWithSomeFile], [RNASeq S (File tags 'Tsv)])
            -> WorkflowConfig config [RNASeq S (File '[GeneQuant] 'Tsv)]
geneId2Name (ori_input, quantifications) = do
    outdir <- asks _rnaseq_output_dir >>= getPath
    anno_fl <- asks (fromJust . _rnaseq_annotation)
    liftIO $ do
        id2Name <- fmap (M.fromList . map (\x -> (geneId x, geneName x))) $
            readGenes anno_fl
        let prefix = outdir ++ "/"
            suffix = "_gene_quant.tsv"
            getExp filtFn = concatMap split $ concatMap split $
                ori_input & mapped.replicates.mapped.files %~
                    map fromSomeFile . filter filtFn . lefts

            convertTaiji out fl = do
                geneId2Name_ [6] True id2Name (fl^.location) out
                return $ location .~ out $ emptyFile
        quantifications' <- mapM (mapFileWithDefName prefix suffix convertTaiji) quantifications

        let encodeExpr = getExp $ \x -> x `hasTag` GeneQuant && x `hasTag` ENCODE
        encodeExpr' <- mapM (mapFileWithDefName prefix suffix convertTaiji) encodeExpr

        let customizedExpr = getExp (\x -> x `hasTag` GeneQuant && not (x `hasTag` ENCODE))
            convertCustom out fl = do
                geneId2Name_ [2] False id2Name (fl^.location) out
                return $ location .~ out $ emptyFile
        customizedExpr' <- mapM (mapFileWithDefName prefix suffix convertCustom) customizedExpr

        return $ quantifications' ++ encodeExpr' ++ customizedExpr'

geneId2Name_ :: [Int]    -- ^ Fields to be kept in the results. 1-based indexing.
             -> Bool     -- ^ Whether the input contains header.
             -> M.HashMap B.ByteString (CI B.ByteString)
             -> FilePath -- ^ Input
             -> FilePath -- ^ Output
             -> IO ()
geneId2Name_ idx has_header id2Name input output = do
    c <- B.lines <$> B.readFile input
    let content = if has_header then tail c else c
        header = if has_header
            then let xs = B.split '\t' $ head c
                 in B.intercalate "\t" $ "name" : map (\i -> xs!!(i-1)) idx
            else B.intercalate "\t" $ "name" : map (\i -> B.pack $ "value_" ++ show i) idx
    B.writeFile output $ B.unlines $
        header : map (B.intercalate "\t" . convert . B.split '\t') content
  where
    convert xs = getGeneName (head xs) : map (\i -> xs!!(i-1)) idx
    getGeneName x = case M.lookup x id2Name of
        Just x' -> original x'
        Nothing -> if mk x `S.member` geneSet
            then x
            else x <> "(UNKNOWN_Taiji)"
    geneSet = S.fromList $ M.elems id2Name
{-# INLINE geneId2Name_ #-}

readExpr :: [RNASeq S (File '[GeneQuant] 'Tsv)]
         -> IO [(T.Text, [[(CI B.ByteString, Double)]])]
readExpr quantifications = forM quantifications $ \e -> do
    c <- B.readFile $ runIdentity (e^.replicates) ^. files.location
    let result = map (\xs ->
            let fs = B.split '\t' xs in (mk $ head fs, readDouble $ fs!!1)) $
            tail $ B.lines c
    return (fromJust $ e^.groupName, [result])
{-# INLINE readExpr #-}

-- | Combine RNA expression data into a table and output
mkTable :: RNASeqConfig config
        => [RNASeq S (File '[GeneQuant] 'Tsv)]
        -> WorkflowConfig config (Maybe (File '[] 'Tsv))
mkTable es = do
    results <- liftIO $ readExpr es
    if null results
        then return Nothing
        else do
            outdir <- asks _rnaseq_output_dir >>= getPath
            let output = outdir ++ "/" ++ "expression_profile.tsv"
                (expNames, values) = unzip $ M.toList $
                    fmap (map (second average) . combine) $ M.fromListWith (++) $ results
            liftIO $ B.writeFile output $ B.unlines $
                B.pack (T.unpack $ T.intercalate "\t" $ "Name" : expNames) :
                map (\(x,xs) -> B.intercalate "\t" $ original x : map toShortest xs)
                    (combine values)
            return $ Just $ location .~ output $ emptyFile

combine :: [[(CI B.ByteString, Double)]] -> [(CI B.ByteString, [Double])]
combine xs = flip map names $ \nm -> (nm, map (M.lookupDefault 0.01 nm) xs')
  where
    names = nubSort $ concatMap (fst . unzip) xs
    xs' = map (fmap average . M.fromListWith (++) . map (second return)) xs
{-# INLINE combine #-}

average :: [Double] -> Double
average [a,b]   = (a + b) / 2
average [a,b,c] = (a + b + c) / 3
average xs      = foldl1' (+) xs / fromIntegral (length xs)
{-# INLINE average #-}
