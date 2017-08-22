{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE GADTs                 #-}
{-# LANGUAGE OverloadedLists       #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE TemplateHaskell       #-}
module Taiji.Pipeline.RNASeq.Functions
    ( rnaMkIndex
    , rnaDownloadData
    , rnaAlign
    , getFastq
    , quantPrep
    , quantification
    , geneId2Name
    , averageExpr
    , mkTable
    ) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Download
import           Bio.Pipeline.NGS
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Report
import           Bio.Pipeline.Utils
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc                    (readDouble)
import           Control.Lens
import           Control.Monad                     (forM)
import           Control.Monad.IO.Class            (liftIO)
import           Control.Monad.Reader              (asks)
import           Data.Bifunctor                    (bimap, second)
import           Data.Bitraversable                (bitraverse)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Coerce                       (coerce)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Either                       (lefts, rights)
import qualified Data.HashMap.Strict               as M
import           Data.List
import           Data.List.Ordered                 (nubSort)
import           Data.Maybe                        (fromJust, fromMaybe,
                                                    mapMaybe)
import           Data.Promotion.Prelude.List       (Elem)
import           Data.Singletons                   (SingI)
import qualified Data.Text                         as T
import           Scientific.Workflow

import           Taiji.Pipeline.RNASeq.Config

type RNASeqWithSomeFile = RNASeq [Either SomeFile (SomeFile, SomeFile)]

rnaMkIndex :: RNASeqConfig config => () -> WorkflowConfig config (FilePath, FilePath)
rnaMkIndex _ = do
    genome <- asks (fromJust . _rnaseq_genome_fasta)
    dir <- asks (fromJust . _rnaseq_star_index)
    anno <- asks (fromJust . _rnaseq_annotation)
    rsemIndex <- asks (fromJust . _rnaseq_rsem_index)
    liftIO $ (,) <$> starMkIndex "star" dir [genome] anno 100 <*>
        rsemMkIndex rsemIndex anno [genome]

rnaAlign :: RNASeqConfig config
         => ContextData FilePath (MaybePairExp RNASeq '[] '[Pairend] 'Fastq)
         -> WorkflowConfig config (
                Either (RNASeq (File '[] 'Bam, File '[] 'Bam))
                       (RNASeq (File '[Pairend] 'Bam, File '[Pairend] 'Bam)) )
rnaAlign (ContextData idx input) = do
    dir <- asks _rnaseq_output_dir >>= getPath
    liftIO $ starAlign (dir, "bam") idx (return ()) input

rnaDownloadData :: RNASeqConfig config
                => [RNASeqWithSomeFile]
                -> WorkflowConfig config [RNASeqWithSomeFile]
rnaDownloadData dat = do
    dir <- asks _rnaseq_output_dir >>= getPath
    liftIO $ dat & traverse.replicates.traverse.files.traverse %%~ download dir
  where
    download dir input@(Left (SomeFile fl)) = if getFileType fl == SRA
        then bimap SomeFile (bimap SomeFile SomeFile) <$>
                sraToFastq dir (coerce fl :: File '[] 'SRA)
        else return input
    download _ x = return x

getFastq :: ((FilePath, FilePath), [RNASeqWithSomeFile])
         -> ContextData FilePath [ MaybePairExp RNASeq '[] '[Pairend] 'Fastq ]
getFastq ((idx,_), inputs) = ContextData idx $ flip concatMap inputs $ \input ->
    fromMaybe (error "A mix of single and pairend fastq was found") $
        splitExpByFileEither $ input & replicates.mapped.files %~ f
  where
    f :: [Either SomeFile (SomeFile, SomeFile)] -> [MaybePair '[] '[Pairend] 'Fastq]
    f fls = map (bimap fromSomeFile (bimap fromSomeFile fromSomeFile)) $
        filter (either (`someFileIs` Fastq) g) fls
      where
        g (x,y) = x `someFileIs` Fastq && y `someFileIs` Fastq

quantPrep :: ( (FilePath, FilePath)
             , [ Either (RNASeq (File tags1 'Bam, File tags1 'Bam))
                        (RNASeq (File tags2 'Bam, File tags2 'Bam)) ] )
          -> ContextData FilePath [ Either (RNASeq (File tags1 'Bam))
                                           (RNASeq (File tags2 'Bam)) ]
quantPrep ((_, idx), input) = ContextData idx $ map (bimap fun fun) input
  where
    fun x = x & replicates.mapped.files %~ fst

quantification :: (RNASeqConfig config, SingI tags1, SingI tags2)
               => ContextData FilePath ( Either (RNASeq (File tags1 'Bam))
                                           (RNASeq (File tags2 'Bam)) )
               -> WorkflowConfig config (RNASeq (File '[GeneQuant] 'Tsv, File '[TranscriptQuant] 'Tsv))
quantification (ContextData idx input) = do
    dir <- asks _rnaseq_output_dir >>= getPath
    let fun :: SingI tag => RNASeq (File tag 'Bam) -> _
        fun = rsemQuant dir idx (return ())
    liftIO $ either fun fun input

-- | Retrieve gene names
geneId2Name :: (RNASeqConfig config, Elem 'GeneQuant tags ~ 'True)
            => RNASeq (File tags 'Tsv)
            -> WorkflowConfig config (RNASeq (File tags 'Tsv))
geneId2Name experiment = do
    outdir <- asks _rnaseq_output_dir >>= getPath
    anno <- asks (fromJust . _rnaseq_annotation)
    liftIO $ do
        id2Name <- fmap (M.fromList . map (\x -> (geneId x, original $ geneName x))) $
            readGenes' anno
        let fun output fl = do
                c <- B.readFile $ fl^.location
                B.writeFile output $ B.unlines $ map ( (\xs -> B.intercalate "\t"
                    [M.lookupDefault (head xs `B.append` "(inconvertible)") (head xs)
                    id2Name, xs!!4]) . B.split '\t' ) $ tail $ B.lines c
                return $ location .~ output $ emptyFile
        nameWith outdir "_gene_quant_TPM.tsv" fun experiment

-- | Retrieve gene names and compute the average expression of replicates.
averageExpr :: (RNASeqConfig config, Elem 'GeneQuant tags ~ 'True)
            => RNASeq (File tags 'Tsv)
            -> WorkflowConfig config (RNASeq (File tags 'Tsv))
averageExpr experiment = do
    outdir <- asks _rnaseq_output_dir >>= getPath
    let fls = experiment^..replicates.folded.files
        output = outdir ++ "/" ++ T.unpack (experiment^.eid) ++
            "_average_gene_quant.tsv"
        newFile = location .~ output $ emptyFile
        readExpr fl = do
            c <- B.readFile fl
            return $ map (\xs -> let [a,b] = B.split '\t' xs in (mk a, readDouble b)) $
                B.lines c
    liftIO $ do
        expr <- mapM (readExpr . (^.location)) fls
        B.writeFile output $ B.unlines $
            map (\(a,b) -> B.intercalate "\t" [original a, toShortest $ average b]) $
            combine expr
        return $ replicates .~ [Replicate newFile [] 0] $ experiment

-- | Combine RNA expression data into a table and output
mkTable :: RNASeqConfig config
        => [RNASeq (File tags 'Tsv)]
        -> WorkflowConfig config (Maybe (File '[] 'Tsv))
mkTable es
    | null es = return Nothing
    | otherwise = do
        outdir <- asks _rnaseq_output_dir >>= getPath
        let output = outdir ++ "/" ++ "expression_profile.tsv"
        liftIO $ do
            dat <- forM es $ \e -> do
                let [fl] = e^..replicates.folded.files
                expr <- readExpr $ fl^.location
                return (fromJust $ e^.groupName, expr)
            let (expNames, values) = unzip dat
            B.writeFile output $ B.unlines $
                B.pack (T.unpack $ T.intercalate "\t" $ "Name" : expNames) :
                map (\(x,xs) -> B.intercalate "\t" $ original x : map toShortest xs)
                    (combine values)
            return $ Just $ location .~ output $ emptyFile
  where
    readExpr fl = do
        c <- B.readFile fl
        return $ map (\xs -> let [a,b] = B.split '\t' xs in (mk a, readDouble b)) $
            B.lines c

combine :: [[(CI B.ByteString, Double)]] -> [(CI B.ByteString, [Double])]
combine xs = flip map names $ \nm -> (nm, map (M.lookupDefault 0.01 nm) xs')
  where
    names = nubSort $ concatMap (fst . unzip) xs
    xs' = map (fmap average . M.fromListWith (++) . map (second return)) xs
{-# INLINE combine #-}

average :: [Double] -> Double
average [a]     = a
average [a,b]   = (a + b) / 2
average [a,b,c] = (a + b + c) / 3
average xs      = foldl1' (+) xs / fromIntegral (length xs)
{-# INLINE average #-}
