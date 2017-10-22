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
    , rnaGetFastq
    , rnaGetExpression
    , rnaGetExpressionENCODE
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
import           Data.Monoid                       ((<>))
import           Data.Promotion.Prelude.List       (Elem)
import           Data.Singletons                   (SingI)
import qualified Data.Text                         as T
import           Scientific.Workflow

import           Taiji.Pipeline.RNASeq.Config

type RNASeqWithSomeFile = RNASeq N [Either SomeFile (SomeFile, SomeFile)]

type RNASeqMaybePair tag1 tag2 filetype =
    Either (RNASeq S (File tag1 filetype))
           (RNASeq S (File tag2 filetype, File tag2 filetype))

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
         => Either (RNASeq S (SomeTags 'Fastq))
                   (RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq))
         -> WorkflowConfig config (
                Either (RNASeq S (File '[] 'Bam, File '[] 'Bam))
                       (RNASeq S (File '[Pairend] 'Bam, File '[Pairend] 'Bam)) )
rnaAlign input = do
    dir <- asks _rnaseq_output_dir >>= getPath
    idx <- asks (fromJust . _rnaseq_star_index)
    liftIO $ case input of
        Left e -> if (runIdentity (e^.replicates) ^. files) `hasTag` Gzip
            then starAlign (dir, ".bam") idx (starCores .= 4) (Left $
                e & replicates.traverse.files %~ fromSomeTags
                    :: Either (RNASeq S (File '[Gzip] 'Fastq))
                              (RNASeq S (File '[] 'Fastq, File '[] 'Fastq)) )
            else starAlign (dir, ".bam") idx (starCores .= 4) (Left $
                e & replicates.traverse.files %~ fromSomeTags
                    :: Either (RNASeq S (File '[] 'Fastq))
                              (RNASeq S (File '[] 'Fastq, File '[] 'Fastq)) )
        Right e -> do
            let (f_a, f_b) = runIdentity (e^.replicates) ^. files
            if f_a `hasTag` Gzip && f_b `hasTag` Gzip
                then starAlign (dir, ".bam") idx (starCores .= 4) ( Right $
                    e & replicates.traverse.files %~ bimap fromSomeTags fromSomeTags
                        :: Either (RNASeq S (File '[] 'Fastq))
                                  (RNASeq S (File '[Gzip] 'Fastq, File '[Gzip] 'Fastq)) )
                else starAlign (dir, ".bam") idx (starCores .= 4) ( Right $
                    e & replicates.traverse.files %~ bimap fromSomeTags fromSomeTags
                        :: Either (RNASeq S (File '[] 'Fastq))
                                  (RNASeq S (File '[] 'Fastq, File '[] 'Fastq)) )

rnaDownloadData :: RNASeqConfig config
                => [RNASeqWithSomeFile]
                -> WorkflowConfig config [RNASeqWithSomeFile]
rnaDownloadData dat = do
    dir <- asks _rnaseq_output_dir >>= getPath . (<> (asDir "/Download"))
    liftIO $ dat & traverse.replicates.traverse.files.traverse %%~ downloadFiles dir

rnaGetFastq :: [RNASeqWithSomeFile]
            -> [ Either (RNASeq S (SomeTags 'Fastq))
                        (RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq))
               ]
rnaGetFastq inputs = concatMap split $ concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap castFile (bimap castFile castFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq

rnaGetExpression :: [RNASeqWithSomeFile]
                 -> [RNASeq S (File '[GeneQuant] 'Tsv)]
rnaGetExpression inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map fromSomeFile $ filter (`hasTag` GeneQuant) $ lefts fls

rnaGetExpressionENCODE :: [RNASeqWithSomeFile]
                       -> [RNASeq S (File '[GeneQuant] 'Tsv)]
rnaGetExpressionENCODE inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map fromSomeFile $ filter
        (\x -> x `hasTag` GeneQuant && x `hasTag` ENCODE) $ lefts fls

quantification :: (RNASeqConfig config, SingI tags1, SingI tags2)
               => Either (RNASeq S (File tags1 'Bam))
                         (RNASeq S (File tags2 'Bam))
               -> WorkflowConfig config (RNASeq S
                    (File '[GeneQuant] 'Tsv, File '[TranscriptQuant] 'Tsv))
quantification input = do
    dir <- asks _rnaseq_output_dir >>= getPath
    idx <- asks (fromJust . _rnaseq_rsem_index)
    let fun :: SingI tag => RNASeq S (File tag 'Bam) -> _
        fun = rsemQuant dir idx (rsemCores .= 4)
    liftIO $ either fun fun input

-- | Retrieve gene names
geneId2Name :: (RNASeqConfig config, Elem 'GeneQuant tags ~ 'True)
            => [RNASeq S (File tags 'Tsv)]
            -> WorkflowConfig config [RNASeq S (File tags 'Tsv)]
geneId2Name experiments = do
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
        mapM (mapFileWithDefName (outdir ++ "/") "_gene_quant_TPM.tsv" fun) experiments

-- | Retrieve gene names and compute the average expression of replicates.
averageExpr :: (RNASeqConfig config, Elem 'GeneQuant tags ~ 'True)
            => RNASeq N (File tags 'Tsv)
            -> WorkflowConfig config (RNASeq S (File tags 'Tsv))
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
        return $ replicates .~ return (Replicate newFile [] 0) $ experiment

-- | Combine RNA expression data into a table and output
mkTable :: RNASeqConfig config
        => [RNASeq S (File tags 'Tsv)]
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
average [a,b]   = (a + b) / 2
average [a,b,c] = (a + b + c) / 3
average xs      = foldl1' (+) xs / fromIntegral (length xs)
{-# INLINE average #-}
