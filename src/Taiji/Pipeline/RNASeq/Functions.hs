{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.RNASeq.Functions where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Download
import           Bio.Pipeline.Report
import           Bio.Pipeline.NGS
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import qualified Data.Text as T
import           Control.Lens
import           Control.Monad.IO.Class        (liftIO)
import           Control.Monad.Reader          (asks)
import           Data.Bifunctor                (bimap)
import           Data.Bitraversable            (bitraverse)
import           Data.Coerce                   (coerce)
import           Data.Either                   (lefts, rights)
import           Data.List                     (nub)
import           Data.Maybe                    (fromJust, fromMaybe, mapMaybe)
import           Data.Singletons               (SingI)
import           Scientific.Workflow

import Taiji.Pipeline.RNASeq.Config

rnaMkIndex :: RNASeqConfig config => () -> WorkflowConfig config (FilePath, FilePath)
rnaMkIndex _ = do
    genome <- asks (fromJust . _rnaseq_genome_fasta)
    dir <- asks (fromJust . _rnaseq_star_index)
    anno <- asks (fromJust . _rnaseq_annotation)
    rsemIndex <- asks (fromJust . _rnaseq_rsem_index)
    liftIO $ (,) <$> starMkIndex "star" dir [genome] anno 100 <*>
        rsemMkIndex rsemIndex anno [genome]

downloadData :: FilePath
             -> [RNASeq [Either SomeFile (SomeFile, SomeFile)]]
             -> IO [RNASeq [Either SomeFile (SomeFile, SomeFile)]]
downloadData dir = traverse.replicates.traverse.files.traverse %%~ download
  where
    download :: Either SomeFile (SomeFile, SomeFile)
             -> IO (Either SomeFile (SomeFile, SomeFile))
    download input@(Left (SomeFile fl)) = if getFileType fl == SRA
        then bimap SomeFile (bimap SomeFile SomeFile) <$>
                sraToFastq dir (coerce fl :: File '[] 'SRA)
        else return input
    download x = return x

getFastq :: ((FilePath, FilePath), [RNASeq [Either SomeFile (SomeFile, SomeFile)]])
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
