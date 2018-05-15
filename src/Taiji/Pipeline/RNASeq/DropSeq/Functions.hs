{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Pipeline.RNASeq.DropSeq.Functions
    ( dropSeqMkIndex
    , dropSeqAlign
    , filterSortBam
    ) where

import           Bio.Data.Experiment
import           Bio.Pipeline.NGS.STAR
import           Bio.Pipeline.NGS.Utils
import           Bio.Pipeline.Utils
import           Control.Lens
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import           Data.Either                          (fromLeft)
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
            (runIdentity (input^.replicates) ^. number)
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
            (runIdentity (input^.replicates) ^. number)
        f fl = withTempFile "./" "tmp_file." $ \tmp _ -> do
            filterBam "./" tmp fl >>= sortBamByName "./" output
    input & replicates.traverse.files %%~ liftIO . f
