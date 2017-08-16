{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.RNASeq (builder) where

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
import           Data.Either                   (lefts, rights)
import           Data.List                     (nub)
import           Data.Maybe                    (fromJust, fromMaybe, mapMaybe)
import           Data.Singletons               (SingI)
import           Scientific.Workflow

import           Taiji.Pipeline.RNASeq.Config
import           Taiji.Pipeline.RNASeq.Functions

builder :: Builder ()
builder = do
    nodeS "Make_Index" 'rnaMkIndex $ return ()
    nodeS "Read_Input" [| \_ -> do
        input <- asks _rnaseq_input
        liftIO $ readRNASeq input "RNA-seq"
        |] $ submitToRemote .= Just False
    nodeS "Download_Data" [| \input -> do
        dir <- asks _rnaseq_output_dir >>= getPath
        liftIO $ downloadData dir input
        |] $ submitToRemote .= Just False
    node' "Get_Fastq" 'getFastq $ submitToRemote .= Just False
    nodeSharedPS 1 "Align" [| \(ContextData idx input) -> do
        dir <- asks _rnaseq_output_dir >>= getPath
        liftIO $ starAlign (dir, "bam") idx (return ()) input
        |] $ return ()
    node' "Quant_Prepare" 'quantPrep $ submitToRemote .= Just False

    path ["Make_Index", "Read_Input", "Download_Data"]
    ["Make_Index", "Download_Data"] ~> "Get_Fastq"
    path ["Get_Fastq", "Align"]
    ["Make_Index", "Align"] ~> "Quant_Prepare"
