{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.RNASeq (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.NGS
import           Data.Bifunctor                    (bimap)
import           Control.Lens
import           Control.Monad.IO.Class          (liftIO)
import           Control.Monad.Reader            (asks)
import           Scientific.Workflow

import           Taiji.Pipeline.RNASeq.Config
import           Taiji.Pipeline.RNASeq.Functions

builder :: Builder ()
builder = do
    nodeS "Read_Input" [| \_ -> do
        input <- asks _rnaseq_input
        liftIO $ readRNASeq input "RNA-seq"
        |] $ submitToRemote .= Just False
    nodeS "Download_Data" 'rnaDownloadData $ submitToRemote .= Just False
    node' "Get_Fastq" 'getFastq $ submitToRemote .= Just False
    nodeS "Make_Index" 'rnaMkIndex $ return ()
    nodePS 1 "Align" 'rnaAlign $ return ()
    nodePS 1 "Quant" [| \input -> do
        let fun x = x & replicates.mapped.files %~ fst
        quantification $ bimap fun fun input
        |] $ return ()
    nodePS 1 "Convert_ID_To_Name" [| \input ->
        geneId2Name (input & replicates.mapped.files %~ fst)
        |] $ return ()
    node' "Average_Prep" [| \input ->
        let fun [x] = x
            fun _ = error "Found multiple files"
        in (mergeExps input) & mapped.replicates.mapped.files %~ fun
        |] $ submitToRemote .= Just False
    nodePS 1 "Average" 'averageExpr $ return ()
    nodeS "Make_Expr_Table" 'mkTable $ return ()

    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index", "Align",
        "Quant", "Convert_ID_To_Name", "Average_Prep", "Average", "Make_Expr_Table"]
