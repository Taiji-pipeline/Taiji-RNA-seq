{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.RNASeq (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Bio.Pipeline.NGS
import           Control.Lens
import           Control.Monad.IO.Class          (liftIO)
import           Control.Monad.Reader            (asks)
import           Data.Bifunctor                  (bimap)
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
    node' "Get_Fastq" 'rnaGetFastq $ submitToRemote .= Just False
    nodeS "Make_Index" 'rnaMkIndex $ remoteParam .= "--mem=40000 -p gpu"
    nodePS 1 "Align" 'rnaAlign $
        remoteParam .= "--ntasks-per-node=4 --mem=40000 -p gpu"
    nodePS 1 "Quant" [| \input -> do
        let fun x = x & replicates.mapped.files %~ snd
        quantification $ bimap fun fun input
        |] $
        remoteParam .= "--ntasks-per-node=4 --mem=40000 -p gpu"
    nodeS "Convert_ID_To_Name" [| \input -> geneId2Name $
        input & mapped.replicates.mapped.files %~ fst
        |] $ return ()
    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index", "Align",
        "Quant", "Convert_ID_To_Name"]
    nodeS "Make_Expr_Table" 'mkTable $ return ()
    ["Download_Data", "Convert_ID_To_Name"] ~> "Make_Expr_Table"
