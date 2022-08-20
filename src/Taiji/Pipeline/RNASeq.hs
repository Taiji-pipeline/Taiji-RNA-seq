{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Pipeline.RNASeq
    ( builder
    , RNASeqConfig(..)
    ) where

import           Bio.Data.Experiment.Parser
import           Data.Bifunctor                          (bimap)
import           Control.Workflow
import Bio.Data.Experiment.Types

import           Taiji.Pipeline.RNASeq.Types
import           Taiji.Pipeline.RNASeq.Functions
import Taiji.Prelude

builder :: Builder ()
builder = do
    node "Read_Input" [| \() -> do
        input <- asks _rnaseq_input
        liftIO $ mkInputReader input "RNA-seq" (\_ x -> RNASeq x)
        |] $ doc .= "Read input data information."
    nodePar "Download_Data" [| rnaDownloadData |] $ return ()
    uNode "Get_Fastq" [| return . rnaGetFastq |]
    node "Make_Index" [| rnaMkIndex |] $ memory .= 40
    nodePar "Align" [| rnaAlign |] $ do
        nCore .= 4
        memory .= 40
    nodePar "Quant" [| \input -> do
        quantification $ input & replicates.mapped.files %~
            bimap (fromJust . snd) (fromJust . snd)
        |] $ do
        nCore .= 4
        memory .= 40
    path ["Read_Input", "Download_Data", "Get_Fastq", "Make_Index", "Align",
        "Quant"]

    node "Convert_ID_To_Name" [| \(ori, input) -> do
        let input' = input & mapped.replicates.mapped.files %~ fst
        geneId2Name (ori, input')
        |] $ return ()
    ["Download_Data", "Quant"] ~> "Convert_ID_To_Name"
    node "Make_Expr_Table" [| mkTable |] $ return ()
    ["Convert_ID_To_Name"] ~> "Make_Expr_Table"