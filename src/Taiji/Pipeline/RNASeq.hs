{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Pipeline.RNASeq
    ( inputReader
    , builder
    , RNASeqConfig(..)
    ) where

import           Bio.Data.Experiment.Parser
import           Data.Bifunctor                          (bimap)
import qualified Data.Text                               as T
import           Control.Workflow
import qualified Data.IntMap.Strict as M

import           Taiji.Pipeline.RNASeq.Types
import           Taiji.Pipeline.RNASeq.Functions
import Taiji.Prelude

inputReader :: String    -- ^ The key
            -> Builder ()
inputReader key = do
    node "Read_Input" [| \_ -> do
        input <- asks _rnaseq_input
        es <- liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readRNASeqTSV input key
            else readRNASeq input key
        forM_ es $ \e -> forM_ (M.toList $ e^.replicates) $ \(i, rep) ->
            when (length (rep^.files) > 1) $ error $ printf
                "replicate %d in \"%s\" contains more than 2 files"
                i (T.unpack $ e^.eid)
        return es
        |] $ doc .= "Read RNA-seq data information from input file."
    ["Read_Input"] ~> "Download_Data"

builder :: Builder ()
builder = do
    nodePar "Download_Data" 'rnaDownloadData $ return ()
    node "Get_Fastq" [| return . rnaGetFastq |] $ return ()
    node "Make_Index" 'rnaMkIndex $ memory .= 40
    nodePar "Align" 'rnaAlign $ do
        nCore .= 4
        memory .= 40
    nodePar "Quant" [| \input -> do
        quantification $ input & replicates.mapped.files %~
            bimap (fromJust . snd) (fromJust . snd)
        |] $ do
        nCore .= 4
        memory .= 40
    path ["Download_Data", "Get_Fastq", "Make_Index", "Align",
        "Quant"]

    node "Convert_ID_To_Name" [| \(ori, input) -> do
        let input' = input & mapped.replicates.mapped.files %~ fst
        geneId2Name (ori, input')
        |] $ return ()
    ["Download_Data", "Quant"] ~> "Convert_ID_To_Name"
    node "Make_Expr_Table" 'mkTable $ return ()
    ["Convert_ID_To_Name"] ~> "Make_Expr_Table"