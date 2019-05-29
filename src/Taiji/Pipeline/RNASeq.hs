{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Pipeline.RNASeq
    ( inputReader
    , builder
    , RNASeqConfig(..)
    ) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Control.Lens
import           Control.Monad                           (forM_, when)
import           Control.Monad.IO.Class                  (liftIO)
import           Control.Monad.Reader                    (asks)
import           Data.Bifunctor                          (bimap)
import           Data.Maybe                              (fromJust)
import qualified Data.Text                               as T
import           Control.Workflow
import           Text.Printf                             (printf)

import           Taiji.Pipeline.RNASeq.Config
import           Taiji.Pipeline.RNASeq.Functions

inputReader :: String    -- ^ The key
            -> Builder ()
inputReader key = do
    node "Read_Input" [| \_ -> do
        input <- asks _rnaseq_input
        es <- liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readRNASeqTSV input key
            else readRNASeq input key
        forM_ es $ \e -> e & replicates.itraversed<.files %%@~ ( \i fls ->
            when (length fls > 1) $ error $ printf
                "replicate %d in \"%s\" contains more than 2 files"
                i (T.unpack $ e^.eid)
            )
        return es
        |] $ doc .= "Read RNA-seq data information from input file."
    ["Read_Input"] ~> "Download_Data"

builder :: Builder ()
builder = do
    node "Download_Data" 'rnaDownloadData $ return ()
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