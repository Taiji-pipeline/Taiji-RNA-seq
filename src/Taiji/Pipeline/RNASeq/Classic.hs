{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Taiji.Pipeline.RNASeq.Classic
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
import           Scientific.Workflow
import           Text.Printf                             (printf)

import           Taiji.Pipeline.RNASeq.Classic.Config
import           Taiji.Pipeline.RNASeq.Classic.Functions
import           Taiji.Pipeline.RNASeq.Common

inputReader :: String    -- ^ The key
            -> Builder ()
inputReader key = do
    nodeS "Read_Input" [| \_ -> do
        input <- asks _rnaseq_input
        es <- liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readRNASeqTSV input key
            else readRNASeq input key
        forM_ es $ \e -> e & replicates.itraversed<.files %%@~ ( \i fls ->
            when (length fls > 1) $ error $ printf
                "replicate %d in \"%s\" contains more than 2 files"
                (i :: Int) (T.unpack $ e^.eid)
            )
        return es
        |] $ do
            submitToRemote .= Just False
            note .= "Read RNA-seq data information from input file."
    ["Read_Input"] ~> "Download_Data"

builder :: Builder ()
builder = do
    nodeS "Download_Data" 'rnaDownloadData $ submitToRemote .= Just False
    node' "Get_Fastq" 'rnaGetFastq $ submitToRemote .= Just False
    nodeS "Make_Index" 'rnaMkIndex $ do
        remoteParam .= "--mem=40000"  -- slurm
        -- remoteParam .= "-l vmem=40G"  -- sge
    nodePS 1 "Align" 'rnaAlign $
        remoteParam .= "--ntasks-per-node=4 --mem=40000"  -- slurm
        -- remoteParam .= "-l vmem=10G -pe smp 4"  -- sge
    nodePS 1 "Quant" [| \input -> do
        quantification $ input & replicates.mapped.files %~
            bimap (fromJust . snd) (fromJust . snd)
        |] $
        remoteParam .= "--ntasks-per-node=4 --mem=40000"  -- slurm
        -- remoteParam .= "-l vmem=10G -pe smp 4"  -- sge
    path ["Download_Data", "Get_Fastq", "Make_Index", "Align",
        "Quant"]

    nodeS "Convert_ID_To_Name" [| \(ori, input) -> do
        let input' = input & mapped.replicates.mapped.files %~ fst
        geneId2Name (ori, input')
        |] $ return ()
    ["Download_Data", "Quant"] ~> "Convert_ID_To_Name"
    nodeS "Make_Expr_Table" 'mkTable $ return ()
    ["Convert_ID_To_Name"] ~> "Make_Expr_Table"
