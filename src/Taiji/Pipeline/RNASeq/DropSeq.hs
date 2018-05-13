{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.RNASeq.DropSeq where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Control.Lens
import           Control.Monad.IO.Class                  (liftIO)
import           Control.Monad.Reader                    (asks)
import           Data.Either
import           Scientific.Workflow

import           Taiji.Pipeline.RNASeq.Common
import           Taiji.Pipeline.RNASeq.DropSeq.DeBarcode

builder :: Builder ()
builder = do
    nodeS "Read_Input" [| \_ -> do
        input <- asks _dropSeq_input
        liftIO $ if ".tsv" == reverse (take 4 $ reverse input)
            then readRNASeqTSV input "Drop-seq"
            else readRNASeq input "Drop-seq"
        |] $ submitToRemote .= Just False

    node' "Get_Fastq" [|
        map (\x -> x & replicates.mapped.files %~ fromRight undefined) .
        filter (\x -> isRight $ runIdentity (x^.replicates) ^. files ) .
        rnaGetFastq
        |] $ submitToRemote .= Just False
    nodePS 1 "Tag_Fastq" 'tagAndFilter $
        remoteParam .= "--mem=20000 -p gpu"
    nodeS "Barcode_Stat" [| \xs -> do
        let xs' = flip map xs $ \x -> x & replicates.mapped.files %~
                (\(_,b,c) -> (b,c))
        mapM countBarcodeBaseFreq xs'
        |] $ return ()
    path ["Read_Input", "Get_Fastq", "Tag_Fastq", "Barcode_Stat"]
