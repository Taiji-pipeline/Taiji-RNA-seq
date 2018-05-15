{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.RNASeq.DropSeq
    ( builder
    , DropSeqConfig(..)
    ) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser
import           Control.Lens
import           Control.Monad.IO.Class                  (liftIO)
import           Control.Monad.Reader                    (asks)
import           Data.Either
import           Scientific.Workflow

import           Taiji.Pipeline.RNASeq.Common
import           Taiji.Pipeline.RNASeq.DropSeq.Config
import           Taiji.Pipeline.RNASeq.DropSeq.DeBarcode
import           Taiji.Pipeline.RNASeq.DropSeq.Functions
import           Taiji.Pipeline.RNASeq.DropSeq.QC

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
    nodePS 1 "Tag_Fastq" 'deBarcode $ remoteParam .= "--mem=20000 -p gpu"
    nodeS "Barcode_Stat" [| \xs -> do
        let xs' = flip map xs $ \x -> x & replicates.mapped.files %~
                (\(_,b,c) -> (b,c))
        mapM countBarcodeBaseFreq xs'
        |] $ return ()
    path ["Read_Input", "Get_Fastq", "Tag_Fastq", "Barcode_Stat"]

    nodeS "Make_Index" 'dropSeqMkIndex $ remoteParam .= "--mem=40000 -p gpu"
    nodePS 1 "Align" [| \x -> dropSeqAlign $ x & replicates.mapped.files %~
        (\(a,_,_) -> a)
        |] $ remoteParam .= "--ntasks-per-node=4 --mem=40000 -p gpu"
    nodePS 1 "Filter_Bam" 'filterSortBam $ return ()
    nodePS 1 "Barcode_Stat_Aligned" 'barcodeStat $ return ()
    path ["Tag_Fastq", "Make_Index", "Align", "Filter_Bam", "Barcode_Stat_Aligned"]

    node' "Quantification_Prep" [| \(x, y) -> zipExp x y |] $ submitToRemote .= Just False
    nodePS 1 "Quantification" 'dropSeqQuantification $ return ()
    ["Filter_Bam", "Barcode_Stat_Aligned"] ~> "Quantification_Prep"
    path ["Quantification_Prep", "Quantification"]
