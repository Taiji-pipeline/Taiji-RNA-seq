{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.Utils
import           Control.Lens                  ((.=))
import           Data.Aeson                    (FromJSON, ToJSON)
import           Data.Default
import           Data.Maybe                    (fromJust)
import           GHC.Generics                  (Generic)
import           Scientific.Workflow

import Taiji.Pipeline.RNASeq

data RNASeqOpts = RNASeqOpts
    { output_dir     :: Directory
    , star_index     :: Maybe FilePath
    , rsem_index     :: Maybe FilePath
    , genome         :: Maybe FilePath
    , input          :: FilePath
    , annotation     :: Maybe FilePath
    , cellBarcodeLen :: Int
    , molBarcodeLen  :: Int
    } deriving (Generic)

instance FromJSON RNASeqOpts
instance ToJSON RNASeqOpts

instance Default RNASeqOpts where
    def = RNASeqOpts
        { output_dir = asDir "output"
        , star_index = Nothing
        , rsem_index = Nothing
        , genome = Nothing
        , input = "input.yml"
        , annotation = Nothing
        , cellBarcodeLen = 12
        , molBarcodeLen = 8
        }

instance RNASeqConfig RNASeqOpts where
    _rnaseq_output_dir = output_dir
    _rnaseq_star_index = star_index
    _rnaseq_rsem_index = rsem_index
    _rnaseq_genome_fasta = genome
    _rnaseq_input = input
    _rnaseq_annotation = annotation

mainWith defaultMainOpts
    { programHeader = "Taiji-RNA-Seq"
    , workflowConfigType = Just ''RNASeqOpts } $ do
        inputReader "RNA-seq"
        builder
