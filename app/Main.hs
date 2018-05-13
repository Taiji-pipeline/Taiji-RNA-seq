{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE QuasiQuotes #-}

module Main where

import           Bio.Pipeline.Utils
import           Control.Lens                         ((.=))
import           Data.Aeson                           (FromJSON, ToJSON)
import           Data.Default
import           GHC.Generics                         (Generic)
import           Scientific.Workflow

import           Taiji.Pipeline.RNASeq                (builder)
import           Taiji.Pipeline.RNASeq.Classic.Config
import           Taiji.Pipeline.RNASeq.DropSeq.Config

data RNASeqOpts = RNASeqOpts
    { outputDir      :: Directory
    , starIndex      :: Maybe FilePath
    , rsemIndex      :: Maybe FilePath
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
        { outputDir = asDir "output"
        , starIndex = Nothing
        , rsemIndex = Nothing
        , genome = Nothing
        , input = "input.yml"
        , annotation = Nothing
        , cellBarcodeLen = 12
        , molBarcodeLen = 8
        }

instance RNASeqConfig RNASeqOpts where
    _rnaseq_output_dir = outputDir
    _rnaseq_star_index = starIndex
    _rnaseq_rsem_index = rsemIndex
    _rnaseq_genome_fasta = genome
    _rnaseq_input = input
    _rnaseq_annotation = annotation

instance DropSeqConfig RNASeqOpts where
    _dropSeq_input = input
    _dropSeq_output_dir = outputDir
    _dropSeq_cell_barcode_length = cellBarcodeLen
    _dropSeq_molecular_barcode_length = molBarcodeLen

-- | Instantiate the "ATACSeqConfig".
initialization :: () -> WorkflowConfig RNASeqOpts ()
initialization _ = return ()

mainWith defaultMainOpts
    { programHeader = "Taiji-RNA-Seq", workflowConfigType = Just ''RNASeqOpts }
    builder
