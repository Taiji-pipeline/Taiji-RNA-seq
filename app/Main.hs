{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.Utils
import           Control.Lens                  ((.=))
import           Data.Aeson                    (FromJSON, ToJSON)
import           Data.Default
import           GHC.Generics                  (Generic)
import           Scientific.Workflow
import           Taiji.Pipeline.RNASeq (builder)
import qualified Taiji.Pipeline.RNASeq.Config as C

data RNASeqOpts = RNASeqOpts
    { outputDir :: Directory
    , starIndex  :: Maybe FilePath
    , rsemIndex  :: Maybe FilePath
    , genome    :: Maybe FilePath
    , input     :: FilePath
    , annotation :: Maybe FilePath
    } deriving (Generic)

instance C.RNASeqConfig RNASeqOpts where
    _rnaseq_output_dir = outputDir
    _rnaseq_star_index = starIndex
    _rnaseq_rsem_index = rsemIndex
    _rnaseq_genome_fasta = genome
    _rnaseq_input = input
    _rnaseq_annotation = annotation

instance Default RNASeqOpts where
    def = RNASeqOpts
        { outputDir = asDir "output"
        , starIndex = Nothing
        , rsemIndex = Nothing
        , genome = Nothing
        , input = "input.yml"
        , annotation = Nothing
        }

instance FromJSON RNASeqOpts
instance ToJSON RNASeqOpts

-- | Instantiate the "ATACSeqConfig".
initialization :: () -> WorkflowConfig RNASeqOpts ()
initialization _ = return ()

mainWith defaultMainOpts { programHeader = "Taiji-RNA-Seq" } $ do
    nodeS "Initialization" 'initialization $ submitToRemote .= Just False
    ["Initialization"] ~> "Read_Input"
    builder
