{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}

module Main where

import           Bio.Pipeline.Utils
import           Control.Lens                  ((.=))
import           Data.Aeson                    (FromJSON)
import           Data.Maybe                    (fromJust)
import Data.Binary (Binary)
import           GHC.Generics                  (Generic)

import           Control.Workflow
import qualified Control.Workflow.Coordinator.Drmaa as D
import Control.Workflow.Main
import Data.Proxy (Proxy(..))

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
instance Binary RNASeqOpts

instance RNASeqConfig RNASeqOpts where
    _rnaseq_output_dir = output_dir
    _rnaseq_star_index = star_index
    _rnaseq_rsem_index = rsem_index
    _rnaseq_genome_fasta = genome
    _rnaseq_input = input
    _rnaseq_annotation = annotation

decodeDrmaa :: String -> Int -> FilePath -> IO D.DrmaaConfig
decodeDrmaa ip port _ = D.getDefaultDrmaaConfig
    ["remote", "--ip", ip, "--port", show port]

build "wf" [t| SciFlow RNASeqOpts |] $ inputReader "RNA-seq" >> builder

main :: IO ()
main = defaultMain "" cmd wf
  where
    cmd = [ runParser decodeDrmaa
          , viewParser
          , remoteParser (Proxy :: Proxy D.Drmaa) ]