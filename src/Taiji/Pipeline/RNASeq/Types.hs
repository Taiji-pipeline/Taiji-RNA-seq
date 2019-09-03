{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE LambdaCase #-}
module Taiji.Pipeline.RNASeq.Types where

import qualified Data.Text as T
import Shelly hiding (FilePath)

import Taiji.Prelude

class RNASeqConfig config where
    _rnaseq_assembly :: config -> Maybe String
    _rnaseq_genome_fasta :: config -> Maybe FilePath
    _rnaseq_star_index :: config -> Maybe FilePath
    _rnaseq_annotation :: config -> Maybe FilePath
    _rnaseq_rsem_index :: config -> Maybe FilePath
    _rnaseq_input :: config -> FilePath
    _rnaseq_output_dir :: config -> Directory

getGenomeFasta :: RNASeqConfig config => ReaderT config IO FilePath
getGenomeFasta = asks _rnaseq_genome_fasta >>= \case
    Nothing -> error "genome fasta is missing"
    Just fasta -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack fasta 
       if exist
           then return fasta
           else asks _rnaseq_assembly >>= \case
               Nothing -> error "genome fasta is missing"
               Just assembly -> do
                   liftIO $ fetchGenome fasta assembly
                   return fasta

getAnnotation :: RNASeqConfig config => ReaderT config IO FilePath
getAnnotation = asks _rnaseq_annotation >>= \case
    Nothing -> error "annotation is missing"
    Just anno -> do
       exist <- liftIO $ shelly $ test_f $ fromText $ T.pack anno
       if exist
           then return anno
           else asks _rnaseq_assembly >>= \case
               Nothing -> error "annotation is missing"
               Just assembly -> do
                   liftIO $ fetchAnnotation anno assembly
                   return anno