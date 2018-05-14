{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}

module Taiji.Pipeline.RNASeq.DropSeq.QC
    (barcodeStat) where

import           Bio.Data.Bam
import           Bio.Data.Experiment
import           Bio.HTS                              (qName)
import           Bio.Pipeline.Utils
import           Bio.Utils.Misc                       (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import Data.Ord (comparing)
import qualified Data.ByteString.Char8                as B
import qualified Data.HashMap.Strict                  as M
import           Data.List
import           Data.Monoid                          ((<>))
import qualified Data.Text                            as T
import           Scientific.Workflow
import           Text.Printf                          (printf)

import           Taiji.Pipeline.RNASeq.DropSeq.Config

barcodeStat :: DropSeqConfig config
            => RNASeq S (File tags 'Bam)
            -> WorkflowConfig config (RNASeq S (File '[] 'Tsv))
barcodeStat input = do
    dir <- asks _dropSeq_output_dir >>= getPath
    let f x = do
            B.writeFile output $ B.unlines $
                map (\(a,b) -> B.pack (show a) <> "\t" <> B.pack (show b)) $
                sortBy (flip (comparing snd)) $ M.toList x
            return $ location .~ output $ emptyFile
        output = printf "%s/%s_rep%d_barcode_stat.txt" dir (T.unpack $ input^.eid)
            (runIdentity (input^.replicates) ^. number)
    input & replicates.traverse.files %%~ (\x -> liftIO $ fun x >>= f)
  where
    fun fl = withBamFile (fl^.location) $ \h -> runConduit $
        readBam h .| foldlC f M.empty
      where
        f m x = let [cellBc, _, _] = B.split ':' $ qName x
                in M.insertWith (+) (readInt cellBc) (1::Int) m
