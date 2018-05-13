{-# LANGUAGE DataKinds #-}
module Taiji.Pipeline.RNASeq.Common
    ( rnaGetFastq
    ) where

import           Bio.Data.Experiment
import           Control.Lens

rnaGetFastq :: [ RNASeq N [Either SomeFile (SomeFile, SomeFile)] ]
            -> [ RNASeq S ( Either (SomeTags 'Fastq)
                                   (SomeTags 'Fastq, SomeTags 'Fastq) )
               ]
rnaGetFastq inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map (bimap castFile (bimap castFile castFile)) $
        filter (either (\x -> getFileType x == Fastq) g) fls
      where
        g (x,y) = getFileType x == Fastq && getFileType y == Fastq
