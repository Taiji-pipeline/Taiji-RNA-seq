module Taiji.Pipeline.RNASeq.Config where

import           Bio.Pipeline.Utils

class RNASeqConfig config where
    _rnaseq_genome_fasta :: config -> Maybe FilePath
    _rnaseq_star_index :: config -> Maybe FilePath
    _rnaseq_annotation :: config -> Maybe FilePath
    _rnaseq_rsem_index :: config -> Maybe FilePath
    _rnaseq_input :: config -> FilePath
    _rnaseq_output_dir :: config -> Directory
