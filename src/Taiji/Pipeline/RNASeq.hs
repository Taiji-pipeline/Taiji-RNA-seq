{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Pipeline.RNASeq (builder) where

import           Scientific.Workflow

import qualified Taiji.Pipeline.RNASeq.Classic as Classic
import qualified Taiji.Pipeline.RNASeq.DropSeq as DropSeq

builder :: Builder ()
builder = do
    namespace "RNASeq" Classic.builder
    namespace "DropSeq" DropSeq.builder
