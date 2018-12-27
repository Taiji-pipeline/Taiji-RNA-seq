{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}

module Taiji.Pipeline.RNASeq.DropSeq.QC
    ( barcodeStat
    , barcodeCorrect
    , hasNMismatch
    ) where

import           Bio.Data.Bam
import           Bio.Data.Experiment
import           Bio.HTS                              (qName)
import           Bio.Pipeline.Utils
import           Bio.Utils.Misc                       (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import qualified Data.ByteString.Char8                as B
import           Data.Either                          (either)
import qualified Data.HashMap.Strict                  as M
import           Data.List
import           Data.Ord                             (comparing)
import           Data.Serialize                       (decode, encode)
import qualified Data.Text                            as T
import           IGraph
import           Scientific.Workflow
import           Text.Printf                          (printf)
import IGraph.Exporter.GEXF

import           Taiji.Pipeline.RNASeq.DropSeq.Config

type Barcode = Int

barcodeCorrect :: RNASeq S (File '[] 'Other)
               -> IO ()
barcodeCorrect input = do
    m <- fmap (either error id . decode) $ B.readFile $ input^.replicates._2.files.location :: IO ([(Int, Int)])
    let gr = mkBarCodeGraph m
    writeGEXF "out.gexf" $ emap (\_ -> defaultEdgeAttributes) $ nmap (\(_,x) -> defaultNodeAttributes{ _nodeLabel = show $ snd x }) gr

mkBarCodeGraph :: [(Int, Int)] -> Graph 'U (Int, Int) ()
mkBarCodeGraph bcs = fromLabeledEdges $ zip es $ repeat ()
  where
    es = filter (\(a, b) -> not $ hasNMismatch 2 (fst a) (fst b)) $ comb bcs
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb _      = []

{-
collapseBarcode :: Graph 'U (Barcode, Int) -> [(Barcode, [Barcode])]
collapseBarcode gr = do
    let x = maximumBy (length . neighbors) $ nodes gr
-}

-- | Whether two sequences have at least N mismatches.
hasNMismatch :: Int
             -> Int   -- seq 1
             -> Int   -- seq 2
             -> Bool
hasNMismatch n x y = go 0 x y
  where
    go acc a b
        | a == 0 && b == 0 = False
        | otherwise = if a `mod` base == b `mod` base
            then go acc (a `div` base) (b `div` base)
            else if acc + 1 >= n
                then True
                else go (acc+1) (a `div` base) (b `div` base)
    base = 5

barcodeStat :: DropSeqConfig config
            => RNASeq S (File tags 'Bam)
            -> WorkflowConfig config (RNASeq S (File '[] 'Other))
barcodeStat input = do
    dir <- asks _dropSeq_output_dir >>= getPath
    let output = printf "%s/%s_rep%d_barcode_stat.bin" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ ( \x -> liftIO $ do
        r <- fun x
        B.writeFile output $ encode $ sortBy (flip (comparing snd)) $ M.toList r
        return $ location .~ output $ emptyFile )
  where
    fun fl = withBamFile (fl^.location) $ \h -> runConduit $
        readBam h .| foldlC f M.empty
      where
        f m x = let [cellBc, _, _] = B.split ':' $ qName x
                in M.insertWith (+) (readInt cellBc) (1::Int) m
