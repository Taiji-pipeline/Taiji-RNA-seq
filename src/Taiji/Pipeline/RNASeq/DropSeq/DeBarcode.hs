{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Pipeline.RNASeq.DropSeq.DeBarcode
    ( deBarcode
    , countBarcodeBaseFreq
    , dnaToInt
    , intToDna
    ) where

import           Bio.Data.Experiment
import           Bio.Pipeline.Utils
import           Bio.Utils.Misc                       (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad
import           Control.Monad.IO.Class               (liftIO)
import           Control.Monad.Reader                 (asks)
import qualified Data.ByteString                      as BS
import qualified Data.ByteString.Char8                as B
import           Data.Conduit.Internal                (zipSources)
import           Data.Conduit.Zlib                    (gzip, multiple, ungzip)
import qualified Data.HashMap.Strict                  as M
import           Data.List                            (sortBy)
import qualified Data.Matrix.Unboxed                  as MU
import qualified Data.Matrix.Unboxed.Mutable          as MUM
import           Data.Monoid                          ((<>))
import           Data.Ord                             (comparing)
import           Data.Singletons
import qualified Data.Text                            as T
import           Scientific.Workflow
import           Text.Printf

import           Taiji.Pipeline.RNASeq.DropSeq.Config

countBarcodeBaseFreq :: DropSeqConfig config
                     => RNASeq S (File tags1 'Other, File tags2 'Other)
                     -> WorkflowConfig config (T.Text, ([[Int]], [[Int]]))
countBarcodeBaseFreq e = do
    lenCellBc <- asks _dropSeq_cell_barcode_length
    lenMolBc <- asks _dropSeq_molecular_barcode_length
    liftIO $ do
        mat1 <- runResourceT $ runConduit $ sourceFile (cellBc^.location) .|
            linesUnboundedAsciiC .| mapC f .| posBaseCount lenCellBc
        mat2 <- runResourceT $ runConduit $ sourceFile (molBc^.location) .|
            linesUnboundedAsciiC .| mapC f .| posBaseCount lenMolBc
        return (e^.eid, (MU.toLists mat1, MU.toLists mat2))
  where
    (cellBc, molBc) = runIdentity (e^.replicates) ^. files
    f x = let [a, b] = B.split '\t' x
          in (a, readInt b)

-- | Nucleotide count at each position. Assuming possible nts are:
-- A, C, G, T, N. Output is a 5 X N matrix.
posBaseCount :: PrimMonad m
             => Int  -- ^ Length of input dna
             -> ConduitT (B.ByteString, Int) o m (MU.Matrix Int)
posBaseCount n = do
    mat <- MUM.replicate (5, n) 0
    mapM_C $ \(dna, copy) -> forM_ [0 .. n-1] $ \j -> do
        let x = dna `B.index` j
            i = case x of
                'A' -> 0
                'C' -> 1
                'G' -> 2
                'T' -> 3
                'N' -> 4
                _   -> error "Unrecognized base"
        MUM.unsafeRead mat (i,j) >>= MUM.unsafeWrite mat (i,j) . (+copy)
    MU.unsafeFreeze mat
{-# INLINE posBaseCount #-}

deBarcode :: DropSeqConfig config
          => RNASeq S (SomeTags 'Fastq, SomeTags 'Fastq)
          -> WorkflowConfig config ( RNASeq S
                    ( File '[Gzip] 'Fastq
                    , File '[] 'Other
                    , File '[] 'Other ) )
deBarcode dropseq = do
    outdir <- asks _dropSeq_output_dir >>= getPath
    lenCellBc <- asks _dropSeq_cell_barcode_length
    lenMolBc <- asks _dropSeq_molecular_barcode_length
    let output = printf "%s/%s_debarcode.fastq.gz" outdir
            (T.unpack $ dropseq^.eid)
        flCellBc = printf "%s/%s_cell_barcode_stats.tsv" outdir
            (T.unpack $ dropseq^.eid)
        flMolBc = printf "%s/%s_molecule_barcode_stats.tsv" outdir
            (T.unpack $ dropseq^.eid)
        fl = ( location .~ output $ emptyFile
             , location .~ flCellBc $ emptyFile
             , location .~ flMolBc $ emptyFile )
    if flRead1 `hasTag` Gzip && flRead2 `hasTag` Gzip
        then let f1 = fromSomeTags flRead1 :: File '[Gzip] 'Fastq
                 f2 = fromSomeTags flRead2 :: File '[Gzip] 'Fastq
             in liftIO $ deBarcodeHelp f1 f2 lenCellBc lenMolBc output flCellBc flMolBc
        else let f1 = fromSomeTags flRead1 :: File '[] 'Fastq
                 f2 = fromSomeTags flRead2 :: File '[] 'Fastq
             in liftIO $ deBarcodeHelp f1 f2 lenCellBc lenMolBc output flCellBc flMolBc
    return $ dropseq & replicates.mapped.files .~ fl
  where
    (flRead1, flRead2) = runIdentity (dropseq^.replicates) ^. files

deBarcodeHelp :: (SingI tags1, SingI tags2)
              => File tags1 'Fastq
              -> File tags2 'Fastq
              -> Int      -- ^ Cell barcode
              -> Int      -- ^ Molecular barcode
              -> FilePath  -- ^ output Fastq
              -> FilePath  -- ^ cell map
              -> FilePath  -- ^ mol map
              -> IO ()
deBarcodeHelp flRead1 flRead2 lenCellBc lenMolBc output flCellBc flMolBc = do
    let taggedFastq = zipSources (sourceFastq flRead1) (sourceFastq flRead2) .|
            tagging lenCellBc lenMolBc
    runResourceT $ runConduit $ zipSources (iterateC (+1) (0::Int)) taggedFastq .|
        passthroughSink (foldlC mkBarCodeMap (M.empty, M.empty)) writeBarCodeStats .|
        concatMapC changeName .| unlinesAsciiC .| gzip .| sinkFile output
  where
    sourceFastq fastq = if fastq `hasTag` Gzip
        then sourceFile (fastq^.location) .| multiple ungzip .| linesUnboundedAsciiC
        else sourceFile (fastq^.location) .| linesUnboundedAsciiC
    changeName (idx, x) = case x of
        Nothing -> []
        Just ((cellBc, molBc), (_,b,c,d)) ->
            let name = "@" <> B.intercalate ":"
                    (map (B.pack . show) [cellBc, molBc, idx])
            in [name, b, c, d]
    writeBarCodeStats (cellBcMap, molBcMap) = liftIO $ do
        B.writeFile flCellBc $ B.unlines $
            map (\(a,b) -> intToDna lenCellBc a <> "\t" <> B.pack (show b)) $
            sortBy (flip (comparing snd)) $ M.toList cellBcMap
        B.writeFile flMolBc $ B.unlines $
            map (\(a,b) -> intToDna lenMolBc a <> "\t" <> B.pack (show b)) $
            sortBy (flip (comparing snd)) $ M.toList molBcMap
    mkBarCodeMap m (_, Nothing) = m
    mkBarCodeMap (!cellBcMap, !molBcMap) (_, Just ((cellBc, molBc), _)) =
        ( M.insertWith (+) cellBc (1::Int) cellBcMap
        , M.insertWith (+) molBc (1::Int) molBcMap )
{-# INLINE deBarcodeHelp #-}

dnaToInt :: B.ByteString -> Int
dnaToInt = fst . B.foldl' f (0, 1)
  where
    f (acc, i) 'A' = (acc, i * 5)
    f (acc, i) 'C' = (acc + 1 * i, i * 5)
    f (acc, i) 'G' = (acc + 2 * i, i * 5)
    f (acc, i) 'T' = (acc + 3 * i, i * 5)
    f (acc, i) 'N' = (acc + 4 * i, i * 5)
    f _ _          = error "Unexpected character!"
{-# INLINE dnaToInt #-}

intToDna :: Int   -- ^ length of the resulting bytestring
         -> Int
         -> B.ByteString
intToDna n = B.pack . reverse . go 1 []
  where
    go !i !acc !x
        | m == 0 && i >= n = c : acc
        | otherwise = go (i+1) (c:acc) m
      where
        c = case x `mod` 5 of
            0 -> 'A'
            1 -> 'C'
            2 -> 'G'
            3 -> 'T'
            4 -> 'N'
            _ -> undefined
        m = x `div` 5
{-# INLINE intToDna #-}

type TaggedFastq = Maybe
    ( (Int, Int)
    , (B.ByteString, B.ByteString, B.ByteString, B.ByteString) )

tagging :: Monad m
        => Int
        -> Int
        -> ConduitT (B.ByteString, B.ByteString) TaggedFastq m ()
tagging lenCellBc lenMolBc = do
    l1 <- await
    case l1 of
        Just (_, name) -> do
            Just (sequ, l2) <- await
            Just (_, l3) <- await
            Just (qual, l4) <- await
            case getBarCodes lenCellBc lenMolBc sequ qual of
                Just bc -> yield $ Just (bc, (name, l2, l3, l4))
                Nothing -> yield Nothing
            tagging lenCellBc lenMolBc
        Nothing -> return ()
{-# INLINE tagging #-}

getBarCodes :: Int   -- ^ Length of cell barcode
            -> Int   -- ^ Length of molecular barcode
            -> B.ByteString
            -> B.ByteString
            -> Maybe (Int, Int)
getBarCodes n1 n2 fastqSeq fastqSeqQual =
    let (cb, mb) = (B.take n1 fastqSeq, B.take n2 $ B.drop n1 fastqSeq)
        (cb_q, mb_q) = (B.take n1 fastqSeqQual, B.take n2 $ B.drop n1 fastqSeqQual)
    in if nBelowQ 10 cb_q > 1 || nBelowQ 10 mb_q > 1
        then Nothing
        else Just (dnaToInt cb, dnaToInt mb)
{-# INLINE getBarCodes #-}

nBelowQ :: Int -> BS.ByteString -> Int
nBelowQ q = BS.foldl' f 0
  where
    f acc x = if fromIntegral x - 33 < q then acc + 1 else acc
{-# INLINE nBelowQ #-}
