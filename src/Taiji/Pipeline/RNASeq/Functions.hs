{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE GADTs                 #-}
{-# LANGUAGE OverloadedLists       #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE PartialTypeSignatures #-}
module Taiji.Pipeline.RNASeq.Functions
    ( rnaGetFastq
    , rnaMkIndex
    , rnaDownloadData
    , rnaAlign
    , quantification
    , geneId2Name
    , mkTable
    ) where

import           Bio.RealWorld.GENCODE
import           Data.Bifunctor                       (second, bimap)
import qualified Data.ByteString.Char8                as B
import           Data.CaseInsensitive                 (CI, mk, original)
import           Data.Either                          (lefts)
import qualified Data.HashMap.Strict                  as M
import qualified Data.HashSet                         as S
import           Data.List.Ordered                    (nubSort)
import           Data.List.Singletons (Elem)
import           Data.Singletons                      (SingI)
import qualified Data.Text                            as T

import           Taiji.Pipeline.RNASeq.Types
import Taiji.Prelude
import Taiji.Utils (readGenesValidated)
import qualified Taiji.Utils.DataFrame as DF

type RNASeqWithSomeFile = RNASeq N [Either SomeFile (SomeFile, SomeFile)]

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

rnaMkIndex :: RNASeqConfig config => [a] -> ReaderT config IO [a]
rnaMkIndex input
    | null input = return input
    | otherwise = do
        genome <- getGenomeFasta 
        starIndex <- asks (fromJust . _rnaseq_star_index)
        anno <- getAnnotation 
        rsemIndex <- asks (fromJust . _rnaseq_rsem_index)
        liftIO $ do
            _ <- starMkIndex "STAR" starIndex [genome] anno 100
            _ <- rsemMkIndex rsemIndex anno [genome]
            return input

rnaAlign :: RNASeqConfig config
         => RNASeq S ( Either (SomeTags 'Fastq)
                              (SomeTags 'Fastq, SomeTags 'Fastq) )
         -> ReaderT config IO (
                RNASeq S ( Either (File '[] 'Bam, Maybe (File '[] 'Bam))
                    (File '[PairedEnd] 'Bam, Maybe (File '[PairedEnd] 'Bam))))
rnaAlign input = do
    dir <- asks _rnaseq_output_dir >>= getPath
    idx <- asks (fromJust . _rnaseq_star_index)
    let outputGenome = printf "%s/%s_rep%d_genome.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
        outputTranscriptome = printf "%s/%s_rep%d_transcriptome.bam" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
        opt = defaultSTAROpts & starCores .~ 4
                              & starTranscriptome .~ Just outputTranscriptome
        f (Left fl) = if fl `hasTag` Gzip
            then let fl' = Left $ fromSomeTags fl :: _ (_ '[Gzip] _) _
                 in starAlign outputGenome idx fl' opt
            else let fl' = Left $ fromSomeTags fl :: _ (_ '[] _) _
                 in starAlign outputGenome idx fl' opt
        f (Right (f1, f2)) = if f1 `hasTag` Gzip && f2 `hasTag` Gzip
            then let fl' = Right (fromSomeTags f1, fromSomeTags f2) ::
                        _ _ (_ '[Gzip] _, _)
                 in starAlign outputGenome idx fl' opt
            else let fl' = Right (fromSomeTags f1, fromSomeTags f2) ::
                        _ _ (_ '[] _, _)
                 in starAlign outputGenome idx fl' opt
    input & replicates.traverse.files %%~ liftIO . f

rnaDownloadData :: RNASeqConfig config
                => RNASeqWithSomeFile
                -> ReaderT config IO RNASeqWithSomeFile
rnaDownloadData dat = do
    tmp <- fromMaybe "./" <$> asks _rnaseq_tmp_dir
    dat & replicates.traverse.files.traverse %%~ (\fl -> do
        dir <- asks _rnaseq_output_dir >>= getPath . (<> (asDir "/Download"))
        liftIO $ downloadFiles dir tmp fl )

quantification :: (RNASeqConfig config, SingI tags1, SingI tags2)
               => RNASeq S ( Either (File tags1 'Bam)
                                    (File tags2 'Bam) )
               -> ReaderT config IO (RNASeq S
                    (File '[GeneQuant] 'Tsv, File '[TranscriptQuant] 'Tsv))
quantification input = do
    dir <- asks _rnaseq_output_dir >>= getPath
    idx <- asks (fromJust . _rnaseq_rsem_index)
    let output = printf "%s/%s_rep%d.bam" dir (T.unpack $ input^.eid)
            (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . either
        (\x -> rsemQuant output idx x opt)
        (\x -> rsemQuant output idx x opt)
  where
    opt = defaultRSEMOpts & rsemCores .~ 4

-- | Retrieve gene names
geneId2Name :: (RNASeqConfig config, Elem 'GeneQuant tags ~ 'True)
            => ([RNASeqWithSomeFile], [RNASeq S (File tags 'Tsv)])
            -> ReaderT config IO [RNASeq S (File '[GeneQuant] 'Tsv)]
geneId2Name ([], []) = return []
geneId2Name (ori_input, quantifications) = do
    outdir <- asks _rnaseq_output_dir >>= getPath
    anno_fl <- getAnnotation
    liftIO $ do
        id2Name <- fmap (M.fromList . map (\x -> (geneId x, geneName x))) $
            readGenesValidated anno_fl
        let getExp filtFn = concatMap split $ concatMap split $
                ori_input & mapped.replicates.mapped.files %~
                    map fromSomeFile . filter filtFn . lefts
            process fun input = do
                let output = printf "%s/%s_rep%d_gene_quant.tsv" outdir
                        (T.unpack $ input^.eid) (input^.replicates._1)
                input & replicates.traverse.files %%~ fun output
            convertTaiji out fl = do
                geneId2Name_ [6] True id2Name (fl^.location) out
                return $ location .~ out $ emptyFile
        quantifications' <- mapM (process convertTaiji) quantifications

        let encodeExpr = getExp $ \x -> x `hasTag` GeneQuant && x `hasTag` ENCODE
        encodeExpr' <- mapM (process convertTaiji) encodeExpr

        let customizedExpr = getExp (\x -> x `hasTag` GeneQuant && not (x `hasTag` ENCODE))
            convertCustom out fl = do
                geneId2Name_ [2] False id2Name (fl^.location) out
                return $ location .~ out $ emptyFile
        customizedExpr' <- mapM (process convertCustom) customizedExpr

        return $ quantifications' ++ encodeExpr' ++ customizedExpr'

geneId2Name_ :: [Int]    -- ^ Fields to be kept in the results. 1-based indexing.
             -> Bool     -- ^ Whether the input contains header.
             -> M.HashMap B.ByteString (CI B.ByteString)
             -> FilePath -- ^ Input
             -> FilePath -- ^ Output
             -> IO ()
geneId2Name_ idx has_header id2Name input output = do
    c <- B.lines <$> B.readFile input
    let content = if has_header then tail c else c
        header = if has_header
            then let xs = B.split '\t' $ head c
                 in B.intercalate "\t" $ "name" : map (\i -> xs!!(i-1)) idx
            else B.intercalate "\t" $ "name" : map (\i -> B.pack $ "value_" ++ show i) idx
    B.writeFile output $ B.unlines $
        header : map (B.intercalate "\t" . convert . B.split '\t') content
  where
    convert xs = getGeneName (head xs) : map (\i -> xs!!(i-1)) idx
    getGeneName x = case M.lookup x id2Name of
        Just x' -> original x'
        Nothing -> if mk x `S.member` geneSet
            then x
            else x <> "(UNKNOWN_Taiji)"
    geneSet = S.fromList $ M.elems id2Name
{-# INLINE geneId2Name_ #-}

-- | Combine RNA expression data into a table and output
mkTable :: RNASeqConfig config
        => [RNASeq S (File '[GeneQuant] 'Tsv)]
        -> ReaderT config IO (Maybe (File '[] 'Tsv))
mkTable [] = return Nothing
mkTable input = do
    outdir <- asks _rnaseq_output_dir >>= getPath
    let output = outdir ++ "/" ++ "expression_profile.tsv"
    liftIO $ do
        geneNames <- getGeneName input
        (sampleNames, vals) <- unzip . combine <$> mapM (f geneNames) input
        DF.writeTable output (T.pack . show) $ DF.mkDataFrame
            (map (T.pack . B.unpack . original) geneNames) sampleNames $
            transpose vals
        return $ Just $ location .~ output $ emptyFile
  where
    f genes fl = do
        (sampleName, vals) <- readExpr fl
        return $! (sampleName, map (\g -> M.lookupDefault 0.01 g vals) genes)
    combine = map (\x -> (fst $ head x, map average $ transpose $ map snd x)) .
        groupBy ((==) `on` fst) . sortBy (comparing fst)
    average [a,b]   = (a + b) / 2
    average [a,b,c] = (a + b + c) / 3
    average xs      = foldl1' (+) xs / fromIntegral (length xs)

readExpr :: RNASeq S (File '[GeneQuant] 'Tsv)
         -> IO (T.Text, M.HashMap (CI B.ByteString) Double)
readExpr input = do
    c <- B.readFile $ input^.replicates._2.files.location
    return ( fromJust $ input^.groupName
           , M.fromListWith max $ map f $ tail $ B.lines c )
  where
    f xs = let fs = B.split '\t' xs in (mk $ head fs, readDouble $ fs!!1)
{-# INLINE readExpr #-}

getGeneName :: [RNASeq S (File '[GeneQuant] 'Tsv)]
            -> IO [CI B.ByteString]
getGeneName = fmap S.toList . foldM f S.empty
  where
    f names x = foldl' (flip S.insert) names . map (mk . head . B.split '\t') .
        tail . B.lines <$> B.readFile (x^.replicates._2.files.location)
{-# INLINE getGeneName #-}