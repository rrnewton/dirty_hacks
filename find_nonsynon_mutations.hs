#!/usr/bin/env runghc

-- | USAGE: You should be able to just run this script.
-- 
--   HOWEVER, it requires a full Haskell platform install (namely the 'text'
--   package).

{-# LANGUAGE OverloadedStrings #-}
import Control.Monad
import qualified Data.ByteString.Char8 as B
import Data.Text as T
import Data.Text.Encoding as E
import Data.Text.Read 
-- import qualified Data.List as L
import Data.Char
import Data.Either
import System.Environment
import Prelude as P
import Debug.Trace
import Data.Map (fromListWith, toList)
import Data.Maybe (catMaybes)

-- | Look for three-line sequences of sequence text, alignment, sequence text.  For
-- example:
--
-- 57_3               1 ATGAATAACATTACCATTAATTACGCAAAGGACAGCAAATTAACCGACTT     50
--                      ||||||||||||||||||||||||||||||||||||||||||||||||||
-- contig00193_2      1 ATGAATAACATTACCATTAATTACGCAAAGGACAGCAAATTAACCGACTT     50
findChunks :: [Text] -> [(Text,Text,Text)]
findChunks [] = []
findChunks [_] = []
findChunks [_,_] = []
findChunks (ln1:ln2:ln3:rst) 
  | [_nam1, num1a, seq1, num1b] <- T.words ln1
  , [syms] <- T.words ln2
  , [_nam2, num2a, seq2, num2b] <- T.words ln3
  , posNum num1a, posNum num1b
  , posNum num2a, posNum num2b
  = 
    if not (T.all validAlignChar syms) 
    then error $ "findChunks: Got what seems to be a matching entry, but middle line contains invalid alignment char: "++show syms
    else 
--     trace ("Passed: "++ show (ln1)) $ 
     (seq1,syms,seq2) : findChunks rst

findChunks (hd:rst) = 
--  trace ("FAILED: "++show (P.take 3 $ hd:rst)) $ 
  findChunks rst

main :: IO ()
main = do
  args <- getArgs
  case args of 
   [] -> error "Command should take two args: <.water file> <.aln file> "
   [wat, aln] -> do 
     P.putStrLn$  "Processing files: "++ show args
     f1 <- B.readFile wat
     f2 <- B.readFile aln
     let t1 = E.decodeUtf8 f1
         t2 = E.decodeUtf8 f2
     let goodLines = -- P.filter (not . P.null . T.words) $ 
                     P.filter (not . isComment) $ 
                     T.lines t1
         chunks = findChunks goodLines
         (seq1_ls,aligns_ls,seq2_ls) = unzip3 chunks
         seq1   = T.concat seq1_ls
         aligns = T.concat aligns_ls
         seq2   = T.concat seq2_ls
         len1   = T.length seq1
         len2   = T.length aligns
         len3   = T.length seq2

     -- The second file has the amino acids:
     let coding = P.filter (T.any (== '*')) $ 
                  T.lines t2
         whitespaceDistro = P.map (T.length . T.takeWhile isSpace) coding
         whitespaceMax    = P.minimum whitespaceDistro
         codingChanges0   = T.concat $ P.map (T.drop whitespaceMax) coding
         codingChanges    = T.dropWhile (== ' ') codingChanges0

--     mapM_ print chunks
     P.putStrLn$ "Parsed "++show(P.length chunks)++" valid chunks of sequence alignment data."
     P.putStrLn$ "Length of seq1, aligns, seq2: "++show(len1,len2,len3)
     unless (len1 == len2 && len2 == len3) $ 
        error "Cannot proceed! Length mismatch."
     P.putStrLn$ "Number of total mutations: "++ show(T.count "." aligns)
     P.putStrLn$ "Number of total insertions: "++ show(T.count "-" aligns)
     P.putStrLn$ "Number of amino acids: "++show (T.length codingChanges)
     P.putStrLn$ " (Dropped "++show (T.length codingChanges - T.length codingChanges0)++" spaces from beginning of amino acid allignment.)"
     P.putStrLn$ " (Leading whitespace length distribution for coding patterns: "++show whitespaceDistro++")"
     unless (T.length codingChanges == len1 `quot` 3) $
       error ("Expected number of coding changes to be "++show(len1 `quot` 3))
--     P.putStrLn$ "Coding or non-coding changes?: "++show codingChanges

     codingMuts <- forM [0.. len1-1] $ \ ix -> do 
       let st = seq1 `index` ix
           en = seq2 `index` ix
       case aligns `index` ix of 
         '.' -> do P.putStr$ "  Mutation at base-pair position "++show ix++", "++[st]++" -> "++[en]
                   case codingChanges `index` (ix `quot` 3) of
                     '*' -> do P.putStrLn ", non-coding."
                               return (Just (Left (st,en)))
                     _   -> do P.putStrLn ", CODING."
                               return (Just (Right (st,en)))
         _   -> return Nothing
     P.putStrLn "\nFinal tables: FROM, TO, COUNT"
     let printTable changes = do 
          P.putStrLn "----------------------------"
          forM_ (frequency changes) $ \ ((st,en),cnt) -> do
            P.putStrLn$ "  "++[st]++" "++[en]++" "++show cnt
     P.putStrLn "\nNON-synonomous changes:"
     printTable $ rights $ catMaybes codingMuts
     P.putStrLn "\nSynonomous changes:"
     printTable $ lefts $ catMaybes codingMuts
     return ()

-- Is a given line a comment line:
isComment :: Text -> Bool
-- isComment :: ByteString -> Bool
isComment t = f (strip t)
 where
  f x | T.null x        = False
      | T.head x == '#' = True
      | otherwise       = False

frequency :: (Ord a) => [a] -> [(a, Int)]
frequency xs = toList (fromListWith (+) [(x, 1) | x <- xs])

posNum :: Text -> Bool
posNum t = 
  case decimal t of
    Left _  -> False
    Right _ -> True

validAlignChar :: Char -> Bool
validAlignChar '.' = True
validAlignChar '-' = True
validAlignChar '|' = True
validAlignChar _   = False
