#!/usr/bin/env runghc

import Data.ByteString.Char8 as B
import Data.Text as T
import System.Environment
import Prelude as P

findChunks :: [ByteString] -> ()
findChunks _ = 
  undefined

main :: IO ()
main = do
  args <- getArgs
  case args of 
   [] -> error "Command should take two args: <.water file> <.aln file> "
   [wat, aln] -> do 
     P.putStrLn$  "Processing files: "++ show args
     f1 <- B.readFile wat
     f2 <- B.readFile aln
     let goodLines = P.filter (not . P.null . B.words) $ 
                     P.filter (not . isComment) $ 
                     B.lines f1
--     findChunks (B.lines f1)
     return ()

-- Is a given line a comment line:
isComment :: ByteString -> Bool
isComment _ = False
