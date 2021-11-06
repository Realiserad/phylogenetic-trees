-- IO i Haskell
-- Bastian Fredriksson 12/02-13

{-
Läs in genetiska sekvenser i det s.k Fasta-formatet,
antingen från fil eller från stdin. En indatafil är
på formatet:

>id1 Beskrivning
ACTCGTACGTACGT
>id2
ACTCGTACGTACGT

Där teckent ">" betecknar en kommentar. Detta tecken
åtföljs av ett obligatorisk id, och en valfri
kommentar. Nästa rad (eller rader) innehåller själva 
sekvensen (protein eller dna). I FASTA-formatet finns 
ingen begränsning för sekvensernas längd, men alla 
sekvenser skall vara lika långa.

Notera att endast ASCII-tecken är tillåtna i indata,
annars riskerar man att få felet hGetContents: invalid 
argument (invalid byte sequence). Det är också viktigt
att alla radbrytningar är på UNIX-formatet '\n' och
inte på formatet '\r\n' som används i Windows.
-}

import System.IO
import Data.List
-- Import labwork from F3.
import F3
import NJ
-- Required to get arguments passed as input on
-- execution.
import System.Environment

-- Do some IO stuff. Usage:
-- F4.exe [inFile] [outFile]
-- If no argument is given, the program will read from
-- stdin and print to stdout.
-- If one argument is given, the program will read from
-- file and write to stdout.
-- If two arguments are given, the program will read
-- from inFile and write to outFile.
main = do
    args <- getArgs
    let argsLen = length args
    if argsLen == 0 then do
        -- Read from stdin and write to stdout
        contents <- getContents
        putStrLn $ getTree contents
    else if argsLen == 1 then do
        let inFile = args !! 0
        -- Read from file and write to stdout
        contents <- readFile inFile
        putStrLn $ getTree contents
    else if argsLen == 2 then do
        let inFile = args !! 0
        let outFile = args !! 1
        -- Read and write to file
        contents <- readFile inFile  
        writeFile outFile $ getTree contents
    else
        error "Too many arguments"

-- Parse a string with genetic sequences on FASTA format using
-- recursive descent.
parse :: String -> ([String], [String])
parse (startSymbol:content) = 
    if startSymbol /= '>'
      then error "File must begin with comment"
      else parseId content [] [] ""

parseId [] _ _ _ = error "File cannot end with identifier"
parseId (h:t) ids seqs id
    | h == ' ' = parseComment t ((reverse id):ids) seqs
    | h == '\n' = parseSequence t ((reverse id):ids) seqs ""
    | otherwise = parseId t ids seqs (h:id)

parseComment [] _ _ = error "File cannot end with comment"
parseComment (h:t) ids seqs
    | h == '\n' = parseSequence t ids seqs ""
    | otherwise = parseComment t ids seqs

parseSequence [] ids seqs seq = (ids, ((reverse seq):seqs))
parseSequence (h:t) ids seqs seq
    | h == '\n' && null t = parseSequence t ids seqs seq
    | h == '\n' && isPrefixOf ">" t = parseId (drop 1 t) ids 
        ((reverse seq):seqs) ""
    | h == '\n' = parseSequence t ids seqs seq
    | otherwise = parseSequence t ids seqs (h:seq)

-- Takes an input string, and transforms it into a tree.
getTree :: String -> String
getTree str = res
  where
    (ids, seqs) = parse str
    res = neighbor ids $ makeDistanceMatrix $ 
        map fromSeqString seqs
