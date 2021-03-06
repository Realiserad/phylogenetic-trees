module F3 where -- used by F4

import Data.List
import NJ

-- Evolutionära träd i Haskell
-- Bastian Fredriksson 16/11-13

-- A type Evol that contains genetic data that can be compared.
class Evol a 
  where
    distance :: a -> a -> Float
    -- Create a distance matrix, with the evolutionary distance
    -- between all pair of species, defined in a list of Profiles
    -- or in a list of MolSeqs
    makeDistanceMatrix :: [a] -> [[Float]]
    makeDistanceMatrix pl = res
      where
        -- Solve this problem using pattern matching. We create the
        -- matrix by creating rows. Each row is a sublist with elements
        -- were element in row R and column C is the distance between
        -- the R:th and C:th profile in the list. The result will be
        -- a symmetric matrix with zeros on the diagonal. 
        res = createMatrix [] pl pl
        createMatrix res [] _ = res
        createMatrix res (h:t) pl = (createRow h pl):
            (createMatrix res t pl)
        createRow c [] = []
        createRow c (h:t) = (distance c h):(createRow c t)

instance Evol MolSeq
  where
    distance seq1 seq2 = seqDistance seq1 seq2

instance Evol Profile 
  where
    distance p1 p2 = profileDistance p1 p2

-- A datatype Profile containing a summary of a
-- set with dna or protein sequences.
data Profile = Profile { 
    -- A profile matrix containing a fingerprint
    -- of the set with sequences it was created from.
    profile :: [[(Char, Float)]],
    -- The number of sequences used to create the
    -- profile matrix.
    numOfSequences :: Int,
    -- True if the profile matrix was made from
    -- strings of dna.
    isDna :: Bool,
    -- A unique id used to identify the profile.
    -- Will be set to the the first molecule
    -- sequence as default.
    profileName :: String 
} deriving (Show)

-- Wrapper for a molecular sequence.
data MolSeq = MolSeq {
    molSequence :: String
} deriving (Show)

-- Calculate the evolutionary distance between two species
-- using a pair of profile matrices.
profileDistance :: Profile -> Profile -> Float
profileDistance x y = res
  where
    -- Extract the profile matrices
    p1 = profile x
    p2 = profile y
    -- Flatten the matrices down to lists of tuples.
    a = concat p1; b = concat p2
    -- Create a list of indexes.
    ind = [0.. (length a) - 1]
    -- Count the distance between the lists by summing
    -- up the absolute value of the elementwise difference.
    res = foldl (\sum i -> sum + (abs (snd (a !! i) - 
      snd (b !! i)))) 0.0 ind

-- Calculate the evolutionary distance between two species
-- using a pair of molecule sequences.
seqDistance :: MolSeq -> MolSeq -> Float
seqDistance seq1 seq2
    | length s1 /= length s2 = 
        error "Sequences differ in length"
    | isDnaSeq s1 /= isDnaSeq s2 = 
        error "Both sequences must be of the same type"
    -- Sequence is DNA, use Jukes-Cantor's formula
    | isDnaSeq s1 = jukesCantor
    -- Sequence is protein, use Poisson's ratio
    | otherwise = poisson
    where
      s1 = molSequence seq1
      s2 = molSequence seq2
      a = (mut s1 s2 0.0) / (fromIntegral $ length s1)
      mut [] [] n = n
      mut (h1:t1) (h2:t2) n 
          | h1 == h2 = mut t1 t2 n 
          | otherwise = mut t1 t2 $ n + 1.0
      poisson
          | a > 0.94 = 3.7 
          | otherwise = -19/20 * log (1 - 20/19 * a)
      jukesCantor
          | a > 0.74 = 3.3
          | otherwise = -3/4 * log (1 - 4/3 * a)

-- Create a profile from a list of sequences.
fromSeqStrings :: [String] -> Profile
fromSeqStrings strs = Profile {
    profile = makeProfileMatrix strs,
    numOfSequences = length strs,
    isDna = isDnaSeq $ head strs,
    profileName = head strs
}

-- Create a MolSeq from a sequence string.
fromSeqString :: String -> MolSeq
fromSeqString str = MolSeq {
    molSequence = str
}

-- Determine whether the sequence is a dna sequence by
-- examining if it contains other molecules than A, T, C
-- or G. If no other molecules are found, it is considered
-- to be a dna sequence, otherwise it is probably a
-- protein sequence.
isDnaSeq :: String -> Bool
isDnaSeq str = res
  where
    res = scan str
    scan [] = True
    scan (h:t) | elem h nucleotides = scan t | otherwise = False

-- Create a profile matrix from a list of sequences.
nucleotides = "ACGT"
aminoacids = sort "ARNDCEQGHILKMFPSTWYVX"
makeProfileMatrix :: [String] -> [[(Char, Float)]]
makeProfileMatrix strs = res
  where
    -- Count the number of sequences.
    n = length strs
    -- Define our molecules to be either aminoacids or
    -- nucletides depending on the type of genes we are
    -- dealing with.
    dna = isDnaSeq (head strs)
    molecules = if dna then nucleotides else aminoacids
    -- Create tuples (A,0), (C,0), (G,0), (T,0) or
    -- (A,0), (R,0), (N,0), (D,0), (C,0) ect.
    tuples = zip molecules $ replicate (length molecules) 0.0
    -- Create a list of sublists with tuples. Each tuple
    -- contains a letter, and the corresponding number of
    -- occurrences for that letter in the sequences for a
    -- specific index position, divided with the number of
    -- sequences. For example, if we have dna sequences 
    -- ["ATCG","CTCC"] our first sublist will contain the tuples 
    -- (A,1/2) and (C,1/2) because A and C occurs one time as the 
    -- first letter in the sequences. So, the list tmp1 would be
    -- [[(A,1/2),(C,1/2)],[(T,1)],[(C,1)],[(C,1/2),(G,1/2)]]
    tmp1 = map (map (\x -> ((head x), fromIntegral (length x)/
      fromIntegral n)) . group . sort)
           $ transpose strs
    -- Introduce an equality test for tuples. Tuples 
    -- (A, B) and (C, D) are equal if A=C holds. This 
    -- function is used by unionBy to exclude copies.
    equalFst a b = (fst a) == (fst b)
    -- Fill the matrix with zeros and sort the result, so
    -- the first column contains tuples with the letter A.
    -- From our previous example, we would obtain the result
    --        [(A,1/2) (C,1/2) (G,0)   (T,0)]
    -- res =  [(A,0)   (C,0)   (G,0)   (T,1)]
    --        [(A,0)   (C,1)   (G,0)   (T,0)]
    --        [(A,0)   (C,1/2) (G,1/2) (T,0)]
    res = map sort $ map (\l -> unionBy equalFst l tuples) tmp1
    
{-   _            _    __      _ _ 
    | |          | |  / _|    | | |
    | |_ ___  ___| |_| |_ __ _| | |
    | __/ _ \/ __| __|  _/ _` | | |
    | ||  __/\__ \ |_| || (_| | | |
     \__\___||___/\__|_| \__,_|_|_|
    Testfall inför redovisningen
-}
--       |  fly  | mouse | human |
--       |       |       |       |
-- fly   |   0   |       |       |
-- mouse |  26   |   0   |       |
-- human |  28   |  17   |   0   |
t1 = res
  where
    (species, sl) = unzip dnastr
    dm = makeDistanceMatrix $ map (\x -> fromSeqStrings x) $ take 3 sl
    res = neighbor (take 3 species) dm

--         | fly |mouse|human|fish|lizard|penguin|eagle|snake|horse|
--         |     |     |     |    |      |       |     |     |     |
-- fly     |  0  |     |     |    |      |       |     |     |     |
-- mouse   | 25  |  0  |     |    |      |       |     |     |     |
-- human   | 28  | 17  |  0  |    |      |       |     |     |     |
-- fish    | 32  | 23  | 22  | 0  |      |       |     |     |     |
-- lizard  | 36  | 30  | 25  | 32 |  0   |       |     |     |     |
-- penguin | 25  | 16  | 14  | 23 |  26  |  0    |     |     |     |
-- eagle   | 33  | 26  | 22  | 29 |  27  | 20    |  0  |     |     |
-- snake   | 30  | 25  | 20  | 28 |  31  | 21    | 21  |  0  |     |
-- horse   | 24  | 18  | 22  | 29 |  30  | 19    | 28  |  27 |  0  |
t2 = res
  where
    (species, sl) = unzip dnastr
    dm = makeDistanceMatrix $ map (\x -> fromSeqStrings x) sl
    res = neighbor species dm

-- FOXP4 from lab F2
t3 = res
  where
    (species, sl) = unzip foxp4
    dm = makeDistanceMatrix $ map (\x -> fromSeqString x) sl
    res = neighbor species dm

{-   _            _      _       _        
    | |          | |    | |     | |       
    | |_ ___  ___| |_ __| | __ _| |_ __ _ 
    | __/ _ \/ __| __/ _` |/ _` | __/ _` |
    | ||  __/\__ \ || (_| | (_| | || (_| |
     \__\___||___/\__\__,_|\__,_|\__\__,_|
    DNA-strängar från några olika arter
    Källa: http://www.ncbi.nlm.nih.gov/
-}
-- Format [("NAME" [[..], [..]]), ..]
    -- Bananfluga
dnastr = [("FLY", [
    "ATGTCGTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGA",
    "GGCATTGGTCTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTG",
    "AAGAACCTGGTGATCCTCGACCGCATTGAGAACCCGGCTGCCATT",
    "GCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACCTTCTAC",
    "CCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTG",
    "AAGACCATCTTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAAC",
    "GGAGCTGGTATCCTGGACGATCACCAGATCGAGCGCACCATTGCC",
    "GTCAACTACACTGGCCTGGTCAACACCACGACGGCCATTCTGGAC",
    "TTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAAC",
    "ATTGGATCCGTCACTGGATTCAATGCCATCTACCAGGTGCCCGTC",
    "TACTCCGGCACCAAGGCCGCCGTGGTCAACTTCACCAGCTCCCTG",
    "GCGAAACTGGCCCCCATTACCGGCGTGACGGCTTACACTGTGAAC",
    "CCCGGCATCACCCGCACCACCCTGGTGCACACGTTCAACTCCTGG"]
    -- Skogsmus
    ), ("MOUSE", [
    "ATGGCTGTAGCCGTGGCCGCCGCCGGGGTCCTAATGGGCTCGGAG",
    "AACTTGAGTACCTGTCTTTGGTGTCGAAGGTTTGCACGGAGCTGG",
    "GGATCTGGCTGAATTTGTCATCAGTCTTGCTGAGAAAAATACTAC",
    "GTCAAAAATGGTGCAGAATTCACAGATTCTCTTATTAGCAACCTC",
    "ACTCTTCCCCGTCCTCTGCCAGCCGGACAACCCTTCAGCTCGGAC",
    "GTCGCTGTGGACGTCTTGAAAGAGCTGGAAGCTTTAATGCCCAGC",
    "ACCCTGAGCACAGGGACAGGACAAAGAAGAAGAAGAGGAGTCGGA",
    "CCGAGACCGGGACCGGGACCGGGACCGGGACAGGGACCGTGACAA",
    "GAACGAGACAGAGAGCGAGATAGAGAGAGAGACCATAAGCGGAGA",
    "AGATCGGAAAGACCGAGAAAAGTATGGGGAGAGGAACTTGGACAG",
    "CCTCCCCCGGAAGAGCCTGCCATTGGAGATATTTACAATGGCAAA",
    "GCTTTGTGCAGCTGGAAGGCCTAAGGAAGCGGTGGGAAGGCCTGG",
    "AGGTCGAGTGGCCAATGTGGCTGATGTTGTGAGTAAAGGCCAGAG",
    "ACTGGAACCAAGACCAGCCTGAGCATGAAGGATGTGGACCAAGAG",
    "GACGGAGAAATCTGGTTGGCGAGACCAATGAAGAGACGTCAATGC",
    "CTCCCTTGTCAGCGCTCCTGAAGTAGAGGATGACTCGCTGGAGCG",
    "CCAGAGAAGTGGGAGATCAAACAGATGATTGCTGCCAATGTGCTT",
    "ACGAAGAGACTGGCATTCTCCCCAAGGTGGATGATGAAGAAGATG",
    "AGAGGAGCCTCCGTTCCTGAGAGGGCATACCAAGCAGAGCATGGA",
    "AACCCAGATGGCTCACTCTCCCAAGCAGCCATGATGCAGAGTGCC",
    "AGCAGGCCCAGCGGGAAGCAGAGATGGATTCTATACCCATGGGGC",
    "GCCCGATGCGGAAGGCAGGCAGATCGCTGCCAACATGAGGGGGAT",
    "GAGTGGAAGAAGCATGCCTTTGGGGGCAACAAAGCTTCTTATGGG",
    "AGCAGAGAGAGAGCCTGCCCATCTACAAACTGAAGGAGCAGCTGG",
    "CCTAATTGTAATTGGAGAGACAGGGTCTGGGAAGACAACACAGAT",
    "TATACTTCCAGGGGCAAGATTGGGTGTACCCAGCCCAGAAGAGTG",
    "TCTCAGAAGAGTTTGGTTGTTGCCTTGGTCAAGAGGTGGGCTACA",
    "CCCAGAAACAGTCATCAAGTACATGACAGATGGCATGCTGCTAAG",
    "ACTCAGTATGCCATCATTATGTTGGATGAAGCACACGAGAGGACC",
    "TGTTGAAAAAGACAGTTCAAAAACGACAAGACATGAAGCTGATTG",
    "GAAGTTTTCTCAGTACTTCTATGAAGCACCCATCTTCACAATCCC",
    "CTGTACACAAAGGAGCCGGAAACAGACTATTTGGACGCTAGCCTG",
    "CAGAGCCGCCAGGTGACATCTTGGTTTTCTTGACTGGTCAGGAAG",
    "GTATGAAAGAATGAAGTCCCTGGGACCTGATGTCCCAGAGTTGAT",
    "CCGAGTGAGATGCAGACCCGAATCTTTGATCCAGCTCCCCCTGGT",
    "ACATTGCAGAAACATCTCTGACCATCGACGGCATCTACTATGTGG",
    "AGTGTACAATTCTAAAACAGGGATCGATCAGCTTGTGGTGACACC",
    "AGAGCTGGCCGGGCTGGGAGAACTGGCCCAGGAAAGTGTTACAGA",
    "ATGAAATGCTGACCACCAATGTACCGGAAATCCAGAGAACCAACC",
    "GGCCATGGGTATTAATGACCTGCTGTCCTTTGACTTCATGGATGC",
    "GCCATGGAGCAGCTGTATACACTGGGGGCCCTGGACGACGAGGGC",
    "TGGCAGAGTTCCCTCTGGAGCCAATGCTGTGCAAAATGCTCATCA",
    "AGAGATGCTAACTATCGTGTCCATGCTGTCTGTGCAGAACGTCTT",
    "CTTGCAGATCAGAAGAAAGCCAAATTCCACCAGACAGAAGGGGAC",
    "ACTCCTGGAAGAACAACAAATTCTCCAACCCGTGGTGCTATGAGA",
    "CCGTGCCCAGGACATTCGGAAGCAGATGTTAGGCATAATGGACAG",
    "GGCAAATCCACAGTCCGAGTGCAGAAGGCCATCTGCAGTGGGTTC",
    "CACAGGAAGGCTACCGGACTCTGATAGACCAGCAGGTGGTGTATA"]
    -- Människa
    ), ("HUMAN", [
    "GATCTTATCTATCATGTTCACCTCCCAAGAGGTGAACATATCCCC",
    "CATTAATATTTAATGCATGACCATGTGCAGACTTGGGAGGAAAAA",
    "TCCTTAATAAACAAGGATGTTTCTGCATCATTTCCCCACAACACC",
    "TTTAAGCAAATGCATTGTTTTTCCAGTTATATATCTGGTAGAGAT",
    "CGATCTCCTTTTATTTTGATGACCCAGCATGGCTGAACACTCAGT",
    "TTCAGCATTAGAGATGCCAGCCCTGTAGGATATAAAACAGGAACA",
    "TACTCAAGTCTTAGAAGCACCACTTGTCTTTTTTCAAGGGAGAGA",
    "AAGGGAGGGAGTCACTCACTTGAACGGTTCCCTTAGGCTGTGTGG",
    "CTGACAGTGGGAAATGCACTGGAGACGATGACTGGCAAAGCCCTC",
    "CTGACAGCAAAGGGTTTGTCACAATGACAACTATACACTCCCAAT",
    "GGTATATTATGAGTGACTGAAGTTTAGAATAAATTAATAAATATT",
    "AAGGTCTAGTAAGGCTAAGGATATAACAAGAAAATAATATGAATA",
    "GAGTAAGTTACAAATGGCTTCAGGAAGGGGAGAGAGGAAGAAGAG",
    "AGGGCTAATTTTATGAAAGCTTTGGGAAGTTTTAAGAAAAAGAAA",
    "TATGCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
    "CCTAAGAAGACTATGAGACACTAAGAGAAAAATTAAGGTAAAAAA",
    "AGGGAGGAGGGAGGAGGTTAAGACATTTTACTATGTGCTGTGAAT",
    "ATGCAATATATATACATATATACACACATATACATATGTATTTAA",
    "TTTAGAGATATGGTTTCACTATGTCACTCTGCCCAGGCTGCAGTA",
    "TAGCACATTATAGCCTTGAACTCCTGGGCTCAAGCAACCCTCCTG",
    "TACTAGCATATGCCACCATGTCCACCTTTATGCTTTTTAAAGTGA",
    "TCAACTTAATAATAAAAACATTTCAAATGTAAAGAAATTTACAAA",
    "TTGGGCAAAGGGAATGAACAGACACTTTTCAAAAGAATACATGCA",
    "AAGTTCAACATCACTGATCATTAGAGAAATGCAAATCAAAACCAT",
    "AGAATAGCTATCATTAAAAAGTCAAAAAATAACAGATGCTAGTGA",
    "TACACTGTTGTTGGGTGTGCAAATCAGTTCAATCATTGTGCAAGG",
    "GCAGAGCTACCATTCGACCCAGTAATCCCACTACTGGGTATATAC",
    "ATAAAGACACATGCATACAAATGTTCATTGCAGCACTGTTCACAA",
    "ATGCCCATCAATGACAGATTGGATAAAGAAAATGTGGTACATATA",
    "AAAAAATGATATCATGTCTTTTGCTGGAATATGGATGGACCTTCT",
    "AACAGAAAACCAAATATAGCATACTCTCAGTTATAAGTGGGAGCT",
    "AGAATAAAACAGACACTGGGGTCTACTTGAGGGTGGAGGGTGAGA",
    "CTATTGGGTACTAGGTTTAATACCTGGGTGATGAAATGATCTGTA",
    "ACCTATGTAACAAATGCCCCTAAACTTAAAATAAAAGTTAAAAAA",
    "ATCTACCTGGTAATATGAAAAACACATATCTTTCATTCATTCCTT",
    "TGGGAGTTAGTAAAAGTCCACATTGAGATATGAGACCCACCACTG",
    "AATCCCAGCACTTTGGGAGGCCGATGCTGGTGGATCACCTAAGGT",
    "ACATGGTGAAACCCCCATCTCTACTAAAAATACAAAAATTAGCTG",
    "CCAGCTACTAGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCAGG",
    "TCATGCCATTGCACTCCAGCCTGGGCAACAAGAGCAAGACTCTGT",
    "CCACCACCATCATTTTGCAAGTGTTACCACTATTGTGTGTTAATA",
    "TTCTTTGTATTCCTAATTGTAATAGCTTTGTATTTGAAAAATTAT",
    "TTTGTATGCGATGACAACAGAATATATTATCATGCTCCTTTTGTG",
    "AATTTGTGATTTTGCTTTAATTTGAAATATTAATTTCAAATATGT",
    "ACAGTAAATCTGTGGATTAAGTAATGTCTTAGTAGGTATTGGGAA",
    "ATATTGTCATTGTTTATTTCAAAGCCAGTTAAAATTCTGCAAAGC",
    "AATTTATAAAATACCGAGATTACGGTGTATAAACAACTTTAGATT",
    "TGTAATATATGCTTCATTCAAAGTAGCTAAGGGCTGTACCTGGCT",
    "AAAAGGAATACTGAGTAGCTGGGACCTCCTGAGTAGCTGGGACCA",
    "ATTACTGTTTAGAGAATAACATTTGATGGAATCATGCTTTTACTT",
    "CTGACATTAACATCCCAAATCCTTAGCATGGCCTACAAGGCCCTG",
    "GCTGCCTCATTTAATAACTCTTTGTCTCTTTCCCAGATCCAGCCA",
    "CAAGACAAGCTCTTCCCAGAACCTGACCTTTGTACCTGTTCTTTA",
    "AATTACTTATCATCTATCATAATTCAGGTTACATGGCACTAACTC",
    "CTTCTCCAACCAAATTAGGAACAATTATATGGCCACATAGTATCG",
    "TTGGGAGATTTTGTTGTTTAACACTTGTTTTCACTATAAGACTGT",
    "CTGTTTGTTCACTCCTGCCACAGTCAGAATAGTGCCTGGAATATG",
    "TAAATGAACAAGTGAATAAATGGATATTGTCTCATTTTTAGAACA",
    "ATCTGGTATGTCACGTAGGTAATTTACAAGGGCTACAATTTCAGC",
    "GGTCTTGATAGGTCTCTTGATGTCATTTCACTTCAGATTCTTCTT",
    "CTGTCTTGTCCAAATTGTTACTAAGAATCAAGAGAGATATCTGAC",
    "ACACGATTGAAATAATGCTAGCCAATATGGTTATTATTAGAAACC",
    "TAATACTTATTGCAGACTCAAATGTGCTTATTCTAAAACAAGTAA",
    "AATCCACGGAGTTCATTCTAATCCACATTCAACACTATCATGTAC",
    "CCTGTGATTTTTCAGGTTCACTTTTCTAAACTTGTGAATTAAATA",
    "AAAAAACTCTTGTAATTGTTGCCCATTTCAGGAGAAATCTTGCAT",
    "AACTGAGGGCTGTGGTTTAAACAAAATCTTGAGAATGTTTTTTGA",
    "ACAAAATGATATAGACAAAGGTAACTTTTAATAGAACCAGTCACT",
    "TGCTTAGCTAAGCAACAGAGAAGGTAAAATACTAATTCAATTCAT",
    "CCAAGTATGTGCTCACTGAATAAGCTGCTAAGGTTTGGTGGTTAC",
    "CATCACAGTCCAACATTCACAGAGTTTAAAAGCCTACCAAGAATC"]
    -- Kvastfening
    ), ("FISH", [
    "TCTGAAGGAAAACAAAGGAAAGGAGGAGTCAAACAGTCTGCAGTG",
    "CAAAAAACAATAGAAAAAGTGCTTCTTGTGACCACTCTGAATGTG",
    "CTTAAAAGGTCAACAGCACTCACATAGAAGAAATAAACCATGTCA",
    "CCTTTACTGAAACACCCTAATGAACAGAGAAAAATTCACTCTGAA",
    "GTGGAAAAGATTTGAAACTGTCTTTAGAGCCTGACCACTCTCAAC",
    "CAAATGTTCTGATTGTGGGAAGTGTCATAGATGGTCTGGTCAATT",
    "GAAGAGAAATCACATAAATGTACGGAATGTGGAAAAACTTTTAGA",
    "AACACATCCACACAGAAGATAAACCACACAAATGTACAGAATGCA",
    "CCTTAATAAACATCAACGTATCCATACAGGAGAGAAACCCTACAA",
    "AGATACTCTGAACACCTTAACAGACACCAACGAACTCACACAGGA",
    "GTGGGAAGAGTTTCAGAGATTCAGCAGCCCTGAACTCACACCAAC",
    "CAAATGCACTAAATGTGGAAAGAGCTTCAGAGATTCTTCAGCCAT",
    "GGAGAGAAACCATATAGATGTACAGAATGTGGGAAGAATTTCAGA",
    "AACGTATTCACACAGGAGAGAAACCATTCAGATGTACAGAATGTG",
    "CCTTGACAGACACCAACAAATTCACGCAGGCGAGAAGCTATACAA"]
    -- Geckoödla
    ), ("LIZARD", [
    "GAAAAACCACCCAGTACTAAAAATTGTTAATAGCTCATTTATCGA",
    "CCTCCCCTCCCCATCAAATATTTCAATCTGATGGAACTTTGGATC",
    "ACTACTAGGAATCTGCCTTATTATACAAATTATAACCGGCCTATT",
    "CTTAGCAATACACTACACCGCCGACACCTCCATCGCATTCTCTTC",
    "AATCGCCCATATTTGCCGGAACGTTAACTACGGCTGACTAATTCG",
    "TAACCTCCACGCTAATGGTGCCTCTTTATTTTTCATCTGCCTTTA",
    "CTTACATATCGGACGAGGCCTGTACTATAGCTCCTACTTATATAA",
    "AGAAACATGAAACATCGGAGTCCTCCTCCTCCTCCTAACCATAGC",
    "AACCGCCTTCATGGGTTATGTCCTACCTTGAGGACAAATATCATT",
    "GAAAAACCACCCAGTACTAAAAATTGTTAATAGCTCATTTATCGA",
    "CCTCCCCTCCCCATCAAATATTTCAATCTGATGGAACTTTGGATC",
    "ACTACTAGGAATCTGCCTTATTATACAAATTATAACCGGCCTATT",
    "CTTAGCAATACACTACACCGCCGACACCTCCATCGCATTCTCTTC",
    "AATCGCCCATATTTGCCGGAACGTTAACTACGGCTGACTAATTCG",
    "TAACCTCCACGCTAATGGTGCCTCTTTATTTTTCATCTGCCTTTA",
    "CTTACATATCGGACGAGGCCTGTACTATAGCTCCTACTTATATAA"]
    -- Pingvin
    ), ("PENGUIN", [
    "GTCCCTGTAGCTTATAACCAAAGCATGGCACTGAAGATGCCAAGA",
    "TGATGCTATATGCACCCAAGGACAAAAGACTTAGTCCTAACCTTA",
    "CGTTAATTTCTGCTAAACATATACATGCAAGTATCTGCGCCCCAG",
    "GTAAATGCCCCAGAGCTCTTATTGACTAGACAAAGGGAGCGGGCA",
    "CAGGCTCGCCCATTGCTGCAGCCCAAGACGCCTTGCTTAGCCACA",
    "CCCCACGGGTATTCAGCAGTAATTAGTATTAAGCAATAAGTGAAA",
    "CTTGACTTAGTTATAGCAGCATTCAGGGTTGGTAAATCTTGTGCC",
    "GCCACCGCGGTCACACAAGAGACCCAAATTAACTGTCACACGGCG",
    "AAAGAGTGGTACCATGCTATCCTAACAACTAAGATCAAAATGCAG",
    "TAAGCTGTCATAAGCCCAAGATGTACCTAAAAACACCCCCAAGAC",
    "ATCTTAGCGCCCCCGACTGATTAAACTCCACGAAAGCTAGGACAC",
    "AACTGGGATTAGATACCCCACTATGCCTAGCCCTAAATCTTGATA",
    "CTTCTATCACCAAAGTATCCGCCTGAGAACTACGAGCACAAACGC",
    "TAAAACTCTAAGGACTTGGCGGTGCCCCAGACCCACCTAGAGGAG",
    "CTGTTCTATAATCGATAACCCACGATGCACCCGACCGCCCCTTGC",
    "AAAACAGCCTATATACCGCCGTCGCCAGCCCACCTCCCCTGAAAG",
    "ACAACAGTGGACATAATAGCCCCAACCCCGCTAACAAGACAGGTC",
    "AGGTATAGCCTATGGGGCGGAAGAAATGGGCTACATTTTCTAAGA",
    "AGATAACTTCACGGAAGGGGACATGAAATTGCCCCCGGAAGGCGG",
    "TTTAGCAGTAAGGTGGGACAATAAGGCCCACTTTAAGCTGGCCCT",
    "AGGCACGTACATACCGCCCGTCACCCTCCTCACAAGCTACGCACA",
    "CACATAACTAATACCCCCCCCCCGCTGAAGATGAGGTAAGTCGTA",
    "CAAGGTAAGTGTACCGGAAGGTGCACTTAGCACATCAAGACGTAG",
    "TATAACATAAAGCACTCAGCTTACACCTGAAAGATATCTGCTACC",
    "ACCAGATCGTCTTGAAGCCAAACTCTAGCCCAACCACGACACAAC",
    "AAACAATTAAAAATTTACTCCACTATAAATTAAACTAAAACATTC",
    "TCTGCCCCAGTATGGGCGACAGAAAAGGCCTACTGGCGCGATAGA",
    "ATCCGTACCGTAAGGGAAAGATGAAATAGCAATGAAAACCCGAGC",
    "AAAAACAGCAAAGATAAACCCTTGTACCTCTTGCATCATGATTTA",
    "CAAGAACAACCAAGCAAAACGAATTTAAGCTTGCCACCCCGAAAT",
    "CGAGCGAGCTACTTACAAGCAGCTACCCCCTGAGCGAACCCGTCT",
    "TGTTGCAAAAGAGTGGGATGACTTGTTAGTAGAGGTGAAAAGCCT",
    "CCGAGCCGGGTGATAGCTGGTTGCCCGTGAAACGAATCTAAGTTC",
    "CTCTTAATTTTCTCTGCCCCCCGGACATCTAACCACAACCACCAT",
    "TGGCAAATCAAGAGCAATTTAAAGGAGGTACAGCTCCTTTAAAAG",
    "GAATACAACCTCCCCTAGCGGATAACTACCCAACTGACCCCAAAC",
    "GTAGGCCCTCAAGCAGCCATCAACAAAGAGTGCGTCAAAGCTCTA",
    "CCCAAAAAATCCAAAAACAATAGGACTCCCTTATCTATAGCGGGC",
    "AACCTATAATAATAGGAGAATTAATGCTAAAATGAGTAATCGGGG",
    "CCACCCTCTCGAGCGCAAACTTACATCACCACATTATTAACAAGC",
    "TTAAACTAATACTCCAACCTCAAACAAGCTAAGTATTAAACTTGC",
    "CTGTTAGCCCGACTCAGGAGCGCCCATTAGAAAGATTCAAATCTG",
    "AAAAGGAACTAGGCAACCCCAAGGCCCGACTGTTTACCAAAAACA",
    "AGCCTTCAGCCAACCAAGTATTGAAGGTGATGCCTGCCCAGTGAC",
    "CTCTGTTCAACGGCCGCGGTATCCTAACCGTGCGAAGGTAGCGCA",
    "TCAATTGTCCCATAAATCGAGACTTGTATGAATGGCTAAACGAGG",
    "CTTAACTGTCTCTTGCAGACAATCAATGAAATTGATCTTCCTGTG",
    "AAAAGCAGGAATAAACCCATAAGACGAGAAGACCCTGTGGAACTT",
    "AAAATCAGCGGCCACCACATACAACCCAAAACCTACCAGGCCTAC",
    "GCCCCAAGTAAAACACTGGCCCGCATTTTTCGGTTGGGGCGACCT",
    "GGAGAAAAACAAACCCTCCAAAAACAAGACCACACATCTTGACCA",
    "GAGCAACCTCTCAACGTACTAACAGCAACCAGACCCAATATAATT",
    "ACCAATGGACCAAGCTACCCCAGGGATAACAGCGCAATCCCCTCC",
    "AGAGCCCATATCGACAAGGGGGTTTACGACCTCGATGTTGGATCA",
    "GACACCCTAATGGTGCAGCCGCTATTAAGGGTTCGTTTGTTCAAC"]
    -- Kungsörn
    ), ("EAGLE", [
    "AAATACTCCCACCCAGCCTCCCCCTTAAACTCTAATCTGCGCACT",
    "AAATCCCACTATTCACCCCACCTTCTATCCTAAGCAACACATTCA",
    "ACACCCGCACCACATCCACAACAAACATACCCCCATTTCACCTCC",
    "AAGCCCCCACTTTCAAGACCAGGCTTCACAACACATACAATCTAA",
    "ACCAAGAACACTCCTACCTTCAATCAACACCCAATCCCAAAATCC",
    "ACTAACGAAACGGTTTCAAGTATCCAAAAACACCAAAACACTAAC",
    "CACGTTTTTCACGTATTTTTTCACATTTTTCACTTTCACATTTTT",
    "TTGTTCACCAAAAGCACTGGAGTCTCATTAATTTTTTCAAATCAC",
    "ATATCTCATATTTTATTTGTGTGTACACTAAACACACCAAAATAT",
    "TAAAGAGACCCCCCTACTGAATACACCAAACACATAACACCACCA",
    "AAATACTCCCACCCAGCCTCCCCCTTAAACTCTAAATAATCTGCG",
    "CACTAAATCCCACCATTCACCCCACCTTCTATCCTAAGCAACACA",
    "TTCAACACCCGCACCACATCCACAACAAACATACCCCCATTTCAC",
    "CTCCAAGCCCCCACTTTCAAGACCAGGCTTCACAACACATACAAT",
    "CTAAACCAAGAACACTCCTACCTTCAATCAACACCCAATCCCAAA",
    "ATCCACTAACGAAACGGTTTCAAGTATCCAAAAACACCAAAACAC",
    "TAACCACGTTTTTCACGTATTTTTTCACATTTTTCACTTTCACAT",
    "TTTTTTGTTCACCAAAAGCACTGGAGTCTCATTAATTTTTTCAAA",
    "TCACATATCTCATATTTTATTTGTGTGTACACTAAACACACCAAA",
    "ATATTAAAGAGACCCCCCTACTGAATACACCAAACACATAACACC"]
    -- Snok
    ), ("SNAKE", [
    "TTATTTAACCTACTACCAGTAGGCCTAAATATTTCAACCTGATGA",
    "AACTTCGGATCTATATTACTAACTTGTCTGGGACTACAAACCATC",
    "ACAGGTTTCTTCCTAGCAATTCATTATACAGCCAATATCAACTTA",
    "GCCTTCTCATCCATTATTCACATTACACGCGATGTGCCCTACGGA",
    "TGAATAATACAAAATACCCACGCAATTGGCGCATCAATATTTTTT",
    "ATCTGCATCTACACCCATATTGCACGTGGACTTTACTATGGCTCC",
    "TACCTTAACAAAGAAGTATGACTATCAGGAACAACCCTACTAATT",
    "ATCCTTATAGCCACAGCATTCTTCGGCTATGTCCTCCCATGAGGA",
    "CAAATATCCTTTTGAGCAGCGACAGTAATTACAAACCTCCTAACT",
    "GCCGTACCCTACCTAGGAAACACCCTCACAACCTGACTCTGAGGG",
    "GGATTCTCAATTAATGACCCAACCTTAACCCGATTCTTTGCCCTA",
    "CACTTTCTCCTCCCTTTTATCATTATTTCACTATCTTCAATCCAC",
    "ATTATACTTCTTCACACTGAAGGGTCAAGCAACCCCCTTGGAACA",
    "AACTCGGATATCGATAAAATCCCATTCCACCCCTACCACTCACAC",
    "AAAGACATACTATTATTTACCATAATAATTACCATATTATTCATT",
    "ATCCTATCATTCTCCCCCGATGTATTTAATGACCCAGAAAACTTT",
    "TCAAAAGCAAACCCTCTAGTAACCCCCCAACACATTAAGCCAGAA",
    "TGATATTTCCTATTTGCCTACGGCATTCTCCGATCTATCCCTAAT",
    "AAACTCGGGGGAACAATTGCCCTGGTTCTATCAATCGTTATCTTA"]
    -- Häst
    ), ("HORSE", [
    "ATGGACATCGCTATCCAGCACCCCTGGTTCAAACGCGCCCTGGGC",
    "CCCTTCTACCCCAGCCGGCTCTTCGACCAGTTTTTCGGCGAGGGC",
    "CTCTTTGAGTATGACCTGCTGCCCTTCCTGTCCTCCACCATCAGC",
    "CCCTACTACCGCCAGTCCCTCTTCCGCACCGTGCTGGACTCCGGC",
    "ATCTCTGAGGTCCGGTCTGACAGGGACAAGTTCGTCATCTTCCTG",
    "GATGTGAAGCACTTCTCTCCCGAGGACCTCACCGTGAAGGTGCAG",
    "GAGGATTTTGTGGAGATCCACGGCAAACACAACGAGAGGCAGGAC",
    "GACCATGGCTACATTTCCCGTGAGTTCCACCGTCGCTACCGCCTG",
    "CCTTCCAATGTGGACCAGACGGCGCTCTCCTGTTCCGTGTCCGCG",
    "GATGGCATGCTCACCTTCTCTGGCCCCAAGATCCCGTCCGGCATG",
    "GACGCCGGCCACAGCGAGCGAGCCATCCCCGTGTCTCGGGGTTGA",
    "CGCGATGTGATTTCTGCCCAGTGCTCTGAATGTCAAAGTGAAGAA",
    "AGTCAATGAAGCGCGAGTAGATGGCGGGAGTACCTATGACTCTCT",
    "TAAGGTAGCCAAATGCCTCGTCATCTAATTAGTGACGCGCATGAA",
    "TGGATGAAGGAGATTCCCACTGTCCCTAACTGTTATCCAGCGAAA",
    "CCACAGCCAAGGCAACGGGCTTGGCGGAATCAGCGGGGAAAGAAG",
    "ACCCTGTTGAGCTTGACTCTAGTCTGGCATGGTGAAGAGACATGA",
    "GAGGTGTAGAACAAGTGGAAGGCCCCCTTGCGGGCCGCCGGTGAA",
    "ATACCACTACTCTGATCGTTTTTTCACTGACCGGGTGAGGTGGGG",
    "GGCGAGCCCCAAGGGGCTCTCGCTTCTGGCGCCAAGCGCCCGGCC",
    "CGGCCGCCGGGCAAGACCCGCTCCGGGGACAGTGCCAGATGGGGA",
    "GTTTGACTGGGGCGGTACACCTGTCAAACGGTAACGCAGGTGTCC",
    "TAAGGCGAGCTCAGGGAGGACAGAAACCTCCCGTGGAGCAGAAGG",
    "GCAAAAGCTCGCTTGATCTTGATTTTCAGTACGAATACAGACGGT",
    "GAAAGCGGGGCCTCACGATCCTTCTGACCTTAGGGGTTTTAAGCA"]
    )]
-- Proteinsekvenser från människa, ko, hund, 
-- råtta, mus och groda.
-- Format [("NAME", sequence), ..]
foxp4 = [("FOXP4_HUMAN", 
    "EMSPAELLHFQQQQALQVARQFLLQQASGLSSPGNNDSKQSAVQVP\
    \VSVAMMSPQMLTPQQMQQILSPPQLQALLQQQQALMLQQLQEYYKK\
    \QQEQLHLQLLTQQQAGKPQPKEALGNKQLAFQQQLLQMQQLQQQHL\
    \LNLQRQGLVSLQPNQASGPLQTLPQAVCPTDLPQLWKGEGAPAEDS\
    \VKQEGLDLTGTAATSFAAPKVSPPLSHHTLPNGQPTRRDSSSHEET\
    \SPLYGHGECKWPGCETLCEDLGQFIKHLNTEHALDDRSTAQCRVQM\
    \QVVQQLEIQLRLQAMMAHLHMRPSEPKPFSQPVTVSADSFPDGLVH\
    \PPTSAAAPVTPLRPPGLGSASLHGGGPARRRSSDKFCSPISSELAQ\
    \NHEFYKNADVRPPFTYASLIRQAILETPDRQLTLNEIYNWFTRMFA\
    \YFRRNTATWKNAVRHNLSLHKCFVRVENVKGAVWTVDEREYQKRRP\
    \PKMTGSPTLVKNMISGLSYGALNASYQAALAESSFPLLNSPGMLNS\
    \ASSLLPLSHDDVGAPVEPLPSNGSSPRLSPQYSHQVQVKEEPAEED\
    \RQPGPLGAPNPSASGPPEDRDLEEELPGEEL"),
    ("FOXP4_COW",
    "EMSPAELLHFQQQQALQVARQFLLQQASGLSSPGNNDSKQSAVQVP\
    \VSVAMMSPQMLTPQQMQQILSPPQLQALLQQQQALMIQQLQEYYKK\
    \QQEQLHLQLLTQQQAGKQQPKEALGNKQLAFQQQLLQMQQLQQQHL\
    \LNLQRQGLVSLQPSQASGPLQTLPQAVCPTDLPQLWKGEGAPAEDS\
    \VKQEGLDLTGTATTSFAAPKVSPPLSHHTLPNGQPTRRDSWGLALT\
    \AAIVGAGGLLLPGHTKLCSACGEPVRHLNTEHALDDRSTAQCRVQM\
    \QVVQQLEIQVWPQAVGAGRGGAPARPKPFSQPVTVSADSFPDGLAH\
    \PPTSAAAPVTPLRPPGLGSASLHSGGPARRRSSDKFCSPISSELAQ\
    \NHEFYKNADVRPPFTYASLIRQAILETPDRQLTLNEIYNWFTRMFA\
    \YFRRNTATWKNAVRHNLSLHKCFVRVENVKGAVWTVDEREYQKRRP\
    \PKMTGSPTLVKNMISGLSYGTLNASYQAALAESSFPLLNSPGMLNS\
    \ASSLLPLGHDDAGAPVEPLPSNGSSPRLSPQYSHQVQVKEEPAEED\
    \RRPGPMGPPNPSTAGPPEDRDLEEELPGEEL"),
    ("FOXP4_DOG", 
    "EMSPAELLHFQQQQALQVARQLLLQQASGLSSPGNNDSKQSAVQVP\
    \VSVAMMSPQMLTPQQMQQILSPPQLQALLQQQQALMLQQLQEYYKK\
    \QQEQLHLQLLTQQQAGKQQPKEGVGRADCTFRDALLPTWSSPQQHR\
    \RNQQRQGLVSLQPSQASGPLQTLPQAVCPTDLPQLWKGEGAPAEDS\
    \VKQEGLDLTGSATTSFAAPKVSPPLSHHTLPNGQPTRRDSSSHEET\
    \SPLYGHGECKWPGCETLCEDLGQFIKHLNTEHALDDRSTAQCRVQM\
    \QVVQQLEIQLRLQAMMAHLHMRPSEPKPFSQPVTVSADSFPDGLVH\
    \PPTSAAAPVTPLRPPGLSSASLHSGGPARRRSSDKFCSPISSELAQ\
    \NHEFYKNADVRPPFTYASLIRQAILETPDRQLTLNEIYNWFTRMFA\
    \YFRRNTATWKNAVRHNLSLHKCFVRVENVKGAVWTVDEREYQKRRP\
    \PKMTGSPTLVKNMISGLSYGALNASYQAALAESSFPLLNSPGMLNS\
    \ASSLLPLSHDEVGAPVEPLPSNGSSPRLSPQYSHQVQVKEEPAEED\
    \RRPGPLGPPNPSTAGPPEDRDLEEELPGEEL"),
    ("FOXP4_RAT", 
    "EMSPAELLHFQQQQALQVARQFLLQQASSLNSPGNNDSKQSAVQVP\
    \VSVAMMSQQMLTPQQMQQILSPPQLQALLQQQQALMLQQLQEYYKK\
    \QQEQLHLQLLSQQQAGKQQPKEALGNKQLAFQQQLLQMQQLQQQHL\
    \LNLQRQGLVSLQPSQASGPLQALPQAVCPTDLPQLWKGEGAPAEDG\
    \GRQEGLDLASPAATSFASPKVSPPLSHHPLPNGQPTRRDSSSHEET\
    \SPLYGHGECKWPGCETLCEDLGQFIKHLNTEHALDDRSTAQCRVQM\
    \QVVQQLEIQLRLQAMMAHLHMRPSEPKPFSQPVTVSADPFPDGLVH\
    \PPTSAAAPVTPLRPPGLGSASLHGGGPARRRSNDKFCSPISSELAQ\
    \NHEFYKNADVRPPFTYASLIRQAILETPDRQLTLNEIYNWFTRMFA\
    \YFRRNTATWKNAVRHNLSLHKCFVRVENVKGAVWTVDEREYQKRRP\
    \PKMTGSPTLVKNMISGLSYGALNASYQAALAESSFPLLSSPGMLNS\
    \ASSLLPLSQDDMGAPGEPLPSNGSSPRLSPQYSHQIQVKEEPAEED\
    \RRPGPLGAPNPSTVGPPEDRDLEEDLAGEDI"),
    ("FOXP4_MOUSE", 
    "EMSPAELLHFQQQQALQVARQFLLQQASSLNSPGNNDSKQSAVQVP\
    \VSVAMMSQQMLTPQQMQQILSPPQLQALLQQQQALMLQQLQEYYKK\
    \QQEQLHLQLLTQQQAGKQQPKEALGNKQLAFQQQLLQMQQLQQQHL\
    \LNLQRQGLVSLQPSQASGPLQALPQAVCPTDLPQLWKGEGAPAEDS\
    \GRQEGLDLASTAVTSFASPKVSPPLSHHPLPNGQPTRRDSSSHEET\
    \SPLYGHGECKWPGCETLCEDLGQFIKHLNTEHALDDRSTAQCRVQM\
    \QVVQQLEIQLRLQAMMAHLHMRPSEPKPFSQPVTVSADPFPDGLVH\
    \PPTSAAAPVTPLRPPGLGSASLHSGGPARRRSNDKFCSPISSELAQ\
    \NHEFYKNADVRPPFTYASLIRQAILETPDRQLTLNEIYNWFTRMFA\
    \YFRRNTATWKNAVRHNLSLHKCFVRVENVKGAVWTVDEREYQKRRP\
    \PKMTGSPTLVKNMISGLSYGALNASYQAALAESSFPLLSNPGMLNS\
    \ASSLLPLSQEDLGVPGEPLPSNGSSPRLSPQYSHQIQVKEEPAEED\
    \RRPGPLGAPNPSTVGPPEDRDLEEDLGGEDI"),
    ("FOXP4_FROG", 
    "ELSPAELLHFQQQQALQMARQLLLQQATGLSSPSSTDNKQPSVQVP\
    \VSVAMMSPGMITPQQMQQILSPTQLQAVLQQQQALMLQQLQEYYKK\
    \QQEQLHLQLLSQQQAGKQQPKESLGNKQLAFQQQLLQMQQLQQQHL\
    \LNLQRQNLVGLQSGQGPLPIQSLPQAVSPSDLQQLLKEMSSNQEES\
    \SKQDTVDLTTSITTSFPNSKVSLPTIHPSLPNGQNTRRDSMSHYES\
    \SPLYGHGECRWPGCEALCEDMGQFIKHLNTEHALDDRSTAQCRVQM\
    \QVVQQLEIQLRLQAMMTHLHMRPSEPKPFSQPNKMSPDTFPDGLPQ\
    \PPTSATAPITPLRTSVISSSSLPSVGPVRRRIVDKFSTPISSELAQ\
    \NHEFYKNAEVRPPFTYASLIRQAILDTPDRQLTLNEIYNWFTRMFA\
    \YFRRNTATWKNAVRHNLSLHKCFVRVENVKGAVWTVDELEYQKRRP\
    \PKMTGSPTLVKNMISGLGYSALNASYQAALAESSFPLLNSPPLHNS\
    \SGSVLHGGHDDVTSTGEPGNSNGSSPRLSPQYSQSIHVKEEPAEDD\
    \VRPASLSAPTNQTTVLPEDRDIEPETPMEDL")]