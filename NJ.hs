-- Bastian Fredriksson 11/29-13

-- A module NJ that implements the Neighbour Joining algorithm
-- used to create phylogenetic trees.

module NJ(neighbor) where

-- Import dependencies, declare Data.Set and Data.Map as 
-- qualified to avoid name clashes with the default Prelude 
-- package.
import qualified Data.Set as Set
import qualified Data.Map as Map
import Data.List
import Data.Maybe

-- A corresponding distance map, where each pair of species
-- are mapped to their respective evolutionary distance.
type DistanceMap = Map.Map (String, String) Float

-- A map with all edges in a graph. An edge is a connection between
-- two nodes in a graph. This type maps the contents of one node 
-- to the contents of the adjacent node.
type Edges = Map.Map String String

-- A set of nodes in a graph. Each node contains the name of the
-- species if the node is a leaf, or the name of the subtree, if
-- the node has two children. 
type Vertices = Set.Set String

-- A kit used to parse a Newick tree.
type Tree = (Vertices, Edges)

-- The label of a node in the tree.
nodeLabel = "v"

-- The label of the topmost vertex in the tree.
topLabel = ""

-- Convert a distance matrix to a distance map. For example given the
-- distance matrix          | fish  | human |  bird |
--                   fish   | 0     |       |       |
--                   human  | 23    | 0     |       |
--                   bird   | 12    | 27    |  0    |
-- and the list of species [fish, human, bird]
-- the corresponding distance map will be (fish, human) -> 23,
-- (fish, bird) -> 12, (bird, human) -> 27. Note that the 
-- number of elements in the distance map is exactly
-- (n^2-n) / 2 where n is the number of species.
makeDistanceMap :: [String] -> [[Float]] -> DistanceMap
makeDistanceMap s dm = res
  where
    -- Zip values in distance matrix with a pair of species
    -- and remove any equal pairs.
    tmp = filter (\(a, b, c) -> a /= b) $ concat
        $ map (\i -> zip3 (replicate (length s) $ s !! i) s $ dm !! i)
        [0.. (length s) - 1]
    -- Create a map from the list. Sort the tuple (s1, s2) first,
    -- to make sure that the map doesn't contain any duplicates.
    res = Map.fromList $ map (\(s1, s2, v) -> 
        ((sortTuple (s1, s2)), v)) tmp
        
-- Sort a tuple in ascending alphabetic order.
sortTuple :: (String, String) -> (String, String)
sortTuple (a1, a2) = if a1 > a2 then (a2, a1) else (a1, a2)

-- The selection function. Variable v is a set of vertices in the 
-- Newick tree we are building, dm is a distance map created from 
-- a distance matrix, and s1, s2 are species.
select :: Vertices -> DistanceMap -> String -> String -> Float
select v dm s1 s2 = (fromIntegral $ Set.size v - 2) * 
    (lookupValue s1 s2 dm) -
    -- Use list comprehension to sum up the difference
    -- D(s1,z) + D(s2,z) where D is the distance map and
    -- z is a member in the set of vertices.    
    sum [(lookupValue s1 z dm) + (lookupValue s2 z dm) 
	    | z <- Set.elems v]

-- Lookup a value in the distance map dm using the tuple (s1, s2)
-- as key. Returns the corresponding value if the key was found or
-- the default value zero if the key was not found.
lookupValue :: String -> String -> DistanceMap -> Float
lookupValue s1 s2 dm = fromMaybe 0 $ 
    Map.lookup (sortTuple (s1, s2)) dm

-- Render a phylogenetic tree using a list of species and
-- a distance matrix describing the evolutionary distance
-- between each pair of species.
neighbor :: [String] -> [[Float]] -> String
neighbor species dm = tree
  where
    -- Create a distance map to use
    di = makeDistanceMap species dm
    -- Create vertices to reduce as a set of strings from
    -- our list of species.
    fi = Set.fromList species
    -- Vertices in our Newick tree (will be a bunch of species
    -- to begin with)
    vi = fi
    -- Edges in our Newick tree (will be empty in the beginning
    -- because all species are together in the root of the tree)
    ei = Map.empty
    -- Build the tree using neighbour joining and parse it to
    -- Newick format.
    tree = parse $ build di fi vi ei 1

-- Build the tree using neighbour joining. This function only
-- implements the recursive part, which is made up by the
-- following steps:
-- If F_i > 3
--  (a, b) := min Select(x, y) where x and y are members of F_i
--  v_i = new vertex
--  F_i+1 := V_i union {v_i}
--  E_i+1 := E_i union {(a, v_i), (b, v_i)}
--  D_i+1 := 
--   such as
--    D_i+1(x, y) = D_i(x, y) 
--     where x and y are members of F_i \ {a, b}
--    D_i+1(x, v_i) = (D_i(x, a) + D_i(x, b) / 2
--     where x is a member of V_i \ {a, b}
--  increment i and recurse
-- Else
--  create top of tree and link edges to vertices in F_i
build :: DistanceMap -> Vertices -> Vertices -> Edges -> Int -> Tree
build di fi vi ei i
    | Set.size fi > 3 =
      -- The tree is not ready yet, and we need to reduce some of
      -- its nodes.
      let
        -- Pick the pair of species s1 and s2 that minimizes the
        -- selection function.
        (a, b) = snd (minimum [(select fi di s1 s2, (s1, s2)) | 
            s1 <- Set.elems fi, s2 <- Set.elems fi, s1 /= s2])
		-- Create reduced set of species F_i & V_i \ {a, b}
        fiRed = Set.difference fi $ Set.fromList [a, b]
        viRed = Set.difference vi $ Set.fromList [a, b]
        -- Subtract a and b found above from our set of
        -- vertices, and add a new branch.
        fin = Set.union fiRed $ Set.singleton stamp
        -- Add branch to our Newick tree holding a and b.
        vin = Set.union viRed $ Set.singleton stamp
        -- Add two new edges to the Newick tree, where the
        -- first edge is leading to a and the other edge to b.
        ein = Map.union ei $ Map.fromList [(a, stamp), (b, stamp)]
        -- Create a new distance map with the key (a, b) removed,
		-- and the new "row" of keys (x, v_i) added, where x are
		-- the set of elements from V_i \ {a, b}.
        din = Map.fromList
	    -- Create a new map with all remaining pairs of species
		-- in fiRed. Remove pairs where species are equal.
		-- Remember to store the key tuple sorted, since
		-- lookupValue will sort before making the lookup.
		-- For example, insertion of the key (b, a) is invalid
		-- because lookupValue will reorder the elements as (a, b)
		-- and then make the lookup, resulting in Nothing.
            ([((sortTuple (x, y)), (lookupValue x y di)) |
            x <- Set.elems fiRed, y <- Set.elems fiRed, x /= y]++
            -- Add the new set of key-value pairs describing
            -- the distance between each branch or vertex in the
            -- tree and the newly added branch.
            [((sortTuple (x, stamp)), ((lookupValue x a di 
            + lookupValue x b di) / 2)) | x <- Set.elems viRed])
      in
        build din fin vin ein $ i + 1
    | otherwise =
        -- Link the top of the tree to the three remaining nodes.
        (vi, Map.union ei $ Map.fromList $ zip (Set.toList fi) 
            $ replicate 3 topLabel)
    where
      -- Branches does not contain any species, but should have
      -- some kind of stamp telling the name of the subtree.
      stamp = nodeLabel++show i
      
-- Parse a tree kit to Newick format. This method is probably slow
-- as molasses in January.
parse :: Tree -> String
parse (v, e) = res
  where
    -- A leaf should be parsed as the content of the leaf.
	-- A branch should be parsed as (..) where .. is the content of
	-- the branch. Begin to parse the top nodes of the tree defined
	-- in the list of vertices. Parse each of the subtrees recursively
    -- using the list of edges until the whole tree is parsed.
    nodes = Set.toList v; edges = Map.toList e
    res = "("++(parseNodes nodes edges)++");"
    parseNodes [] _ = "" -- Base case
    parseNodes (h:t) e
        -- Don't append separator if tail is null
        | isPrefixOf nodeLabel h && null t =
        "("++(parseNodes (getChildren h e) e)++")"++parseNodes t e
        | isPrefixOf nodeLabel h =
        -- The node is a branch, extract its children using the
        -- list of edges, evaluate the result recursively and parse 
        -- the rest of the tree.
        "("++(parseNodes (getChildren h e) e)++")"++
            ", "++parseNodes t e
        | null t = h++parseNodes t e
        | otherwise = h++", "++parseNodes t e
    getChildren key edges =
        -- Find children that has key as branch
        [child | (child, branch) <- edges, branch == key]
