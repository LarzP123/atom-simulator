{-# LANGUAGE RankNTypes #-}
module AtomSimulator.MathHelper where

import AtomSimulator.Units
import Data.Metrology.Poly
import Numeric.LinearAlgebra
import Data.Ord
import Data.List
import Data.Array
import GHC.Natural

{-| An array represnting a 3D space. Each point contains a 3D coordinate identifier for its corresponding space on
the grid and a value contained at that point -} 
type Grid3D a = Array Idx3 a

-- | A point representing a given coordinate on a 3D grid
type Idx3 = (Natural, Natural, Natural)

-- | Gets the difference between 2 natural numbers (subtraction and absolute value all in one)
(|--|) :: Natural -> Natural -> Natural
(|--|) a b
  | a >= b    = a - b
  | otherwise = b - a

-- | Determines number of moves along an axis need to be taken one step in order to get between points
manhattan :: Idx3 -> Idx3 -> Natural
manhattan (i1,j1,k1) (i2,j2,k2) = (i1|--|i2) + (j1|--|j2) + (k1|--|k2)

{-| A grid representing physical space. Contains a length saying the spacing between each item in the grid
  Also contains a 3d grid/triple array containing all of the points in the grid with arbitrary units.-}
data PhysicalGrid a = PhysicalGrid { spacing :: Length, samples :: Grid3D a } deriving (Show)

{-| Given 1. a function that takes a 3d point in space and outputs a value,
2. a integer specifying the spacing between each point in a grid,
3. a physical distance between each point in a grid
outputs a physical grid containning values calculated from the function at each spot in space -}
createGridFromFunc :: (Length -> Length -> Length -> a) -> Natural -> Length -> PhysicalGrid a
createGridFromFunc f n s = PhysicalGrid s (listArray ((0,0,0),(n-1,n-1,n-1)) values)
  where
    values = [ f (fromIntegral x |*| s) (fromIntegral y |*| s) (fromIntegral z |*| s)
      | x <- [0..n-1], y <- [0..n-1], z <- [0..n-1] ]

{-| This will find the eigen values and eigen vectors of a matrix and outputs them in order
from lowest eigenvalue to highest eigenvalue. This handles Matrixes of dimensioned units.-}
solveMatrix :: forall dim lcsu unit. ValidDLU dim lcsu unit => unit ->
  [[Qu dim lcsu Double]] -> ([Qu dim lcsu Double], [[Double]])
solveMatrix unit matrix = (vals', vecs')
  where
    n            = length matrix
    bareMatrix   = (n >< n) (foldMap (map (# unit)) matrix)
    (vals, vecs) = eigSH (trustSym bareMatrix)
    pairs = sortBy (comparing fst) $ zip (toList vals) (map toList (toColumns vecs))
    vals' = map ((% unit) . fst) pairs
    vecs' = map snd pairs
