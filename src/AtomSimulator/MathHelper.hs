{-# LANGUAGE RankNTypes #-}
module AtomSimulator.MathHelper where

import AtomSimulator.Units
import Data.Metrology.Poly
import Numeric.LinearAlgebra
import Data.Ord
import Data.List
import Data.Array (Array)
import Data.Array.IArray
import GHC.Natural
import Data.Array.Base (amap, UArray)
import Data.Kind (Type)

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

-- | Distance in cartesian space between 2 points
cartDist :: Idx3 -> Idx3 -> Length -> Length
cartDist (i1,j1,k1) (i2,j2,k2) spacing = spacing |*| sqrt(fromIntegral $ (i1|--|i2)^2 + (j1|--|j2)^2 + (k1|--|k2)^2)

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

type Vec = UArray Idx3 Double

dotA :: Vec -> Vec -> Double
dotA a b = sum (zipWith (*) (elems a) (elems b))

normA :: Vec -> Double
normA a = sqrt (dotA a a)

scaleA :: Double -> Vec -> Vec
scaleA c = amap (* c)

subA, addA :: Vec -> Vec -> Vec
subA a b = listArray (bounds a) (zipWith (-) (elems a) (elems b))
addA a b = listArray (bounds a) (zipWith (+) (elems a) (elems b))

startingVector :: Idx3 -> Vec
startingVector bnd =
  let unNormalized = listArray ((0,0,0),bnd) [ 0.5 + 0.3 * sin (fromIntegral i) | i <- [1 .. rangeSize ((0,0,0),bnd)] ]
  in amap (/ normA unNormalized) unNormalized 

go :: Int -> Vec -> Vec -> Double -> [Vec] -> (Vec -> Vec) -> ([Double], [Double], [Vec])
go 0 _ _ _ _ _ = ([], [], [])
go j q qPrev betaPrev qsRev matvec =
  let w0    = matvec q
      alpha = dotA w0 q
      w1    = w0 `subA` scaleA alpha q `subA` scaleA betaPrev qPrev
      w2    = foldl' (\acc q -> acc `subA` scaleA (dotA acc q) q) w1 (q : qsRev)
      beta  = normA w2
      qNext = if beta < 1e-10 then zeroVec q else scaleA (1 / beta) w2
      (alphas, betas, qs) = go (j - 1) qNext q beta (q : qsRev) matvec
  in (alpha : alphas, beta : betas, q : qs)

zeroVec :: Vec -> UArray Idx3 Double
zeroVec = amap (const 0)

-- | Lowest eigenpair (ground state) via Lanczos.
lanczosLowest :: forall dim lcsu unit. (ValidDLU dim lcsu unit)
  => (Vec -> Vec, unit) -> Idx3 -> (Qu dim lcsu Double, Vec)
lanczosLowest (matvec, unit) bnd =
  let
    (alphas, betas, qs) = go m v0 (zeroVec v0) 0 [] matvec

    tridiag = assoc (m, m) 0 $
      [ ((i, i), a) | (i, a) <- zip [0 ..] alphas ] ++
      concat [ [((i, i+1), b), ((i+1, i), b)] | (i, b) <- zip [0 ..] (init betas) ]

    (thetaAll, sVecs) = eigSH (trustSym tridiag)
    theta = last (toList thetaAll)
    s = last (toColumns sVecs)

    y = foldl1 addA (zipWith scaleA (toList s) qs)
  in (theta % unit, y)
  where
    m = 30
    v0 = startingVector bnd

neighborsOf :: Idx3 -> Idx3 -> [Idx3]
neighborsOf (i, j, k) (iMax, jMax, kMax) = concat
      [ [(i - 1, j, k) | i > 0], [(i + 1, j, k) | i < iMax]
      , [(i, j - 1, k) | j > 0], [(i, j + 1, k) | j < jMax]
      , [(i, j, k - 1) | k > 0], [(i, j, k + 1) | k < kMax]
      ]
