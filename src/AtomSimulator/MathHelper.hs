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
data PhysicalGrid a = PhysicalGrid { spacing :: Length, samples :: Grid3D a } deriving (Show, Eq)

-- | Like createGridFromFunc but takes a 0 centered function and makes it centered at middle of grid we are creating
createGridFromFuncShifted :: (Length -> Length -> Length -> a) -> Natural -> Length -> PhysicalGrid a
createGridFromFuncShifted f n s = createGridFromFunc shifted n s
  where
    offset = (fromIntegral n / 2) |*| s
    shifted x y z = f (x |-| offset) (y |-| offset) (z |-| offset)

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

-- | Given an existing Vector, creates a Zero vector matching the same samples
zeroVec :: Vec -> UArray Idx3 Double
zeroVec = amap (const 0)

{- | Projects a new trial vector onto the set of existing known vectors. Then it subtracts those vectors from the trial vector
to get back a new trial vector. This is to try to eliminate orbital sharing via paul-exclusion principle. -}
projectOut :: Vec -> [Vec] -> Vec
projectOut = foldl' (\acc q -> acc `subA` scaleA (dotA acc q) q)

-- | Runs a lanczos iteration up to j steps. The more steps, the closer to actually solving a matrix this is
lanczosIter :: [Vec] -> Int -> Vec -> Vec -> Double -> [Vec] -> (Vec -> Vec) -> [(Double, Double, Vec)]
lanczosIter _ 0 _ _ _ _ _ = []
lanczosIter occupied j q qPrev betaPrev qsRev matvec =
  let 
    w0 :: Vec
    w0 = matvec q
    alpha :: Double
    alpha = dotA w0 q
    w1 :: Vec
    w1 = w0 `subA` scaleA alpha q `subA` scaleA betaPrev qPrev
    w2 :: Vec
    w2 = projectOut (foldl' (\acc qq -> acc `subA` scaleA (dotA acc qq) qq) w1 (q : qsRev)) occupied
    beta :: Double
    beta = normA w2
    qNext :: Vec
    qNext = if beta < 1e-10 then zeroVec q else scaleA (1 / beta) w2
  in (alpha, beta, q) : lanczosIter occupied (j - 1) qNext q beta (q : qsRev) matvec

-- | Lowest eigenpair (ground state) via Lanczos.
lanczosLowestExcluding :: forall dim lcsu unit. (ValidDLU dim lcsu unit)
  => [Vec] -> (Vec -> Vec, unit) -> Idx3 -> (Qu dim lcsu Double, Vec)
lanczosLowestExcluding occupied (matvec, unit) bnd =
  let
    m'   = min m (rangeSize (bounds v0) - length occupied)
    v0'  = orthonormalStart (startingVector bnd) occupied
    (alphas, betas, qs) = unzip3 $ lanczosIter occupied m' v0' (zeroVec v0') 0 [] matvec

    tridiag = assoc (m', m') 0 $
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

-- | Normalizes a vector to unit length, returning the zero vector instead if its norm is too small
normVec :: Vec -> Vec
normVec v =
  let n = normA v
  in if n < 1e-8 then zeroVec v else scaleA (1 / n) v

-- | Projects a vector out of the subspace spanned by a list of occupied vectors (making a new orthogonal vector) and normalizes
orthonormalStart :: Vec -> [Vec] -> Vec
orthonormalStart = (normVec .) . projectOut

-- | Lowest eigenpair (ground state) via Lanczos.
lanczosLowest :: forall dim lcsu unit. (ValidDLU dim lcsu unit)
  => (Vec -> Vec, unit) -> Idx3 -> (Qu dim lcsu Double, Vec)
lanczosLowest = lanczosLowestExcluding []

neighborsOf :: Idx3 -> Idx3 -> [Idx3]
neighborsOf (i, j, k) (iMax, jMax, kMax) = concat
      [ [(i - 1, j, k) | i > 0], [(i + 1, j, k) | i < iMax]
      , [(i, j - 1, k) | j > 0], [(i, j + 1, k) | j < jMax]
      , [(i, j, k - 1) | k > 0], [(i, j, k + 1) | k < kMax]
      ]

-- | Adds 2 grids with a unit type together Currently assumes same spacing too
addGrids :: PhysicalGrid (Qu dim lcsu Double) -> PhysicalGrid (Qu dim lcsu Double) -> PhysicalGrid (Qu dim lcsu Double)
addGrids (PhysicalGrid s a) (PhysicalGrid _ b) =
  PhysicalGrid s (listArray (bounds a) (zipWith (|+|) (elems a) (elems b)))

-- | A function to find the potential energy at a given point in a harmonic system
harmonicFunc :: SpringConstant -> Length -> Length -> Length -> Energy
harmonicFunc springConstant x y z =
  let
    r2 :: DistSq
    r2 = (x |*| x) |+| (y |*| y) |+| (z |*| z)
  in 0.5 |*| springConstant |*| r2

-- | Given two category labels and a count, generate every way to split that many items avoiding mirror-image duplicates
binaryPartitions :: a -> a -> Natural -> [[a]]
binaryPartitions labelA labelB total =
  [ replicate (fromIntegral numA) labelA ++ replicate (fromIntegral (total - numA)) labelB
  | numA <- [0 .. total `div` 2]
  ]
