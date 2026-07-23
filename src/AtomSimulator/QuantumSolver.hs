{-# LANGUAGE FlexibleContexts #-}
module AtomSimulator.QuantumSolver where

import AtomSimulator.Units
import AtomSimulator.MathHelper
import Data.Array
import Data.Metrology.Poly
import Data.List
import Data.Ord
import GHC.Natural

{- Creates a 3D potential energy grid from a function that says a given potential energy at a point in space -}
createPotentialGrid :: (Length -> Length -> Length -> Energy) -> Natural -> Length -> PhysicalGrid Energy
createPotentialGrid = createGridFromFunc

{- Given a potential energy 3D space, and any particle specific properties (like mass) this will return an energy
matrix for solving for possible solutions to this system -}
hamiltonian3D :: PhysicalGrid Energy -> [[Energy]]
hamiltonian3D (PhysicalGrid s v) = [ [ f p q | q <- pts ] | p <- pts ]
  where
    pts = range ((0,0,0), snd $ bounds v)
    kineticCoeff = (0.0381 % ElectronVolt) |*| qSq (1 % Nanometer)
    k :: Energy
    k = kineticCoeff |/| qSq s
    f :: Idx3 -> Idx3 -> Energy
    f p q
      | p == q             = (6 |*| k) |+| (v ! p)
      | manhattan p q == 1 = (-1) |*| k
      | otherwise          = 0 % ElectronVolt

{- Given a grid saying the potential energy in the system, and any particle specific properties (like mass) this will
solve for possible eigen energies the particle can hold and also for probability densities of the particle at that given
energy (the eigen vectors) -}
solve3D :: PhysicalGrid Energy -> [(Energy, PhysicalGrid Density)]
solve3D (PhysicalGrid s v) = sortBy (comparing fst) (zip energies densityGrids)
  where
    (energies, waveFuncs) = solveMatrix ElectronVolt (hamiltonian3D (PhysicalGrid s v))
    waveFuncGrids :: [Grid3D Double]
    waveFuncGrids = map (listArray ((0,0,0),snd $ bounds v)) waveFuncs
    densitiesGrid3D :: [Grid3D Density]
    densitiesGrid3D = map (fmap ((|/| (s |*| s |*| s)) . (% Number) . (^2))) waveFuncGrids
    densityGrids :: [PhysicalGrid Density]
    densityGrids = map (PhysicalGrid s) densitiesGrid3D