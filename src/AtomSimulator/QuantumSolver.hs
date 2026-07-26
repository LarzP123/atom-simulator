{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeOperators #-}
module AtomSimulator.QuantumSolver where

import AtomSimulator.Units
import AtomSimulator.MathHelper
import Data.Array (Array, range)
import Data.Array.IArray (listArray, elems, (!), IArray (bounds))
import Data.Metrology.Poly
import Data.List
import Data.Ord
import GHC.Natural
import GHC.RTS.Flags (ProfFlags(ccsLength))
import Data.Metrology.Show

{- Creates a 3D potential energy grid from a function that says a given potential energy at a point in space -}
createPotentialGrid :: (Length -> Length -> Length -> Energy) -> Natural -> Length -> PhysicalGrid Energy
createPotentialGrid = createGridFromFunc

{- Given a potential energy 3D space, and any particle specific properties (like mass) this will return an energy
matrix for solving for possible solutions to this system -}
hamiltonianMatrix :: PhysicalGrid Energy -> [[Energy]]
hamiltonianMatrix (PhysicalGrid s v) = [ [ f p q | q <- pts ] | p <- pts ]
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

hamiltonianVector :: forall dim lcsu unit. (Energy ~ Qu dim lcsu Double, ValidDLU dim lcsu unit)
  => unit -> PhysicalGrid Energy -> (Vec -> Vec, unit)
hamiltonianVector unit (PhysicalGrid s v) = (matvec, unit)
  where
    matvec psi = listArray (bounds v)
      [ (6 |*| k |*| psiP |+| (v ! p |*| psiP) |-| (k |*| (neighborSum psi p % Number))) # unit
      | p <- range $ bounds v, let psiP = psi ! p % Number ]
    kineticCoeff = (0.0381 % ElectronVolt) |*| qSq (1 % Nanometer)
    k :: Energy
    k = kineticCoeff |/| qSq s
    neighborSum psi p = sum [ psi ! q | q <- neighborsOf p (snd $ bounds v) ]

type Orbital = (Energy, PhysicalGrid Density)

waveToProb :: Length -> Double -> Density
waveToProb s = (|/| (s |*| s |*| s)) . (% Number) . (^2)

{- Given a grid saying the potential energy in the system, and any particle specific properties (like mass) this will
solve for possible eigen energies the particle can hold and also for probability densities of the particle at that given
energy (the eigen vectors) -}
solveAllOrbitals :: PhysicalGrid Energy -> [Orbital]
solveAllOrbitals (PhysicalGrid s v) = sortBy (comparing fst) (zip energies densityGrids)
  where
    (energies, waveFuncs) = solveMatrix ElectronVolt (hamiltonianMatrix (PhysicalGrid s v))
    waveFuncGrids :: [Grid3D Double]
    waveFuncGrids = map (listArray (bounds v)) waveFuncs
    densitiesGrid3D :: [Grid3D Density]
    densitiesGrid3D = map (fmap $ waveToProb s) waveFuncGrids
    densityGrids :: [PhysicalGrid Density]
    densityGrids = map (PhysicalGrid s) densitiesGrid3D

solveLowestOribtal :: PhysicalGrid Energy -> Orbital
solveLowestOribtal grid@(PhysicalGrid s v) = (energy, densityGrid)
  where
    (energy, waveFunc) = lanczosLowest (hamiltonianVector ElectronVolt grid) (snd (bounds v))
    waveFuncGrid  :: Grid3D Double
    waveFuncGrid = listArray (bounds v) (elems waveFunc)
    densityGrid3D :: Grid3D Density
    densityGrid3D = fmap (waveToProb s) waveFuncGrid
    densityGrid :: PhysicalGrid Density
    densityGrid = PhysicalGrid s densityGrid3D

coulombPotentialFromDensity :: Charge -> Charge -> PhysicalGrid Density -> PhysicalGrid Energy
coulombPotentialFromDensity q1 q2 (PhysicalGrid spacing grid) = 
  let
    potentialGrid :: Grid3D Energy
    potentialGrid = listArray (bounds grid) [ potentialAt p | p <- allPoints ]
  in PhysicalGrid spacing potentialGrid
  where
    allPoints = range $ bounds grid
    amtOfParticleAtPoint :: Idx3 -> Double
    amtOfParticleAtPoint point = (grid ! point) |*| (spacing |*| spacing |*| spacing) # Number
    sumInverseDists :: Idx3 -> InverseLength
    sumInverseDists p = foldr1 (|+|) [ (amtOfParticleAtPoint q % Number) |/| cartDist p q spacing | q <- allPoints, q /= p ]
    relCoulombConst :: EnergyDist
    relCoulombConst = q1 |*| q2 |/| (4 |*| pi |*| vacuumPermittivity)
    potentialAt :: Idx3 -> Energy
    potentialAt p = relCoulombConst |*| sumInverseDists p