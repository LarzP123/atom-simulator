{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeOperators #-}
{-# OPTIONS_GHC -Wmissing-signatures -Wmissing-local-signatures #-}
{-# LANGUAGE InstanceSigs #-}
module AtomSimulator.QuantumSolver where

import AtomSimulator.Units
import AtomSimulator.MathHelper
import Data.Array (Array, range)
import Data.Array.IArray (listArray, elems, (!), IArray (bounds))
import Data.Metrology.Poly
import Data.List
import Data.Ord
import GHC.RTS.Flags (ProfFlags(ccsLength))
import Data.Metrology.Show
import Data.Function (on)
import Text.Printf (printf)
import Numeric.LinearAlgebra (Convert(double))

-- ===REGION Creating Potential Energy Grids
-- | Creates a 3D potential energy grid from a function that says a given potential energy at a point in space
createPotentialGrid :: (Length -> Length -> Length -> Energy) -> Int -> Length -> PhysicalGrid Energy
createPotentialGrid = createGridFromFuncShifted

-- | Finds the harmonic potential at a position given a spring constant
harmonicPotential :: SpringConstant -> Length -> Length -> Length -> Energy
harmonicPotential = harmonicFunc

-- | Given an electron density, this will find the potential energy for another electron to exist nearby
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
-- ===ENDREGION

-- ===REGION Wave Solution Definitions
-- | The solution for an electron cloud. This contains the eigen-energy of the particle, a probability density, and the spin of the particle if from a multi-particle solver
newtype WaveSolution = WaveSolution (Energy, PhysicalGrid Density, Maybe Spin) deriving (Eq)
{- | The wave function of an electron (psi). This contains the eigen-energy of the particle, a vector containning indexes and results of the
wavefunction at various points in space, and the particles spin -}
type WaveState = (Energy, Vec, Spin)

instance Show WaveSolution where
  show :: WaveSolution -> String
  show (WaveSolution (e, _, spin)) = case spin of
    Just s  -> show s ++ " " ++ show e
    Nothing -> show e

-- | PLEASE DELETE ME
showInHartree :: WaveSolution -> String
showInHartree (WaveSolution (e, _, spin)) = case spin of
  Just s  -> show s ++ " " ++ printf "%.4f Hartree" (e # Hartree)
  Nothing -> printf "%.4f Hartree" (e # Hartree)

instance Ord WaveSolution where
  compare = comparing (\(WaveSolution (e, _, _)) -> e)

-- | Given a wave function/result, convert it to a probabilty. (ie. psi -> |psi|^2)
waveToProb :: Length -> Double -> Density
waveToProb s = (|/| (s |*| s |*| s)) . (% Number) . (^2)
-- ===ENDREGION

-- ===REGION One State Solution Matrix Solving
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

{- Given a grid saying the potential energy in the system, and any particle specific properties (like mass) this will
solve for possible eigen energies the particle can hold and also for probability densities of the particle at that given
energy (the eigen vectors) -}
solveAllOrbitals :: PhysicalGrid Energy -> [WaveSolution]
solveAllOrbitals (PhysicalGrid s v) = sort (map WaveSolution (zip3 energies densityGrids (repeat Nothing)))
  where
    (energies, waveFuncs) = solveMatrix ElectronVolt (hamiltonianMatrix (PhysicalGrid s v))
    waveFuncGrids :: [Grid3D Double]
    waveFuncGrids = map (listArray (bounds v)) waveFuncs
    densitiesGrid3D :: [Grid3D Density]
    densitiesGrid3D = map (fmap $ waveToProb s) waveFuncGrids
    densityGrids :: [PhysicalGrid Density]
    densityGrids = map (PhysicalGrid s) densitiesGrid3D
-- ===ENDREGION

-- ===REGION One State Solution Vector Solving
-- | Computes the hamiltonian given a potential energy grid. This is done using vector storage.
hamiltonianVector :: forall dim lcsu unit. (Energy ~ Qu dim lcsu Double, ValidDLU dim lcsu unit)
  => unit -> PhysicalGrid Energy -> (Vec -> Vec, unit)
hamiltonianVector unit (PhysicalGrid s v) = (matvec, unit)
  where
    matvec :: Vec -> Vec
    matvec psi = listArray (bounds v)
      [ (6 |*| k |*| (psi ! p % Number) |+| (v ! p |*| (psi ! p % Number)) |-| (k |*| (neighborSum psi p % Number))) # unit
      | p <- range $ bounds v ]
    kineticCoeff = (0.0381 % ElectronVolt) |*| qSq (1 % Nanometer)
    k :: Energy
    k = kineticCoeff |/| qSq s
    neighborSum :: Vec -> Idx3 -> Double
    neighborSum psi p = sum [ psi ! q | q <- neighborsOf p (snd $ bounds v) ]

-- | Gets the lowest orbital solution in a 1 particle system using lanczos quick matrix calculations
solveLowestOribtal :: PhysicalGrid Energy -> WaveSolution
solveLowestOribtal grid@(PhysicalGrid s v) = WaveSolution (energy, densityGrid, Nothing)
  where
    (energy, waveFunc) = lanczosLowest (hamiltonianVector ElectronVolt grid) (snd (bounds v))
    waveFuncGrid  :: Grid3D Double
    waveFuncGrid = listArray (bounds v) (elems waveFunc)
    densityGrid3D :: Grid3D Density
    densityGrid3D = fmap (waveToProb s) waveFuncGrid
    densityGrid :: PhysicalGrid Density
    densityGrid = PhysicalGrid s densityGrid3D
-- ===ENDREGION

-- ===REGION Orthogonality Forcing
-- | Generates every combination of spin up/spin down unique assignments for a given number of electrons
candidateSpinAssignments :: Int -> [[Spin]]
candidateSpinAssignments = binaryPartitions SpinUp SpinDown

-- | The bare one-electron energy <psi|T+V_ext|psi> — kinetic + external potential only
oneElectronEnergy :: PhysicalGrid Energy -> Vec -> Energy
oneElectronEnergy extPotential psi = (dotA psi (h0 psi) / dotA psi psi) % ElectronVolt
  where
    h0 :: Vec -> Vec
    (h0, _) = hamiltonianVector ElectronVolt extPotential

-- | Correct total Hartree-Fock energy: E = 1/2 * sum(h_i + epsilon_i)
totalEnergyHF :: PhysicalGrid Energy -> [WaveState] -> Energy
totalEnergyHF extPotential states =
  (0.5 % Number) |*| foldl1' (|+|)
    [ h_i |+| epsilon_i | (epsilon_i, psi, _) <- states, let h_i = oneElectronEnergy extPotential psi ]
-- ===ENDREGION

-- ===REGION Hartree Fock
{-| Hartree-Fock exchange contribution. Given the orbitals of electrons already placed into the
system, this is what gets subtracted from the plain Hamiltonian's action on a trial wavefunction. -}
exchangeTerm :: Length -> [Vec] -> Vec -> Vec
exchangeTerm _ [] psi = zeroVec psi
exchangeTerm s previousOrbitals psi = foldl1' addA (map contributionFrom previousOrbitals)
  where
    pts :: [Idx3]
    pts = range (bounds psi)
    exchangeConst :: EnergyDist
    exchangeConst = electronCharge |*| electronCharge |/| (4 |*| pi |*| vacuumPermittivity)
    contributionFrom :: Vec -> Vec
    contributionFrom phi = listArray (bounds psi)
      [ (exchangeConst |*| overlapPotentialAt p |*| (phi ! p % Number)) # ElectronVolt | p <- pts ]
      where
        overlapPotentialAt :: Idx3 -> InverseLength
        overlapPotentialAt p = foldr1 (|+|)
          [ ((phi ! q * psi ! q) % Number) |/| cartDist p q s | q <- pts, q /= p ]

{- | Builds a matrix-vector operator representing the Fock operator acting on a test wavefunction, by taking the
plain single-particle Hamiltonian (kinetic + external potential) and subtracting the Hartree-Fock exchange contribution. -}
hamiltonianVectorHF :: forall dim lcsu unit. (Energy ~ Qu dim lcsu Double, ValidDLU dim lcsu unit)
  => unit -> PhysicalGrid Energy -> [Vec] -> (Vec -> Vec, unit)
hamiltonianVectorHF unit grid previousOrbitals =
  let
    localMatvec :: Vec -> Vec
    (localMatvec, _) = hamiltonianVector unit grid
    fockMatvec :: Vec -> Vec
    fockMatvec psi = localMatvec psi `subA` exchangeTerm (spacing grid) previousOrbitals psi
  in (fockMatvec, unit)

-- | An pass on adjusting and finding new states when solving for a self consitent field of multiple electrons
scfPass :: [WaveState] -> PhysicalGrid Energy -> [WaveState]
scfPass states extPotential = map resolve [0 .. length states - 1]
  where
    othersDensity :: Int -> PhysicalGrid Density
    othersDensity i = case [ mapVecToGrid extPotential psi waveToProb | (j, (_, psi, _)) <- zip [0 ..] states, j /= i ] of
      [] -> toZeroDensity extPotential
      densities -> foldl1' addGrids densities
    sameSpinOthers :: Int -> Spin -> [Vec]
    sameSpinOthers i mySpin =
      [ psi | (j, (_, psi, sp)) <- zip [0 ..] states, j /= i, sp == mySpin ]
    resolve :: Int -> WaveState
    resolve i =
      let
        mySpin :: Spin
        (_, _, mySpin) = states !! i
        repulsion :: PhysicalGrid Energy
        repulsion = coulombPotentialFromDensity electronCharge electronCharge (othersDensity i)
        effectivePotential :: PhysicalGrid Energy
        effectivePotential = addGrids extPotential repulsion
        occupied :: [Vec]
        occupied = sameSpinOthers i mySpin
        energy :: Energy
        waveFunc :: Vec
        (energy, waveFunc) = lanczosLowestExcluding occupied
          (hamiltonianVectorHF ElectronVolt effectivePotential occupied) (bndGrid extPotential)
      in (energy, waveFunc, mySpin)

-- | An initial first greedy pass when doing a self-consitent field setup.
greedyPass :: [Spin] -> PhysicalGrid Density -> [(Vec, Spin)] -> PhysicalGrid Energy -> [WaveState]
greedyPass [] _ _ _ = []
greedyPass (mySpin : restSpins) cumulativeDensity previousElectrons extPotential =
  let 
    repulsion :: PhysicalGrid Energy
    repulsion = coulombPotentialFromDensity electronCharge electronCharge cumulativeDensity
    effectivePotential :: PhysicalGrid Energy
    effectivePotential = addGrids extPotential repulsion
    sameSpinPrevious :: [Vec]
    sameSpinPrevious = [ psi | (psi, sp) <- previousElectrons, sp == mySpin ]
    energy :: Energy
    waveFunc :: Vec
    (energy, waveFunc)  = lanczosLowestExcluding sameSpinPrevious
      (hamiltonianVectorHF ElectronVolt effectivePotential sameSpinPrevious) (bndGrid extPotential)
    newCumulativeDensity :: PhysicalGrid Density
    newCumulativeDensity = addGrids cumulativeDensity (mapVecToGrid extPotential waveFunc waveToProb)
  in (energy, waveFunc, mySpin) : greedyPass restSpins newCumulativeDensity (previousElectrons ++ [(waveFunc, mySpin)]) extPotential

-- | Does a self consistent field calculation of repeated electron orbital position checking and potential energy generating on repeat
runSCF :: Int -> PhysicalGrid Energy -> [Spin] -> [WaveState]
runSCF maxIterations extPotential spins = iterateSCF 0 initialStates
  where
    initialStates :: [WaveState]
    initialStates = greedyPass spins (toZeroDensity extPotential) [] extPotential
    converged :: [WaveState] -> [WaveState] -> Bool
    converged old new = all (< 1.0e-6) (zipWith energyDiff old new)
      where
        energyDiff :: WaveState -> WaveState -> Double
        energyDiff (e1, _, _) (e2, _, _) = abs ((e1 |-| e2) # ElectronVolt)
    iterateSCF :: Int -> [WaveState] -> [WaveState]
    iterateSCF itersDone states
      | itersDone >= maxIterations = states
      | converged states next      = next
      | otherwise                  = iterateSCF (itersDone + 1) next
      where
        next :: [WaveState]
        next = scfPass states extPotential

-- | Finds hartree fock solutions for total energy and electron orbital solutions given a number of electrons, a max iteration, and a starting external potential
solveInteractingElectronsHF :: Int -> Int -> PhysicalGrid Energy -> (Energy, [WaveSolution])
solveInteractingElectronsHF numElectrons maxIterations extPotential =
  (totalEnergyHF extPotential bestStates, toSpinOrbitals bestStates)
  where
    candidates :: [[WaveState]]
    candidates = map (runSCF maxIterations extPotential) (candidateSpinAssignments numElectrons)
    bestStates :: [WaveState]
    bestStates = minimumBy (comparing scoreOf) candidates
    scoreOf :: [WaveState] -> Double
    scoreOf states = sum [ e # ElectronVolt | (e, _, _) <- states ]
    toSpinOrbitals :: [WaveState] -> [WaveSolution]
    toSpinOrbitals states = [ WaveSolution (energy, mapVecToGrid extPotential psi waveToProb, Just sp) | (energy, psi, sp) <- states ]
-- ===ENDREGION
