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
import Data.Function (on)
import Text.Printf (printf)

{- Creates a 3D potential energy grid from a function that says a given potential energy at a point in space -}
createPotentialGrid :: (Length -> Length -> Length -> Energy) -> Natural -> Length -> PhysicalGrid Energy
createPotentialGrid = createGridFromFuncShifted

-- | Finds the harmonic potential at a position given a spring constant
harmonicPotential :: SpringConstant -> Length -> Length -> Length -> Energy
harmonicPotential = harmonicFunc

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

-- | Spin is an intrinsic form of angular momentum carried by particles. There are 2 orthogonal states of spin for spin 1/2 spin particles
data Spin = SpinUp | SpinDown deriving (Eq, Show)
-- | The solution for an electron cloud. This contains the eigen-energy of the particle, a probability density, and the spin of the particle if from a multi-particle solver
newtype WaveSolution = WaveSolution (Energy, PhysicalGrid Density, Maybe Spin) deriving (Eq)
{- | The wave function of an electron (psi). This contains the eigen-energy of the particle, a vector containning indexes and results of the
wavefunction at various points in space, and the particles spin -}
type WaveState = (Energy, Vec, Spin)

instance Show WaveSolution where
  show (WaveSolution (e, _, spin)) = case spin of
    Just s  -> show s ++ " " ++ show e
    Nothing -> show e

showInHartree :: WaveSolution -> String
showInHartree (WaveSolution (e, _, spin)) = case spin of
  Just s  -> show s ++ " " ++ printf "%.4f Hartree" (e # Hartree)
  Nothing -> printf "%.4f Hartree" (e # Hartree)

instance Ord WaveSolution where
  compare = comparing (\(WaveSolution (e, _, _)) -> e)

waveToProb :: Length -> Double -> Density
waveToProb s = (|/| (s |*| s |*| s)) . (% Number) . (^2)

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
    (localMatvec, _) = hamiltonianVector unit grid
    fockMatvec psi = localMatvec psi `subA` exchangeTerm (spacing grid) previousOrbitals psi
  in (fockMatvec, unit)

-- | Generates every combination of spin up/spin down unique assignments for a given number of electrons
candidateSpinAssignments :: Natural -> [[Spin]]
candidateSpinAssignments = binaryPartitions SpinUp SpinDown

runSCF :: Natural -> PhysicalGrid Energy -> [Spin] -> [WaveState]
runSCF maxIterations extPotential@(PhysicalGrid s v) spins = iterateSCF 0 initialStates
  where
    bnd = snd (bounds v)
    zeroDensity = PhysicalGrid s (fmap (const (waveToProb s 0)) v)
    densityOf psi = PhysicalGrid s (fmap (waveToProb s) (listArray (bounds v) (elems psi)))

    initialStates :: [WaveState]
    initialStates = greedyPass spins zeroDensity []

    greedyPass [] _ _ = []
    greedyPass (mySpin : restSpins) cumulativeDensity previousElectrons =
      let repulsion          = coulombPotentialFromDensity electronCharge electronCharge cumulativeDensity
          effectivePotential  = addGrids extPotential repulsion
          sameSpinPrevious    = [ psi | (psi, sp) <- previousElectrons, sp == mySpin ]
          (energy, waveFunc)  = lanczosLowestExcluding sameSpinPrevious
                                  (hamiltonianVectorHF ElectronVolt effectivePotential sameSpinPrevious) bnd
          newCumulativeDensity = addGrids cumulativeDensity (densityOf waveFunc)
      in (energy, waveFunc, mySpin)
          : greedyPass restSpins newCumulativeDensity (previousElectrons ++ [(waveFunc, mySpin)])

    scfPass states = map resolve [0 .. length states - 1]
      where
        othersDensity i = case [ densityOf psi | (j, (_, psi, _)) <- zip [0 ..] states, j /= i ] of
          []        -> zeroDensity
          densities -> foldl1' addGrids densities

        sameSpinOthers i mySpin =
          [ psi | (j, (_, psi, sp)) <- zip [0 ..] states, j /= i, sp == mySpin ]

        resolve i =
          let (_, _, mySpin)    = states !! i
              repulsion          = coulombPotentialFromDensity electronCharge electronCharge (othersDensity i)
              effectivePotential = addGrids extPotential repulsion
              occupied           = sameSpinOthers i mySpin
              (energy, waveFunc) = lanczosLowestExcluding occupied
                                      (hamiltonianVectorHF ElectronVolt effectivePotential occupied) bnd
          in (energy, waveFunc, mySpin)

    converged old new = all (< 1.0e-6) (zipWith energyDiff old new)
      where energyDiff (e1, _, _) (e2, _, _) = abs ((e1 |-| e2) # ElectronVolt)

    iterateSCF itersDone states
      | itersDone >= maxIterations = states
      | converged states next      = next
      | otherwise                  = iterateSCF (itersDone + 1) next
      where next = scfPass states

solveInteractingElectronsHF :: Natural -> Natural -> PhysicalGrid Energy -> (Energy, [WaveSolution])
solveInteractingElectronsHF numElectrons maxIterations extPotential =
  (totalEnergyHF extPotential bestStates, toSpinOrbitals bestStates)
  where
    (PhysicalGrid s v) = extPotential
    densityOf psi = PhysicalGrid s (fmap (waveToProb s) (listArray (bounds v) (elems psi)))
    candidates :: [[WaveState]]
    candidates = map (runSCF maxIterations extPotential) (candidateSpinAssignments numElectrons)
    bestStates :: [WaveState]
    bestStates = minimumBy (comparing scoreOf) candidates
    scoreOf :: [WaveState] -> Double
    scoreOf states = sum [ e # ElectronVolt | (e, _, _) <- states ]
    toSpinOrbitals :: [WaveState] -> [WaveSolution]
    toSpinOrbitals states = [ WaveSolution (energy, densityOf psi, Just sp) | (energy, psi, sp) <- states ]

-- | The bare one-electron energy <psi|T+V_ext|psi> — kinetic + external potential only
oneElectronEnergy :: PhysicalGrid Energy -> Vec -> Energy
oneElectronEnergy extPotential psi = (dotA psi (h0 psi) / dotA psi psi) % ElectronVolt
  where
    (h0, _) = hamiltonianVector ElectronVolt extPotential

-- | Correct total Hartree-Fock energy: E = 1/2 * sum(h_i + epsilon_i)
totalEnergyHF :: PhysicalGrid Energy -> [WaveState] -> Energy
totalEnergyHF extPotential states =
  (0.5 % Number) |*| foldl1' (|+|)
    [ h_i |+| epsilon_i | (epsilon_i, psi, _) <- states, let h_i = oneElectronEnergy extPotential psi ]
