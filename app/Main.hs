{-# LANGUAGE FlexibleContexts #-}
import AtomSimulator.Units
import AtomSimulator.QuantumSolver
import Data.Metrology.Poly
import Data.List (sortBy)
import Data.Ord (comparing)
import Text.Printf

{-| Currently this just returns output for a solution for a quantum well as a demonstration-}
main :: IO ()
main = do
  let 
    k :: SpringConstant
    k = (0.25 % Hartree) |/| (1 % Bohr) |/| (1 % Bohr)
    grid = createPotentialGrid (harmonicPotential k) 12 (0.03 % Nanometer)

  putStrLn "Single-particle harmonic spectrum (no e-e interaction)"
  putStrLn "Expect energies 1x0.75, 3x1.25, 6x1.75"
  mapM_ (putStrLn . showWaveSolutionInUnit Hartree) (take 10 (solveAllOrbitals grid))

  putStrLn "\nHooke's atom, 2 interacting electrons (Hartree-Fock)"
  let (systemEnergy, solutions) = solveInteractingElectronsHF 2 5 grid
  mapM_ (putStrLn . showWaveSolutionInUnit Hartree) solutions
  printf "Total: %.4f Hartree  (must be >= 2 Hartree -- variational bound)\n" (systemEnergy # Hartree)