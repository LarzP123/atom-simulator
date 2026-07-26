{-# LANGUAGE FlexibleContexts #-}
import AtomSimulator.Units
import AtomSimulator.QuantumSolver
import Data.Metrology.Poly
import Data.List (sortBy)
import Data.Ord (comparing)
import GHC.Natural (Natural)
import Text.Printf

{-| Currently this just returns output for a solution for a quantum well as a demonstration-}
main :: IO ()
main = do
  let n     = 9
      s     = 0.1 % Nanometer
      zeroFunc _ _ _ = 0 % ElectronVolt
      grid    = createPotentialGrid zeroFunc n s

  putStrLn "Dense solver (solve3D)"
  let 
    denseEnergies :: [Energy]
    denseEnergies = take 5 (map fst (solveAllOrbitals grid))
  mapM_ (\e -> printf "%.5f eV\n" (e # ElectronVolt)) denseEnergies

  putStrLn "\nSparse Lanczos solver"
  let 
    sparseEnergies :: Energy
    sparseEnergies = fst (solveLowestOribtal grid)
  (\e -> printf "%.5f eV\n" (e # ElectronVolt)) sparseEnergies
