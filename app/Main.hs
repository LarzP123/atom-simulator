import AtomSimulator.Units
import AtomSimulator.QuantumSolver
import Data.Metrology

{-| Currently this just returns output for a solution for a quantum well as a demonstration-}
main :: IO ()
main = do
  let n     = 9
      s     = 0.1 % Nanometer
      zeroFunc _ _ _ = 0 % ElectronVolt
      well3D = createPotentialGrid zeroFunc n s
  putStrLn "\n3D infinite cubic well (10x10x10 grid, 0.1 nm spacing), lowest 4 levels:"
  mapM_ (putStrLn . ("  E = " ++) . show . fst) (take 4 (solve3D well3D))