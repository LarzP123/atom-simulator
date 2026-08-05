module AtomSimulator.OutputFormatter where
import AtomSimulator.MathHelper
import AtomSimulator.QuantumSolver
import AtomSimulator.Units
import Data.Array
import Data.List
import Data.Metrology

-- | Converts a list of strings to an array that is parsable valid json
jsonArr :: [String] -> String
jsonArr xs = "[" ++ intercalate "," xs ++ "]\n"

-- | Converts a PhysicalGrid to a readable json array of those numbers. You supply a unit to convert the grid out of that unit before writing as json
gridToJSON :: forall dim unit. (ValidDLU dim DefaultLCSU unit, Show unit)
  => unit -> PhysicalGrid (Qu dim DefaultLCSU Double) -> String
gridToJSON unit (PhysicalGrid _ samplesArr) =
  jsonArr
    [ jsonArr
        [ jsonArr [ show (probAt (x, y, z)) | z <- [0 .. izMax] ]
        | y <- [0 .. iyMax]
        ]
    | x <- [0 .. ixMax]
    ]
  where
    (ixMax, iyMax, izMax) = snd $ bounds samplesArr
    probAt p = (samplesArr ! p) # unit

-- | Serializes one solved electron orbital to a JSON object.
electronToJSON :: forall dim lcsu unit. (Density ~ Qu dim lcsu Double, ValidDLU dim lcsu unit, Show unit, Unit unit)
  => unit -> Int -> WaveSolution -> String
electronToJSON densityUnit idx (WaveSolution (energy, densityGrid, spin)) =
  "{\n"
    ++ intercalate
      ",\n"
      [ "\"index\":" ++ show idx
      , "\"energyHartree\":" ++ show (energy # Hartree)
      , "\"spin\":\"" ++ show spin ++ "\""
      , "\"dims\":" ++ jsonArr [show (ixMax + 1), show (iyMax + 1), show (izMax + 1)]
      , "\"spacingNm\":" ++ show (spacing densityGrid # Nanometer)
      , "\"density\":" ++ gridToJSON densityUnit densityGrid
      ]
    ++ "\n}"
  where
    (ixMax, iyMax, izMax) = bndGrid densityGrid