{-# LANGUAGE FlexibleContexts #-}
import Test.Hspec
import Data.Metrology.Poly
import Data.Array.IArray (bounds, elems, (!), range)
import AtomSimulator.QuantumSolver
import AtomSimulator.Units
import Data.List (sort)
import Data.Maybe (listToMaybe, isJust)

{-| Gets the relative error between 2 values. Allows us to test numerical solutions mostly line
up with analytical solutions -}
relError :: Double -> Double -> Double
relError ref val = abs (val - ref) / abs ref

-- | Test input
main :: IO ()
main = hspec $ do
  describe "Potential Grid Creation" $ do
    let
      k :: SpringConstant
      k = (0.5 % Hartree) |/| ((1 % Bohr) |*| (1 % Bohr))
      gridSpacingIn :: Length
      gridSpacingIn = 0.03 % Nanometer
      mid :: Int
      mid = 5
      gridSpacingOut :: Length
      grid3D :: Grid3D Energy
      PhysicalGrid gridSpacingOut grid3D = createPotentialGrid (harmonicPotential k) (2*mid+1) gridSpacingIn
    it "has correct spacing" $ do
      gridSpacingOut `shouldBe` gridSpacingIn
    it "has zero potential at center" $ do
      (grid3D ! (mid,mid,mid)) # Hartree `shouldSatisfy` (< 0.01)
      (grid3D ! (mid,mid,mid)) # Hartree `shouldSatisfy` (< 0.01)
    it "has potential increases with X distance" $ do
      (grid3D ! (mid,mid,mid) # Hartree) `shouldSatisfy` (<= (grid3D ! (0,mid,mid)) # Hartree)
      (grid3D ! (mid,mid,mid) # Hartree) `shouldSatisfy` (<= (grid3D ! (2*mid,mid,mid)) # Hartree)
    it "has potential increases with Y distance" $ do
      (grid3D ! (mid,mid,mid) # Hartree) `shouldSatisfy` (<= (grid3D ! (mid,0,mid)) # Hartree)
      (grid3D ! (mid,mid,mid) # Hartree) `shouldSatisfy` (<= (grid3D ! (mid,2*mid,mid)) # Hartree)
    it "has potential increases with Z distance" $ do
      (grid3D ! (mid,mid,mid) # Hartree) `shouldSatisfy` (<= (grid3D ! (mid,mid,0)) # Hartree)
      (grid3D ! (mid,mid,mid) # Hartree) `shouldSatisfy` (<= (grid3D ! (mid,mid,2*mid)) # Hartree)
  describe "Harmonic Potential Function" $ do
    let
      k :: SpringConstant
      k = (1.0 % Hartree) |/| ((1 % Bohr) |*| (1 % Bohr))
    it "has zero potential at the origin" $ do
      let 
        v :: Energy
        v = harmonicPotential k (0 % Bohr) (0 % Bohr) (0 % Bohr)
      v `shouldSatisfy` (< 0.001 % Hartree)
    it "has symmetric potential" $ do
      let
        v_pos, v_neg :: Energy
        v_pos = harmonicPotential k (1 % Bohr) (0 % Bohr) (0 % Bohr)
        v_neg = harmonicPotential k ((-1) % Bohr) (0 % Bohr) (0 % Bohr)
      abs ((v_pos # Hartree) - (v_neg # Hartree)) `shouldSatisfy` (< 0.001)
    it "has potential scale with the spring constant" $ do
      let 
        k1, k2 :: SpringConstant
        k1 = (1.0 % Hartree) |/| ((1 % Bohr) |*| (1 % Bohr))
        k2 = 2 |*| k1
        v1, v2 :: Energy
        v1 = harmonicPotential k1 (1 % Bohr) (0 % Bohr) (0 % Bohr)
        v2 = harmonicPotential k2 (1 % Bohr) (0 % Bohr) (0 % Bohr)
      relError (v2 |/| v1 # Number) 2  `shouldSatisfy` (< 0.05)
  describe "Wave Function to Probability Conversion" $ do
    let
      spacing :: Length
      spacing = 0.05 % Nanometer
    it "has zero wavefunction gives zero probability" $ do
      (waveToProb spacing 0.0 # InverseCubicNanometer) `shouldBe` 0.0
    it "has always non-negative probability" $ do
      let
        probs :: [Double]
        probs = [waveToProb spacing x # InverseCubicNanometer | x <- [-2,-1,0,1,2]]
      all (>= 0.0) probs `shouldBe` True
    it "is symetric" $ do
      let
        p_pos, p_neg :: Double
        p_pos = waveToProb spacing 2.5 # InverseCubicNanometer
        p_neg = waveToProb spacing (-2.5) # InverseCubicNanometer
      p_pos `shouldBe` p_neg
    it "is homogenous of degree 2" $ do
      let
        p1, p2, p3, p4 :: Double
        p1 = waveToProb spacing 1.0 # InverseCubicNanometer
        p2 = waveToProb spacing 2.0 # InverseCubicNanometer
        p3 = waveToProb spacing 3.0 # InverseCubicNanometer
        p4 = waveToProb spacing 4.0 # InverseCubicNanometer
      abs (p2 - 4.0 * p1) `shouldSatisfy` (< 1e-10)
      abs (p3 - 9.0 * p1) `shouldSatisfy` (< 1e-10)
      abs (p4 - 16.0 * p1) `shouldSatisfy` (< 1e-10)
    it "density scales as 1/spacing^3" $ do
      let 
        s1, s2 :: Length
        s1 = 0.1 % Nanometer
        s2 = 0.05 % Nanometer
        prob1, prob2, expected_ratio :: Double
        prob1 = waveToProb s1 1.0 # InverseCubicNanometer
        prob2 = waveToProb s2 1.0 # InverseCubicNanometer
        expected_ratio = (0.1 / 0.05) ^ (3::Int)
      abs ((prob2 / prob1) - expected_ratio) `shouldSatisfy` (< 0.01)
  describe "WaveSolution Definition" $ do
    let
      k :: InverseLengthSq
      k = 1 |/| ((4 % Nanometer) |*| (4 % Nanometer))
      amp :: Energy
      amp = 20 % ElectronVolt
      grid :: PhysicalGrid Energy
      grid = createPotentialGrid (sinePotential amp k) 12 (0.03 % Nanometer)
      sol :: WaveSolution
      sol = solveLowestOribtal grid
    it "can be shown" $ do
      let
        output :: String
        output = show sol
      length output `shouldSatisfy` (> 0)
    it "WaveSolutions can be compared by spin" $ do
      let
        withSpinUp :: WaveSolution -> WaveSolution
        withSpinUp (WaveSolution (e, grid, _)) = WaveSolution (e, grid, Just SpinUp)
        flipSpinSol :: WaveSolution -> WaveSolution
        flipSpinSol (WaveSolution (e, grid, spin)) = WaveSolution (e, grid, flipSpin spin)
        solSpinUp, solFlippedTwice :: WaveSolution
        solSpinUp = withSpinUp sol
        solFlippedTwice = (flipSpinSol . flipSpinSol) solSpinUp
      sol `shouldNotBe` solSpinUp
      solSpinUp `shouldBe` solFlippedTwice
    it "WaveSolutions can be ordered" $ do
      let
        k1, k2 :: SpringConstant
        k1 = (0.3 % Hartree) |/| ((1 % Bohr) |*| (1 % Bohr))
        k2 = (0.7 % Hartree) |/| ((1 % Bohr) |*| (1 % Bohr))
        grid1, grid2 :: PhysicalGrid Energy
        grid1 = createPotentialGrid (harmonicPotential k1) 12 (0.03 % Nanometer)
        grid2 = createPotentialGrid (harmonicPotential k2) 12 (0.03 % Nanometer)
        sol1, sol2 :: WaveSolution
        sol1 = solveLowestOribtal grid1
        sol2 = solveLowestOribtal grid2
        sorted :: [WaveSolution]
        sorted = sort [sol2, sol1]
      listToMaybe sorted `shouldBe` Just sol1
    it "spin is Nothing for single particle" $ do
      let
        spin :: Maybe Spin
        WaveSolution (_, _, spin) = sol
      spin `shouldBe` Nothing
  describe "Single Particle Infinite Spherical Well Ground Energy Correct" $ do
    let 
      numericalSol :: Length -> Energy
      numericalSol radius = 
          let numericalEnergy :: Energy
              WaveSolution (numericalEnergy, _, _) = solveLowestOribtal (createPotentialGrid (sphericalWellPotential (100000 % MegaElectronVolt) radius) 40 (radius |/| 20))
          in numericalEnergy
      analyticalSol :: Length -> Energy
      analyticalSol radius = hbar2 |*| pi^2 |/| (2 |*| electronMass |*| radius |*| radius)
    it "solves rad=0.1nm" $ do
      relError (numericalSol (0.1 % Nanometer) # ElectronVolt)
               (analyticalSol (0.1 % Nanometer) # ElectronVolt) `shouldSatisfy` (< 0.1)
    it "solves rad=1nm" $ do
      relError (numericalSol (1 % Nanometer) # ElectronVolt)
               (analyticalSol (1 % Nanometer) # ElectronVolt) `shouldSatisfy` (< 0.1)
    it "solves rad=2nm" $ do
      relError (numericalSol (2 % Nanometer) # ElectronVolt)
               (analyticalSol (2 % Nanometer) # ElectronVolt) `shouldSatisfy` (< 0.1)
    it "solves rad=5nm" $ do
      relError (numericalSol (5 % Nanometer) # ElectronVolt)
               (analyticalSol (5 % Nanometer) # ElectronVolt) `shouldSatisfy` (< 0.1)
  describe "Hartree-Fock Solver for Multiple Electrons" $ do
    let 
      k :: SpringConstant
      k = (0.25 % Hartree) |/| ((1 % Bohr) |*| (1 % Bohr))
      grid :: PhysicalGrid Energy
      grid = createPotentialGrid (harmonicPotential k) 10 (0.04 % Nanometer)
      totalE :: Energy
      solutions :: [WaveSolution]
      (totalE, solutions) = solveInteractingElectronsHF 2 5 grid
    it "returns info on 2 electrons" $ do
      length solutions `shouldBe` 2
    it "gets the correct total energy" $ do
      relError (totalE # Hartree) 2 `shouldSatisfy` (< 0.1)
    it "follows Pauli exclusion principle" $ do
      let
        spin1, spin2 :: Maybe Spin
        WaveSolution (_, _, spin1) : WaveSolution (_, _, spin2) : _ = solutions
      flipSpin spin1 `shouldBe` spin2
    it "energies should be similar" $ do
      let
        energy1, energy2 :: Energy
        WaveSolution (energy1, _, _) : WaveSolution (energy2, _, _) : _ = solutions
      relError (energy1 # ElectronVolt) (energy2 # ElectronVolt) `shouldSatisfy` (< 0.1)
    it "particle energy should be correct" $ do
      let
        particleEnergy :: Energy
        WaveSolution (particleEnergy, _, _) : _ = solutions
      relError (particleEnergy # Hartree) 1.25 `shouldSatisfy` (< 0.1)
  describe "Correct EigenEnergies Infinite Quantum Dot" $ do
    let
      numericalSol :: Length -> [Energy]
      numericalSol length =
        -- we intentionally have the spacing and length divider off by 1 because of how well boundaries done
        let
          solutions :: [WaveSolution]
          solutions = solveAllOrbitals (createPotentialGrid zeroPotential 10 (length |/| 11))
        in map (\(WaveSolution (e, _, _)) -> e) solutions
      analyticalSol :: Length -> Int -> Energy
      analyticalSol length nXnYnZSq = h2 |/| (8 |*| electronMass |*| length |*| length) |*| (fromIntegral nXnYnZSq % Number)
    describe "1nm" $ do
      let
        numericalSols :: [Energy]
        numericalSols = numericalSol (1 % Nanometer)
        analyticalSols :: Int -> Energy
        analyticalSols = analyticalSol (1 % Nanometer)
      it "has the correct ground energy" $ do
        let
          numericalSol0 :: Energy
          numericalSol0 : _ = numericalSols
        relError (numericalSol0 # ElectronVolt) (analyticalSols 3 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 1st energy" $ do
        relError (numericalSols !! 1 # ElectronVolt) (analyticalSols 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 2nd energy" $ do
        relError (numericalSols !! 2 # ElectronVolt) (analyticalSols 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 3rd energy" $ do
        relError (numericalSols !! 3 # ElectronVolt) (analyticalSols 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 4th energy" $ do
        relError (numericalSols !! 4 # ElectronVolt) (analyticalSols 9 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 5th energy" $ do
        relError (numericalSols !! 5 # ElectronVolt) (analyticalSols 9 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 6th energy" $ do
        relError (numericalSols !! 6 # ElectronVolt) (analyticalSols 9 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 7th energy" $ do
        relError (numericalSols !! 7 # ElectronVolt) (analyticalSols 11 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 8th energy" $ do
        relError (numericalSols !! 8 # ElectronVolt) (analyticalSols 11 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 9th energy" $ do
        relError (numericalSols !! 9 # ElectronVolt) (analyticalSols 11 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 10th energy" $ do
        relError (numericalSols !! 10 # ElectronVolt) (analyticalSols 12 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has degeneracies at excitations 1-3" $ do
        relError (numericalSols !! 1 # ElectronVolt) (numericalSols !! 3 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has degeneracies at excitations 4-6" $ do
        relError (numericalSols !! 4 # ElectronVolt) (numericalSols !! 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has degeneracies at excitations 7-9" $ do
        relError (numericalSols !! 7 # ElectronVolt) (numericalSols !! 9 # ElectronVolt) `shouldSatisfy` (<0.1)
    describe "10nm" $ do
      let
        numericalSols :: [Energy]
        numericalSols = numericalSol (10 % Nanometer)
        analyticalSols :: Int -> Energy
        analyticalSols = analyticalSol (10 % Nanometer)
      it "has the correct ground energy" $ do
        let
          numericalSol0 :: Energy
          numericalSol0 : _ = numericalSols
        relError (numericalSol0 # ElectronVolt) (analyticalSols 3 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 1st energy" $ do
        relError (numericalSols !! 1 # ElectronVolt) (analyticalSols 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 2nd energy" $ do
        relError (numericalSols !! 2 # ElectronVolt) (analyticalSols 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 3rd energy" $ do
        relError (numericalSols !! 3 # ElectronVolt) (analyticalSols 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 4th energy" $ do
        relError (numericalSols !! 4 # ElectronVolt) (analyticalSols 9 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 5th energy" $ do
        relError (numericalSols !! 5 # ElectronVolt) (analyticalSols 9 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 6th energy" $ do
        relError (numericalSols !! 6 # ElectronVolt) (analyticalSols 9 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 7th energy" $ do
        relError (numericalSols !! 7 # ElectronVolt) (analyticalSols 11 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 8th energy" $ do
        relError (numericalSols !! 8 # ElectronVolt) (analyticalSols 11 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 9th energy" $ do
        relError (numericalSols !! 9 # ElectronVolt) (analyticalSols 11 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has the correct 10th energy" $ do
        relError (numericalSols !! 10 # ElectronVolt) (analyticalSols 12 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has degeneracies at excitations 1-3" $ do
        relError (numericalSols !! 1 # ElectronVolt) (numericalSols !! 3 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has degeneracies at excitations 4-6" $ do
        relError (numericalSols !! 4 # ElectronVolt) (numericalSols !! 6 # ElectronVolt) `shouldSatisfy` (<0.1)
      it "has degeneracies at excitations 7-9" $ do
        relError (numericalSols !! 7 # ElectronVolt) (numericalSols !! 9 # ElectronVolt) `shouldSatisfy` (<0.1)