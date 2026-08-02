{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE TypeFamilies      #-}
{-# LANGUAGE TypeOperators #-}
{-# OPTIONS_GHC -Wmissing-signatures -Wmissing-local-signatures #-}
module AtomSimulator.Units where
import Data.Metrology
import Data.Metrology.TH
import Data.Metrology.Show ()

-- ===REGION Base Units
declareMonoUnit "Nanometer" (Just "nm")
-- | Length is a physical dimension that quantifies the space between points
type Length = MkQu_D Nanometer
declareDerivedUnit "MicroMeter" [t| Nanometer |] 1e3 (Just "um")
declareDerivedUnit "Meter" [t| Nanometer |] 1e9 (Just "m")
declareDerivedUnit "Bohr" [t| Nanometer |] 0.0529177249 (Just "m")

declareMonoUnit "ElectronVolt" (Just "eV")
{-| Energy is a conserved property in physical systems. Particles can have different potential energies at
different points in space and the energy that a particle in a system have can be quantized.-}
type Energy = MkQu_D ElectronVolt
declareDerivedUnit "Joule" [t| ElectronVolt |] 6.241509074e18 (Just "J")
declareDerivedUnit "MegaElectronVolt" [t| ElectronVolt |] 1e6 (Just "MeV")
declareDerivedUnit "Hartree" [t| ElectronVolt |] 27.2114 (Just "E h")

declareMonoUnit "AtomicUnit" (Just "Au")
-- | Mass is a property of particles that decreases the wavelength of particles
type Mass = MkQu_D AtomicUnit
declareDerivedUnit "Kilogram" [t| AtomicUnit |] 6.02214076e26 (Just "kg")

declareMonoUnit "ElementaryCharge" (Just "e")
-- | Charge is a property of particles that attracts particles with opposite sign charge and repels ones with same sign charge
type Charge = MkQu_D ElementaryCharge

-- | Spin is an intrinsic form of angular momentum carried by particles. There are 2 orthogonal states of spin for spin 1/2 spin particles
data Spin = SpinUp | SpinDown deriving (Eq, Show)
-- ===ENDREGION

-- ===REGION Derived Unit Types
-- | Density is the amount of something in a given volume of 3D space.
type Density  = Length %^ MThree
{-| Speed Squared is an abstract concept that represents the amount of distance covered by an object in a certain
amount of time... squared -}
type SpeedSq = MkQu_D (ElectronVolt :/ AtomicUnit)
-- | Measure of electric polarizability of a dielectric material
type Permittivity = MkQu_D (ElementaryCharge :* ElementaryCharge :/ ElectronVolt :/ Nanometer)
-- | Energy times distance
type EnergyDist = MkQu_D (ElectronVolt :* Nanometer)
-- | Length to negative one power
type InverseLength = Length %^ MOne
-- | This is a constant for the slope of a force that scales linearly with distance. Like the force on a spring
type SpringConstant = MkQu_D (ElectronVolt :/ Nanometer :/ Nanometer)
-- | A distance times a distances, or an area
type DistSq = MkQu_D (Nanometer :* Nanometer)
-- ===ENDREGION

-- ===REGION Physical Constants
-- | Light travels a constant amount of distance in a given amount of time
lightSpeedSq :: SpeedSq
lightSpeedSq = (2.99792458*10^8)^2 |*| (1 % Joule) |/| (1 % Kilogram)

-- | This is the mass of an electron particle at rest. 
electronMass :: Mass
electronMass = (0.510998 % MegaElectronVolt) |/| lightSpeedSq
-- | This is the mass of a proton particle at rest. 
protonMass :: Mass
protonMass = (938.272 % MegaElectronVolt) |/| lightSpeedSq
-- | Permittivity of free space
vacuumPermittivity :: Permittivity
vacuumPermittivity = (55.26349406 % ElementaryCharge) |*| (1 % ElementaryCharge) |/| (1 % ElectronVolt) |/| (1 % MicroMeter)
-- ===ENDREGION

-- | The charge of an electron
electronCharge :: Charge
electronCharge = (-1) % ElementaryCharge
-- | The charge of a proton
protonCharge :: Charge
protonCharge = 1 % ElementaryCharge
