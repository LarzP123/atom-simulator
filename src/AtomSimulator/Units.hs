{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE TypeFamilies      #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE TypeOperators #-}

module AtomSimulator.Units where

import Data.Metrology
import Data.Metrology.TH
import Data.Metrology.Show ()

declareMonoUnit "Nanometer" (Just "nm")

declareDerivedUnit "Meter" [t| Nanometer |] 1e9 (Just "m")

declareMonoUnit "ElectronVolt" (Just "eV")
declareDerivedUnit "Joule" [t| ElectronVolt |] 6.241509074e18 (Just "J")
declareDerivedUnit "MegaElectronVolt" [t| ElectronVolt |] 1e6 (Just "MeV")

declareMonoUnit "AtomicUnit" (Just "Au")
declareDerivedUnit "Kilogram" [t| AtomicUnit |] 6.02214076e26 (Just "kg")

{-| Energy is a conserved property in physical systems. Particles can have different potential energies at
different points in space and the energy that a particle in a system have can be quantized.-}
type Energy   = MkQu_D ElectronVolt
-- | Length is a physical dimension that quantifies the space between points
type Length   = MkQu_D Nanometer
-- | Mass is a property of particles that decreases the wavelength of a particle 
type Mass     = MkQu_D AtomicUnit

-- | Density is the amount of something in a given volume of 3D space.
type Density  = Length %^ MThree
{-| Speed Squared is an abstract concept that represents the amount of distance covered by an object in a certain
amount of time... squared -}
type SpeedSq = MkQu_D (ElectronVolt :/ AtomicUnit)

-- | A meter per second is a standard SI unit for velocity. This value is that number squared
meterPerSecondSq :: SpeedSq
meterPerSecondSq = (1 % Joule) |/| (1 % Kilogram)

-- | Light travels a constant amount of distance in a given amount of time
lightSpeedSq :: SpeedSq
lightSpeedSq = (2.99792458*10^8)^2 |*| meterPerSecondSq

-- | This is the mass of an electron particle at rest. 
electronMass :: Mass
electronMass = (0.510998 % MegaElectronVolt) |/| lightSpeedSq
-- | This is the mass of a proton particle at rest. 
protonMass :: Mass
protonMass = (938.272 % MegaElectronVolt) |/| lightSpeedSq
