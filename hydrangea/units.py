#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Provides functions for unit transformation."""

from pdb import set_trace

# ----------------------------------
# Define constants (all in SI units)
# ----------------------------------

# Lengths
CENTIMETRE = CM = 1e-2
METRE = M = 1
KILOMETRE = KM = 1e3
PARSEC = PC = 3.08567758e16
KILOPARSEC = KPC = 1e3 * PARSEC
MEGAPARSEC = MPC = 1e6 * PARSEC

# Masses
GRAM = 1e-3
KILOGRAM = KG = 1
#SOLAR_MASS = MSUN = 1.98847e30
SOLAR_MASS = MSUN = 1.989e30

PROTON_MASS = MP = 1.6726219e-27

# Times
SECOND = S = 1
HOUR = 3600.0
DAY = 24 * HOUR
YEAR = 365.24 * DAY
MEGAYEAR = MYR = 1e6 * YEAR
GIGAYEAR = GYR = 1e9 * YEAR

# ENERGIES
ELECTRON_VOLT = EV = 1.60217662e-19
KILO_ELECTRON_VOLT = KEV = 1e3 * EV

# TEMPERATURES
KELVIN = 1

# NATURAL CONSTANTS
NEWTONS_CONSTANT = GNEWTON = 6.67428e-11
SPEED_OF_LIGHT = C = 299792458
GAMMA = 5/3


DIM_MASS = 'Mass'
DIM_LENGTH = 'Length'
DIM_TIME = 'Time'
DIM_VELOCITY = 'Velocity'
DIM_ENERGY = 'Energy'
DIM_SPECIFIC_ENERGY = 'SpecificEnergy'
DIM_MASS_RATE = 'MassRate'
DIM_SPECIFIC_ANGULAR_MOMENTUM = 'SpecificAngularMomentum'
DIM_MOMENT_OF_INERTIA = 'MomentOfInertia'
DIM_TEMPERATURE = 'Temperature'
DIM_DENSITY = 'Density'
DIM_ENTROPY = 'Entropy'


# -------------------
# Define target units
# -------------------

CGS_UNITS = {
    DIM_MASS: GRAM,
    DIM_LENGTH: CENTIMETRE,
    DIM_TIME: SECOND,
    DIM_VELOCITY: CM/S,
    DIM_ENERGY: GRAM * CM**2 / S**2,
    DIM_SPECIFIC_ENERGY: CM**2 / S**2,
    DIM_MASS_RATE: GRAM / S,
    DIM_SPECIFIC_ANGULAR_MOMENTUM: CM**2 / S,
    DIM_MOMENT_OF_INERTIA: GRAM * CM**2,
    DIM_TEMPERATURE: KELVIN,
    DIM_DENSITY: GRAM / CM**3,
    DIM_ENTROPY: GRAM**(1-GAMMA) * CM**(3*GAMMA-1) * S**(-2),
    }

ASTRO_UNITS = {
    DIM_MASS: SOLAR_MASS,
    DIM_LENGTH: MEGAPARSEC,
    DIM_TIME: GIGAYEAR,
    DIM_VELOCITY: KM / S,
    DIM_ENERGY: SOLAR_MASS * (KM / S)**2,
    DIM_SPECIFIC_ENERGY: (KM / S)**2,
    DIM_MASS_RATE: SOLAR_MASS / YEAR,
    DIM_SPECIFIC_ANGULAR_MOMENTUM: MEGAPARSEC * KM / S,
    DIM_MOMENT_OF_INERTIA: SOLAR_MASS * MEGAPARSEC**2,
    DIM_TEMPERATURE: KELVIN,
    DIM_DENSITY: PROTON_MASS / CM**3,
    DIM_ENTROPY: GRAM**(1-GAMMA) * CM**(3*GAMMA-1) * S**(-2),
    }


# --------------------------
# Define data set dimensions
# --------------------------

SUBHALO_DIMENSIONS = {
    'BlackHoleMass': DIM_MASS,
    'BlackHoleMassAccretionRate': DIM_MASS_RATE,
    'CentreOfMass': DIM_LENGTH,
    'CentreOfPotential': DIM_LENGTH,
    'GasSpin': DIM_SPECIFIC_ANGULAR_MOMENTUM,
    'GroupNumber': None,
    'HalfMAssProjRad': DIM_LENGTH,
    'HalfMassRad': DIM_LENGTH,
    'IDMostBound': None,
    'InertiaTensor': DIM_MOMENT_OF_INERTIA,
    'InitialMassWeightedBirthZ': None,
    'InitialMassWeightedStellarAge': None,
    'IronFromSNIa': None,
    'IronFromSNIaSmoothed': None,
    'KineticEnergy': DIM_ENERGY,
    'Mass': DIM_MASS,
    'MassFromAGB': DIM_MASS,
    'MassFromSNII': DIM_MASS,
    'MassFromSNIa': DIM_MASS,
    'MassTwiceHalfMassRad': DIM_MASS,
    'MassType': DIM_MASS,
    'MassWeightedEntropy': DIM_ENTROPY,
    'MassWeightedPotential': DIM_SPECIFIC_ENERGY,
    'MassWeightedTemperature': DIM_TEMPERATURE,
    'Metallicity': None,
    'MetalsFromAGB': DIM_MASS,
    'MetalsFromSNII': DIM_MASS,
    'MetalsFromSNIa': DIM_MASS,
    'Parent': None,
    'SFR': DIM_MASS_RATE,
    'SmoothedMetallicity': None,
    'Spin': DIM_SPECIFIC_ANGULAR_MOMENTUM,
    'StarFormationRate': DIM_MASS_RATE,
    'StellarInitialMass': DIM_MASS,
    'StellarVelDisp': DIM_VELOCITY,
    'StellarVelDisp_HalfMassRad': DIM_VELOCITY,
    'SubGroupNumber': None,
    'SubLength': None,
    'SubLengthType': None,
    'SubOffset': None,
    'ThermalEnergy': DIM_ENERGY,
    'TotalEnergy': DIM_ENERGY,
    'Velocity': DIM_VELOCITY,
    'VelDisp': DIM_VELOCITY,
    'Vmax': DIM_VELOCITY,
    'VmaxRadius': DIM_LENGTH
    }

FOF_DIMENSIONS = {
    'ContaminationCount': None,
    'ContaminationMass': DIM_MASS,
    'FirstSubhaloID': None,
    'GroupCentreOfPotential': DIM_LENGTH,
    'GroupLenth': None,
    'GroupMass': DIM_MASS,
    'GroupOffset': None,
    'Group_M_Crit200': DIM_MASS,
    'Group_M_Crit2500': DIM_MASS,
    'Group_M_Crit500': DIM_MASS,
    'Group_M_Mean200': DIM_MASS,
    'Group_M_Mean2500': DIM_MASS,
    'Group_M_Mean500': DIM_MASS,
    'Group_M_TopHat200': DIM_MASS,
    'Group_R_Crit200': DIM_LENGTH,
    'Group_R_Crit2500': DIM_LENGTH,
    'Group_R_Crit500': DIM_LENGTH,
    'Group_R_Mean200': DIM_LENGTH,
    'Group_R_Mean2500': DIM_LENGTH,
    'Group_R_Mean500': DIM_LENGTH,
    'Group_R_TopHat200': DIM_LENGTH,
    'NumOfSubhalos': None
    }

PARTICLE_DIMENSIONS = {
    'AExpMaximumTemperature': None,
    'Coordinates': DIM_LENGTH,
    'Density': DIM_DENSITY,
    'Entroyp': DIM_ENTROPY,
    'GroupNumber': None,
    'HostHalo_TVir_Mass': DIM_TEMPERATURE,
    'InternalEnergy': DIM_ENERGY,
    'IronMassFracFromSNIa': None,
    'Mass': DIM_MASS,
    'MaximumTemperature': DIM_TEMPERATURE,
    'MetalMassFracFromAGB': None,
    'MetalMassFracFromSNII': None,
    'MetalMassFracFromSNIa': None,
    'MetalMassWeightedRedshift': None,
    'Metallicity': None,
    'OnEquationOfState': None,
    'ParticleIDs': None,
    'SmoothedIronMassFracFromSNIa': None,
    'SmoothedMetallicity': None,
    'SmoothingLength': DIM_LENGTH,
    'StarFormationRate': DIM_MASS_RATE,
    'Temperature': DIM_TEMPERATURE,
    'TotalMassFromAGB': DIM_MASS,
    'TotalMassFromSNII': DIM_MASS,
    'TotalMassFromSNIa': DIM_MASS,
    'Velocity': DIM_VELOCITY
    }


def get_dimensions(base_group, dataset_name):
    """Get the dimensions of a specified data set."""
    real_name = dataset_name.split('/')[-1]

    # Shortcut for any kind of element abundance:
    if real_name in ['Carbon', 'Helium', 'Hydrogen', 'Iron',
                     'Magnesium', 'Neon', 'Nitrogen', 'Oxygen', 'Silicon']:
        return None

    if base_group == 'Subhalo':
        try:
            dimension = SUBHALO_DIMENSIONS[real_name]
        except KeyError:
            print(f"Unknown data set '{real_name}'!")
            return None
        return dimension

    elif base_group == 'FOF':
        try:
            dimension = FOF_DIMENSIONS[real_name]
        except KeyError:
            print(f"Unknown data set '{real_name}'!")
            return None
        return dimension

    elif base_group.startswith('PartType'):
        try:
            dimension = PARTICLE_DIMENSIONS[real_name]
        except KeyError:
            print(f"Unknown data set '{real_name}'!")
            return None
        return dimension

    else:
        print(f"No dimension information for group '{base_group}'...")
        return None


def get_unit_conversion_factor(dimension, system1, system2):
    """Get the conversion factor for a dimension from unit system 1 to 2.

    This is the number that values in system 1 must be multiplied with
    to convert them to system 2.
    """
    if dimension is None:
        return None

    try:
        if system1.lower() == 'si':
            system1_to_si = 1
        elif system1.lower() == 'cgs':
            system1_to_si = CGS_UNITS[dimension]
        elif system1.lower() == 'astro':
            system1_to_si = ASTRO_UNITS[dimension]
        else:
            print(f"Unsupported unit system {system1}!")
            set_trace()

        if system2.lower() == 'si':
            system2_to_si = 1
        elif system2.lower() == 'cgs':
            system2_to_si = CGS_UNITS[dimension]
        elif system2.lower() == 'astro':
            system2_to_si = ASTRO_UNITS[dimension]
        else:
            print(f"Unsupported unit system {system2}!")
            set_trace()
    except KeyError:
        print(f"Unknown dimension {dimension}!")
        set_trace()

    return system1_to_si / system2_to_si


def si_to_astro_factor(dimension):
    """Find the conversion factor from CGS to astronomical."""
    try:
        return 1/ASTRO_UNITS[dimension]
    except KeyError:
        print(f"Unknown dimension '{dimension}'!")
        set_trace()


def si_to_cgs_factor(dimension):
    """Find the conversion factor from CGS to SI."""
    try:
        return 1/CGS_UNITS[dimension]
    except KeyError:
        print(f"Unknown dimension '{dimension}'!")
        set_trace()

