import numpy as np

# Fundamental constants
h              = 6.6260715e-34 # Planck constant in m^2 kg / s
hbar           = h/(2*np.pi)   #     h / 2pi
C              = 2.99792458e8  # Speed of light in m/s
kb             = 1.380649e-23  # Boltzmann constant in JK^-1
ge             = 2.002319      # Electron g-factor
meeV           = .511e6        # mass of electron
mue            = 9.27e-24      # magnetic moment of electron
game           = -28024.9*1e6  # electron gamma factor in Hz T^-1
finestructure  = 1/137

# Unit conversions
BohrInAng      = .529177210544
HartreeInEV    = 27.211386245981
hcInEVAngstrom = 1973.2698044
BohrInEV       = BohrInAng/hcInEVAngstrom
K2eV           = 8.617333262e-5
cmm1_in_eV     = 1.23984198E-4
J2eV           = 6.242e18
DaltonInEV     = 931.49410372*1e6

# H and T constants
delWH          = 1.4204057517667e9 # Hyperfine splitting of hydrogen in Hz
delWT          = 1.516701396e9     # Hyperfine splitting of tritium in Hz
gIH            = 5.585694702  # Hydrogen nuclear g-factor
gIT            = 5.95792492   # Tritium nuclear g-factor
mH             = 1.00784      # Mass of hydrogen in Dalton
mT             = 3.01604928   # Mass of tritium in Dalton
muT            = mT*DaltonInEV/2
muH            = mH*DaltonInEV/2
gamH           = 42.577*1e6   # Gyromagnetic constant for H In Hz T^-1
gamT            = 45.415*1e6  # Gyromagnetic constant for T In Hz T^-1