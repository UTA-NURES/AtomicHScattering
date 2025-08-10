# ================================================================================
# the various potentials that can be used in the subsequent calculations.
# ================================================================================

import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
from constants import *

path=os.path.dirname(os.path.abspath(__file__))


#=====================================
# Silvera potentials in analytic form
#=====================================

# From European Scientific Journal. 2012 Oct 1;8(24).
def Silvera_Triplet(Rho):
    R = Rho * hcInEVAngstrom / 4.16
    D = 1.28
    F=(R>D)+(R<D)*np.exp(-(D/R-1)**2)
    return 6.46 * K2eV * (4.889e4*np.exp(0.0968-8.6403*R-2.427*R**2)-(1.365/R**6+0.425/R**8+0.183/R**10)*F)

# From Reviews of Modern Physics, 52(2), p.393. and Progress in Low Temperature Physics. 1986 Jan 1;10:139-370.
def Silvera_Triplet2(Rho):
    R = Rho * hcInEVAngstrom / BohrInAng
    P = np.exp(0.09678-1.10173*R-0.03945*R**2)+np.exp(-(10.04/R-1)**2)*(-6.5/R**6-124/R**8-3285/R**10)
    return P * HartreeInEV

def Silvera_Triplet2NoR8(Rho):
    R = Rho * hcInEVAngstrom / BohrInAng
    P = np.exp(0.09678-1.10173*R-0.03945*R**2)+np.exp(-(10.04/R-1)**2)*(-6.5/R**6)
    return P * HartreeInEV

# From Reviews of Modern Physics, 52(2), p.393.
def Silvera_J(Rho):
    R = Rho * hcInEVAngstrom / BohrInAng
    P = np.exp(-.288-.275*R-.176*R**2+.0068*R**3)
    return P * HartreeInEV

# From Reviews of Modern Physics, 52(2), p.393
# Note - Silvera says this is not an accurate representation, use with extreme caution
def Silvera_Singlet(R):
    return Silvera_Triplet(R) - Silvera_J(R)

#==========================
# The adiabatic correction
#==========================

# Adiabatic correction from  Journal of Molecular Spectroscopy. 1990 Oct 1;143(2):237-50.
dat_Adiabatic=pd.read_excel(path+"/InputData/AdiabaticCorrection_Kolos.xlsx")
Hvals=np.array(dat_Adiabatic.Hprime)
DVals=np.array(dat_Adiabatic.D)
Corr=((Hvals[0:-1]-Hvals[-1])/(DVals[0:-1]-DVals[-1]))
CorrInterp=interp1d(dat_Adiabatic.R[0:-1]*BohrInAng/hcInEVAngstrom,Corr,bounds_error=False,kind='linear',fill_value=(Corr[0],Corr[-1]))

def UnCorrectedSilvera_Triplet(R):
    return Silvera_Triplet(R)*(1-CorrInterp(R))

def TCorrectedSilvera_Triplet(R):
    return Silvera_Triplet(R)*(1-(2/3)*CorrInterp(R))


#===========================
# Kolos Potentials
#===========================

#The Journal of Chemical Physics. 1965 Oct 1;43(7):2429-41.
dat_KolosSinglet1 = np.genfromtxt(path+"/InputData/Singlet_Kolos1965.csv", delimiter=',', skip_header=1)

#Kolos 1974, Chemical Physics Letters. 1974 Feb 15;24(4):457-60.
dat_KolosSinglet2 = np.genfromtxt(path+"/InputData/Singlet_Kolos1974.csv", delimiter=',', skip_header=1)

interp_KolosSinglet1 = interp1d(dat_KolosSinglet1[:,0], dat_KolosSinglet1[:,1], kind = 'cubic',bounds_error=False, fill_value='extrapolate')
interp_KolosSinglet2 = interp1d(dat_KolosSinglet2[:,0], -dat_KolosSinglet2[:,1], kind = 'cubic',bounds_error=False,fill_value='extrapolate')

# From Journal of molecular spectroscopy, 54(2), pp.303-311
def Kolos_Singlet1(rho):
    return (interp_KolosSinglet1(rho/BohrInAng * hcInEVAngstrom) + 1) * HartreeInEV

#Kolos 1974, Chemical Physics Letters. 1974 Feb 15;24(4):457-60.
def Kolos_Singlet2(rho):
    return (interp_KolosSinglet2(rho/ BohrInAng * hcInEVAngstrom) + 1) * HartreeInEV

def Kolos_SingletCombo(rho):
    Rp = rho / BohrInAng * hcInEVAngstrom
    split1_2=6.3  # Make the split within the range where function has support (between last 2 bins)
    split2_3=11.5 #  [As above]
    Decide1=(Rp<split1_2)
    Decide2=(Rp>=split1_2)*(Rp<split2_3)
    Decide3=(Rp>=split2_3)
    return (Decide1*(interp_KolosSinglet1(Rp) + 1)  + Decide2*(interp_KolosSinglet2(Rp) + 1) + Decide3*(interp_KolosSinglet2(split2_3) + 1)*(Rp/split2_3)**-6) * HartreeInEV

dat_KolosTriplet = np.genfromtxt(path+"/InputData/Triplet_Kolos1974.csv", delimiter=',', skip_header=1)
interp_KolosTriplet = interp1d(dat_KolosTriplet[:,0], dat_KolosTriplet[:,1], kind = 'cubic', bounds_error = False, fill_value = 'extrapolate')

# From Journal of molecular spectroscopy, 54(2), pp.303-311
def Kolos_Triplet(rho):
    Rp = rho / BohrInAng * hcInEVAngstrom 
    maxx = np.max(interp_KolosTriplet.x)
    Decide = (Rp < maxx)
    return (Decide * (interp_KolosTriplet(Rp) + 1) + (1 - Decide) * (interp_KolosTriplet(maxx) + 1) * (Rp / maxx) ** -6) * HartreeInEV


#==========================
# Jamieson Potentials
#==========================
# All from Physical Review A. 2000 Mar 6;61(4):042705.

dat_Jamieson = pd.read_excel(path+"/InputData/Potential_Jamieson.xlsx")

interp_JamiesonSinglet = interp1d(dat_Jamieson.R, dat_Jamieson.Singlet, kind = 'cubic',bounds_error=False, fill_value='extrapolate')
interp_JamiesonTriplet = interp1d(dat_Jamieson.R, dat_Jamieson.Triplet, kind = 'cubic',bounds_error=False,fill_value='extrapolate')

def Jamieson_Singlet(rho):
    return (interp_JamiesonSinglet(rho/BohrInAng * hcInEVAngstrom) + 1) * HartreeInEV

def Jamieson_Triplet(rho):
    return (interp_JamiesonTriplet(rho/ BohrInAng * hcInEVAngstrom) + 1) * HartreeInEV

def Jamieson_TripletCombo(rho):
    Rp=rho/BohrInAng * hcInEVAngstrom
    maxx=np.max(interp_JamiesonTriplet.x)
    Decide = (Rp < maxx)
    return (Decide * (interp_JamiesonTriplet(Rp)+1) + (1-Decide)*(interp_JamiesonTriplet(maxx)+1)*(Rp/maxx)**-6 ) * HartreeInEV

def Jamieson_SingletCombo(rho):
    Rp = rho / BohrInAng * hcInEVAngstrom
    maxx = np.max(interp_JamiesonSinglet.x)
    Decide = (Rp < maxx)
    return (Decide * (interp_JamiesonSinglet(Rp) + 1) + (1 - Decide) * (interp_JamiesonSinglet(maxx) + 1) * (Rp / maxx) ** -6 )* HartreeInEV

def Jamieson_Kolos_Mixed(rho):
    Rp = rho / BohrInAng * hcInEVAngstrom
    maxx = np.max(interp_JamiesonSinglet.x)
    minx = np.min(interp_JamiesonSinglet.x)
    Decide1=(Rp<minx)
    Decide2 = (Rp>=minx)*(Rp<maxx)
    Decide3=(Rp>=maxx)
    V1 = (interp_func(Rp) + 1) * Decide1
    V2 = (interp_func_sing(Rp) + 1) * Decide2
    V3 = (interp_func_sing(split2_3) + 1) * (Rp / split2_3) ** -6 * Decide3

    return (V1 + V2 + V3) * HartreeInEV

# ==========================================
#For illustration only -
# radial dependence of the dipole potential
#===========================================

# From Stoof et al, Physical Review B 38.7 (1988): 4688.
def DipoleRadialPart(rho):
    muel = np.sqrt(4 * np.pi * finestructure) / (2 * meeV)
    return muel**2/(4*np.pi*rho**3)*(4*np.pi/5)


# ==========================================
# Compilation of useful potentials
#===========================================

Triplets = {"Silvera":Silvera_Triplet,
            "Silvera2":Silvera_Triplet2,
            "Uncorrected Silvera":UnCorrectedSilvera_Triplet,
            "SilveraNoR8":Silvera_Triplet2NoR8,
            "Kolos":Kolos_Triplet,
            "Jamieson":Jamieson_TripletCombo}

Singlets = {"Kolos":Kolos_SingletCombo,
            "Pure Jamieson":Jamieson_SingletCombo,
            "Kolos Jamieson Mixed":Jamieson_Kolos_Mixed}
