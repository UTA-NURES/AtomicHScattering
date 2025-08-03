import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
from constants import *

def Silvera_Triplet(Rho):
    R = Rho * hcInEVAngstrom / 4.16
    D = 1.28
    F=(R>D)+(R<D)*np.exp(-(D/R-1)**2)
    return 6.46 * K2eV * (4.889e4*np.exp(0.0968-8.6403*R-2.427*R**2)-(1.365/R**6+0.425/R**8+0.183/R**10)*F)

def Silvera_Triplet2(Rho):
    R = Rho * hcInEVAngstrom / BohrInAng
    P = np.exp(0.09678-1.10173*R-0.03945*R**2)+np.exp(-(10.04/R-1)**2)*(-6.5/R**6-124/R**8-3285/R**10)
    return P * HartreeInEV

def Silvera_J(Rho):
    R = Rho * hcInEVAngstrom / BohrInAng
    P = np.exp(-.288-.275*R-.176*R**2+.0068*R**3)
    return P * HartreeInEV

def Silvera_Singlet(R):
    return Silvera_Triplet(R) - Silvera_J(R)

path=os.path.dirname(os.path.abspath(__file__))
dat_Adiabatic=pd.read_excel(path+"/InputData/AdiabaticCorrection_Kolos.xlsx")
Hvals=np.array(dat_Adiabatic.Hprime)
DVals=np.array(dat_Adiabatic.D)
Corr=((Hvals[0:-1]-Hvals[-1])/(DVals[0:-1]-DVals[-1]))
CorrInterp=interp1d(dat_Adiabatic.R[0:-1]*BohrInAng/hcInEVAngstrom,Corr,bounds_error=False,kind='linear',fill_value=(Corr[0],Corr[-1]))

def UnCorrectedSilvera_Triplet(R):
    return Silvera_Triplet(R)*(1-CorrInterp(R))

def TCorrectedSilvera_Triplet(R):
    return Silvera_Triplet(R)*(1-(2/3)*CorrInterp(R))
