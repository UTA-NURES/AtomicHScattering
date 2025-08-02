import numpy as np
from constants import *

def Silvera_Triplet(Rho):
    R = Rho * hcInEVAngstrom / 4.16
    D = 1.28
    F=(R>D)+(R<D)*np.exp(-(D/R-1)**2)
    return 6.46 * K2eV * (4.889e4*np.exp(0.0968-8.6403*R-2.427*R**2)-(1.365/R**6+0.425/R**8+0.183/R**10)*F)

def SilveraTriplet2(Rho):
    R = Rho * hcInEVAngstrom / BohrInAng
    P = np.exp(0.09678-1.10173*R-0.03945*R**2)+np.exp(-(10.04/R-1)**2)*(-6.5/R**6-124/R**8-3285/R**10)
    return P * HartreeInEV

def SilveraJ(Rho):
    R = Rho * hcInEVAngstrom / BohrInAng
    P = np.exp(-.288-.275*R-.176*R**2+.0068*R**3)
    return P * HartreeInEV

def SilveraSinglet(R):
    return Silvera_Triplet(R) - SilveraJ(R)

dat_Adiabatic=pd.read_excel("InputData/AdiabaticCorrection_Kolos.xlsx")
Hvals=np.array(dat_Adiabatic.Hprime)
DVals=np.array(dat_Adiabatic.D)
Corr=((Hvals-Hvals[-1])/(DVals-DVals[-1]))[0:-1]
CorrInterp=interp1d(dat.R[0:-1]*BohrInAng/hcInEVAngstrom,Corr,bounds_error=False,kind='linear',fill_value=(Corr[0],Corr[-1]))

def UnCorrectedSilveraTripletH(R):
    return Silvera_Triplet(R)*(1-CorrInterp(R))

def TCorrectedSilveraTriplet(R):
    return Silvera_Triplet(R)*(1-(2/3)*CorrInterp(R))
