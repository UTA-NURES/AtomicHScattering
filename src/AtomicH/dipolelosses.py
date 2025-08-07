import constants
import elastic
import potentials
import hyperfine
import spinbasis
import sympy
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad


DipoleChannels=[]
DipoleChannels.append({'alpha':'d','beta':'d','alphaprime':'a','betaprime':'a'})
DipoleChannels.append({'alpha':'d','beta':'d','alphaprime':'a','betaprime':'c'})
DipoleChannels.append({'alpha':'d','beta':'d','alphaprime':'a','betaprime':'d'})
DipoleChannels.append({'alpha':'d','beta':'d','alphaprime':'c','betaprime':'c'})
DipoleChannels.append({'alpha':'d','beta':'d','alphaprime':'c','betaprime':'d'})

# A useful function for comparing to Stoof et al
def B_Naught(B_Values):
    return (1 + B_Values/3.17e-3)

def p_of_temp(mu, T):
    return np.sqrt(2 * mu * constants.kb * constants.J2eV * T)

def pprime(pin, epsa, epsb, epsprimea, epsprimeb, mu):
    E = pin ** 2 / (2 * mu)
    Eprime = E + epsa + epsb - epsprimea - epsprimeb
    pprime = np.sqrt(2 * mu * Eprime)
    return pprime

def p_abs(mu, pin, epsa, epsb, epsprimea, epsprimeb):
    psquared = pin**2 + mu * (epsa + epsb - epsprimea - epsprimeb)
    return np.sqrt(psquared)

def GetIntegral(rhos,alphain, betain, alphaout, betaout, mu, temp, potential, how_to_int, lin=0, lout=2):

    P1 = p_of_temp(mu, temp)
    P2 = pprime(P1, alphain, betain, alphaout, betaout, mu)

    InState =  np.array(elastic.Wave_Function(rhos, P1, lin, mu, potential, how_to_int)[0])
    OutState = np.array(elastic.Wave_Function(rhos, P2, lout, mu, potential, how_to_int)[0])

    Integrand = interp1d(rhos, InState * OutState / rhos ** 3, kind='quadratic')
    Integral = quad(Integrand, rhos[0], rhos[-1])[0] / (P1 * P2)
    return Integral


def GetSpatialPart(channel=DipoleChannels[0], B_value=1e-5, consts=constants.HydrogenConstants, Temperature=5e-4, potential=potentials.Silvera_Triplet,rhos=np.linspace(1e-9,0.75,2000),lin=0,lout=2,how_to_int='Radau'):
    HFLevels=hyperfine.AllHFLevels(B_value, consts)

    aHf =  HFLevels[channel['alpha']]
    bHf =  HFLevels[channel['beta']]
    apHf = HFLevels[channel['alphaprime']]
    bpHf = HFLevels[channel['betaprime']]

    Pin = p_of_temp(consts.mu, Temperature)
    Pout = pprime(Pin, aHf, bHf, apHf, bpHf, consts.mu)

    Integral    = GetIntegral(rhos,aHf, bHf, apHf, bpHf, consts.mu, Temperature, potential, how_to_int,lin,lout)
    SpatialPart = Pout * consts.mu * Integral**2

    return(SpatialPart)


def GetSpinPart(channel=DipoleChannels[0], B_value=1e-5, consts=constants.HydrogenConstants):
    NormDiff = 4*np.sqrt(6)
    Rets=spinbasis.GetRotatedElements()
    theta = hyperfine.Theta(2*consts.delW, B_value, consts.gam)
    value = 0
    for m in Rets.keys():
        El = ( spinbasis.GetElement(Rets[m], channel['alpha'], channel['beta'], 1, channel['alphaprime'], channel['betaprime'], 1)) **2
        try:
            value += El.subs(spinbasis.sr2, np.sqrt(2)) \
                    .subs(spinbasis.sr3, np.sqrt(3)) \
                    .subs(spinbasis.c, np.cos(theta)) \
                    .subs(spinbasis.s, np.sin(theta))
        except:
            value += 0
    SpinPart = value*NormDiff**2
    return(SpinPart)


def GetGFactor( channel=DipoleChannels[0],  B_value=1e-5, consts=constants.HydrogenConstants, Temperature=5e4, potential=potentials.Silvera_Triplet,rhos=np.linspace(1e-9,0.75,2000),lin=0,lout=2):
    mue = np.sqrt(4 * np.pi * constants.finestructure) / (2 * constants.meeV)
    Pre_Factor = 1 / (5 * np.pi) * mue ** 4 * constants.NatUnits_cm3sm1

    SpatialMatrixElementSq = GetSpatialPart( channel, B_value, consts, Temperature, potential, rhos, lin, lout,'Radau')
    SpinMatrixElementSq    = GetSpinPart(    channel, B_value, consts)

    return(Pre_Factor * SpatialMatrixElementSq * SpinMatrixElementSq)


def GetSummedGFactor( channel=DipoleChannels[0],  B_value=1e-5, consts=constants.HydrogenConstants, Temperature=5e4, potential=potentials.Silvera_Triplet,rhos=np.linspace(1e-9,0.75,2000)):

    degeneracies =  [1,      1,      3,      5,      5,      5,      5]
    PWaves =       [[0, 2], [2, 0], [2, 2], [2, 4], [4, 2], [4, 4], [4, 6]]

    G=0
    for pi in range(0,len(PWaves)):
        G+=GetGFactor(channel,  B_value, consts, Temperature, potential,rhos,PWaves[pi][0],PWaves[pi][1])*degeneracies[pi]
    return G
