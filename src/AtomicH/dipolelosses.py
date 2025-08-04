import constants
import elastic
import potentials
import hyperfine
import spinbasis
import sympy
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

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



def GetSpatialPart(mu, alpha,beta,alphaprime,betaprime, B_values,HFLevels, potential,temp=5e-4,rhos=np.linspace(1e-9,0.75,20000),lin=0,lout=2):
    SpatialPart = []
    for i in range(0,len(HFLevels[alpha])):
        aHf = HFLevels[alpha][i]
        bHf = HFLevels[beta][i]
        apHf = HFLevels[alphaprime][i]
        bpHf = HFLevels[betaprime][i]

        Pin = p_of_temp(mu, temp)
        Pout = pprime(Pin, aHf, bHf, apHf, bpHf, mu)

        if B_values[i]<= 8.5:
            Integral = GetIntegral(rhos,aHf, bHf, apHf, bpHf, mu, temp, potential, 'Radau',lin,lout)
        else:
            Integral = GetIntegral(rhos,aHf, bHf, apHf, bpHf, mu, temp, potential, 'RK45',lin,lout)

        SpatialPart.append(Pout * mu * Integral**2)
        i = i + 1

    SpatialPart = np.array(SpatialPart)
    return(SpatialPart)


def GetSpinPart(delW, alpha,beta,alphaprime,betaprime,B_values, gam):


    NormDiff = 4*np.sqrt(6)
    rotator = []
    Rets=spinbasis.GetRotatedElements()
    for b in B_values:
        theta = hyperfine.Theta(delW, b, gam)
        value = 0
        for m in Rets.keys():
            El = ( spinbasis.GetElement(Rets[m], alpha, beta, 1, alphaprime, betaprime, 1)) **2
            try:
                value += El.subs(spinbasis.sr2, np.sqrt(2)) \
                        .subs(spinbasis.sr3, np.sqrt(3)) \
                        .subs(spinbasis.c, np.cos(theta)) \
                        .subs(spinbasis.s, np.sin(theta))
            except:
                value += 0
        rotator.append(value)

    SpinPart = np.array(rotator)*NormDiff**2

    return(SpinPart)


def GetGFactor(alpha='d',beta='d',alphaprime='a',betaprime='a',which='T', B_values=np.logspace(-3,1,50),potential=potentials.Silvera_Triplet,temp=5e-4,rhos=np.linspace(1e-9,0.75,20000)):

    if(which=='T'):
        delW = constants.delWT
        gam  = constants.gamT
        mu   = constants.muT
        gI   = constants.gIT

    elif(which=='H'):
        delW = constants.delWH
        gam  = constants.gamH
        mu   = constants.muH
        gI   = constants.gIH

    mue = np.sqrt(4 * np.pi * constants.finestructure) / (2 * constants.meeV)
    Pre_Factor = 1 / (5 * np.pi) * mue ** 4 * constants.NatUnits_cm3sm1

    HFLevels = hyperfine.AllHFLevels(B_values,delW, mu, gI)
    Spatials = GetSpatialPart(mu, alpha,beta,alphaprime,betaprime, B_values,HFLevels, potential,temp,rhos)
    Spins    = GetSpinPart(delW, alpha,beta,alphaprime,betaprime,B_values, gam)

    return(Pre_Factor * Spatials * Spins)

