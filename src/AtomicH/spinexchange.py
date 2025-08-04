import elastic
import dipolelosses
import hyperfine
import constants
import spinbasis
import potentials
import numpy as np


def GetSpatialPart(mu, alpha,beta,alphaprime,betaprime, HFLevels, triplet_potential, singlet_potential,temp=5e-4,rhos=np.linspace(1e-9,100*constants.BohrInAng/constants.hcInEVAngstrom,100),how_to_int='Radau'):
    SpatialPart = []

    for i in range(0,len(HFLevels[alpha])):
        aHf = HFLevels[alpha][i]
        bHf = HFLevels[beta][i]
        apHf = HFLevels[alphaprime][i]
        bpHf = HFLevels[betaprime][i]

        Pin    = dipolelosses.p_of_temp(mu, temp)
        Pout   = dipolelosses.pprime(Pin, aHf, bHf, apHf, bpHf, mu)
        Pabs   = dipolelosses.p_abs(mu, Pin, aHf, bHf, apHf, bpHf)

        const = np.pi / (mu * Pin)*constants.NatUnits_cm3sm1
        tdeltaaa = elastic.GetPhaseShift(rhos, Pabs, 0, mu, triplet_potential, how_to_int)[-1]
        sdeltaaa = elastic.GetPhaseShift(rhos, Pabs, 0, mu, singlet_potential, how_to_int)[-1]

        SpatialPart.append(const * (Pin * Pout / Pabs ** 2) * (np.sin(tdeltaaa - sdeltaaa) ** 2))

    SpatialPart = np.array(SpatialPart)

    return (SpatialPart)


def GetSpinPart(alpha,beta,alphaprime, betaprime,delW, B_values, gamN):
    SpinPart = []
    for Bs in B_values:
        th = hyperfine.Theta(delW, Bs, gamN)

        trans = spinbasis.TransformMatrix(spinbasis.TripletProj - spinbasis.SingletProj, spinbasis.Rotator)
        El = (spinbasis.GetElement(trans, alpha, beta, 1, alphaprime, betaprime, 1)) ** 2
        Value = El.subs(spinbasis.sr2, np.sqrt(2)).subs(spinbasis.sr3, np.sqrt(3)).subs(spinbasis.c, np.cos(th)).subs(spinbasis.s, np.sin(th))

        SpinPart.append(Value)
    SpinPart = np.array(SpinPart)
    return SpinPart


def GetGFactor(alpha='d',beta='d',alphaprime='a',betaprime='a',which='T', B_values=np.logspace(-3,1,50),triplet_potential=potentials.Silvera_Triplet,singlet_potential=potentials.Silvera_Singlet,temp=5e-4,rhos=np.linspace(1e-9,0.75,20000)):
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

    HFLevels = hyperfine.AllHFLevels(B_values,delW, mu, gI)

    SpinPart    = GetSpinPart(alpha,beta,alphaprime, betaprime,delW, B_values, gam)
    SpatialPart = GetSpatialPart(mu, alpha,beta,alphaprime,betaprime, HFLevels, triplet_potential, singlet_potential,temp,rhos)

    return (SpinPart*SpatialPart)

