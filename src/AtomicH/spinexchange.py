import elastic
import dipolelosses
import hyperfine
import constants
import spinbasis
import potentials
import numpy as np

SpinExChannels=[]
SpinExChannels.append({'alpha':'c','beta':'c','alphaprime':'a','betaprime':'a'})
SpinExChannels.append({'alpha':'c','beta':'c','alphaprime':'a','betaprime':'c'})
SpinExChannels.append({'alpha':'c','beta':'c','alphaprime':'b','betaprime':'d'})


def GetSpatialPart(channel=SpinExChannels[0], B_value=1e-5, consts=constants.HydrogenConstants, Temperature=5e-4, triplet_potential=potentials.Silvera_Triplet,singlet_potential=potentials.Kolos_SingletCombo,rhos=np.linspace(1e-9,0.75,2000),how_to_int='Radau'):
    HFLevels = hyperfine.AllHFLevels(B_value, consts)

    aHf =  HFLevels[channel['alpha']]
    bHf =  HFLevels[channel['beta']]
    apHf = HFLevels[channel['alphaprime']]
    bpHf = HFLevels[channel['betaprime']]

    Pin  = dipolelosses.p_of_temp(consts.mu, Temperature)
    Pout = dipolelosses.pprime(Pin, aHf, bHf, apHf, bpHf, consts.mu)
    Pabs = dipolelosses.p_abs(consts.mu, Pin, aHf, bHf, apHf, bpHf)

    const = np.pi / (consts.mu * Pin)*constants.NatUnits_cm3sm1
    tdeltaaa = elastic.GetPhaseShift(rhos, Pabs, 0, consts.mu, triplet_potential, how_to_int)[-1]
    sdeltaaa = elastic.GetPhaseShift(rhos, Pabs, 0, consts.mu, singlet_potential, how_to_int)[-1]

    SpatialPart = (const * (Pin * Pout / Pabs ** 2) * (np.sin(tdeltaaa - sdeltaaa) ** 2))

    return (SpatialPart)


def GetSpinPart(channel=SpinExChannels[0], B_value=1e-5, consts=constants.HydrogenConstants):

    th = hyperfine.Theta(2*consts.delW, B_value, consts.gam)
    trans = spinbasis.TransformMatrix(spinbasis.TripletProj - spinbasis.SingletProj, spinbasis.Rotator)
    El = (spinbasis.GetElement(trans, channel['alpha'], channel['beta'], 1, channel['alphaprime'], channel['betaprime'], 1)) ** 2
    SpinPart = El.subs(spinbasis.sr2, np.sqrt(2)).subs(spinbasis.sr3, np.sqrt(3)).subs(spinbasis.c, np.cos(th)).subs(spinbasis.s, np.sin(th))

    return (SpinPart)


def GetGFactor(channel=SpinExChannels[0],  B_value=1e-5, consts=constants.HydrogenConstants, Temperature=5e4, triplet_potential=potentials.Silvera_Triplet,singlet_potential=potentials.Kolos_SingletCombo,rhos=np.linspace(1e-9,0.75,20000),how_to_int='Radau'):

    SpinPart    = GetSpinPart(    channel, B_value, consts)
    SpatialPart = GetSpatialPart( channel, B_value, consts, Temperature, triplet_potential, singlet_potential, rhos, how_to_int)

    return (SpinPart*SpatialPart)

