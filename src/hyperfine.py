import numpy as np
from constants import *


def GetHyperFineLevels(pm, mf, delW, mN, gI):
    gL = 1  # Orbital g-factor
    gS = 2  # Electron spin g-factor
    L = 0  # Orbital Angular Momentum
    S = .5  # Electron spin
    I = .5  # Nuclear spin
    J = .5  # Total Angular Momentum

    muN = mue * meeV / (mN * 1e9)  # magnetic moment of nucleus

    gJ = gL * (J * (J + 1) + L * (L + 1) - S * (S + 1)) / (2 * J * (J + 1)) + ge * (
                J * (J + 1) - L * (L + 1) + S * (S + 1)) / (2 * J * (J + 1))

    x = B_Values * (gJ * mue - gI * muN) / (h * delW)
    Term1 = -h * delW / (2 * (2 * I + 1)) * np.ones_like(B_Values)
    Term2 = muN * gI * mf * B_Values

    if (abs(mf) == abs(I + .5)):
        sgn = mf / (I + .5)
        Term3 = h * delW / 2 * (1 + sgn * x)
    else:
        Term3 = pm * h * delW / 2 * np.sqrt(1 + 2 * mf * x / (I + .5) + x ** 2)

    delE = (Term1 + Term2 + Term3) / h

    return delE * h * J2eV


def AllHFLevels(delW, mN, gI):
    delEs = []
    for pm in [-1, 1]:
        F = .5 + pm / 2
        for mF in np.arange(-F, F + 1, 1):
            delEs.append(GetHyperFineLevels(pm, mF, delW, mN, gI))
    delEs = np.array(delEs)
    delEs = np.sort(delEs, axis=0)
    delEDict = {}
    for i in range(0, 4):
        letter = chr(97 + i)
        delEDict[letter] = delEs[i]
    return delEDict