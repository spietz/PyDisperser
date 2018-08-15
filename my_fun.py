#############################################
# Expressions from the exercise formulation #
#############################################

import numpy as np
from scipy.integrate import simps


def ffun(Re):
    return (2.5*np.log(Re)+5)


def ufun(y, Re):
    return (2.5*np.log(y+np.exp(-2.0)/Re)+ffun(Re))


def vfun(y, Re):
    return (-0.017*y+0.04)*ffun(Re)


def Um(Re, n_int=1000):
    Y = np.linspace(0, 1, n_int)
    return np.trapz(ufun(Y, Re), x=Y)
    # return simps(ufun(Y, Re), x=Y)

