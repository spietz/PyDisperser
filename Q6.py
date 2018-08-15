# -*- coding: utf-8 -*-
#!/usr/bin/python2

#################
# PLOT PROFILES #
#################

from mymod import PyDisperser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import my_fun

plt.close("all")

# Customizing plots
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 14

## Inputs #############################
N = 15000
T = np.linspace(0.0, 10.0, 3)
# Remember that beta<1. beta = W/kappa, kappa = 0.41 (eq. 4.102)
WVec = np.linspace(0, 0.25, 15)
y_l = 0.0
y_u = 1.0
max_step = 1
Re = 1000

# Init stuff
# Time between saves
dT = np.diff(T)
Xmean = np.zeros((len(dT), len(WVec)))
Xvar = np.zeros((len(dT), len(WVec)))
Um = np.zeros(len(WVec))
D1 = np.zeros(len(WVec))
TL_Xvar = np.zeros((len(dT), len(WVec)))
TL_Xmean = np.zeros((len(dT), len(WVec)))
Y = np.zeros((len(dT), len(WVec), N))
X = np.zeros((len(dT), len(WVec), N))

# ## Functions for plots
# def init():
#     """initialize animation"""
#     im.set_data(np.random.random((5,5)))
#     return [im]

# def animate(i):
#     """perform animation step"""

# Loop Re numbers
for r in range(0, len(WVec)):
    W = WVec[r]

    # declare/redeclare new Disperser object from wrapped cpp, mummm cpp...
    D = PyDisperser(Re, y_l, y_u)
    D.gaussian_u = False
    D.u_f = False
    D.gaussian_v = False
    D.bursting_process = False
    D.bouncing_particles = True
    D.direction_likelyhood = True
    D.fall_vel = W
    D.displacement = 0.0

    # set particle distribution at x=0
    Y[0, r, :] = np.linspace(y_l, y_u, N)

    # Get results for each T by updating with dT
    for t in range(0, len(dT)):

        D.update(N, X[t, r, :], Y[t, r, :], dT[t], max_step)

        if (t+1 != len(dT)):
            X[t+1, r, :] = X[t, r, :]
            Y[t+1, r, :] = Y[t, r, :]

        # Calculate mean and variance
        Xmean[t, r] = np.mean(X[t, r, :])
        Xvar[t, r] = np.var(X[t, r, :])

    # Linear regressions
    p = np.polyfit(T[1:len(T)], Xmean[:, r], 1)
    TL_Xmean[:, r] = np.polyval(p, T[1:len(T)])
    Um[r] = p[0]

    p = np.polyfit(T[1:len(T)], Xvar[:, r], 1)
    TL_Xvar[:, r] = np.polyval(p, T[1:len(T)])
    D1[r] = 0.5*p[0]


sti = 'plots/heavy'
fname = '%s/Wrun.txt' % (sti)
## Write to file #######################
f = open(fname, 'w')
f.write('W\tRe\tUm_ana\tD1_ana\tUm_num\tD1_num\tUm_err\tD1_err\n')
for r in range(0, len(WVec)):
    f.write('%0.3f\t%i\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n'
            % (WVec[r],Re, my_fun.Um(Re), 5.93,
               Um[r], D1[r],
               100*np.abs(my_fun.Um(Re)-Um[r])/my_fun.Um(Re),
               100*np.abs(5.93-D1[r])/5.93))
f.close()

## CONCENTRATION PLOTS ##
# Using imshow

