#################
# PLOT PROFILES #
#################

from mymod import PyDisperser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import my_fun                   # expressions from the exercise formulation

plt.close("all")

# Customizing plots
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['font.size'] = 14

## Inputs #############################
Ny = 100
Re = 1000.0
y_l = 0.0
y_u = 1.0

# New Disperser object from wrapped cpp, mummm cpp...
D = PyDisperser(1000, y_l, y_u)
D.gaussian_u = True
D.gaussian_v = True
D.bursting_process = True
D.bouncing_particles = False
D.direction_likelyhood = False

# Init stuff
Y = np.linspace(y_l, y_u, Ny)
dY = Y.copy()
dX = np.zeros(Ny)
dT = np.zeros(Ny)
U = np.zeros(Ny)
V = np.zeros(Ny)
dT = np.zeros(Ny)
L = np.zeros(Ny)

# Get results from object by taking a single step for all y
for i in range(0, len(Y)):
    U[i]  = D.get_u(Y[i])
    V[i]  = D.get_v(Y[i])
    dT[i] = D.get_l(Y[i])/abs(V[i])
    L[i]  = D.get_l(Y[i])

D.update(Ny, dX, dY, 30.0, 1)  # tF skal død og pine være en double
dY = dY - Y

## Plotting ###########################

# from expressions:
fig1, axarr1 = plt.subplots(1, 2, sharey=True, num=1)
fig1.clf()
fig1.suptitle(r'$Re_f=$%i' % (Re))

## Horizontal particle velocity #######
axarr1[0] = fig1.add_subplot(1, 2, 1)
axarr1[0].plot(my_fun.ufun(Y, Re), Y, '-b', label=r'$\bar u(y)$')
axarr1[0].plot(my_fun.Um(Re)*np.ones(Ny), Y, '--r', label=r'$u_m$')
axarr1[0].set_xlabel(r'$\frac{u}{U_f}$')
axarr1[0].set_ylabel(r'$\frac{y}{h}$')
axarr1[0].grid(True)
axarr1[0].legend(loc=0)

## Vertical particale velocity ########
axarr1[1] = fig1.add_subplot(1, 2, 2)
axarr1[1].plot(my_fun.vfun(Y, Re), Y, '-')
axarr1[1].set_xlabel(r'$\frac{v}{U_f}$')
axarr1[1].grid(True)
for label in axarr1[1].get_yticklabels():
    label.set_visible(False)

fig1.show()

fname = 'plots/profiles1.pdf'
plt.savefig(fname)

# from disperser:

fig2, axarr2 = plt.subplots(1, 2, sharey=True, num=2)
fig2.clf()
fig2.suptitle(r'$Re_f=$%i' % (Re))

## Horizontal particle velocity #######
axarr2[0] = fig2.add_subplot(1, 2, 1)
axarr2[0].plot(U, Y, 'x')
axarr2[0].set_xlabel(r'$\frac{u}{U_f}$')
axarr2[0].set_ylabel(r'$\frac{y}{h}$')
axarr2[0].grid(True)

## Vertical particale velocity ########
axarr2[1] = fig2.add_subplot(1, 2, 2)
axarr2[1].plot(V, Y, 'x')
axarr2[1].set_xlabel(r'$\frac{v}{U_f}$')
axarr2[1].grid(True)
for label in axarr2[1].get_yticklabels():
    label.set_visible(False)

fig2.show()

fname = 'plots/profiles2.pdf'
plt.savefig(fname)

fig3, axarr3 = plt.subplots(1, 2, sharey=True, num=3)
fig3.clf()
fig3.suptitle(r'$Re_f=$%i' % (Re))

## Horizontal particle movement #######
axarr3[0] = fig3.add_subplot(1, 2, 1)
axarr3[0].plot(dX, Y, 'x')
axarr3[0].set_xlabel(r'$\frac{\Delta x}{h}$')
axarr3[0].set_ylabel(r'$\frac{y}{h}$')
axarr3[0].grid(True)

## Vertical particale movement ########
axarr3[1] = fig3.add_subplot(1, 2, 2)
axarr3[1].plot(dY, Y, 'x')
axarr3[1].set_xlabel(r'$\frac{\Delta y}{h}$')
axarr3[1].grid(True)
for label in axarr3[1].get_yticklabels():
    label.set_visible(False)

fig3.show()

fname = 'plots/profiles3.pdf'
plt.savefig(fname)

fig4, axarr4 = plt.subplots(1, 2, sharey=True, num=4)
fig4.clf()
fig4.suptitle(r'$Re_f=$%i' % (Re))

## Time step size #####################
axarr4[0] = fig4.add_subplot(1, 2, 1)
axarr4[0].plot(dT, Y, 'x')
axarr4[0].set_xlabel(r'$\frac{dt}{T}$')
axarr4[0].grid(True)

## lenght scale of turbulence  ########
axarr4[1] = fig4.add_subplot(1, 2, 2)
axarr4[1].plot(L, Y, 'x')
axarr4[1].set_xlabel(r'$\frac{l}{h}$')
axarr4[1].grid(True)
for label in axarr4[1].get_yticklabels():
    label.set_visible(False)

fig4.show()

fname = 'plots/profiles4.pdf'
plt.savefig(fname)
