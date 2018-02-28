#!/usr/bin/python

from os import listdir
import re
import argparse
import numpy as np
from scipy.integrate import simps
from scipy.interpolate import splrep, splev
from matplotlib import pyplot as plt
#regexp = r'density_{0}_[0-9]+_{1}.txt'.format(name, h)
#files = [f for f in os.listdir('.') if re.match(regexp, f)]
#mags = [int(mag.split('_')[2]) for mag in flies]
#mags.sort()
PRECISION = 0.01
parser = argparse.ArgumentParser(description='Getting radius of cluster and mean background stellar density from radial density profile')
parser.add_argument("ifile", help='input file name with columns: r, density, mean density, density-sigma, density+sigma, (density.txt)')
args = parser.parse_args()
#a = np.loadtxt('density_{0}_{1}_{2}.txt'.format(NAME, jmaglim, h), unpack=True, usecols=(0,1,3,4))
r, dens, densms, densps = np.loadtxt(args.ifile, unpack=True, usecols = (0, 1, 3, 4))

rmax = np.amax(r)
densmax = int(np.amax(densps))+1
density_spline = splrep(r, dens)
density_msigma_spline = splrep(r, densms)
density_psigma_spline = splrep(r, densps)

radius = np.arange(0, rmax+0.01, 0.01)

density = splev(radius, density_spline)
density_low = splev(radius, density_msigma_spline)
density_up = splev(radius, density_psigma_spline)

rect_deviation = np.zeros_like(radius)

for i in range(len(radius)):
    rect = density[i]*(rmax-radius[i])
    integral = simps(density[i:], radius[i:])
    rect_deviation[i] = abs(rect-integral)
    
#    plt.title("{0}, h={1}', Jlim={2}m".format(NAME, h, jmaglim), fontsize=26)
plt.title("{}".format(args.ifile), fontsize=26)
plt.plot(radius, density, 'k', linewidth=1)
plt.plot(radius, density_low, 'k', linestyle='dashed', linewidth=1)
plt.plot(radius, density_up, 'k', linestyle='dashed', linewidth=1)

for i in range(len(radius)-50):
    if rect_deviation[i] < 0.01:
        R, Fb = radius[i], density[i]
        dFb = max(density_up[i]-Fb, Fb-density_low[i])
        ibound = i
        break

for r in radius[:ibound:-1]:
    rho_l = abs(splev(r, density_msigma_spline)-Fb)
    if rho_l < PRECISION:
        dR1 = R - r
        break
    rho_u = abs(splev(r, density_psigma_spline)-Fb)
    if rho_u < PRECISION:
        dR1 = R - r
        break

for r in radius[ibound:]:
    rho_l = abs(splev(r, density_msigma_spline)-Fb)
    if rho_l < PRECISION:
        dR2 = r - R
        break
    rho_u = abs(splev(r, density_psigma_spline)-Fb)
    if rho_u < PRECISION:
        dR2 = r - R
        break

dR = max(dR1, dR2)
if dR > R:
    dR = R
if dFb > Fb:
    dFb = Fb
#x1, y1 = [0, rmax], [Fb, Fb]
#x2, y2 = [R, R], [0, densmax]
print("R = {0:5.2f}, dR = {1:5.2f}, Fb = {2:5.2f}, dFb = {3:5.2f}".format(R, dR, Fb, dFb))
#plt.plot(x1, y1, x2, y2, color='k', linewidth=1)
ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(26)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(26)
ax.axhline(Fb, c='k', alpha=0.5, lw=0.5, label = "Fb = {0:5.2f}±{1:5.2f}".format(Fb, dFb))
ax.axvline(R, c='k', alpha=0.5, lw=0.5, label = "R = {0:5.2f}±{1:5.2f}".format(R, dFb))
ax.tick_params(width=2)
ax.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
ax.set_xlim(0, rmax)
ax.set_ylim(0, densmax)
ax.legend()
ax.set_xlabel(r'distance from center, $arcmin$', fontsize=26)
ax.set_ylabel(r'radial density, $1/arcmin^2$', fontsize=26)
plt.show()
plt.savefig('profile_parameters.png', dpi=300)
