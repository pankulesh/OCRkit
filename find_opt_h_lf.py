#!/usr/bin/python
# This program gives luminosity function by the list of stellar magnitudes
# by the fixed kernel estimator.
# Case of initial field from data of 2MASS.
# Hystogram with selectable step is also determined.
# This code select stars in accordance with geometry conditions on x and y coordinates.
# density is an array of density values for cluster, densityb is an array of density values for
# hyst - density hystogram array for cluster, hystb - density hystogram array for comparing fie
# argument - argument array (distance from the beginning of LF),
# argum - an array of "bootstrap" argument values,
# distrib - an array of distribution function for r,
# distrib2 - an array of second derivative of estimated distribution function
#         (it needs for spline approximation of distribution function),
# dens_mag - an array of density values for "bootstrap" estimation,
# s1 - an array of sum of "bootstrap" density values (for mean values
#      determination),
# s2 - an array of sum of squared "bootstrap" density values (for
#      dispersion determination),
# d_mean - an array of mean "bootstrap" density values,
# d_disp -  an array of dispersion of "bootstrap" density values.
# x0 and y0 - coordinates of the cluster centre
# rc - radius of the cluster in arcminutes

import numpy as np
from scipy import interpolate
#from multiprocessing import Pool, cpu_count, freeze_support
from datetime import datetime
import argparse
import matplotlib.pyplot as plt

nboot = 20

class LF_region:
    
    def __init__(self):
        self.density = np.zeros(ibound2)
        self.d_mean = np.zeros(ibound)
        self.d_disp = np.zeros(ibound)
        self.n_tot = 0
    
    def bootstrap(self):
        density_max = np.amax(self.density)
        app_densitytck = interpolate.splrep(argument, self.density)
        dens_mags = [np.zeros(ibound) for _ in range(nboot)]
        begin = datetime.now()
        for k in nboot:
        # Density estimation for "bootstrap" set of radial distances values.
            begin1 = datetime.now()-begin
            for mag in (neiman_method(mag0, maglim+delta-mag0, density_max, app_densitytck) for _ in range(self.n_tot)):
                rstar = mag - mag0
                imin, imax = make_imin_imax(rstar, delta, step, ibound)
                if imin <= ibound and imax >= 0:
                    dens_mags[k][imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
            print('{0} for {1} cycle'.format(datetime.now()-begin1, k))
        print('{0} in total'.format(datetime.now()-begin))
        s1, s2 = (np.zeros(ibound) for _ in range(2))

        print(datetime.now()-begin)

        for i in range(ibound):
            for k in range(nboot):
                s1[i] += dens_mags[k][i]
                s2[i] += dens_mags[k][i]**2

        self.d_mean = s1/nboot
        self.d_disp = np.sqrt((s2-nboot*d_mean**2)/(nboot-1))

def make_imin_imax(rstar, delta, step, ibound):
    if (rstar < delta):
        imin = 0
    else:
        imin = int((rstar-delta)/step)
    if imin < 0:
        imin = 0
    imax = int((rstar+delta)/step)
    if imax > ibound:
        imax = ibound
    return (imin, imax)

#Density estimator function
def de(i, step, rstar, delta):
    ri = i*step
    bracket = 1 - (rstar-ri)**2/(delta**2)
    return (3*bracket/4/delta)

def plot_lf(lf):
    plt.plot(lf.argument, lf.density, label='{}'.format(lf.delta))
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    ax.set_xlabel(r'G magnitude, $^m$', fontsize=26)
    ax.set_ylabel(r'Luminocity function', fontsize=26)
    ax.tick_params(width=2)

infile = 'asu_coord.txt'
step, mag0, maglim, x0, y0, rc = [0.01, 8, 17, -2.3, 0.4, 15]
deltas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]    
lfs = []
for delta in deltas:
    ibound = int((maglim-mag0)/step)
    ibound2 = int((maglim-mag0+delta)/step)
    lfs.append(LF_region())
    lfs[-1].delta = delta
    lfs[-1].argument = np.array([mag0+i*step for i in range(ibound2)])
    lf_cluster = lfs[-1]
    with open(infile, 'r') as f:
        for line in f:
            words = line.split()
            mag, x, y = [float(words[i]) for i in (-1, 0, 1)]
            if mag > maglim+delta or mag < mag0-delta:
                continue
            rstar = mag-mag0
            imin, imax = make_imin_imax(rstar, delta, step, ibound)
            if imin <= ibound and imax >= 0:
                #Cluster field
                r = np.sqrt((x-x0)**2+(y-y0)**2)
                if r < rc:
                    lf_cluster.n_tot+=1
                    lf_cluster.density[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]

    for i in range(ibound2):
        if i >= ibound:
            lf_cluster.density[i] = lf_cluster.density[ibound-1]
    print('{} ready'.format(delta))

plt.title("LFs with different kernel halfwidth", fontsize=26)
for lf in lfs:
    plot_lf(lf)
plt.legend()
plt.show()
