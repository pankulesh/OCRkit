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
from timeit import default_timer as dt
import argparse
import matplotlib.pyplot as plt

nboot = 20

class LF_region:
    
    def __init__(self):
        self.density = np.zeros(ibound2)
        self.d_mean = np.zeros(ibound2)
        self.d_disp = np.zeros(ibound2)
        self.n_tot = 0
    
    def bootstrap(self):
        density_max = np.amax(self.density)
        app_densitytck = interpolate.splrep(argument, self.density)
        dens_mags = [np.zeros(ibound) for _ in range(nboot)]
        for k in range(nboot):
        # Density estimation for "bootstrap" set of radial distances values.
            for mag in (neiman_method(mag0, maglim+delta, density_max, app_densitytck) for _ in range(self.n_tot)):
                rstar = mag - mag0
                imin, imax = make_imin_imax(rstar, delta, step, ibound)
                if imin <= ibound and imax >= 0:
                    dens_mags[k][imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
        s1, s2 = (np.zeros(ibound2) for _ in range(2))

        for i in range(ibound):
            for k in range(nboot):
                s1[i] += dens_mags[k][i]
                s2[i] += dens_mags[k][i]**2

        self.d_mean = s1/nboot
        self.d_disp = np.sqrt((s2-nboot*self.d_mean**2)/(nboot-1))

def neiman_method(a, b, omega_max, tck):
    while True:
        x1 = a+(b-a)*np.random.random()
        x2 = omega_max*np.random.random()
        if x2 < interpolate.splev(x1, tck):
            return x1

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

def plot_lf(kind, argument, density, disp):
    plt.title("LF for {}".format(kind), fontsize=18)
    arg = argument
    d = density
    down = density-disp
    up = density+disp
    plt.plot(arg, d, 'k')
    plt.plot(arg, down, 'k', linestyle='dashed', linewidth=1)
    plt.plot(arg, up, 'k', linestyle='dashed', linewidth=1)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    ax.set_xlabel(r'Apparent {} magnitude, $^m$'.format(mag_type), fontsize=15)
    ax.set_ylabel(r'Luminocity function', fontsize=15)
    ax.set_xlim(mag0, maglim)
    ax.set_ylim(np.amin(down), np.amax(up))
    ax.tick_params(width=2)
    ax.grid(alpha=0.7, linestyle='dashed', linewidth=0.5)
    plt.savefig('{0}_{1}_lf_{2}.png'.format(name, catalog_type, kind), dpi=100, transparent=False)
    plt.close()

def write_and_plot_dens(kind, density, disp):
    with open('{0}_{1}_lf_{2}.txt'.format(name, catalog_type, kind), 'w') as f:
        for i in range(ibound):
            f.write('{0:7.3f}  {1:15.6f}  {2:15.6f}  {3:15.6f}\n'.format(argument[i], density[i], density[i]-disp[i], density[i]+disp[i]))
    plot_lf(kind, argument, density, disp)

parser = argparse.ArgumentParser(description='This program gives luminosity function by the list of stellar magnitudes.')
parser.add_argument("infile", type=str, help='input file name (coord.txt)')
parser.add_argument("parameters", type=str, help='file with parameters')
parser.add_argument("name", type=str, help="cluster's name")
parser.add_argument("-g", "--gaia", action='store_true', help='Input type flag, transformed 2MASS coord.txt by default')
args = parser.parse_args()
param = args.parameters
name = args.name
magxy = (3, 18, 19)
mag_type = 'J'
catalog_type = '2MASS'

if args.gaia:
    magxy = (-1, 0, 1)
    mag_type = 'G'
    catalog_type = 'Gaia'

with open(param, 'r') as f:
    words = f.readline().split()
    delta, step, mag0, maglim, x0, y0, rc, x1, y1 = [float(x) for x in words]
    
x2 = 2*x0-x1
y2 = 2*y0-y1

x3 = -y1+y0+x0
y3 = x1-x0+y0

x4 = 2*x0-x3
y4 = 2*y0-y3

xsys = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]

ibound = int((maglim-mag0)/step)+1
ibound2 = int((maglim-mag0+delta)/step)
# ibound - number of point corresponding boundary of LF plot
# ibound2 - number of point corresponding faintest stars taking part in LF estimation

# n_tot - number of stars used for LF estimate.
# We must know it for bootstrap density estimates.

lf_cluster = LF_region()
lf_ring = LF_region()
lf_circles = [LF_region() for _ in range(4)]
for i in range(4):
    lf_circles[i].xc, lf_circles[i].yc = xsys[i]

with open(args.infile, 'r') as f:
    for line in f:
        words = line.split()
        mag, x, y = [float(words[i]) for i in magxy]
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
            #Comparing field as ring
            elif rc < r < rc*np.sqrt(2):
                lf_ring.n_tot+=1
                lf_ring.density[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
            #Comparing fields as circles
            for j in lf_circles:
                r = np.sqrt((x-j.xc)**2+(y-j.yc)**2)
                if r < rc:
                    j.n_tot+=1
                    j.density[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]

print("Cluster region: {0} stars,\nDisc region: {1} stars,\nCircles: {2}, {3}, {4}, {5} stars,\nibound={1}".format(lf_cluster.n_tot, lf_ring.n_tot, lf_circles[0].n_tot, lf_circles[1].n_tot, lf_circles[2].n_tot, lf_circles[3].n_tot, ibound))
argument = np.zeros(ibound2)

# Star radial distances are distributed with distribution function distrib.
# For points with i>ibound we consider density(i)=density(ibound).

for i in range(ibound2):
    argument[i]=mag0+i*step
    if i >= ibound:
        for j in (lf_cluster, lf_ring, *lf_circles):
            j.density[i] = j.density[ibound-1]

begin0 = dt()
for j in (lf_cluster, lf_ring, *lf_circles):
    j.bootstrap()
print('{0:.1f} s for bootstrapping'.format(dt()-begin0))
densityb_2c, db_disp_2c, densityb_4c, db_disp_4c = [np.zeros(ibound2) for _ in range(4)]
densityd_d, dd_disp_d, densityd_2c, dd_disp_2c, densityd_4c, dd_disp_4c = [np.zeros(ibound2) for _ in range(6)]

#Averaging 2 circles and 4 circles regions
densityb_2c = 0.5*(lf_circles[0].density+lf_circles[1].density)
db_disp_2c = 0.5*np.sqrt(lf_circles[0].d_disp**2+lf_circles[1].d_disp**2)
densityb_4c = 0.25*(lf_circles[0].density+lf_circles[1].density+lf_circles[2].density+lf_circles[3].density)
db_disp_4c = 0.25*np.sqrt(lf_circles[0].d_disp**2+lf_circles[1].d_disp**2+lf_circles[2].d_disp**2+lf_circles[3].d_disp**2)

write_and_plot_dens('cluster_and_back', lf_cluster.density, lf_cluster.d_disp)
write_and_plot_dens('background_ring', lf_ring.density, lf_ring.d_disp)
write_and_plot_dens('background_2c', densityb_2c, db_disp_2c)
write_and_plot_dens('background_4c', densityb_4c, db_disp_4c)

#Counting LF of cluster and its confidence interval
densityd_d = lf_cluster.density - lf_ring.density
dd_disp_d = np.sqrt(lf_cluster.d_disp**2+lf_ring.d_disp**2)

densityd_2c = lf_cluster.density - densityb_2c
dd_disp_2c = np.sqrt(lf_cluster.d_disp**2+db_disp_2c**2)

densityd_4c = lf_cluster.density - densityb_4c
dd_disp_4c = np.sqrt(lf_cluster.d_disp**2+db_disp_4c**2)

write_and_plot_dens('cluster_ring', densityd_d, dd_disp_d)
write_and_plot_dens('cluster_2c', densityd_2c, dd_disp_2c)
write_and_plot_dens('cluster_4c', densityd_4c, dd_disp_4c)
