#!/usr/bin/python
# This program gives radial profile of surface density
# from the list of stellar rectangular coordinates X and Y.
# Confidence interval is determined by smoothed bootstrap algorithm.
# Case of file with rectangular coordinates obtained from downloads
# from 2MASS PSC via Aladin.
# Surface density is determined by kernel estimator.
# Hystogram with selectable step is also determined.

# density is an array of density values,
# hyst - density hystogram array,
# radial - argument array (distance from the cluster center),
# distrib - an array of distribution function for r,
# dens_mag - an array of density values for "bootstrap" estimation,
# s1 - an array of sum of "bootstrap" density values (for mean values
#      determination),
# s2 - an array of sum of squared "bootstrap" density values (for
#      dispersion determination),
# d_mean - an array of mean "bootstrap" density values,
# d_disp -  an array of dispersion of "bootstrap" density values.

import numpy as np
from numpy import pi, f2py
from scipy import interpolate
from multiprocessing import Pool, cpu_count
from datetime import datetime
import argparse
import matplotlib.pyplot as plt

source = """
        subroutine dens_estimator(i, step, rstar, delta, dens)
        integer i
        real*8 step, rstar, rstar2, delta, delta2, delta4, delta6
        real*8 pi, pi2, ri, ri2, bracket, bracket2
        real*8 arg, fimax, dens1, dens2, dens
cf2py   intent(out) dens
        pi = 3.14159265358979
        pi2 = pi*pi
        rstar2 = rstar*rstar
        delta2 = delta*delta
        delta4 = delta2*delta2
        delta6 = delta4*delta2
        ri = i*step
        ri2 = ri*ri
        bracket = 1.0 -(ri2+rstar2)/delta2
        bracket2 = bracket*bracket
        if (ri.lt.(delta-rstar)) then
            dens = 3.0*bracket2/pi/delta2+6.0*ri2*rstar2/pi/delta6
        else
            if (rstar.eq.0.) rstar = 0.01
            arg = (ri*ri+rstar*rstar-delta2)/(2.*ri*rstar)
            fimax=acos(arg)
            dens1 = 3.0*fimax*bracket2/pi2/delta2
            dens2 = 6.0*ri2*rstar2*fimax/pi2/delta6
            dens3 = 12.0*ri*rstar*bracket*sin(fimax)/pi2/delta4
            dens4 = 3.0*ri2*rstar2*sin(2.0*fimax)/pi2/delta6
            dens = dens1 + dens2 + dens3 + dens4
        endif
        return
        end
"""
try:
    import densest
except:
    f2py.compile(source.encode(), modulename='densest')
    import densest
dens_estimator = densest.dens_estimator
nboot = 50

def neiman_method(a, b, omega_max, tck):
    while True:
        x1 = a+(b-a)*np.random.random()
        x2 = (omega_max*1.01)*np.random.random()
        if x2 < interpolate.splev(x1, tck):
            return x1

def make_imin_imax(rstar, delta, step, ibound):
    if (rstar < delta):
        imin = int((rstar-delta)/step)+1
    else:
        imin = int((rstar-delta)/step)+2
    if imin < 0:
        imin = 0
    imax = int((rstar+delta)/step)+1
    if imax > ibound:
        imax = ibound
    return (imin, imax)

#def dens_estimator(i, step, rstar, delta):
#    pi2 = pi**2
#    rstar2 = rstar**2
#    delta2 = delta**2
#    delta4 = delta**4
#    delta6 = delta**6
#    ri = i*step
#    ri2 = ri**2
#    bracket = 1.0-(ri2+rstar2)/delta2
#    bracket2 = bracket**2
#    if ri < delta-rstar:
#        return 3.0*bracket2/pi/delta2+6.0*ri2*rstar2/pi/delta6
#    else:
#        if rstar == 0:
#            rstar = 0.01
#        fimax = np.arccos((ri2+rstar2-delta2)/(2*ri*rstar))
#        return 3*fimax*bracket2/pi2/delta2 + 6*ri2*rstar2*fimax/pi2/delta6 + 12*ri*rstar*bracket*np.sin(fimax)/pi2/delta4 + 3*ri2*rstar2*np.sin(2*fimax)/pi2/delta6


begin = datetime.now()
parser = argparse.ArgumentParser(description='Radial density profile estimator.')
parser.add_argument("-i", "--ifile", default='coord.txt', help='input file name, default coord.txt')
parser.add_argument("-o", "--ofile", default='', help='output file name prefix, no prefix by default')
parser.add_argument("-p", "--parameters", default='pparam.txt', help='file with parameters OR string with them, default pparam.txt')
parser.add_argument("-s", "--hystogram", action='store_true', help='make hystogram.txt')
parser.add_argument("-c", "--cores", default=eval("cpu_count()"), help='CPU cores amount for multiprocessing, counted as maximum available by default')
args = parser.parse_args()
cores = int(args.cores)
param = args.parameters

try:
    with open(param, 'r') as f:
        words = f.readline().split()
    delta, step, hstep, x0, y0, rmax, jmaglim = [float(x) for x in words]

except FileNotFoundError:
    delta, step, hstep, x0, y0, rmax, jmaglim = [float(x) for x in param.split()]

else:
    print("Wrong parameters!")
    exit(2)

ibound = int(rmax/step)+2
ibound2 = int((rmax + delta)/step)+2
hyst, d_mean, d_disp, s1, s2 = (np.zeros(ibound) for _ in range(5))
density = np.zeros(ibound2)
#ibound - number of point corresponding boundary of density plot
#ibound2 - number of point corresponding most remote stars (from the center) taking part in density estimation

rmin = 0.0
rmax = int(rmax/step)*step

# for the case when rmax is not an integer number of steps
# n_tot - number of stars used for density profile estimate.
# We must know it for bootstrap density estimates.

n_tot = 0

with open(args.ifile, 'r') as f:
    for line in f:
        words = line.split()
        jmag, x, y = [float(words[i]) for i in (3, 18, 19)]
        if jmag > jmaglim:
            continue
        rstar = np.sqrt((x-x0)**2+(y-y0)**2)
        k = int(rstar/hstep)
        #k - number of ring for hystogram
        hyst[k] += 1.0/(pi*hstep**2*(2*k+1))
        imin, imax = make_imin_imax(rstar, delta, step, ibound)
        if imin <= ibound and imax >= 0:
            n_tot+=1
        #The star contributes into density in points between imin and imax.
            for i in range(imin, imax):
                density[i] += dens_estimator(i, step, rstar, delta)
print("n_tot={0}, ibound={1}".format(n_tot, ibound))

radial = np.zeros(ibound2)
distrib = np.zeros(ibound2)

# Star radial distances are distributed with distribution function distrib.
# For points with i>ibound we consider density(i)=density(ibound).
for i in range(ibound2):
    radial[i]=i*step
    if i >= ibound:
        density[i] = density[ibound-1]
    #distrib[i] = 2*pi*radial[i]*density[i]
distrib = 2*pi*radial*density

# Distribution approximated by fortran subroutine
# determination of first derivatives of distribution function at boundary points.
#distrib1_1 = (distrib[1]-distrib[0])/(radial[1]-radial[0])
#distrib1_n = (distrib[ibound2-1]-distrib[ibound2-2])/(radial[ibound2-1]-radial[ibound2-2])
#splines.spline(radial, distrib, ibound2, distrib1_1, distrib1_n, distrib2)

# Distribution is approximated by SciPy cubic spline
distrib_max = np.amax(distrib)
app_distribtck = interpolate.splrep(radial, distrib)

# Cycle for "smoothed bootstrap" estimation of confidence interval
# for density.
# nboot is number of secondary density estimates.
# Program uses Neiman method for random points distributing.

def get_dens_mag(x):
    dens_mag = np.zeros(ibound)
    # Density estimation for "bootstrap" set of radial distances values.
    de = dens_estimator
    for rstar in (neiman_method(rmin, rmax+delta, distrib_max, app_distribtck) for _ in range(n_tot)):
        imin, imax = make_imin_imax(rstar, delta, step, ibound)
        if imin <= ibound and imax >= 0:
           dens_mag[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
    return dens_mag

with Pool(cores) as p:
    dens_mags = p.map(get_dens_mag, range(nboot))
print(datetime.now()-begin)

for i in range(ibound):
    for k in range(nboot):
        s1[i] += dens_mags[k][i]
        s2[i] += dens_mags[k][i]**2

# without multiprocessing
#for k in range(nboot):
#    print(k+1)
#    dens_mag = np.zeros(ibound)
# Density estimation for "bootstrap" set of radial distances values.
#    de = dens_estimator
#    for rstar in (neiman_method(rmin, rmax+delta, distrib_max, app_distribtck) for _ in range(n_tot)):
#        imin, imax = make_imin_imax(rstar, delta, step, ibound)
#        if imin <= ibound and imax >= 0:
#            dens_mag[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
#            for i in range(imin, imax):
#                dens_mag[i] += dens_estimator(i, step, rstar, delta)
#    for i in range(ibound):
#        s1[i] += dens_mag[i]
#        s2[i] += dens_mag[i]**2
#print(datetime.now()-begin)
# Cycle for "smoothed bootstrap" estimation of confidence interval ended 

for i in range(ibound):
    d_mean[i] = s1[i]/nboot
    d_disp[i] = np.sqrt((s2[i]-nboot*d_mean[i]**2)/(nboot-1))

#d_mean = s1/nboot
#d_disp = np.sqrt(np.abs(s2/(nboot-1)-d_mean**2*nboot/(nboot-1)))

with open(args.ofile+'density.txt', 'w') as f:
    for i in range(ibound):
        f.write('{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\n'.format(i*step, density[i], d_mean[i], density[i]-d_disp[i], density[i]+d_disp[i], hyst[int(i*step/hstep)]))
if args.hystogram:
    with open(args.ofile+'hystogram.txt', 'w') as f:
        for i in range(int(rmax/hstep)):
            f.write('{0:.3f}\t{1:.6f}\n'.format(i*hstep-hstep/2, hyst[i]))
