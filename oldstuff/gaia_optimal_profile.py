#!/usr/bin/python

import matplotlib.pyplot as plt
from numpy import pi, f2py
from scipy import interpolate
from multiprocessing import Pool, cpu_count
import numpy as np
import densest
dens_estimator = densest.dens_estimator
nboot = 20
cores = 8
NAME = "IC2714"
step = 0.05
maglim = 17
hs = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
rmax = 30
coord_in = 'asu_coord.txt'
#jmaglims = [16, 15, 14, 13, 12, 11]

def neiman_method(a, b, omega_max, tck):
    while True:
        x1 = a+(b-a)*np.random.random()
        x2 = omega_max*np.random.random()
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

for h in hs:
    x0, y0 = 0, 0
    delta = h
    ibound = int(rmax/step)
    ibound2 = int((rmax + delta)/step)
    d_mean, d_disp, s1, s2 = (np.zeros(ibound) for _ in range(4))
    density = np.zeros(ibound2)
    rmin = 0.0
    rmax = int(rmax/step)*step
    n_tot = 0

    with open(coord_in, 'r') as f:
        for line in f:
            words = line.split()
            x = float(words[0])
            y = float(words[1])
            mag = float(words[-1])
            if mag > maglim:
                continue
            rstar = np.sqrt((x-x0)**2+(y-y0)**2)
            imin, imax = make_imin_imax(rstar, delta, step, ibound)
            if imin <= ibound and imax >= 0:
                n_tot+=1
            #The star contributes into density in points between imin and imax.
                for i in range(imin, imax):
                    density[i] += dens_estimator(i, step, rstar, delta)
    print("n_tot={0}, ibound={1}".format(n_tot, ibound))

    radial = np.zeros(ibound2)
    distrib = np.zeros(ibound2)

    for i in range(ibound2):
        radial[i]=i*step
        if i >= ibound:
            density[i] = density[ibound-1]
    distrib = 2*pi*radial*density

    distrib_max = np.amax(distrib)
    app_distribtck = interpolate.splrep(radial, distrib)
    # without multiprocessing
    for k in range(nboot):
       print(k+1)
       dens_mag = np.zeros(ibound)
# Density estimation for "bootstrap" set of radial distances values.
       de = dens_estimator
       for rstar in (neiman_method(rmin, rmax+delta, distrib_max, app_distribtck) for _ in range(n_tot)):
           imin, imax = make_imin_imax(rstar, delta, step, ibound)
           if imin <= ibound and imax >= 0:
               dens_mag[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
       for i in range(ibound):
           s1[i] += dens_mag[i]
           s2[i] += dens_mag[i]**2
    
    for i in range(ibound):
        d_mean[i] = s1[i]/nboot
        d_disp[i] = np.sqrt((s2[i]-nboot*d_mean[i]**2)/(nboot-1))

    with open('density_{0}_{1}_{2}.txt'.format(NAME, maglim, h), 'w') as f:
        for i in range(ibound):
            f.write('{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\n'.format(i*step, density[i], d_mean[i], density[i]-d_disp[i], density[i]+d_disp[i]))
    print(str(maglim)+'_'+str(h)+' is done.\n')
fig = plt.figure(figsize=(10,10))
for h in hs:
    a = np.loadtxt('density_{0}_{1}_{2}.txt'.format(NAME, maglim, h), unpack=True, usecols=(0,1))
    r = a[0]
    plt.plot(r, a[1], linewidth=1, label = "h={}'".format(h))
    
ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(26)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(26)
ax.set_xlabel(r'distance from center, $arcmin$', fontsize=26)
ax.set_ylabel(r'radial density, $1/arcmin^2$', fontsize=26)
ax.set_xlim(0, 10)
ax.set_ylim(0, 20)
ax.tick_params(width=2)
ax.legend()
plt.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
fig.savefig('densities_{0}_{1}.png'.format(NAME, maglim), dpi=100, transparent=False)
plt.show()
plt.close()



