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
rac = 169.365
decc = -62.740
centers = {
    11: [-5, -3, 3, 5, -5, -3, 3, 5],
    12: [-5, -3, 3, 5, -5, -3, 3, 5],
    13: [-5, -3, 3, 5, -5, -3, 3, 5],
    14: [-4, -2, 2, 4, -5, -4, -1, 0],
    15: [5, 6, 9, 10, 0, 1, 4, 5],
    16: [3, 4, 6, 7, 0, 1, 4, 5]
}
h = 2
rmax = 30
hstep = h*2
coord_in = 'coord.txt'
jmaglims = [16, 15, 14, 13, 12, 11]
class LinDens:
    def __init__(self, d):
        self.outer_min, self.inner_min, self.inner_max, self.outer_max = [float(x) for x in d]

    def count_ibdens(self):
        self.ibound = int((self.outer_max-self.outer_min)/self.step)+1
        self.dens = np.zeros((2, self.ibound))

    def count_densstep(self):
        for i in range(self.ibound):
            self.dens[0][i] = self.outer_min + i*self.step

    def interval(self):
        if self.r < self.delta:
            self.imin = 0
        else:
            self.imin = int((self.r - self.delta)/self.step)+1
        if self.imin > self.ibound:
            raise ValueError('Out of bounds')
        self.imax = int((self.r+self.delta)/self.step)+1
        if self.imax > self.ibound:
            self.imax = self.ibound

    def calculate_dens(self):
        for i in range(self.imin, self.imax):
            ri = i*self.step
            bracket = 1.0 - (self.r - ri)**2/self.delta**2
            if bracket < 0.0:
                bracket = 0.0
            self.dens[1][i] += 15.0*bracket**2/16.0/self.delta

    def is_between(self):
        return self.inner_min <= self.dot <= self.inner_max


def makecenterparameters(alpha, delta, m, h, s):
    k = 1/60/np.cos(np.pi*delta/180)
#   x0, x1 = -1*x3, -1*x2
#   y0, y1, y2, y3 = x0, x1, x2, x3
    x0, x1, x2, x3, y0, y1, y2, y3 = centers[m]
    with open('centerparameters_{0}_{1}_{2}.txt'.format(NAME, m, h), 'w') as f:
        f.write('{:.1f}\n'.format(m))
        f.write('{0:.1f}\t{1:.4f}\t{2:.4f}\n'.format(h, h/60, k*h))
        f.write('{0:.2f}\t{1:.6f}\t{2:.6f}\n'.format(s, s/60, k*s))
        f.write('{0:.1f}\t{1:.1f}\t{2:.1f}\t{3:.1f}\n'.format(x0, x1, x2, x3))
        f.write('{0:.1f}\t{1:.1f}\t{2:.1f}\t{3:.1f}\n'.format(y0, y1, y2, y3))
        f.write('{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n'.format(delta+x0/60, delta+x1/60, delta+x2/60, delta+x3/60))
        f.write('{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n'.format(alpha+k*y0, alpha+k*y1, alpha+k*y2, alpha+k*y3))
    params = (m, (h, h/60, k*h), (s, s/60, k*s), (x0, x1, x2, x3), (y0, y1, y2, y3), (delta+x0/60, delta+x1/60, delta+x2/60, delta+x3/60), (alpha+k*y0, alpha+k*y1, alpha+k*y2, alpha+k*y3))
    return params

def findcentres(coord_in, params):
    jlim = float(params[0])
    xld, yld, decld, rald = (LinDens(s) for s in params[3:7])
    xld.title, yld.title, decld.title, rald.title = ['X', 'Y', 'DEC', 'RA']
    xld.delta, decld.delta, rald.delta = [float(x) for x in params[1]]
    xld.step, decld.step, rald.step = [float(x) for x in params[2]]
    yld.delta, yld.step = xld.delta, xld.step
    
    for ld in (xld, yld, rald, decld):
        ld.count_ibdens()
        ld.count_densstep()
    
    with open(coord_in, 'r') as f:
        for line in f:
            words = line.split()
            jmag, rald.dot, decld.dot, xld.dot, yld.dot = [float(words[i]) for i in (3, 0, 1, 18, 19)]
            if jmag > jlim:
                continue
            try:
                for ld in (xld, yld, rald, decld):
                    ld.r = ld.dot - ld.outer_min
                    ld.interval()
            except ValueError:
                continue
            for ld1, ld2 in ((yld, xld), (rald, decld)):
                if ld1.is_between():
                    ld2.calculate_dens()
                if ld2.is_between():
                    ld1.calculate_dens()
    
    for ld in (xld, yld, rald, decld):
        ld.maxim = ld.dens[0][np.argmax(ld.dens[1])]

    with open("centres_{0}_{1}_{2}.txt".format(NAME, jlim, params[1][0]), 'w') as f:
        f.write('{0:.3f}\t{1:.4f}\t#x and y\n'.format(xld.maxim, yld.maxim))
        f.write('{0:.3f}\t{1:.4f}\t#alpha and delta\n'.format(rald.maxim, decld.maxim))
    
    fig = plt.figure(figsize=(12, 8))
    z = 1
    
    for ld in (xld, yld, decld, rald):
        ld.plot = plt.subplot(2, 2, z)
        plt.plot(ld.dens[0], ld.dens[1], 'k')
        ax = plt.gca()
        for xtick in ax.get_xticklabels():
            xtick.set_fontsize(10)
        for ytick in ax.get_yticklabels():
            ytick.set_fontsize(10)
        ld.plot.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.title('{} linear density'.format(ld.title))
        z+=1
    fig.tight_layout()
    fig.savefig('centerprofile_{0}_{1}_{2}.png'.format(NAME, jlim, params[1][0]), dpi=100, frameon=True)
    plt.close()
    maxims = (xld.maxim, yld.maxim)
    return maxims
    
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

for jmaglim in jmaglims:
    params = makecenterparameters(rac, decc, jmaglim, h, step)
    x0, y0 = findcentres(coord_in, params)
    delta = h
    ibound = int(rmax/step)+1
    ibound2 = int((rmax + delta)/step)+1
    hyst, d_mean, d_disp, s1, s2 = (np.zeros(ibound) for _ in range(5))
    density = np.zeros(ibound2)
    rmin = 0.0
    rmax = int(rmax/step)*step
    n_tot = 0

    with open(coord_in, 'r') as f:
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
#            for i in range(imin, imax):
#                dens_mag[i] += dens_estimator(i, step, rstar, delta)
       for i in range(ibound):
           s1[i] += dens_mag[i]
           s2[i] += dens_mag[i]**2

#    def get_dens_mag(x):
#        dens_mag = np.zeros(ibound)
#        de = dens_estimator
#        for rstar in (neiman_method(rmin, rmax+delta, distrib_max, app_distribtck) for _ in range(n_tot)):
#            imin, imax = make_imin_imax(rstar, delta, step, ibound)
#        if imin <= ibound and imax >= 0:
#           dens_mag[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
#        return dens_mag
#
#    with Pool(cores) as p:
#        dens_mags = p.map(get_dens_mag, range(nboot))
#
#    for i in range(ibound):
#        for k in range(nboot):
#            s1[i] += dens_mags[k][i]
#            s2[i] += dens_mags[k][i]**2
    
    for i in range(ibound):
        d_mean[i] = s1[i]/nboot
        d_disp[i] = np.sqrt((s2[i]-nboot*d_mean[i]**2)/(nboot-1))

    with open('density_{0}_{1}_{2}.txt'.format(NAME, jmaglim, params[1][0]), 'w') as f:
        for i in range(ibound):
            f.write('{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\t{5:.6f}\n'.format(i*step, density[i], d_mean[i], density[i]-d_disp[i], density[i]+d_disp[i], hyst[int(i*step/hstep)]))
    print(str(jmaglim)+'_'+str(h)+' is done.\n')
    a = np.loadtxt('density_{0}_{1}_{2}.txt'.format(NAME, jmaglim, h), unpack=True, usecols=(0,1,3,4))
    r = a[0]
    fig = plt.figure(figsize=(10,10))
    plt.title("{0}, h={1}', Jlim={2}m".format(NAME, h, jmaglim), fontsize=26)
    plt.plot(r, a[1], 'k', linewidth=1)
    plt.plot(r, a[2], 'k', linestyle='dashed', linewidth=1)
    plt.plot(r, a[3], 'k', linestyle='dashed', linewidth=1)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    ax.set_xlabel(r'distance from center, $arcmin$', fontsize=26)
    ax.set_ylabel(r'radial density, $1/arcmin^2$', fontsize=26)
    ax.set_xlim(0, np.amax(r))
    ax.set_ylim(0, np.amax(a[3]))
    ax.tick_params(width=2)
    plt.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
    fig.savefig('density_{0}_{1}_{2}.png'.format(NAME, jmaglim, h), dpi=100, transparent=False)
    plt.close()



