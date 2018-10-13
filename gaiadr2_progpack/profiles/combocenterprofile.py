#!/usr/bin/python
from kde_functions import neiman_method, make_imin_imax, dens_estimator

import argparse
import matplotlib.pyplot as plt
from scipy import interpolate
from multiprocessing import Pool, cpu_count
import numpy as np

parser = argparse.ArgumentParser(description='Программа ищет центры скопления и считает радиальные профили плотности для каждой предельной звёздной величины.')
parser.add_argument("infile", help='Имя входного файла')
parser.add_argument("name", help="Имя скопления для заголовков графиков")
parser.add_argument("delta", type=float, help="Полуширина ядра в угловых минутах")
parser.add_argument("-l", "--maglim", default=18, help="Предельная звёздная величина G, по-умолчанию 18")
parser.add_argument("-s", "--step", default=0.05, help="Шаг для построения в угловых минутах, по-умолчанию 0.05'")
parser.add_argument("-r", "--rmax", default=30, help="Максимальное расстояние от центра в угловых минутах, по-умолчанию 30'")
args = parser.parse_args()
infile = args.infile
name = args.name
step = args.step
maglim = args.maglim
rmax = args.rmax
delta = args.delta
cores = cpu_count()

xymag = (0, 1, 13)
magtype = 'G'

nboot = 100
#Этот массив можно и нужно менять, если центры определяются плохо
centers = {
    11: [-5, -3, 3, 5, -5, -3, 3, 5],
    12: [-5, -3, 3, 5, -5, -3, 3, 5],
    13: [-5, -3, 3, 5, -5, -3, 3, 5],
    14: [-5, -3, 3, 5, -5, -3, 3, 5],
    15: [-5, -3, 3, 5, -5, -3, 3, 5],
    16: [-5, -3, 3, 5, -5, -3, 3, 5],
    17: [-5, -3, 3, 5, -5, -3, 3, 5],
    18: [-5, -3, 3, 5, -5, -3, 3, 5]
}
maglims = list(centers.keys())

class LinDens:
    def __init__(self, d):
        self.outer_min, self.inner_min, self.inner_max, self.outer_max = [float(x) for x in d]

    def count_ibdens(self):
        self.ibound = int((self.outer_max-self.outer_min)/self.step)
        self.dens = np.zeros((2, self.ibound))

    def count_densstep(self):
        for i in range(self.ibound):
            self.dens[0][i] = self.outer_min + i*self.step

    def interval(self):
        if self.r < self.delta:
            self.imin = 0
        else:
            self.imin = int((self.r - self.delta)/self.step)
        if self.imin > self.ibound:
            raise ValueError('Out of bounds')
        self.imax = int((self.r+self.delta)/self.step)
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

def findcentres(inarr, maglim, delta, step):
    xld, yld = LinDens(centers[maglim][:4]), LinDens(centers[maglim][4:])
    xld.title, yld.title = ['X', 'Y']
    xld.delta = delta 
    xld.step = step
    yld.delta, yld.step = xld.delta, xld.step
    
    for ld in (xld, yld):
        ld.count_ibdens()
        ld.count_densstep()
    
    for xld.dot, yld.dot, mag in inarr:
            if mag > maglim:
                continue
            try:
                for ld in (xld, yld):
                    ld.r = ld.dot - ld.outer_min
                    ld.interval()
            except ValueError:
                continue
            if yld.is_between():
                xld.calculate_dens()
            if xld.is_between():
                yld.calculate_dens()
    
    for ld in (xld, yld):
        ld.maxim = ld.dens[0][np.argmax(ld.dens[1])]

    with open("centres_{0}_{1}_{2}.txt".format(name, maglim, delta), 'w') as f:
        f.write('{0:.3f}\t{1:.4f}\t#x and y\n'.format(xld.maxim, yld.maxim))
    
    fig, [xld.plot, yld.plot] = plt.subplots(2, sharex=True)
    for ld in (xld, yld):
        ld.plot.plot(ld.dens[0], ld.dens[1], 'k')
        for xtick in ld.plot.get_xticklabels():
            xtick.set_fontsize(10)
        for ytick in ld.plot.get_yticklabels():
            ytick.set_fontsize(10)
        ld.plot.get_xaxis().get_major_formatter().set_useOffset(False)
        ld.plot.set_title('Профиль линейной плотности по {} для скопления {}'.format(ld.title, name))
        ld.plot.set_xlabel('{}, arcmin'.format(ld.title))
        ld.plot.set_ylabel('Линейная плотность, $arcmin^{-1}$')
    fig.tight_layout()
    fig.savefig('centerprofile_{0}_{1}_{2}.png'.format(name, maglim, delta), dpi=100, frameon=True)
    plt.close()
    maxims = (xld.maxim, yld.maxim)
    return maxims


stars = np.loadtxt(infile, usecols=xymag, comments='#')

for maglim in maglims:
    x0, y0 = findcentres(stars, maglim, delta, step)
    ibound = int(rmax/step)
    ibound2 = int((rmax + delta)/step)
    d_mean, d_disp, s1, s2 = (np.zeros(ibound) for _ in range(4))
    density = np.zeros(ibound2)
    rmin = 0.0
    rmax = int(rmax/step)*step
    n_tot = 0

    for x, y, mag in stars:
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
    distrib = 2*np.pi*radial*density

    distrib_max = np.amax(distrib)
    app_distribtck = interpolate.splrep(radial, distrib)

    #with multiprocessing
    
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

    for i in range(ibound):
            for k in range(nboot):
                s1[i] += dens_mags[k][i]
                s2[i] += dens_mags[k][i]**2

#without multiprocessing
#    for k in range(nboot):
#       print(k+1)
#       dens_mag = np.zeros(ibound)
#       # Density estimation for "bootstrap" set of radial distances values.
#       de = dens_estimator
#       for rstar in (neiman_method(rmin, rmax+delta, distrib_max, app_distribtck) for _ in range(n_tot)):
#           imin, imax = make_imin_imax(rstar, delta, step, ibound)
#           if imin <= ibound and imax >= 0:
#               dens_mag[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
#       for i in range(ibound):
#           s1[i] += dens_mag[i]
#           s2[i] += dens_mag[i]**2 

    for i in range(ibound):
        d_mean[i] = s1[i]/nboot
        d_disp[i] = np.sqrt((s2[i]-nboot*d_mean[i]**2)/(nboot-1))

    with open('density_{0}_{1}_{2}.txt'.format(name, maglim, delta), 'w') as f:
        for i in range(ibound):
            f.write('{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\n'.format(i*step, density[i], d_mean[i], density[i]-d_disp[i], density[i]+d_disp[i]))
    print(str(maglim)+'_'+str(delta)+' is done.\n')
    a = np.loadtxt('density_{0}_{1}_{2}.txt'.format(name, maglim, delta), unpack=True, usecols=(0,1,3,4))
    r = a[0]
    fig = plt.figure(figsize=(10,10))
    plt.title("Радиальный профиль плотности {0},\nh={1}', {3}lim={2}m".format(name, delta, maglim, magtype), fontsize=18)
    plt.plot(r, a[1], 'k', linewidth=1)
    plt.plot(r, a[2], 'k', linestyle='dashed', linewidth=1)
    plt.plot(r, a[3], 'k', linestyle='dashed', linewidth=1)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax.set_xlabel(r'Угловое расстояние от центра, $arcmin$', fontsize=16)
    ax.set_ylabel(r'Радиальная плотность, $arcmin^{-2}$', fontsize=16)
    ax.set_xlim(0, np.amax(r))
    ax.set_ylim(0, np.amax(a[3]))
    ax.tick_params(width=2)
    plt.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
    fig.savefig('density_{0}_{1}_{2}.png'.format(name, maglim, delta), dpi=100, transparent=False)
    plt.close()
