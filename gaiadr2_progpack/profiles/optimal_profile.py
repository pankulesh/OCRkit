#!/usr/bin/python
from kde_functions import neiman_method, make_imin_imax, dens_estimator 

import argparse
import matplotlib.pyplot as plt
from scipy import interpolate
from multiprocessing import Pool, cpu_count
import numpy as np

parser = argparse.ArgumentParser(description='Программа строит профили радиальной плотности с разными параметрами KDE для определения оптимального.')
parser.add_argument("infile", help='Имя входного файла')
parser.add_argument("name", help="Имя скопления для заголовков")
parser.add_argument("-l", "--maglim", default=18, help="Предельная звёздная величина, по-умолчанию 18")
parser.add_argument("-s", "--step", default=0.05, help="Шаг для построения в угловых минутах, по-умолчанию 0.05'")
parser.add_argument("-r", "--rmax", default=10, help="Максимальное расстояние от центра, по-умолчанию 10'")
args = parser.parse_args()
infile = args.infile
name = args.name
step = args.step
maglim = args.maglim
rmax = args.rmax
cores = cpu_count()
nboot = 1

hs = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2, 3, 4, 5]

xymag = (0, 1, 13)
magtype = 'G'

stars = np.loadtxt(infile, usecols=xymag, comments='#')

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
    nboot = int(np.sqrt(n_tot))
    radial = np.zeros(ibound2)
    distrib = np.zeros(ibound2)

    for i in range(ibound2):
        radial[i]=i*step
        if i >= ibound:
            density[i] = density[ibound-1]
    distrib = 2*np.pi*radial*density

    distrib_max = np.amax(distrib)
    app_distribtck = interpolate.splrep(radial, distrib)
    # with multiprocessing
    
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

    with open('density_{0}_{1}_{2}.txt'.format(name, maglim, h), 'w') as f:
        for i in range(ibound):
            f.write('{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\t{4:.6f}\n'.format(i*step, density[i], d_mean[i], density[i]-d_disp[i], density[i]+d_disp[i]))
    print(str(maglim)+'_'+str(h)+' is done.\n')

fig = plt.figure(figsize=(10,10))
for h in hs:
    a = np.loadtxt('density_{0}_{1}_{2}.txt'.format(name, maglim, h), unpack=True, usecols=(0,1))
    r = a[0]
    plt.plot(r, a[1], linewidth=1, label = "h={}'".format(h))
    
ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(26)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(26)
ax.set_title("Выбор оптимальной полуширины ядра для {}".format(name))
ax.set_xlabel(r'Угловое расстояние от центра, $arcmin$')
ax.set_ylabel(r'Радиальная плотность, $arcmin^{-2}$')
ax.set_xlim(0, 10)
ax.set_ylim(0, 30)
ax.tick_params(width=2)
ax.legend()
plt.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
fig.savefig('densities_{0}_{1}.png'.format(name, maglim), dpi=100, transparent=False)
plt.show()
plt.close()
