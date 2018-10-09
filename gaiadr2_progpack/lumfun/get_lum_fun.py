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
from kde_functions import neiman_method, make_imin_imax, density_estimator_1D, Region 
import argparse
import matplotlib.pyplot as plt

de = density_estimator_1D

def plot_lf(kind, argument, density, disp):
    plt.title("Функция блеска скопления {}\nдля области сравнения {}".format(clustername, kind), fontsize=12)
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
    ax.set_xlabel(r'Видимая звёздная величина в полосе {}'.format(mag_type), fontsize=15)
    ax.set_ylabel(r'Функция блеска', fontsize=15)
    ax.set_xlim(mag0, maglim)
    ax.invert_xaxis()
    ax.set_ylim(np.amin(down), np.amax(up))
    ax.tick_params(width=2)
    ax.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
    plt.savefig('lf_{}_{}.png'.format(clustername, kind), dpi=100, transparent=False)
    plt.close()

parser = argparse.ArgumentParser(description='Программа, рассчитывающая функции блеска для области скопления с помощью метода KDE')
parser.add_argument("infile", type=str, help='Имя входного файла')
parser.add_argument("clustername", type=str, help='Название скопления')
parser.add_argument("parameters", type=str, help='Имя файла с параметрами в строку: step, mag0, maglim, x0, y0, rc, delta, [x1, y1]')
parser.add_argument('-c', '--circles', action='store_true', help='Укажите этот флаг, если области сравнения - круги вокруг скопления с центром одного из них x1 y1, центры остальных находятся симметрично центру скопления, радиус областей равен rc')

args = parser.parse_args()
clustername = args.clustername
parameters = np.loadtxt(args.parameters, comments='#')

magxy = (13, 0, 1)
mag_type = 'G'
step, mag0, maglim, x0, y0, rc, delta= parameters[:7]

ibound = int((maglim-mag0)/step)+1
ibound2 = int((maglim-mag0+delta)/step)
# ibound - number of point corresponding boundary of LF plot
# ibound2 - number of point corresponding faintest stars taking part in LF estimation

if args.circles:
    x1, y1 = parameters[7:9]
    dist = np.sqrt((x0-x1)**2+(y0-y1)**2)
    defect = 2*rc - dist
    if defect >= 0:
        x1new = x1 + defect/dist*(x1-x0)
        y1new = y1 + defect/dist*(y1-y0)
        x1, y1 = x1new, y1new
        with open(args.parameters, 'w') as f:
            f.write('#step mag0 maglim x0 y0 rc delta x1 y1\n')
            f.write('{:.2f} {:.1f} {:.1f} {:.2f} {:.2f} {:.2f} {:.1f} {:.2f} {:.2f}\n'.format(step, mag0, maglim, x0, y0, rc, delta, x1, y1))

    
    """
         x4y4
    x1y1 x0y0 x2y2
         x3y3

    """

    x2 = 2*x0-x1
    y2 = 2*y0-y1

    x3 = -y1+y0+x0
    y3 = x1-x0+y0

    x4 = 2*x0-x3
    y4 = 2*y0-y3

    xsys = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]
    
    lf_circles = [Region(ibound2) for _ in range(4)]
    for i in range(4):
        lf_circles[i].xc, lf_circles[i].yc = xsys[i]

lf_cluster = Region(ibound2)
lf_ring = Region(ibound2)

input_stars = np.loadtxt(args.infile, comments='#', usecols=magxy)
for mag, x, y in input_stars:
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
        #Ring comparing field 
        elif rc < r < rc*np.sqrt(2):
            lf_ring.n_tot+=1
            lf_ring.density[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
        #Circles comparing fields
        if args.circles:
            for j in lf_circles:
                r = np.sqrt((x-j.xc)**2+(y-j.yc)**2)
                if r < rc:
                    j.n_tot+=1
                    j.density[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]

print("Область скопления: {0} звёзд,\nКольцевая область сравнения: {1} звёзд,\nibound={2}".format(lf_cluster.n_tot, lf_ring.n_tot, ibound))
if args.circles:
    print("Круговые области сравнения: {}, {}, {}, {} звёзд соответственно".format(lf_circles[0].n_tot, lf_circles[1].n_tot, lf_circles[2].n_tot, lf_circles[3].n_tot))
argument = np.zeros(ibound2)

# Star radial distances are distributed with distribution function distrib.
# For points with i>ibound we consider density(i)=density(ibound).

a = [lf_cluster, lf_ring]

if args.circles:
    a += lf_circles

for i in range(ibound2):
    argument[i]=mag0+i*step
    if i >= ibound:
        for j in a:
            j.density[i] = j.density[ibound-1]
for i in a:
    i.bootstrap(argument, ibound, ibound2, mag0, maglim, delta, step)

densityd_d, dd_disp_d = [np.zeros(ibound2) for _ in range(2)]

plot_lf('cluster+back', argument, lf_cluster.density, lf_cluster.d_disp)

def plot_and_write_lf_back(name, density, disp):
    with open('lf_{}_background_{}.txt'.format(clustername, name), 'w') as f:
        for i in range(ibound):
            f.write('{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\n'.format(argument[i], density[i], density[i]-disp[i], density[i]+disp[i])) 
    plot_lf('background_{}'.format(name), argument, density, disp)
    densityd = lf_cluster.density - density
    dd_disp = np.sqrt(lf_cluster.d_disp**2+disp**2)
    with open('lf_{}_{}.txt'.format(clustername, name), 'w') as f:
        for i in range(ibound):
            f.write('{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\n'.format(argument[i], densityd[i], densityd[i]-dd_disp[i], densityd[i]+dd_disp[i]))
    plot_lf(name, argument, densityd, dd_disp)

plot_and_write_lf_back('r', lf_ring.density, lf_ring.d_disp)

if args.circles:
    densityb_2c, db_disp_2c, densityb_4c, db_disp_4c, densityd_2c, dd_disp_2c, densityd_4c, dd_disp_4c = [np.zeros(ibound2) for _ in range(8)]

    #Averaging 2 circles and 4 circles regions
    densityb_2c = 0.5*(lf_circles[0].density+lf_circles[1].density)
    db_disp_2c = 0.5*np.sqrt(lf_circles[0].d_disp**2+lf_circles[1].d_disp**2)
    densityb_4c = 0.25*(lf_circles[0].density+lf_circles[1].density+lf_circles[2].density+lf_circles[3].density)
    db_disp_4c = 0.25*np.sqrt(lf_circles[0].d_disp**2+lf_circles[1].d_disp**2+lf_circles[2].d_disp**2+lf_circles[3].d_disp**2)

    plot_and_write_lf_back('2c', densityb_2c, db_disp_2c)
    plot_and_write_lf_back('4c', densityb_4c, db_disp_4c)
