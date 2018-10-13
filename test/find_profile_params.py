#!/usr/bin/python

from os import listdir
import re
import argparse
import numpy as np
from scipy.integrate import simps
from scipy.interpolate import CubicSpline as spline
from matplotlib import pyplot as plt

PRECISION = 0.01

parser = argparse.ArgumentParser(description='Getting radius of cluster and mean background stellar density from radial density profile')
parser.add_argument("name", help='Name of cluster (as in filenames density_NAME_MAG_H.txt)')
parser.add_argument("h", help='Kernel halfwidth (as in filenames density_NAME_MAG_H.txt)')
args = parser.parse_args()
name = args.name
h = args.h

def count_R_Fb(mag):
    r, dens, densms, densps = np.loadtxt('density_{0}_{1}_{2}.txt'.format(name, mag, h), unpack=True, usecols = (0, 1, 3, 4))
    rmax = np.amax(r)
    densmax = int(np.amax(densps))+1
    density_spline = spline(r, dens)                                                                  
    density_msigma_spline = spline(r, densms)                                                         
    density_psigma_spline = spline(r, densps)                                                         
                                                                                                      
    radius = np.arange(0, rmax+0.01, 0.01)

    density = density_spline(radius)
    density_low = density_msigma_spline(radius)
    density_up = density_psigma_spline(radius)

    rect_deviation = np.zeros_like(radius)

    for i in range(len(radius)):
        rect = density[i]*(rmax-radius[i])
        integral = simps(density[i:], radius[i:])
        rect_deviation[i] = abs(rect-integral)
        
    plt.plot(radius, density, 'k', linewidth=1)
    plt.plot(radius, density_low, 'k', linestyle='dashed', linewidth=1)
    plt.plot(radius, density_up, 'k', linestyle='dashed', linewidth=1)

    for i in range(len(radius)):
        if rect_deviation[i] < 0.01:
            R, Fb = radius[i], density[i]
            dFb = max(density_up[i]-Fb, Fb-density_low[i])
            ibound = i
            break

    dR1, dR2, dR3, dR4 = 0, 0, 0, 0
    roots = (spline(radius, density-Fb, extrapolate=False)).roots()
    roots_l = (spline(radius, density_low-Fb, extrapolate=False)).roots()
    roots_u = (spline(radius, density_up-Fb, extrapolate=False)).roots()
    
    rl_left = np.extract(roots_l < R, roots_l)
    ru_left = np.extract(roots_u < R, roots_u)
    rl_right = np.extract(roots_l > R, roots_l)
    ru_right = np.extract(roots_u > R, roots_u)
    ind_r = (np.argwhere(roots == R))[0][0]
    if ind_r == roots.size-1:
        right_r = rmax - R
    else:
        right_r = roots[ind_r+1]-R
    
    if ind_r == 0:
        left_r = R
    else:
        left_r = R - roots[ind_r-1]
    
    if rl_left.size != 0:
        if np.amin(R-rl_left) < left_r:
            dR1 = np.amin(R-rl_left)

    if ru_left.size != 0:   
        if np.amin(R-ru_left) < left_r:
            dR2 = np.amin(R-ru_left)
    
    if rl_right.size != 0:   
        if np.amin(rl_right-R) < right_r:
            dR3 = np.amin(rl_right-R)

    if ru_right.size != 0:
        if np.amin(ru_right-R) < right_r:
            dR4 = np.amin(ru_right-R)

    dR = max(dR1, dR2, dR3, dR4)

    if dR > R:
        dR = R
    if dFb > Fb:
        dFb = Fb

    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    ax.axhline(Fb, c='b', alpha=0.8, lw=0.5, label = "Fb = {0:5.2f} ± {1:5.2f}".format(Fb, dFb))
    ax.axvline(R, c='r', alpha=0.8, lw=0.5, label = "R = {0:5.2f} ± {1:5.2f}".format(R, dR))
    ax.tick_params(width=2)
    ax.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
    ax.set_xlim(0, rmax)
    ax.set_ylim(0, densmax)
    ax.legend()
    plt.title("Радиальный профиль плотности {0},\nh={1}', maglim={2}m".format(name, h, mag), fontsize=18)
    ax.set_xlabel(r'Угловое расстояние от центра, $arcmin$', fontsize=16)
    ax.set_ylabel(r'Радиальная плотность, $arcmin^{-2}$', fontsize=16)
    plt.savefig('density_{0}_{1}_{2}_RFb.png'.format(name, mag, h), dpi=300)
    plt.close()
    return (R, dR, Fb, dFb)

regexp = r'density_{0}_[0-9]+_{1}.txt'.format(name, h)
files = [f for f in listdir('.') if re.match(regexp, f)]
mags = [int(mag.split('_')[2]) for mag in files]
mags.sort()

with open("parameters_{0}_h={1}.txt".format(name, h), 'w') as f:
    f.write('Mag     R     dR     Fb    dFb\n')
    for mag in mags:
        params = count_R_Fb(mag)
        f.write('{0:2d}  {1:5.2f}  {2:5.2f}  {3:5.2f}  {4:5.2f}\n'.format(mag, *params))
        print("R = {0:5.2f}, dR = {1:5.2f}, Fb = {2:5.2f}, dFb = {3:5.2f} for mag = {4:2d}".format(*params, mag))
