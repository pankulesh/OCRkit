#!/usr/bin/python
import numpy as np
import argparse
from scipy.interpolate import CubicSpline as spline
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Превращает функцию блеска в функцию светимости и функцию масс соответственно с помощью готовой изохроны-соответствия главной последовательности с соответствующими массами')
parser.add_argument("ifile", type=str, help='Входной файл с функцией блеска')
parser.add_argument("isochrone", type=str, help='Входной файл с изохроной: в первой колонке масса, в остальных величины в соответствующих полосах')
parser.add_argument("-j", "--jmag", action='store_true', help='Укажите этот флаг, если нужно соответствие с величиной J (5 колонка вместо 1)')
parser.add_argument("magfix", type=float, help="Укажите соответствующий видимый модуль расстояния")
args = parser.parse_args()
lf_file = args.ifile
isochrone = args.isochrone
magfix = args.magfix
out_file = lf_file[:-4]+f'_mf_magfix={magfix}.txt'

massmag = (0, 1)
magtype = 'G'
if args.jmag:
    massmag = (0, 5)
    magtype = 'J'
mass_absmag = np.loadtxt(isochrone, unpack=True, usecols=massmag)
magmin, magmax = np.amin(mass_absmag[1]), np.amax(mass_absmag[1])

mass_absmag_spline = spline(mass_absmag[1], mass_absmag[0], extrapolate=False)

with open(lf_file, 'r') as f:
    with open(out_file, 'w') as g: 
        for line in f:
            lf_values = [float(i) for i in line.split()]
            mag,lf,lflo,lfhi = lf_values
            magabs = mag - magfix
            if magmin <= magabs <= magmax:
                massarg = float(mass_absmag_spline(magabs))
                masspr = float(mass_absmag_spline(magabs, 1))
                if np.isnan(massarg) or np.isnan(masspr):
                    continue
                mf, mflo, mfhi = [-1/masspr*i for i in lf_values[1:]]
                g.write("{0:10.3f} {1:10.3f} {2:15.6f} {3:15.6f} {4:15.6f} {5:10.6f} {6:15.3f} {7:15.6f} {8:15.6f}\n".format(mag,magabs,lf,lflo,lfhi,massarg,mf,mflo,mfhi))

lumfun = np.loadtxt(out_file, unpack=True, usecols=(1,2,3,4))
massfun = np.loadtxt(out_file, unpack=True, usecols=(5,6,7,8))

f, (lum, mass) = plt.subplots(2)
lum.plot(lumfun[0], lumfun[1], 'k', linewidth=1)
lum.plot(lumfun[0], lumfun[2], 'k', linestyle='dashed', linewidth=1)
lum.plot(lumfun[0], lumfun[3], 'k', linestyle='dashed', linewidth=1)
lum.set_title('Функция светимости')
a = np.arange(np.amin(lumfun[0]), np.amax(lumfun[0])+1, step=1)
lum.set_xticks(a)
lum.set_xlim(np.amin(a), np.amax(a))
lum.set_ylim(0, None)
lum.set_xlabel(f'M{magtype}')
lum.invert_xaxis()

mass.plot(massfun[0], massfun[1], 'k', linewidth=1)
mass.plot(massfun[0], massfun[2], 'k', linestyle='dashed', linewidth=1)
mass.plot(massfun[0], massfun[3], 'k', linestyle='dashed', linewidth=1)
mass.set_title('Функция масс')
mass.set_xlabel('M, Msun')
a = np.arange(np.amin(massfun[0]), np.amax(massfun[0])+0.1, step=0.4)
mass.set_xticks(a)
mass.set_xlim(np.amin(a), np.amax(a))
mass.set_ylim(0, None)

plt.tight_layout()
f.savefig("mf_lf.png", dpi=600)
plt.show()
