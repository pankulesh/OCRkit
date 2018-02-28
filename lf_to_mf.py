#!/usr/bin/python
import numpy as np
import argparse
from os import remove
from os.path import isfile
from scipy.interpolate import CubicSpline as spline
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description='Convert luminocity function into mass function using isochrone.')
parser.add_argument("ifile", type=str, help='input file name with LF')
parser.add_argument("isochrone", type=str, help='isochrone file name with 2 columns: mass and corresponding absolute magnitude')
stroka = "file name like photometry.txt with two rows as this: first row 'r    E(B-V) eE(B-V) distmod edistmod magmin magmax', second row '1168 0.253  0.10    10.294  0.3      10     20'"
parser.add_argument("parameters", type=str, help=stroka)
#parser.add_argument("parameters", type=str, help='string with distant module, lower and upper bound of magnitude, like "10.294 9 16" (in quotes!)')
parser.add_argument("-o", "--ofile", default='lf_mf.txt', help='output file name, default lf_mf.txt')
args = parser.parse_args()
lf_file = args.ifile
isochrone = args.isochrone
out_file = args.ofile

if isfile(out_file):
    remove(out_file)
with open(args.parameters, 'r') as f:
    f.readline()
    _, ebv, _, distmod, _, magmin, magmax = [float(i) for i in f.readline().split()]
#distmod, magmin, magmax = [float(i) for i in args.parameters.split()]
aj = 0.9*3.1*ebv
mass_absmag = np.loadtxt(isochrone, unpack=True)

mass_absmag_spline = spline(mass_absmag[1], mass_absmag[0], extrapolate=False)

with open(lf_file, 'r') as f:
    for line in f:
        lf_values = [float(i) for i in line.split()]
        mag,lf,lflo,lfhi = lf_values
        if magmin <= mag <= magmax:
            magabs = mag - distmod - aj
            massarg = float(mass_absmag_spline(magabs))
            masspr = float(mass_absmag_spline(magabs, 1))
            if np.isnan(massarg) or np.isnan(masspr):
                continue
            mf, mflo, mfhi = [-1/masspr*i for i in lf_values[1:]]
            with open(out_file, 'a') as g:
                g.write("{0:10.3f} {1:10.3f} {2:15.6f} {3:15.6f} {4:15.6f} {5:10.6f} {6:15.3f} {7:15.6f} {8:15.6f}\n".format(mag,magabs,lf,lflo,lfhi,massarg,mf,mflo,mfhi))

lumfun = np.loadtxt(out_file, unpack=True, usecols=(1,2,3,4))
massfun = np.loadtxt(out_file, unpack=True, usecols=(5,6,7,8))

f, (lum, mass) = plt.subplots(2)
lum.plot(lumfun[0], lumfun[1], 'k', linewidth=1)
lum.plot(lumfun[0], lumfun[2], 'k', linestyle='dashed', linewidth=1)
lum.plot(lumfun[0], lumfun[3], 'k', linestyle='dashed', linewidth=1)
lum.set_title('Luminocity function')
lum.invert_xaxis()

mass.plot(massfun[0], massfun[1], 'k', linewidth=1)
mass.plot(massfun[0], massfun[2], 'k', linestyle='dashed', linewidth=1)
mass.plot(massfun[0], massfun[3], 'k', linestyle='dashed', linewidth=1)
mass.set_title('Mass function')
f.savefig("mf_lf.png", dpi=300)
plt.show()
