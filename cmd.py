#!/usr/bin/python
import argparse
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='CMD selection stars from cluster with prepared isochrone corridor.')
parser.add_argument("ifile", help='input file name, (coord.txt)')
parser.add_argument("corridorfile", help='file with isochrone corridor as three columns: left j-h, right j-h, j')
parser.add_argument("jmin", type=float, help='Magnitude less which all stars will be included')
parser.add_argument("-o", "--ofile", default='cmd_coord.txt', help='output file name, cmd_coord.txt by default')
args = parser.parse_args()

coord_in = args.ifile
corridor = args.corridorfile
coord_out = args.ofile
cmin = args.jmin

ci1, ci2, cmag = np.loadtxt(corridor, delimiter='\t', unpack=True)

ci1_spline_tck = interpolate.splrep(cmag, ci1)
ci2_spline_tck = interpolate.splrep(cmag, ci2)

with open(coord_in, 'r') as f:
    with open(coord_out, 'w') as g:
        for line in f:
            words = line.split()
            jmag = float(words[3])
            jhcolor = float(words[15])
            if jmag < cmin:
                g.write(line)
            else:
                
                ci_lo = interpolate.splev(jmag, ci1_spline_tck)
                ci_up = interpolate.splev(jmag, ci2_spline_tck)
                if ci_lo <= jhcolor <= ci_up:
                    g.write(line)

j, jmh  = np.loadtxt(coord_out, unpack=True, usecols=(3, 15))
plt.plot(jmh, j, '.')
plt.show()
