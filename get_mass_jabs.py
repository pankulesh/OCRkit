#!/usr/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt
des = """Making cluster's isochrone from http://stev.oapd.inaf.it/cgi-bin/cmd and photometry parameters distmod and E(B-V);
         distmod and E(B-V) must be in file 'photometry_NAME.txt', isochrone will be saved in 'isochrone_NAME.txt'"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('iso_in', type=str, help="Downloaded isochrone filename")
parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
parser.add_argument('param', type=str, help="E(B-V) and distmod values in quotes")
args = parser.parse_args()
name = args.name
iso_in = args.iso_in
iso_out = "isochrone_{}.txt".format(name)
param = args.param.split()
ebv, distmod = (float(x) for x in param)
ebv2ejh = 0.34        # magic constants
ebv2aj = 2.43*ebv2ejh #
fixmag = distmod + ebv2aj*ebv
dtype = [('mass', float), ('j', float)]
isoarr = np.zeros((500,), dtype=dtype)
arrin = np.genfromtxt(iso_in, skip_header=1, usecols=(3, 13))
i = 0
for m, j in arrin:
    mj = j + fixmag
    isoarr[i] = (m, mj)
    i += 1
iso = np.sort(isoarr[:i], order='j')
np.savetxt(iso_out, iso, fmt='%.3f', delimiter='\t', newline='\n')
isop = np.loadtxt(iso_out, unpack=True)
plt.plot(isop[0], isop[1], ".")
plt.gca().invert_yaxis()
plt.show()

