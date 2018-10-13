#!/usr/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt
des = """Making cluster's isochrone from http://stev.oapd.inaf.it/cgi-bin/cmd and photometry parameters r and E(B-V);
         r and E(B-V) must be in file 'photometry_NAME.txt', isochrone will be saved in 'isochrone_NAME.txt'"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('iso_in', type=str, help="Downloaded isochrone filename")
parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
args = parser.parse_args()
name = args.name
iso_in = args.iso_in
iso_out = "isochrone_{}.txt".format(name)
#param = ['1240', '0.34', '0.2']
pr, ebv = 11.251, 0.393
ebv2ejh = 0.34        # magic constants
ebv2aj = 2.43*ebv2ejh #
fixcolor = ebv2ejh*ebv
fixmag = pr + ebv2aj*ebv
dtype = [('j-h', float), ('j', float)]
isoarr = np.zeros((500,), dtype=dtype)
arrin = np.genfromtxt(iso_in, skip_header=1, delimiter="\t", usecols=(13, 14))
i = 0
for j, h in arrin:
    mj = j + fixmag
#    if mj > 16.0:
#        continue
    jmh = j - h + fixcolor
    isoarr[i] = (jmh, mj)
    i += 1
iso = np.sort(isoarr[:i], order='j')
np.savetxt(iso_out, iso, fmt='%.3f', delimiter='\t', newline='\n')
coord_in = "{}.txt".format(name)
arrin = np.genfromtxt(coord_in, skip_header=2, usecols=(6, 8), unpack=True)
cluster = np.array([arrin[0]-arrin[1], arrin[0]])
isop = np.loadtxt(iso_out, unpack=True)
fig = plt.figure(figsize=(5,5))
plt.plot(cluster[0], cluster[1], ".", color='k')
plt.plot(isop[0], isop[1], "r", linewidth=0.5)
plt.ylim(5, 20)
plt.gca().invert_yaxis()
plt.axes().set_aspect(0.5)
plt.title(name, fontsize=20)
plt.xlabel(r"J-H", fontsize=20)
plt.ylabel(r"J", fontsize=20)
plt.show()
fig.savefig("isochronecluster_{}.png".format(name), dpi=100, figsize=(5,5), transparent=False,bbox_inches='tight', pad_inches=0.1, frameon=False)

