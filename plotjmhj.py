#!/usr/bin/python
import argparse
import numpy as np
from math import pi, log10
from scipy import interpolate
import matplotlib.pyplot as plt
des = """Plotting isochrone with 10' zone in J-H J coordinates,
        10' zone must be named like NAME.txt, isochrone - 'isochrone_NAME.txt'
         """
parser = argparse.ArgumentParser(description=des)
parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
args = parser.parse_args()
name = args.name
coord_in = "{}.txt".format(name)
iso = "isochrone_{}.txt".format(name)
arrin = np.genfromtxt(coord_in, skip_header=2, usecols=(6, 8), unpack=True)
cluster = np.array([arrin[0]-arrin[1], arrin[0]])
isop = np.loadtxt(iso, unpack=True)
fig = plt.figure(figsize=(5,5))
plt.plot(cluster[0], cluster[1], ".", color='k')
plt.plot(isop[0], isop[1], "r", linewidth=0.5)
plt.gca().invert_yaxis()
plt.axes().set_aspect(0.5)
plt.title(name, fontsize=20)
plt.xlabel(r"J-H", fontsize=20)
plt.ylabel(r"J", fontsize=20)
plt.show()
fig.savefig("isochronecluster_{}.png".format(name), dpi=100, figsize=(5,5), transparent=False,bbox_inches='tight', pad_inches=0.1, frameon=False)
