#!/usr/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt
des = """Making cluster's isochrone from http://stev.oapd.inaf.it/cgi-bin/cmd;
         isochrone will be saved in 'isochrone_mass_gabs_NAME.txt'"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('iso_in', type=str, help="Downloaded isochrone filename")
parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
args = parser.parse_args()
name = args.name
iso_in = args.iso_in
iso_out = "isochrone_mass_gabs{}.txt".format(name)

arrin = np.loadtxt(iso_in, skiprows=1, usecols=(3, -3), dtype = [('mass', float), ('g', float)])
iso = np.sort(arrin, order='g')

np.savetxt(iso_out, iso, fmt='%7.3f', newline='\n')
a = np.loadtxt(iso_out, unpack=True)
plt.plot(a[1],a[0], '.')
plt.show()
