#!/usr/bin/python
import numpy as np
import argparse
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description='Plot all profiles, need files "density_NAME_J_h.txt" for all given J, example of usage: python plotprofiles.py NGC7142 2 "11 12 13 14 15 16"')
parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
parser.add_argument('h', type=int, help="Kernel halfwidth (as in filenames)")
parser.add_argument('mags', type=str, help="String with limit magnitudes, i.e. '11 12 13 14 15 16' (quotes are important!)")
args = parser.parse_args()
h = args.h
NAME = args.name
jmaglims = args.mags.split()
for jmaglim in jmaglims:
    a = np.loadtxt('density_{0}_{1}_{2}.txt'.format(NAME, jmaglim, h), unpack=True, usecols=(0,1,3,4))
    r = a[0]
    fig = plt.figure(figsize=(10,10))
    plt.title("{0}, h={1}', Jlim={2}m".format(NAME, h, jmaglim), fontsize=26)
    plt.plot(r, a[1], 'k', linewidth=1)
    plt.plot(r, a[2], 'k', linestyle='dashed', linewidth=1)
    plt.plot(r, a[3], 'k', linestyle='dashed', linewidth=1)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    ax.set_xlabel(r'distance from center, $arcmin$', fontsize=26)
    ax.set_ylabel(r'radial density, $1/arcmin^2$', fontsize=26)
    ax.set_xlim(0, np.amax(r))
    ax.set_ylim(0, np.amax(a[3]))
    ax.tick_params(width=2)
    plt.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
    fig.savefig('density_{0}_{1}_{2}.png'.format(NAME, jmaglim, h), dpi=100, transparent=False)
    print('density_{0}_{1}_{2}.png was saved.'.format(NAME, jmaglim, h))
    plt.close()
