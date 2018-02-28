#!/usr/bin/python
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Specifies center of cluster.')
parser.add_argument("-i", "--ifile", default='coord.txt', help='input file name, default coord.txt')
parser.add_argument("-p", "--parameters", default='parameters.txt', help="parameters filename, default parameters.txt")
parser.add_argument("-c", "--centres", default='centres.txt', help='name of textfile with maximums of linear densities, default centres.txt')
parser.add_argument("-f", "--ldfiles", default='off', help='text linear density files nametail, they are not saved by default')
parser.add_argument("-n", "--no-plots", action='store_true', help='do not show any linear density plots, neither save it as lds.png')
args = parser.parse_args()
coord_in = args.ifile
param_in = args.parameters
centres_out = args.centres

class LinDens:
    def __init__(self, d):
        self.outer_min, self.inner_min, self.inner_max, self.outer_max = [float(x) for x in d.split()]
    
    def count_ibdens(self):
        self.ibound = int((self.outer_max-self.outer_min)/self.step)+1
        self.dens = np.zeros((2, self.ibound))
    
    def count_densstep(self):
        for i in range(self.ibound):
            self.dens[0][i] = self.outer_min + i*self.step
    
    def interval(self):
        if self.r < self.delta:
            self.imin = 0
        else:
            self.imin = int((self.r - self.delta)/self.step)+1
        if self.imin > self.ibound:
            raise ValueError('Out of bounds')
        self.imax = int((self.r+self.delta)/self.step)+1
        if self.imax > self.ibound:
            self.imax = self.ibound
    
    def calculate_dens(self):
        for i in range(self.imin, self.imax):
            ri = i*self.step
            bracket = 1.0 - (self.r - ri)**2/self.delta**2
            if bracket < 0.0:
                bracket = 0.0
            self.dens[1][i] += 15.0*bracket**2/16.0/self.delta 
    
    def is_between(self):
        return self.inner_min <= self.dot <= self.inner_max


with open(param_in, 'r') as f:
    lines = f.readlines()

jlim = float(lines[0])
xld, yld, decld, rald = (LinDens(s) for s in lines[3:7])
xld.title, yld.title, decld.title, rald.title = ['X', 'Y', 'DEC', 'RA']
xld.delta, decld.delta, rald.delta = [float(x) for x in lines[1].split()]
xld.step, decld.step, rald.step = [float(x) for x in lines[2].split()]
yld.delta, yld.step = xld.delta, xld.step

for ld in (xld, yld, rald, decld):
    ld.count_ibdens()
    ld.count_densstep()

with open(coord_in, 'r') as f:
    for line in f:
        words = line.split()
        jmag, rald.dot, decld.dot, xld.dot, yld.dot = [float(words[i]) for i in (3, 0, 1, 18, 19)]
        if jmag > jlim:
            continue
        try:
            for ld in (xld, yld, rald, decld):
                ld.r = ld.dot - ld.outer_min
                ld.interval()
        except ValueError:
            continue
        for ld1, ld2 in ((yld, xld), (rald, decld)):
            if ld1.is_between():
                ld2.calculate_dens()
            if ld2.is_between():
                ld1.calculate_dens()

for ld in (xld, yld, rald, decld):
    if args.ldfiles != 'off':
        np.savetxt('{0}_{1}'.format(ld.title, args.ldfiles), np.transpose(ld.dens), fmt='%.6f', delimiter='\t')
    ld.maxim = ld.dens[0][np.argmax(ld.dens[1])]

with open(centres_out, 'w') as f:
    f.write('{0:.3f}\t{1:.4f}\t#x and alpha\n'.format(xld.maxim, rald.maxim))
    f.write('{0:.3f}\t{1:.4f}\t#y and delta\n'.format(yld.maxim, decld.maxim))

if not args.no_plots:
    fig = plt.figure(figsize=(12, 8))
    z = 1
    for ld in (xld, yld, decld, rald):
        ld.plot = plt.subplot(2, 2, z)
        plt.plot(ld.dens[0], ld.dens[1], 'k')
        ax = plt.gca()
        for xtick in ax.get_xticklabels():
            xtick.set_fontsize(10)
        for ytick in ax.get_yticklabels():
            ytick.set_fontsize(10)
        ld.plot.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.title('{} linear density'.format(ld.title))
        z+=1
    fig.tight_layout()
    fig.savefig('lds.png', dpi=100, frameon=False)
    plt.show()
