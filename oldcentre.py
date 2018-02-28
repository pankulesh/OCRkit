#!/usr/bin/python
import argparse
import matplotlib.pyplot as plt
import numpy as np
size = 1000
parser = argparse.ArgumentParser(description='Specifies center of cluster.')
parser.add_argument("-i", "--ifile", default='coord.txt', help='input file name, default coord.txt')
parser.add_argument("-o", "--ofile", default='centre.txt', help='output files nametail, default centre.txt')
parser.add_argument("-p", "--parameters", default='parameters.txt', help="parameters filename, default parameters.txt")
parser.add_argument("-c", "--centres", default='all_centres.txt', help='centres filename, default all_centres.txt')
args = parser.parse_args()
coord_in = args.ifile
centre_out = args.ofile
param_in = args.parameters
centres_out = args.centres
def interval(r, delta, step, ibound):
    if r<delta:
        imin = 0
    else:
        imin = int((r-delta)/step)+1
    if imin > ibound:
        return 'e', 'e' 
    imax = int((r+delta)/step)
    if imax > ibound:
        imax = ibound
    return imin, imax
def dens_create(imin, imax, r, step, delta, dens):
    for i in range(imin, imax+1):
        ri = i*step
        bracket = 1.0 - (r - ri)**2/delta**2
        if bracket < 0.0:
            bracket = 0.0
        dens[i] = dens[i] + 15.0*bracket**2/16.0/delta
def centre_write(f, ibound, i0, step, dens):
    for i in range(ibound):
        f.write('{0:.6f}\t{1:.6f}\n'.format(i0+i*step, dens[i]))
with open(param_in, 'r') as f:
    jlim = float(f.readline())
    delta, deltadec, deltara = [float(x) for x in f.readline().split()]
    step, stepdec, stepra = [float(x) for x in f.readline().split()]
    x0, x1, x2, xmax = [float(x) for x in f.readline().split()]
    y0, y1, y2, ymax = [float(x) for x in f.readline().split()]
    dec0, dec1, dec2, decmax = [float(x) for x in f.readline().split()]
    ra0, ra1, ra2, ramax = [float(x) for x in f.readline().split()]
ixbound = int((xmax-x0)/step)+1
iybound = int((ymax-y0)/step)+1
irabound = int((ramax-ra0)/stepra)+1
idecbound = int((decmax-dec0)/stepdec)+1
densx = [0.0 for x in range(size)]
densy = [0.0 for x in range(size)]
densra = [0.0 for x in range(size)]
densdec = [0.0 for x in range(size)]
n_tot=0
with open(coord_in, 'r') as f:
    line = f.readline()
    while line:
        words = line.split()
        line = f.readline()
        jmag = float(words[3])
        if jmag > jlim:
            continue
        ra = float(words[0])
        dec = float(words[1])
        x = float(words[18])
        y = float(words[19])
        rx = x-x0
        ry = y-y0
        rra = ra-ra0
        rdec = dec-dec0
        ixmin, ixmax = interval(rx, delta, step, ixbound)
        iymin, iymax = interval(ry, delta, step, iybound)
        iramin, iramax = interval(rra, deltara, stepra, irabound)
        idecmin, idecmax = interval(rdec, deltadec, stepdec, idecbound)
        if 'e' in (ixmin, ixmax, iymin, iymax, iramin, iramax, idecmin, idecmax):
            continue
        n_tot += 1
        if y1 < y < y2:
            dens_create(ixmin, ixmax, rx, step, delta, densx)
        else:
            continue
        if  x1 < x < x2:
            dens_create(iymin, iymax, ry, step, delta, densy)
        else:
            continue
        if dec1 < dec < dec2:
            dens_create(iramin, iramax, rra, stepra, deltara, densra)
        else:
            continue
        if ra1 < ra < ra2:
            dens_create(idecmin, idecmax, rdec, stepdec, deltadec, densdec)
        else:
            continue
with open('{0}_{1}'.format('x', centre_out), 'w') as f:
	centre_write(f, ixbound, x0, step, densx)
with open('{0}_{1}'.format('y', centre_out), 'w') as f:
	centre_write(f, iybound, y0, step, densy)
with open('{0}_{1}'.format('ra', centre_out), 'w') as f:
	centre_write(f, irabound, ra0, stepra, densra)
with open('{0}_{1}'.format('dec', centre_out), 'w') as f:
	centre_write(f, idecbound, dec0, stepdec, densdec)
maximx = x0 + step*max(enumerate(densx), key=lambda x: x[1])[0]
maximy = y0 + step*max(enumerate(densy), key=lambda x: x[1])[0]
maximdec = dec0 + stepdec*max(enumerate(densdec), key=lambda x: x[1])[0]
maximra = ra0 + stepra*max(enumerate(densra), key=lambda x: x[1])[0]
with open(centres_out, 'w') as f:
    f.write('{0:.3f}\t{1:.3f}\t#x and alpha\n'.format(maximx, maximra))
    f.write('{0:.3f}\t{1:.3f}\t#y and delta\n'.format(maximy, maximdec))
