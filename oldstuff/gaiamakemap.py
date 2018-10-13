#!/usr/bin/python
#usage: python/start gaiamakemap.py <coord.txt> <gmax> <halfwidth> -n <clustername> -s <30x30 or 60x60>
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.interpolate import griddata
parser = argparse.ArgumentParser(description='Make density map with maximum magnitude and kernel halfwidth.')
parser.add_argument("infile", help='input file name')
parser.add_argument("gmax", help='maximum magnitude, m')
parser.add_argument("halfwidth", help='kernel halfwidth, arcmin')
parser.add_argument("-s", "--size", default='30x30', help='size of map (30x30 or 60x60), default 30x30')
parser.add_argument("-n", "--cluster-name", default='UNDEFINED', help="set clustername for map titles, clustername isn't specified by default!")
args = parser.parse_args()
coord_in = args.infile
clustername = args.cluster_name
size = args.size
n = 150
magmin = 0
ndim = 2*n + 1
levels = 10
internumber = 100
h = args.halfwidth
gmax = args.gmax
def makemap(xydens, in_file, out_file, name, labels=True, colorbar=True, grid=True):
#    a = np.loadtxt(in_file, delimiter='\t', unpack=True)
    a = np.transpose(xydens)
    X, Y, Z = a[0], a[1], a[2] 
    xi = np.linspace(min(X), max(X), internumber) 
    yi = np.linspace(min(Y), max(Y), internumber)
    zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='linear')  
    fig = plt.figure(figsize=(10,10))
    contour = plt.contour(xi, yi, zi, levels, colors='k')
    if labels:
        plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=10)
    contour_filled = plt.contourf(xi, yi, zi, levels, cmap=plt.cm.Greys)
    plt.xlabel("x, arcmin") 
    plt.ylabel("y, arcmin")
    if colorbar:
        plt.colorbar(contour_filled)
    if grid:
        plt.grid(True)
    plt.axes().set_aspect('equal')
    plt.title(name)
    fig.savefig(out_file, dpi=100, frameon=False)
    plt.close(fig)
    print('Map for {0} saved as {1}\n'.format(in_file, out_file))
    return

def densmap(delta, step, magmax, map_out):
    density = np.zeros((301, 301))
    xs = np.array([(-n+i)*step for i in range(0, ndim)])
    ys = np.array([(-n+i)*step for i in range(0, ndim)])
    xydens = np.zeros((90601, 3))
    k=0
    for star in xyg:
        x,y,mag = star
        if magmin < mag < magmax:
            imin = int((y - delta)/step + n)
            imax = int((y + delta)/step + n)
            if imin < 0:
                imin = 0
            if imax > ndim:
                imax = ndim
            if imax < 0 or imin > ndim:
                continue
            jmin = int((x - delta)/step + n)
            jmax = int((x + delta)/step + n)
            if jmin < 0:
                jmin = 0 
            if jmax > ndim:
                jmax = ndim
            if jmax < 0 or jmin > ndim:
                continue
            k += 1
            for i in range(imin, imax):
                for j in range(jmin, jmax):
                    r2 = (x - xs[j])**2+(y - ys[i])**2
                    bracket = 1. - r2/delta**2
                    density[i][j] += 3*bracket**2/(np.pi*delta**2)
    
    with open(map_out, 'w') as f:
        m = 0
        for i in range(0, ndim):
            for j in range(0, ndim):
                if (xs[j]**2 + ys[i]**2) < n**2:
                    xydens[m] = xs[j], ys[i], density[i][j]
                    m+=1
                    f.write('{0:.3f}\t{1:.3f}\t{2:.6f}\n'.format(xs[j], ys[i], density[i][j]))
    print('{0} strings were handled, {1} is done.'.format(k, map_out))
    return xydens
if size == '30x30':
    step = 0.1
elif size == '60x60':
    step = 0.2
else:
    print('wrong size, stopping...')
    exit(2)
print("""Output files will look like this: map<size of map>_<max magnitude>_<kernel halfwidth>.[txt, png]\n""")
xyg = np.loadtxt(coord_in, skiprows=49, usecols=(0, 1, -1))
name = 'map{0}_{1}_{2}.'.format(size, gmax, h)
if clustername=='UNDEFINED':
    sname = "Densmap: gmax = {0}m, h = {1}'".format(gmax, h)
else:
    sname = "Densmap of {0}: gmax = {1}m, h = {2}'".format(clustername, gmax, h)
xydens = densmap(h, step, gmax, name+'txt')
makemap(xydens, name+'txt', name+'png', sname)
