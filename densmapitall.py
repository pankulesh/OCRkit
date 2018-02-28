#!/usr/bin/python
#usage: python/start densmapitall.py <coord.txt> <clustername> <30x30 or 60x60>
import argparse
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.interpolate import griddata
parser = argparse.ArgumentParser(description='Make 30 density map with different maximum magnitude and kernel halfwidth.')
parser.add_argument("-i", "--ifile", default='coord.txt', help='input file name, default coord.txt')
parser.add_argument("-s", "--size", default='30x30', help='size of map (30x30 or 60x60), default 30x30')
parser.add_argument("-n", "--cluster-name", default='UNDEFINED', help="set clustername for map titles, clustername isn't specified by default!")
args = parser.parse_args()
coord_in = args.ifile
clustername = args.cluster_name
size = args.size
n = 150
jmagmin = 0
ndim = 2*n + 1
levels = 10
internumber = 100
def makemap(in_file, out_file, name, labels=True, colorbar=True, grid=True):
    a = np.loadtxt(in_file, delimiter='\t', unpack=True)
    X,Y,Z=a[0],a[1],a[2] 
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
def densmap(delta, step, jmagmax, map_out):
    density = [[0.0 for q in range(400)] for w in range(400)]
    with open(coord_in, 'r') as f:
        line = f.readline()
        k = 0
        while line:
            k += 1
            words = line.split()
            line = f.readline()
            jmag = float(words[3])
            if (jmagmin > jmag) or (jmag > jmagmax):
                continue
            x = float(words[18])
            y = float(words[19])
            a = x/step + n + 1
            b = y/step + n + 1
            imin = int(b - delta/step) + 1
            if imin > ndim:
                continue
            if imin < 1:
                imin = 1
            imax = int(b + delta/step)
            if imax > ndim:
                imax = ndim
            if imax < 1:
                continue
            jmin = int(a - delta/step) + 1
            if jmin > ndim:
                continue
            if jmin < 1:
                jmin = 1
            jmax = int(a + delta/step) 
            if jmax > ndim:
                jmax = ndim
            if jmax < 1:
                continue
            for i in range(imin, imax+1):
                for j in range(jmin, jmax+1):
                    xsj = (-n+j-1)*step
                    ysi = (-n+i-1)*step
                    r2 = (x - xsj)**2+(y - ysi)**2
                    bracket = 1. - r2/delta**2
                    density[i][j] = density[i][j] + 3*bracket**2/(pi*delta**2)
    with open(map_out, 'w') as f:
        for i in range(1, ndim+1):
            ysi = (-n+i-1)*step
            for j in range(1, ndim+1):
                xsj = (-n+j-1)*step
                if (xsj**2 + ysi**2) < n**2:
                    f.write('{0:.3f}\t{1:.3f}\t{2:.6f}\n'.format(xsj, ysi, density[i][j]))
    print('{0} strings were handled, {1} is done.'.format(k, map_out))
    return
if size == '30x30':
    step = 0.1
elif size == '60x60':
    step = 0.2
else:
    print('wrong size, stopping...')
    exit(2)
maxmags = [11, 12, 13, 14, 15, 16]
deltas = [1, 2, 3, 4, 5]
print("""Output files will look like this: map<size of map>_<max jmagnitude>_<kernel halfwidth>.[txt, png]
    Maps compilation will be in file maps<size>.png
    It will take a time, go make some coffee or tea...\n""")
for j in maxmags:
    for k in deltas:
        name = 'map{0}_{1}_{2}.'.format(size, j, k)
        if clustername=='UNDEFINED':
            sname = "Densmap: jmax = {0}m, h = {1}'".format(j, k)
        else:
            sname = "Densmap of {0}: jmax = {1}m, h = {2}'".format(clustername, j, k)
        densmap(k, step, j, name+'txt')
        makemap(name+'txt', name+'png', sname)
fig, axes = plt.subplots(6, 5, figsize=(25, 30), subplot_kw={'xticks': [], 'yticks': []})
fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
img = []
for i in maxmags:
    for j in deltas:
        img.append(mpimg.imread('map{0}_{1}_{2}.png'.format(size, i, j)))
for ax, image in zip(axes.flat, img):
    ax.axis('off')
    ax.imshow(image)
fig.savefig('maps{}.png'.format(size), dpi=300, frameon=False)
