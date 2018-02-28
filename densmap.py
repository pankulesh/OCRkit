#!/usr/bin/python
import sys
from math import pi
coord_in = sys.argv[1]
map_param_in = sys.argv[2]
map_out = sys.argv[3]
with open(map_param_in, 'r') as f:
    param = f.readline().split()
delta, step, n, jmagmin, jmagmax = (float(x) for x in param) 
n = int(n)
ndim = 2*n + 1
density = [[0.0 for q in range(500)] for w in range(500)]
with open(coord_in, 'r') as f:
    line = f.readline()
    k = 1
    while line:
        k += 1
        print(k)
        words = line.split()
        jmag = float(words[3])
        if (jmagmin > jmag) or (jmag > jmagmax):
            line = f.readline()
            continue
        x = float(words[18])
        y = float(words[19])
        a = x/step + n + 1
        b = y/step + n + 1
        imin = int(b - delta/step) + 1
        if imin > ndim:
            line = f.readline()
            continue
        if imin < 1:
            imin = 1
        imax = int(b + delta/step)
        if imax > ndim:
            imax = ndim
        if imax < 1:
            line = f.readline()
            continue
        jmin = int(a - delta/step) + 1
        if jmin > ndim:
            line = f.readline()
            continue
        if jmin < 1:
            jmin = 1
        jmax = int(a + delta/step) 
        if jmax > ndim:
            jmax = ndim
        if jmax < 1:
            line = f.readline()
            continue
        for i in range(imin, imax+1):
            for j in range(jmin, jmax+1):
                xsj = (-n+j-1)*step
                ysi = (-n+i-1)*step
                r2 = (x - xsj)**2+(y - ysi)**2
                bracket = 1. - r2/delta**2
                density[i][j] = density[i][j] + 3*bracket**2/(pi*delta**2)
        line = f.readline()
with open(map_out, 'w') as f:
    for i in range(1, ndim+1):
        ysi = (-n+i-1)*step
        for j in range(1, ndim+1):
            xsj = (-n+j-1)*step
            if (xsj**2 + ysi**2) < n**2:
                f.write('{0:.3f}\t{1:.3f}\t{2:.6f}\n'.format(xsj, ysi, density[i][j]))
