import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def plot_map(in_file, out_file, name, labels=True, colorbar=True, grid=True, show=False):
    levels = 10
    internumber = 100
    df = pd.read_csv(in_file, header=None, delim_whitespace=True)
    X, Y, Z = df[0], df[1], df[2]
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
    fig.savefig(out_file, dpi=100)
    if show:
        plt.show()
    plt.close(fig)
    print('Карта, построенная по {0}, сохранена как {1}.\n'.format(in_file, out_file))
    return

def make_map(xymag, delta, step, magmax, map_out):
    n = 150
    magmin = 0
    ndim = 2*n + 1
    density = np.zeros((ndim, ndim))
    xs = np.array([(-n+i)*step for i in range(0, ndim)])
    ys = np.array([(-n+i)*step for i in range(0, ndim)])
    k=0
    for x, y, mag in xymag:
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
        f.write("x,y,z\n")
        for i in range(0, ndim):
            for j in range(0, ndim):
                if np.abs(xs[j] + ys[i])+np.abs(xs[j]-ys[j]) <= 2*n*step:
                    m+=1
                    f.write('{0:.3f},{1:.3f},{2:.6f}\n'.format(xs[j], ys[i], density[i][j]))
    print('{0} звёзд обработано, карта {1} готова.'.format(k, map_out))
