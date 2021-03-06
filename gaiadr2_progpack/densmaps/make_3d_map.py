#!/usr/bin/python
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

from scipy.interpolate import griddata

parser = argparse.ArgumentParser(description='Программа для построение карты поверхностной плотности по готовому файлу map_....txt')
parser.add_argument("infile", type=str, help='Имя входного файла с x y z колонками')
parser.add_argument("-c", "--circles", type=str, help='Укажите файл с параметрами и радиусами скопления и области сравнений для рисования этих областей на карте')
args = parser.parse_args()

coord_in = args.infile
spl = coord_in.split('_')
clustername = spl[1]
size = spl[2]
step = int(size.split('x')[0])/300
delta = float(spl[-1].rsplit('.', 1)[0])
maglim = float(spl[3])

name = 'map_{3}_{0}_{1}_{2}.'.format(size, maglim, delta, clustername)
in_file = name+'txt'
out_file = name+'png'
sname = "Карта поверхностной плотности {0}:\nmaglim = {1}, h = {2}'".format(clustername, maglim, delta)
xyz = pd.read_csv(in_file)

def circle(x, y, radius, text):
    from matplotlib.patches import Circle
    circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=1, edgecolor='black', facecolor=(0, 0, 0, .0125))
    ax = plt.gca()
    ax.add_artist(circle)
    ax.text(x, y, text, ha='center', va='top', weight='bold', color='blue')
labels, colorbar, grid = True, True, True
levels = 10
internumber = 100
X, Y, Z = xyz.x, xyz.y, xyz.z 
xi = np.linspace(np.amin(X), np.amax(X), internumber) 
yi = np.linspace(np.amin(Y), np.amax(Y), internumber)
zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='linear')  
ax.plot_surface(xi, yi, zi)
ax.xlabel("x, arcmin") 
ax.ylabel("y, arcmin")
ax.title(sname)

fig.savefig(out_file, dpi=100, frameon=False)
plt.show()
plt.close(fig)
print('Карта, построенная по {0}, сохранена как {1}.\n'.format(in_file, out_file))
