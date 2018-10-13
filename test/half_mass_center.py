#!/usr/bin/python
import argparse
import numpy as np
from multiprocessing import Pool, cpu_count

def get_points_inside(points, center, radius):
    ids = ((points[:,0]-center[0])**2+(points[:,1]-center[1])**2)<=radius**2
    return points[ids]

def get_points_square(points, halfedge):
    id1 = np.abs(points[:,0]) <= halfedge
    newpoints = points[id1]
    id2 = np.abs(newpoints[:,1]) <= halfedge
    return newpoints[id2]
    
def get_half_mass_radii(centerpoint):
    global cluster
    step = 0.1
    hf_rad = 0.
    n_tot = len(cluster)
    n_counted = len(get_points_inside(cluster, centerpoint, hf_rad))
    while n_counted < 0.5*n_tot:
        hf_rad += step
        n_counted = len(get_points_inside(cluster, centerpoint, hf_rad))
#    print(f"Для звезды ({centerpoint[0]}, {centerpoint[1]}) hf_rad={hf_rad:.4f}")
    return hf_rad

parser = argparse.ArgumentParser(description='Скрипт для нахождения центра скопления с помощью вычисления для каждой точки сетки [-5:5]x[-5:5] радиуса половины массы и вычисления минимального из них.')
parser.add_argument("infile", help='Входной файл со звёздами скопления')
#parser.add_argument("radius", type=float, help='Примерный радиус скопления (внутри не должно быть иных максимумов, кроме как центра скопления')
args = parser.parse_args()
inarray = np.loadtxt(args.infile, comments='#', usecols=(0,1))
radius = 5
cluster = get_points_square(inarray, radius)
e = radius
xv, yv = np.meshgrid(np.linspace(-1*e, e, int(2*e/0.1)), np.linspace(-1*e, e, int(2*e/0.1)))
zone = np.column_stack((xv.flatten(), yv.flatten()))
print(len(cluster))
p = Pool(cpu_count())
hfs = list(p.map(get_half_mass_radii, zone))
x0, y0 = zone[np.argmin(hfs)]
print(f"Центр скопления находится в {x0:.2f} {y0:.2f} с hmr={np.amin(hfs):.2f}")

