#!/usr/bin/python
import argparse
import numpy as np
#from kde_map_functions import plot_map, make_map
from sklearn.neighbors import KernelDensity
from multiprocessing import Pool, cpu_count
from sklearn.model_selection import GridSearchCV

parser = argparse.ArgumentParser(description='Программа для построение карты поверхностной плотности конкретного размера, параметра KDE и предельной звёздной величины')
parser.add_argument("infile", type=str, help='Имя входного файла, содержащего координаты и звёздные величины звёзд')
parser.add_argument("name", type=str, help="Название скопления (для заголовков графика)")
parser.add_argument("delta", type=float, help="Величина параметра KDE в угловых минутах")
parser.add_argument("maglim", type=float, help="Предельная звёздная величина звёзд, используемых в построении карты")
parser.add_argument("step", type=float, help="Размер шага в угловых минутах для построения карты, укажите 0.1, 0.2 или 0.4 для размеров карты 30'x30', 60'x60' или 120'x120' соответственно") 
args = parser.parse_args()

coord_in = args.infile
clustername = args.name
step = args.step
delta = args.delta
maglim = args.maglim

xymag_columnsnum = (0, 1, 13)
xymag = np.loadtxt(coord_in, comments='#', usecols=xymag_columnsnum)

size = "{0}x{0}".format(int(300*step))
name = 'map_{3}_{0}_{1}_{2}.'.format(size, maglim, delta, clustername)
sname = "Карта поверхностной плотности {0}:\nmaglim = {1}, h = {2}'".format(clustername, maglim, delta)
xydens = make_map(xymag, delta, step, maglim, name+'txt')
plot_map(xydens, name+'txt', name+'png', sname, show=True)
