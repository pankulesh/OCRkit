#!/usr/bin/python
import numpy as np
from kde_functions import make_imin_imax, density_estimator_1D as de
import argparse
import matplotlib.pyplot as plt

class LF_region:
    
    def __init__(self):
        self.density = np.zeros(ibound2)
        self.d_mean = np.zeros(ibound)
        self.d_disp = np.zeros(ibound)
        self.n_tot = 0 

def plot_lf(lf):
    plt.plot(lf.argument, lf.density, label='{}'.format(lf.delta))
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    ax.set_xlabel(r'Звёздная величина в полосе {}'.format(magtype), fontsize=18)
    ax.set_ylabel(r'Функция блеска', fontsize=18)
    ax.tick_params(width=2)


parser = argparse.ArgumentParser(description='Программа для выбора оптимального параметра для KDE функции блеска.')
parser.add_argument("infile", type=str, help='Имя входного файла')
parser.add_argument("parameters", type=str, help='Файл с параметрами в строку: step, mag0, maglim, x0, y0, rc')

args = parser.parse_args()

magxy = (13, 0, 1)
mag_type = 'G' 

input_region = np.loadtxt(args.infile, comments='#', usecols=magxy)
parameters = np.loadtxt(args.parameters, comments='#')
step, mag0, maglim, x0, y0, rc = parameters[:6]

#Массив с предполагаемыми значениями параметра
deltas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]    

lfs = []

for delta in deltas:
    ibound = int((maglim-mag0)/step)
    ibound2 = int((maglim-mag0+delta)/step)
    lfs.append(LF_region())
    lfs[-1].delta = delta
    lfs[-1].argument = np.array([mag0+i*step for i in range(ibound2)])
    lf_cluster = lfs[-1]
    for mag, x, y in input_region:
        if mag > maglim+delta or mag < mag0-delta:
            continue
        rstar = mag-mag0
        imin, imax = make_imin_imax(rstar, delta, step, ibound)
        if imin <= ibound and imax >= 0:
            r = np.sqrt((x-x0)**2+(y-y0)**2)
            if r < rc:
                lf_cluster.n_tot+=1
                lf_cluster.density[imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]

    for i in range(ibound2):
        if i >= ibound:
            lf_cluster.density[i] = lf_cluster.density[ibound-1]
    print('{} ready'.format(delta))

plt.title("ФБ с разными полуширинами ядра", fontsize=26)

for lf in lfs:
    plot_lf(lf)
plt.legend()
plt.show()
