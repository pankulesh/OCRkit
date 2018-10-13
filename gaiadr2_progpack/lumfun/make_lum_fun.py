#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from multiprocessing import Pool, cpu_count
import argparse
from sys import exit

class Region:
    
    def __init__(self, x, y, r, name, t='circle'):
        self.x0 = x
        self.y0 = y
        self.r = r
        self.name = name
        self.t = t

    def cut_region(self, array_xyz):
        stars = []
        xs, ys, zs = array_xyz.T
        dists = np.sqrt((xs-self.x0)**2+(ys-self.y0)**2)
        for star, dist in zip(array_xyz, dists):
            if self.t == 'circle':
                if dist <= self.r:
                    stars.append(star)
            else:
                if self.r <= dist <= self.r*np.sqrt(2):
                    stars.append(star)
        self.region = np.asarray(stars)
        self.values = self.region.T[2]
        self.n = len(stars)
        return

    def get_optimal_bw(self, mindelta=0.1, maxdelta=1.0, numberofdeltas=10):
        from sklearn.model_selection import GridSearchCV
        grid = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(mindelta, maxdelta, numberofdeltas)}, n_jobs=-1)
        grid.fit(self.values[:, None])
        return grid.best_params_['bandwidth']
    
    def bootstrap(self, argument, kde, bandwidth, nboot=24):
        self.density = get_estimated_density(self.n, argument, kde)
        s1, s2 = (np.zeros_like(argument) for _ in range(2))
        workers = Pool(cpu_count())
        densities = workers.map(bootstrap0, [(self.n, kde, bandwidth, argument)]*nboot)
#без multiprocessing'а
#        for _ in range(nboot):
#            newstars = kde.sample(n_samples=self.n).T[0]
#            newkde = get_kernel_estimator(newstars, bandwidth)
#            newvalues = get_estimated_density(argument, newkde)
#            s1 += newvalues
#            s2 += newvalues*newvalues
        for newvalues in densities:
            s1 += newvalues
            s2 += newvalues*newvalues
        self.d_mean = s1/nboot
        self.d_disp = np.sqrt((s2-nboot*self.d_mean**2)/(nboot-1))
        print(f'Успешно посчитана дисперсия для области {self.name}')


def get_kernel_estimator(values, bandwidth, kerneltype='gaussian'):
    return KernelDensity(kernel=kerneltype, bandwidth=bandwidth).fit(values[:, None])

def get_estimated_density(n, argument, kde):
    return np.exp(kde.score_samples(argument[:, None]))*n

def bootstrap0(t):
    n, kde, bandwidth, argument = t 
    newstars = kde.sample(n_samples=n).T[0]
    newkde = get_kernel_estimator(newstars, bandwidth)
    newvalues = get_estimated_density(n, argument, newkde)
    return newvalues

def plot_lf(kind, argument, density, disp):
    plt.title("Функция блеска скопления {}\nдля области сравнения {}".format(clustername, kind), fontsize=12)
    arg = argument
    d = density
    down = density-disp
    up = density+disp
    plt.plot(arg, d, 'k')
    plt.plot(arg, down, 'k', linestyle='dashed', linewidth=1)
    plt.plot(arg, up, 'k', linestyle='dashed', linewidth=1)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    ax.set_xlabel(r'Видимая звёздная величина в полосе {}'.format(mag_type), fontsize=15)
    ax.set_ylabel(r'Функция блеска', fontsize=15)
    ax.set_xlim(mag0, maglim)
    ax.invert_xaxis()
    ax.set_ylim(np.amin(down), np.amax(up))
    ax.tick_params(width=2)
    ax.grid(alpha=0.2, linestyle='dashed', linewidth=0.5)
    plt.savefig('lf_{}_{}.png'.format(clustername, kind), dpi=100, transparent=False)
    plt.close()

def write_line(i, argument, density, disp):
    return '{0:.3f}\t{1:.6f}\t{2:.6f}\t{3:.6f}\n'.format(argument[i], density[i], density[i]-disp[i], density[i]+disp[i])

def write_file(filename, argument, density, disp):
    with open(filename, 'w') as f:
        for i in range(len(argument)):
            f.write(write_line(i, argument, density, disp))
        print(f'Файл {filename} успешно записан.')

def check_and_fix_overlap(x0, y0, rc, x, y):
    dist = np.sqrt((x0-x)**2+(y0-y)**2)
    defect = 2*rc - dist
    if defect > 0.01:
        xnew = x + defect/dist*(x-x0)
        ynew = y + defect/dist*(y-y0)
        print(f"Круговая область с центром ({x:.2f}, {y:.2f}) слишком близко к области скопления, смещаем в ({xnew:.2f}, {ynew:.2f})...")
        x, y = xnew, ynew
    return (x, y)

parser = argparse.ArgumentParser(description='Программа, рассчитывающая функции блеска для области скопления с помощью метода KDE')
parser.add_argument("infile", type=str, help='Имя входного файла')
parser.add_argument("clustername", type=str, help='Название скопления')
parser.add_argument("parameters", type=str, help='Имя файла с параметрами в строку: step, mag0, maglim, x0, y0, rc, delta, [x1, y1] [x2, y2] [x3, y3] [x4, y4]')
parser.add_argument('-o', '--optimal', action='store_true', help='Укажите этот флаг, если нужно рассчитать оптимальный параметр KDE методом кросс-валидации (bandwidth CV)')
parser.add_argument('-c', '--circles', action='store_true', help='Укажите этот флаг, если нужно рассчитать ФБ с областями сравнения в виде кругов вокруг скопления')
parser.add_argument('-s', '--symm', action='store_true', help='Укажите этот флаг, если в файле с параметрами записаны координаты центра только одной круговой области сравнения — тогда остальные три рассчитаются симметрично относительно центра и запишутся в исходный файл с параметрами.')
args = parser.parse_args()
clustername = args.clustername
parameters = np.loadtxt(args.parameters, comments='#')

xymag = (0, 1, 13)
mag_type = 'G'
step, mag0, maglim, x0, y0, rc, delta = parameters[:7]

if args.circles:
    x1, y1 = parameters[7:9]
    #check circles overlapping and fix it
    x1, y1 = check_and_fix_overlap(x0, y0, rc, x1, y1)
    if args.symm:
        """
             x4y4
        x1y1 x0y0 x2y2
             x3y3

        """

        x2 = 2*x0-x1
        y2 = 2*y0-y1

        x3 = -y1+y0+x0
        y3 = x1-x0+y0

        x4 = 2*x0-x3
        y4 = 2*y0-y3
    else:
        try:
            x2, y2, x3, y3, x4, y4 = parameters[9:]
        except ValueError:
            print(f"Не хватает значений центров круговых областей в {args.parameters}! Запустите скрипт с флагом -s или исправьте файл.")
            exit(2)
        x2, y2 = check_and_fix_overlap(x0, y0, rc, x2, y2)
        x3, y3 = check_and_fix_overlap(x0, y0, rc, x3, y3)
        x4, y4 = check_and_fix_overlap(x0, y0, rc, x4, y4)
    with open(args.parameters, 'w') as f:
        f.write('#step mag0 maglim x0 y0 rc delta x1    y1    x2    y2    x3    y3    x4    y4   \n')
        f.write('{:.2f} {:.1f} {:.1f} {:.2f} {:.2f} {:.2f} {:.1f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'.format(step, mag0, maglim, x0, y0, rc, delta, x1, y1, x2, y2, x3, y3, x4, y4))

    xsys = [(x1, y1, 'Круг 1'), (x2, y2, 'Круг 2'), (x3, y3, 'Круг 3'), (x4, y4, 'Круг 4')]
    
    lf_circles = [Region(x, y, rc, name) for x,y,name in xsys]

lf_cluster = Region(x0, y0, rc, 'Скопление+фон')
lf_ring = Region(x0, y0, rc, 'Кольцо', t='ring')

input_stars = np.loadtxt(args.infile, comments='#', usecols=xymag)

bounds_num = int((maglim-mag0)/step)
argument = np.linspace(mag0, maglim, bounds_num)

a = [lf_cluster, lf_ring]

if args.circles:
    a += lf_circles

for i in a:
    i.cut_region(input_stars)
if args.optimal:
    delta = lf_cluster.get_optimal_bw()
    print("Оптимальная полуширина: ",delta)
for i in a:
    kde = get_kernel_estimator(i.values, delta)
    i.bootstrap(argument, kde, delta)

print("Область скопления: {0} звёзд,\nКольцевая область сравнения: {1} звёзд, число узлов сетки {2}".format(lf_cluster.n, lf_ring.n, bounds_num))

if args.circles:
    print("Круговые области сравнения: {}, {}, {}, {} звёзд соответственно".format(lf_circles[0].n, lf_circles[1].n, lf_circles[2].n, lf_circles[3].n))

plot_lf('cluster+back', argument, lf_cluster.density, lf_cluster.d_disp)
write_file('lf_{}_cluster+back.txt'.format(clustername), argument, lf_cluster.density, lf_cluster.d_disp)

def plot_and_write_lf_back(name, density, disp):
    write_file('lf_{}_background_{}.txt'.format(clustername, name), argument, density, disp)
    plot_lf('background_{}'.format(name), argument, density, disp)
    densityd = lf_cluster.density - density
    dd_disp = np.sqrt(lf_cluster.d_disp**2+disp**2)
    write_file('lf_{}_{}.txt'.format(clustername, name), argument, densityd, dd_disp)
    plot_lf(name, argument, densityd, dd_disp)

plot_and_write_lf_back('r', lf_ring.density, lf_ring.d_disp)

if args.circles:
    densityb_2c, db_disp_2c, densityb_4c, db_disp_4c, densityd_2c, dd_disp_2c, densityd_4c, dd_disp_4c = [np.zeros(bounds_num) for _ in range(8)]

    #Averaging 2 circles and 4 circles regions
    densityb_2c = 0.5*(lf_circles[0].density+lf_circles[1].density)
    db_disp_2c = 0.5*np.sqrt(lf_circles[0].d_disp**2+lf_circles[1].d_disp**2)
    densityb_4c = 0.25*(lf_circles[0].density+lf_circles[1].density+lf_circles[2].density+lf_circles[3].density)
    db_disp_4c = 0.25*np.sqrt(lf_circles[0].d_disp**2+lf_circles[1].d_disp**2+lf_circles[2].d_disp**2+lf_circles[3].d_disp**2)

    plot_and_write_lf_back('2c', densityb_2c, db_disp_2c)
    plot_and_write_lf_back('4c', densityb_4c, db_disp_4c)
