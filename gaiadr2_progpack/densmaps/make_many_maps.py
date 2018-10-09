#!/usr/bin/python
import argparse
import numpy as np
from kde_map_functions import make_map, plot_map
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description='Программа для построения множества карт поверхностной плотности с определенными параметрами KDE и предельными звёздными величинами')
parser.add_argument("infile", help='Название входного файла с координатами и звёздными величинами')
parser.add_argument("name", help="Имя скопления (для заголовков в картах)")
parser.add_argument("step", type=float, help="Размер шага в угловых минутах для построения карты, укажите 0.1, 0.2 или 0.4 для размеров карты 30'x30', 60'x60' или 120'x120' соответственно") 
parser.add_argument("-c", "--combine", action='store_true', help="Укажите этот флаг для объединения всех карт в одну большую картинку")

args = parser.parse_args()
coord_in = args.infile
clustername = args.name
step = args.step
size = "{0}x{0}".format(int(300*step))

#Эти параметры можно задавать вручную, в итоге получится len(maxmags)*len(deltas) карт
maxmags = [14, 15 ,16, 17, 18]
deltas = [4, 5, 6, 7, 8]

#В колонках с этими номерами должны находиться соответсвующие величины
xymag_columns = (0, 1, 13)

print("""Выходные файлы будут выглядеть следующим образом: map_<имя скопления>_<размеры карты>_<предельная величина>_<параметр KDE>.[txt, png]
    Комбинированная карта (если требуется) будет в файле maps_<имя скопления>_<размеры карты>.png
    Построение займёт некоторое время (возможно даже много времени)...\n""")

xymag = np.loadtxt(coord_in, comments='#', usecols=xymag_columns)

for j in maxmags:
    for k in deltas:
        name = 'map_{3}_{0}_{1}_{2}.'.format(size, j, k, clustername)
        sname = "Карта поверхностной плотности {0}:\nmaglim = {1}m, h = {2}'".format(clustername, j, k)
        xydens = make_map(xymag, k, step, j, name+'txt')
        plot_map(xydens, name+'txt', name+'png', sname)

if args.combine:
    fig, axes = plt.subplots(len(maxmags), len(deltas), subplot_kw={'xticks': [], 'yticks': []})
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    img = []
    for i in maxmags:
        for j in deltas:
            img.append(mpimg.imread('map_{3}_{0}_{1}_{2}.png'.format(size, i, j, clustername)))
    for ax, image in zip(axes.flat, img):
        ax.axis('off')
        ax.imshow(image)
    fig.savefig('maps_{}_{}.png'.format(clustername, size), dpi=300, frameon=False)
