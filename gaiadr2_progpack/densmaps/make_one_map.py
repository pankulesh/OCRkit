#!/usr/bin/python
import argparse
import pandas as pd
import numpy as np
from kde_map_functions import plot_map, make_map

def get_xymag(df, namefile='columnsxygmagnames.txt'):
    with open(namefile, 'r') as f:
        columnsnames = f.readline()
    names = dict(zip(['x', 'y', 'Gmag'],columnsnames.split(',')))
    try:
        xymag = np.array([df[names['x']].values, df[names['y']].values, df[names['Gmag']].values]).T
    except:
        columnsnames = input("Введите имена колонок со значениями x, y и G соответственно через запятую (например, x,y,Gmag):\n")
        with open(namefile, 'w') as f:
            f.write(columnsnames)
        names = dict(zip(['x', 'y', 'Gmag'],columnsnames.split(',')))
        xymag = np.array([df[names['x']].values, df[names['y']].values, df[names['Gmag']].values]).T
    return xymag

def make_and_plot_map(xymag, clustername, delta, maglim, step):
    size = "{0}x{0}".format(int(300*step))
    name = 'map_{3}_{0}_{1}_{2}.'.format(size, maglim, delta, clustername)
    sname = "Карта поверхностной плотности {0}:\nmaglim = {1}, h = {2}'".format(clustername, maglim, delta)
    make_map(xymag, delta, step, maglim, name+'txt')
    plot_map(name+'txt', name+'png', sname, show=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Программа для построение карты поверхностной плотности конкретного размера, параметра KDE и предельной звёздной величины')
    parser.add_argument("infile", type=str, help='Имя входного файла, содержащего координаты и звёздные величины звёзд (заголовки колонок без комментирования!)')
    parser.add_argument("name", type=str, help="Название скопления (для заголовков графика)")
    parser.add_argument("delta", type=float, help="Величина параметра KDE в угловых минутах")
    parser.add_argument("maglim", type=float, help="Предельная звёздная величина звёзд, используемых в построении карты")
    parser.add_argument("step", type=float, help="Размер шага в угловых минутах для построения карты, укажите 0.1, 0.2 или 0.4 для размеров карты 30'x30', 60'x60' или 120'x120' соответственно") 
    args = parser.parse_args()

    infile = args.infile
    clustername = args.name
    step = args.step
    delta = args.delta
    maglim = args.maglim
    df = pd.read_csv(infile, delim_whitespace=True)
    xymag = get_xymag(df)
    make_and_plot_map(xymag, clustername, delta, maglim, step)
