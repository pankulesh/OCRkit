#!/usr/bin/python
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from kde_map_functions import plot_map
from scipy.interpolate import griddata

parser = argparse.ArgumentParser(description='Программа для построение карты поверхностной плотности по готовому файлу map_....txt')
parser.add_argument("infile", type=str, help='Имя входного файла с x y z колонками')
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
plot_map(in_file, out_file, sname)
