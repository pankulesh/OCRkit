#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description='Программа по трёхфакторному отбору звёзд скопления.')
parser.add_argument("infile", help='Имя или путь к исходному файлу')
parser.add_argument("outfile",  help="Имя выходного файла")
parser.add_argument("-i","--lighterinput",  default='', help='Введите тут исходный файл для построения графиков, если он не совпадает с infile')
args = parser.parse_args()
infile = args.infile
outfile = args.outfile
if not args.lighterinput:
    initial_array = np.loadtxt(infile, comments='#', usecols=(7,9,11,13))
else:
    initial_array = np.loadtxt(args.lighterinput, comments='#', usecols=(7,9,11,13))
init_stars = np.asarray([[i, *item] for i, item in enumerate(initial_array)])

def clear_bads(stars):
    a = []
    for star in stars:
        if star[1]==99.9999 or star[2]==99.999 or star[3]==99.999:
            continue
        else:
            a.append(star)
    return np.array(a)

def clear_mildly(stars, t, center, s):
    a = []
    for star in stars:
        if abs(star[t]-center)<=s and star[4]<=18:
            a.append(star)
    return np.array(a)

def is_inside(linestar, crit):
    kind, center, halfsize = crit
    words = linestar.split()
    p, ma, md, mag = [float(i) for i in (words[7], words[9], words[11], words[13])]
    if p == 99.9999 or ma == 99.999 or md == 99.999 or mag > 18:
        return False
    if kind == 1:
        return abs(p-center) < halfsize
    elif kind == 2:
        return abs(ma-center) < halfsize
    elif kind == 3:
        return abs(md-center) < halfsize
    else:
        raise ValueError

def write_new_stars(criteria, infile, outfile):
    with open(infile, 'r') as f:
        with open(outfile, 'w') as g:
            counter = 0
            for line in f:
                if line.startswith('#'):
                    g.write(line)
                else:
                    if is_inside(line, criteria):
                        g.write(line)
                        counter += 1
                        print(counter)
    print(f"{outfile} записан успешно, {counter} звёзд отобрано.")

def plot_padm(stars):
    par, mua, mud, mag = stars.T[1:]
    fig, axes = plt.subplots(1,3, figsize=(9,3))
    params = dict(s=0.1, c='k', marker='.')
    axes[0].scatter(mag, par, **params)
    axes[0].set_xlabel('Gmag')
    axes[0].set_title('π')
    axes[1].scatter(mag, mua, **params)
    axes[1].set_xlabel('Gmag')
    axes[1].set_title('μ_α')
    axes[2].scatter(mag, mud, **params)
    axes[2].set_xlabel('Gmag')
    axes[2].set_title('μ_δ')
    plt.show()

newstars = clear_bads(init_stars)
plot_padm(newstars)
#write_new_stars(a, infile, outfile)
wannaex = False
filters = []
while not wannaex:
    t = int(input("По какому графику отбирать? 1, 2, 3?\n"))
    c = float(input("А центральная линия скопления на какой высоте?\n"))
    s = float(input("А полуширина диапазона какая?\n"))
    newstars = clear_mildly(newstars, t, c, s)
    filters.append((t,c,s))
    plot_padm(newstars)
    a = input("Ещё раз? y/n\n")
    if a == 'n':
        wannaex = True
    else:
        print('Продолжаем')
print(filters)
n = int(input("По какому отбору записать новый файл? 0, 1, 2…\n"))
write_new_stars(filters[n], infile, outfile)
