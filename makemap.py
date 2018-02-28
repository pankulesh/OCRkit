#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import argparse
from math import cos,pi
from scipy.interpolate import griddata
parser = argparse.ArgumentParser(description='Make a 2D contour map from XYZ textfile.')
parser.add_argument("-i", "--ifile", default='densgrid.txt', help='input file name, default densgrid.txt')
parser.add_argument("-o", "--ofile", default='densmap.png', help='output file name, default densmap.png')
parser.add_argument("-n", "--name", default='off', help="set map's title, no title by default")
parser.add_argument("-c", "--no-clabel", action='store_true', help='switch off labels on isolines')
parser.add_argument("-b", "--no-colorbar", action='store_true', help='switch off colorbar')
parser.add_argument("-g", "--no-grid", action='store_true', help='switch off additional grid')
parser.add_argument("-z", "--saved-size", default='6', help='size of saved map in inches, default 6x6')
parser.add_argument("-s", "--show-map", action='store_true', help='show map on screen')
parser.add_argument("-e", "--equatorial-center", default='off', help='make map in equatorial coordinates, needs file with string with center in alpha and delta, like "169.365 -62.740"')
parser.add_argument("-m", "--specify-new-center-file", default='off', help='read file name with clarified centre and show it on map, disabled by default')
args = parser.parse_args()
filename = args.ifile
outname = args.ofile
sizesaved = float(args.saved_size)
internumber = 100
levels = 10
a = np.loadtxt(filename, unpack=True)
clcenter = []
if args.specify_new_center_file!='off':
    try:
        clcenter = np.loadtxt(args.specify_new_center_file, unpack=True, usecols=(0,1))
    except:
        print("Wrong clarified centre file input!")
        exit(2)
X,Y,Z=a[0],a[1],a[2]
if args.equatorial_center == 'off':
    xi = np.linspace(min(X), max(X), internumber)                       #
    yi = np.linspace(min(Y), max(Y), internumber)
    zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='linear')
    fig = plt.figure(figsize=(sizesaved,sizesaved))
    contour = plt.contour(xi, yi, zi, levels, colors='k', zorder=2)
    contour_filled = plt.contourf(xi, yi, zi, levels, zorder=1, cmap=plt.cm.Greys)
    ax = plt.gca()
    if clcenter!=[]:
        dot = plt.scatter(clcenter[0][0], clcenter[0][1], s=15, color='w', zorder=10)
        text = plt.annotate('({0:.1f}, {1:.1f})'.format(clcenter[0][0], clcenter[0][1]), xy=(clcenter[0]), fontsize=10, xytext=(clcenter[0][0]+0.5, clcenter[0][1]+0.5))
        text.set_bbox(dict(facecolor='white', edgecolor='None', alpha=0.6))
    ax.set_xlim([min(X), max(X)])
    ax.set_ylim([min(Y), max(Y)])
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    plt.xlabel("x, arcmin")
    plt.ylabel("y, arcmin")
else:
    try: 
        center = np.loadtxt(args.equatorial_center)
    except:
        print("Wrong equatorial centre input!")
        exit(2)
    DEC = np.zeros((len(X)))
    RA = np.zeros((len(Y)))
    for i in range(len(X)):
        DEC[i] = X[i]/60+center[1]
#        RA[i] = -Y[i]/60/cos(DEC[i]*pi/180)+center[0]
        RA[i] = Y[i]/60 + center[0]
    deci = np.linspace(min(DEC), max(DEC), internumber)
    rai = np.linspace(min(RA), max(RA), internumber)
#    overoff = (max(rai)-min(rai)-0.5)/2
    dticks = np.linspace(min(deci), max(deci), 7)
#    aticks = np.linspace(min(rai)+overoff, max(rai)-overoff, 7)
    aticks = np.linspace(min(rai), max(rai), 7)
    zi = griddata((RA, DEC), Z, (rai[None,:], deci[:,None]), method='linear')
    fig = plt.figure(figsize=(sizesaved,sizesaved)) # создать фигуру размером sizesavedxsizesaved дюймов
    contour = plt.contour(rai, deci, zi, levels, colors='k', zorder=2) # построить контурную карту
    contour_filled = plt.contourf(rai, deci, zi, levels, zorder=1, cmap=plt.cm.Greys) # заполнить контуры с цветовой палитрой cmap
    ax = plt.gca()
    if clcenter!=[]:
        dot = plt.scatter(clcenter[1][0], clcenter[1][1], s=15, color='w', zorder=10)
        text = plt.annotate('({0:.3f}, {1:.3f})'.format(clcenter[1][0], clcenter[1][1]), xy=(clcenter[1]), fontsize=10, xytext=(clcenter[1][0]+0.01, clcenter[1][1]+0.01), zorder=10)
        text.set_bbox(dict(facecolor='white', edgecolor='None', alpha=0.6))    
    ax.set_xticks(aticks)
    ax.set_xlim([min(aticks), max(aticks)])
    ax.set_yticks(dticks)
    ax.set_ylim([min(dticks), max(dticks)])
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    plt.xlabel("RA, degrees")
    plt.ylabel("DEC, degrees")
if not args.no_clabel:
    clabels = plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=8)
    for clabel in clabels:
        clabel.set_bbox(dict(facecolor='white', edgecolor='None', alpha=0.3))    
if not args.no_colorbar:
    plt.colorbar(contour_filled) # включить цветовую легенду (полоска сбоку)
if not args.no_grid:
    plt.grid(True) # включить вспомогательную сетку
plt.axes().set_aspect('equal') # сделать равный масштаб по осям
if args.name!='off':
    plt.title(args.name) # установить заголовок графика
fig.savefig(outname, dpi=100, frameon=False) # сохранить график под именем outname, разрешением 200 и без рамки
if args.show_map:
    plt.show() # показать график на экране
