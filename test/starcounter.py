#!/usr/bin/python
import numpy as np
import argparse
import time
from scipy import interpolate, integrate
parser = argparse.ArgumentParser(description='Integrate all profiles, need files "parameters_NAME_h=H.txt" and all "density_NAME_J_H.txt" to be in place, i.e. : "python starcounter.py NGC7142 2"')
parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
parser.add_argument('h', type=str, help="Kernel halfwidth (as in filenames)")
args = parser.parse_args()
h = args.h
name = args.name
mags = np.loadtxt("parameters_{0}_h={1}.txt".format(name, h), dtype='int', unpack=True, skiprows=1, usecols=(0,))
params = np.loadtxt("parameters_{0}_h={1}.txt".format(name, h), skiprows=1, usecols=(0, 1, 3, 4))
def integrateprofile(m):    
    profile = np.loadtxt("density_{0}_{1}_{2}.txt".format(name, m, h), usecols=(0, 1, 3, 4), unpack=True)
    for i in params:
        if int(i[0]) == m:
            fb = i[2]
            dfb = i[3]
            radius = i[1]
            break
#    for i in range(len(profile)):
#        if profile[i][0] > radius:
#            bord = i
#    profile = np.transpose(profile[:bord])
    x = profile[0]
    y = 2*np.pi*x*profile[1]
    ymin = 2*np.pi*x*profile[2]
    ymax = 2*np.pi*x*profile[3]
    yspl = interpolate.CubicSpline(x, y)
    yminspl = interpolate.CubicSpline(x, ymin)
    ymaxspl = interpolate.CubicSpline(x, ymax)
    bkg = np.pi*fb*radius**2
    sigmabkg = np.pi*dfb*radius**2
    n1 = yspl.integrate(0, radius)
#    print('n1=', n1-bkg)
    n2 = yminspl.integrate(0, radius)
#    print('n2=', n2-bkg+sigmabkg)
    n3 = ymaxspl.integrate(0, radius)
#    print('n3=', n3-bkg-sigmabkg)
    sigmaint = max(abs(n2-n1), abs(n3-n1))
    sigmanc = np.sqrt(sigmaint**2+sigmabkg**2)
    return (n1-bkg, sigmanc, sigmaint, sigmabkg)
with open("numstars_{0}_h={1}.txt".format(name, h), 'w') as f:
    f.write("{:>3} {:>10} {:>10} {:>10} {:>10}\n".format("mag", "n_stars", "sigma_n", "sigma_int", "sigma_back"))
for mag in mags:
    num, senum, seint, sebkg = integrateprofile(mag)
    with open("numstars_{0}_h={1}.txt".format(name, h), 'a') as f:
        f.write("{0:>3d} {1:10.2f} {2:10.2f} {3:10.2f} {4:10.2f}\n".format(mag, num, senum, seint, sebkg))
    print("density_{0}_{1}_{2}.txt was counted.".format(name, mag, h))
