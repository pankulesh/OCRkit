#!/usr/bin/python
import numpy as np
import argparse
from scipy.integrate import trapz, simps

parser = argparse.ArgumentParser(description='Count giants addition in number of stars and mass of open cluster with its LF.')
parser.add_argument('in_file', type=str, help='Input file name (lf.txt) with only red giants LF section')
parser.add_argument('mass', type=float, help="Mean giant's mass in solar mass")
args = parser.parse_args()

in_file = args.in_file
mass = args.mass
lumfun = np.loadtxt(in_file, unpack=True, usecols=(0,1,2,3))
mag, lf, lf_lo, lf_up = lumfun
nums, numt = simps(lf, mag), trapz(lf, mag)
Ms, Mt = mass*nums, mass*numt
nums_l, numt_l = simps(lf_lo, mag), trapz(lf_lo, mag)
Ms_l, Mt_l = mass*nums_l, mass*numt_l
nums_u, numt_u = simps(lf_up, mag), trapz(lf_up, mag)
Ms_u, Mt_u = mass*nums_u, mass*numt_u

print("                         COUNTED WITH SIMPSON'S METHOD COUNTED WITH TRAPEZOIDS")
print('Number of stars:         {0:29.2f} {1:23.2f}'.format(nums, numt))
print('Lower bound:             {0:29.2f} {1:23.2f}'.format(nums_l, numt_l))
print('Upper bound:             {0:29.2f} {1:23.2f}'.format(nums_u, numt_u))
print('Mass in solar mass:      {0:29.2f} {1:23.2f}'.format(Ms, Mt))
print('Lower bound:             {0:29.2f} {1:23.2f}'.format(Ms_l, Mt_l))
print('Upper bound:             {0:29.2f} {1:23.2f}'.format(Ms_u, Mt_u))
