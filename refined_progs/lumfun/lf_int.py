#!/usr/bin/python
import numpy as np
import argparse
from scipy.integrate import trapz, simps

parser = argparse.ArgumentParser(description='Integrate luminocity function and ')
parser.add_argument('in_file', type=str, help='Input file name (lf.txt)')
parser.add_argument("-o", "--outfile", default="lf_nums_mass.txt", help='Output filename, default lf_nums_mass.txt')
#parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
#parser.add_argument('h', type=str, help="Kernel halfwidth (as in filenames)")
args = parser.parse_args()

in_file = args.in_file
out_file = args.outfile

lumfun = np.loadtxt(in_file, unpack=True)
mag, lf, lf_lo, lf_up = lumfun
numbers = ["Number of stars from LF: ", 
           "Lower bound:             ", 
           "Upper bound:             ",
           ]
z = 0

with open(out_file, 'w') as f:
    f.write("                         COUNTED WITH SIMPSON'S METHOD COUNTED WITH TRAPEZOIDS\n")

    for x, y in ([mag, lf], [mag, lf_lo], [mag, lf_up]):
        I1s = simps(y, x)
        I1t = trapz(y, x)
        f.write('{0}{I1s:29.2f} {I1t:23.2f}\n'.format(numbers[z], **locals()))
        z +=1

with open(out_file, 'r') as f:
    print(f.read())
