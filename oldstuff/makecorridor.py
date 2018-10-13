#!/usr/bin/python
import sys
import numpy as np
from math import log10
iso_in = sys.argv[1]
cor_out = sys.argv[2]
#with open(param_in, 'r') as f:
#    f.readline()
#    param = f.readline().split()
#param = input("Enter r E(B-V) halfwidth in one string: ").split()
param = ['1240', '0.34', '0.2']
r, ebv, halfwidth = (float(x) for x in param)
skipupper = 1
ebv2ejh = 0.34        # magic constants
ebv2aj = 2.43*ebv2ejh #
pr = 5*log10(r) - 5
fixcolor = ebv2ejh*ebv
fixmag = pr + ebv2aj*ebv
dtype = [('left', float), ('right', float), ('j', float)]
corarr = np.zeros((500,), dtype=dtype)
arrin = np.genfromtxt(iso_in, skip_header=skipupper, delimiter="\t", usecols=(13, 14))
i = 0
for j, h in arrin:
    mj = j + fixmag
    if mj > 16.0:
        continue
    jmh = j - h + fixcolor
    corarr[i] = (jmh - halfwidth, jmh + halfwidth, mj)
    i += 1
cor = np.sort(corarr[:i], order='j')
np.savetxt(cor_out, cor, fmt='%.3f', delimiter='\t', newline='\n') 
