#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
in_file = sys.argv[1]
a = np.loadtxt(in_file, unpack=True, usecols=(-5, 3))
jmh = a[0]
j = a[1]
#fig = plt.figure(figsize=(10,10))
plt.scatter(jmh, j, s=2, color='k')
ax = plt.gca()
ax.set_xlabel('J-H')
ax.set_ylabel('J')
ax.invert_yaxis()
plt.grid(True)
plt.show()
#fig.savefig('j-h_j.png', dpi=100)
