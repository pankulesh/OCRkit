#!/usr/bin/python
from math import cos, pi
help = """alpha delta - coordinates of cluster's centre,
    m - max magnitude, 
    h - kernel halfwidth for x/y,
    s - step for plot in x/y,
    if centre's area is symmetric relative to (0, 0), enter first variant coordinates, string x/y will be like '-s0  -s1 s1  s0'
    if area is assymetric, enter second variant (in parameters.txt coordinates will be like 'x0 x1 x2 x3' and 'y0 y1 y2 y3')"""
print(help)
ins = input("Enter 'alpha delta m h s s0 s1' or 'alpha delta m h s x0 x1 x2 x3 y0 y1 y2 y3': \n").split()
if len(ins) == 7:
    alpha, delta, m, h, s, x3, x2 = (float(x) for x in ins)
    x0, x1 = -1*x3, -1*x2
    y0, y1, y2, y3 = x0, x1, x2, x3
elif len(ins) == 13:
    alpha, delta, m, h, s, x0, x1, x2, x3, y0, y1, y2, y3 = (float(x) for x in ins)
else:
    print("wrong input, stopping...")
    exit(2)
k=1/60/cos(pi*delta/180)
with open('parameters.txt', 'w') as f:
    f.write('{:.1f}\n'.format(m))
    f.write('{0:.1f}\t{1:.4f}\t{2:.4f}\n'.format(h, h/60, k*h))
    f.write('{0:.2f}\t{1:.6f}\t{2:.6f}\n'.format(s, s/60, k*s))
    f.write('{0:.1f}\t{1:.1f}\t{2:.1f}\t{3:.1f}\n'.format(x0, x1, x2, x3))
    f.write('{0:.1f}\t{1:.1f}\t{2:.1f}\t{3:.1f}\n'.format(y0, y1, y2, y3))
    f.write('{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n'.format(delta+x0/60, delta+x1/60, delta+x2/60, delta+x3/60))
    f.write('{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n'.format(alpha+k*y0, alpha+k*y1, alpha+k*y2, alpha+k*y3))
