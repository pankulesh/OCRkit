#!/usr/bin/python
import sys
coord_in = sys.argv[1]
coord_out = sys.argv[2]
with open(coord_in, 'r') as f:
    line = f.readline()
    with open(coord_out, 'w') as g:
        while line:
            p = line.split()
            if float(p[4]) < 9. and float(p[6]) < 9. and float(p[8]) < 9.:
                g.write(line)
        line = f.readline()
