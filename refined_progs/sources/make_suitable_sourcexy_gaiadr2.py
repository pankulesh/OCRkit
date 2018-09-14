#!/usr/bin/python
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Программа для подсчёта x, y, r из исходных файлов cdsclient. Пустые поля в BP и RP заполняются 99.9999.')
parser.add_argument("infile", help='Имя входного файла')
parser.add_argument("outfile",  help="Имя выходного файла")
args = parser.parse_args()
infile = args.infile
outfile = args.outfile
def countxy(alpha0, delta0, alpha, delta):
    #http://cdsarc.u-strasbg.fr/vizier/vizHelp/xy.htx
    #rotational matrix:
    a0, d0 = alpha0*np.pi/180, delta0*np.pi/180
    a, d = alpha*np.pi/180, delta*np.pi/180
    A = [[np.cos(d0)*np.cos(a0), np.cos(d0)*np.sin(a0), np.sin(d0)],
        [-1*np.sin(a0), np.cos(a0), 0],
        [-1*np.sin(d0)*np.cos(a0), -1*np.sin(d0)*np.sin(a0), np.cos(d0)]]
    #create (u, v, w) vector:
    uvw_unrot = [np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)]
    u, v, w = np.matmul(A, uvw_unrot)
    r = np.arccos(u)*180/np.pi*60
    x = r*v/np.sqrt(v**2+w**2)
    y = r*w/np.sqrt(v**2+w**2)
    return (x, y)


with open(outfile, 'w') as g:
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith("#-c="):
                alpha0, delta0 = line[4:].split()
                g.write(line) 
                g.write("#{:>9} {:>9} {:>15} {:>9} {:>15} {:>9} {:>19} {:>6} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>7} {:>7} {:>7}\n".format(*list(range(19))))
                g.write("#{:>9} {:>9} {:>15} {:>9} {:>15} {:>9} {:>19} {:>6} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>7} {:>7} {:>7}\n".format("x", "y", "RA_ICRS", "e_RA_ICRS", "DE_ICRS", "e_DE_ICRS", "Source",  "Plx", "e_Plx", "pmRA", "e_pmRA",  "pmDE", "e_pmDE",  "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag"))
                g.write("#{:>9} {:>9} {:>15} {:>9} {:>15} {:>9} {:>19} {:>6} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>7} {:>7} {:>7}\n".format("arcmin", "arcmin", "deg", "mas", "deg", "mas", " ",  "mas", "mas", "mas/yr", "mas/yr",  "mas/yr", "mas/yr",  "mag", "mag", "mag", "mag", "mag", "mag"))
            if line[0].isdigit():
                words = line.split()
                alpha, delta = words[0], words[2]
                x, y = countxy(float(alpha0), float(delta0), float(alpha), float(delta))
                ealpha, edelta = words[1], words[3]
                source = words[4]
                prlx, eprlx = words[5], words[6]
                mua, emua, mud, emud = words[7:11]
                gmag, egmag = words[14], words[15]
                try:
                    bpmag, ebpmag = words[18], words[19]
                    rpmag, erpmag = words[22], words[23]
                except IndexError:
                    bpmag, rpmag = ["99.9999" for _ in range(2)]
                    ebpmag, erpmag = ["9.9999" for _ in range(2)]
                g.write(" {:>+9.4f} {:>+9.4f} {:>15} {:>9} {:>15} {:>9} {:>19} {:>6} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>6} {:>7} {:>7} {:>7} {:>7}\n".format(x, y, alpha, ealpha, delta, edelta, source, prlx, eprlx, mua, emua, mud, emud, gmag, egmag, bpmag, ebpmag, rpmag, erpmag))

