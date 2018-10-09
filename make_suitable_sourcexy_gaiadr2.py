#!/usr/bin/python
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Программа для подсчёта x, y из исходных файлов cdsclient. Пустые поля в BP и RP заполняются 99.9999.')
parser.add_argument("infile", help='Имя или путь ко входному файлу')
parser.add_argument("outfile",  help="Имя выходного файла")
parser.add_argument("-p", "--parallax", action='store_true', help="Укажите этот флаг, если нужно убрать звёзды, значения параллакса и собственных движений которых не определены")
parser.add_argument("-b", "--bprpmag", action='store_true', help="Укажите этот флаг, если нужно убрать звёзды, значения G_BP и G_RP не определены")
parser.add_argument("-t", "--title", action='store_true', help="Укажите этот флаг, если у заголовков нужно убрать # в начале")
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

def split_at(line, pattern):
    a = []
    a.append(line[:pattern[0]])
    for i in range(len(pattern)-1):
        a.append(line[pattern[i]:pattern[i+1]])
    return a
title = "#{:>9} {:>9} {:>15} {:>10} {:>19} {:>11} {:>23} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8} {:>7} {:>6} {:>7} {:>7} {:>7} {:>7} {:>7} {:>10}\n"
titlenohash = " {:>9} {:>9} {:>15} {:>10} {:>19} {:>11} {:>23} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8} {:>7} {:>6} {:>7} {:>7} {:>7} {:>7} {:>7} {:>10}\n"
star = " {:>+9.4f} {:>+9.4f} {:>15} {:>10} {:>19} {:>11} {:>23} {:>7.4f} {:>7.4f} {:>+8.3f} {:>+8.3f} {:>+8.3f} {:>+8.3f} {:>7.4f} {:>6.4f} {:>7.4f} {:>7.4f} {:>7.4f} {:>7.4f} {:>7.4f} {:>10.4f}\n"
names = ("x", "y", "RA_ICRS", "e_RA_ICRS", "DE_ICRS", "e_DE_ICRS", "Source",  
        "Plx", "e_Plx", "pmRA", "e_pmRA",  "pmDE", "e_pmDE",  "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag", "BP-RP", "e_BP-RPmag")
measures = ("arcmin", "arcmin", "deg", "mas", "deg", "mas", " ",  
        "mas", "mas", "mas/yr", "mas/yr",  "mas/yr", "mas/yr",  "mag", "mag", "mag", "mag", "mag", "mag", "mag", "mag")
with open(outfile, 'w') as g:
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith("#-c="):
                alpha0, delta0 = line[4:].split()
                g.write(line) 
                g.write(title.format(*list(range(21))))
                if args.title:
                    g.write(titlenohash.format(*names))
                    g.write(titlenohash.format(*measures))
                else:
                    g.write(title.format(*names))
                    g.write(title.format(*measures))
            if line.startswith('-'):
                pat = [pos for pos, char in enumerate(line) if char == '\t']	
#                print(pat)
            if line[0].isdigit():
#                words = line.split()
                words = split_at(line, pat)
                alpha, delta = words[0], words[2]
                x, y = countxy(float(alpha0), float(delta0), float(alpha), float(delta))
                ealpha, edelta = words[1], words[3]
                source = words[4]
                try:
                    prlx, eprlx = [float(i) for i in (words[5], words[6])]
                    mua, emua, mud, emud = [float(i) for i in words[7:11]]
                except ValueError:
                    if args.parallax:
                        continue
                    else:
                        prlx, eprlx = [99.9999]*2
                        mua, emua, mud, emud = [99.999]*4 
                gmag, egmag = [float(i) for i in words[14:16]]
                try:
                    bpmag, ebpmag = float(words[18]), float(words[19])
                    rpmag, erpmag = float(words[22]), float(words[23])
                    bpmrp = bpmag - rpmag
                    ebrmrp = np.sqrt(ebpmag*ebpmag+erpmag*erpmag)
                except ValueError:
                    if args.bprpmag:            
                        continue
                    else:
                        bpmag, rpmag, bpmrp = [99.9999 for _ in range(3)]
                        ebpmag, erpmag, ebrmrp = [9.9999 for _ in range(3)]

                values = (x, y, alpha, ealpha, delta, edelta, source, prlx, eprlx, mua, emua, mud, emud, gmag, egmag, bpmag, ebpmag, rpmag, erpmag, bpmrp, ebrmrp)
                g.write(star.format(*values))

