#!/usr/bin/python
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Программа для подсчёта x, y из исходных файлов cdsclient. Пустые поля для ошибок звёздных величин заполняются 9.999.')
parser.add_argument("infile", help='Имя или путь ко входному файлу')
parser.add_argument("outfile",  help="Имя выходного файла")
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

title =       "#{:>9} {:>9} {:>10} {:>10} {:>17} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}\n"
titlenohash = " {:>9} {:>9} {:>10} {:>10} {:>17} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}\n"
star = " {:>+9.4f} {:>+9.4f} {:>10} {:>10} {:>17} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f}\n"
names = ("x", "y", "RAJ2000", "DEJ2000", "2MASS", "J", "e_J", "H", "e_H",  "K", "e_K", "J-H", "e_J-H", "H-K", "e_H-K")    
measures = ("arcmin", "arcmin", "deg", "deg", " ", "mag", "mag", "mag", "mag", "mag", "mag", "mag", "mag", "mag", "mag")
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
            if line[0].isdigit():
                words = split_at(line, pat)
                alpha, delta = words[0].strip(), words[1].strip()
                x, y = countxy(float(alpha0), float(delta0), float(alpha), float(delta))
                source = words[2].strip()
                j, h, k = [float(i) for i in (words[3], words[5], words[7])]
                try:
                    ej = float(words[4])
                except ValueError:
                    ej = 9.999
                try:
                    eh = float(words[6])
                except ValueError:
                    eh = 9.999
                try:
                    ek = float(words[8])
                except ValueError:
                    ek = 9.999
                jmh = j - h
                if ej==9.999 or eh==9.999:
                    ejmh = 9.999
                else:
                    ejmh = np.sqrt(ej*ej+eh*eh)
                hmk = h - k
                if ek==9.999 or eh==9.999:
                    ehmk = 9.999
                else:
                    ehmk = np.sqrt(ek*ek+eh*eh)
                values = (x, y, alpha, delta, source,j, ej, h, eh, k, ek, jmh, ejmh, hmk, ehmk)
                g.write(star.format(*values))

