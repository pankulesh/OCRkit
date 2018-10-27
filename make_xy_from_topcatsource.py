#!/usr/bin/python
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Программа для подсчёта x, y из импортированных из TOPCAT файлов.')
parser.add_argument("infile", help='Имя или путь ко входному файлу')
parser.add_argument("outfile",  help="Имя выходного файла")
parser.add_argument("name", help="Имя скопления для определения центра по готовой встроенной таблице, простите за такой костыль")
#parser.add_argument("-p", "--parallax", action='store_true', help="Укажите этот флаг, если нужно убрать звёзды, значения параллакса и собственных движений которых не определены")
#parser.add_argument("-b", "--bprpmag", action='store_true', help="Укажите этот флаг, если нужно убрать звёзды, значения G_BP и G_RP не определены")
#parser.add_argument("-t", "--title", action='store_true', help="Укажите этот флаг, если у заголовков нужно убрать # в начале")
args = parser.parse_args()
infile = args.infile
outfile = args.outfile

centers = { "NGC6834":(298.050, +29.410),
            "NGC4052":(180.315, -63.215),
            "NGC5316":(208.485, -61.870),
            "NGC5715":(220.875, -57.567),
            "NGC6268":(255.532, -39.715),
            "NGC2099":(088.087, +32.570),
            "Cz38"   :(282.447, +04.960),
            "IC2714" :(169.365, -62.740),
            "NGC1912":(082.215, +35.800),
            "NGC7142":(326.302, +65.792),
}

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
alpha0, delta0 = centers[args.name]
with open(outfile, 'w') as g:
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == '#':
                g.write("#         x           y   "+line[2:])
                continue
            alpha, delta, other = line.split(None, 2)
            x, y = countxy(float(alpha0), float(delta0), float(alpha), float(delta))
            g.write(f"  {x:9.4f}   {y:9.4f}   {line[2:]}")
