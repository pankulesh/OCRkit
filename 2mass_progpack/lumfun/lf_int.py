#!/usr/bin/python
import numpy as np
import argparse
from scipy.integrate import simps

parser = argparse.ArgumentParser(description='Программа интегрирует все функции блеска по имени скопления. Функции блеска должны находиться в файлах типа lf_{имяскопления}_{r, 2c, 4c}.txt.')
parser.add_argument('name', type=str, help='Имя скопления')
parser.add_argument("-o", "--outfile", default="lf_nums_mass.txt", help='Имя выходного файла, по-умолчанию lf_{имяскопления}_numstars.txt')
args = parser.parse_args()

name = args.name
out_file = args.outfile

ts = ('r', '2c', '4c')

def integrate(lumfun):
    mag, lf, lf_lo, lf_up = lumfun
    I = [simps(y, x) for x, y in ([mag, lf], [mag, lf_lo], [mag, lf_up])]
    sigma = max(I[2]-I[0], I[0]-I[1])
    return I, sigma

numbers = ["Тип области сравнения", 
           "Число звёзд", 
           "Погрешность",
           "Нижняя граница", 
           "Верхняя граница",
           ]

with open(out_file, 'w') as f:
    f.write("{:>21} {:>15} {:>15} {:>15} {:>15}\n".format(*numbers))
    for t in ts:
        lumfun = np.loadtxt(f"lf_{name}_{t}.txt", unpack=True)
        I, sigma = integrate(lumfun)
        f.write("{:>21} {:>15.2f} {:>15.2f} {:>15.2f} {:>15.2f}\n".format(t, I[0], sigma, I[1], I[2]))

with open(out_file, 'r') as f:
    print(f.read())
