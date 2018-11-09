#!/usr/bin/python
import numpy as np
import argparse
from scipy.integrate import simps

parser = argparse.ArgumentParser(description='Программа интегрирует функцию светимости и функцию масс скопления, находя число звёзд и массу скопления в солнечных массах.')
parser.add_argument("in_file", type=str, help='Путь к входному файлу с функцией светимости и масс')
parser.add_argument('name', type=str, help='Имя скопления')
args = parser.parse_args()

name = args.name
in_file = args.in_file
out_file = f"lf_mf_{name}_numstars.txt"

lumfun = np.loadtxt(in_file, unpack=True, usecols=(0,2,3,4))
massfun = np.loadtxt(in_file, unpack=True, usecols=(5,6,7,8))

massfun = np.fliplr(massfun)
mag, lf, lf_lo, lf_up = lumfun
massarg, mf, mf_lo, mf_up = massfun
mfm = massarg*mf
mfm_lo = massarg*mf_lo
mfm_up = massarg*mf_up
numbers = ["Число звёзд из ФМ: ",
           "Ошибка числа звёзд:",
           "Нижняя граница:    ",             
           "Верхняя граница:   ",
           "Масса:             ",
           "Ошибка массы:      ",
           "Нижняя граница:    ",
           "Верхняя граница:   "
           ]
N = simps(mf, massarg)
Nlo = simps(mf_lo, massarg)
Nup = simps(mf_up, massarg)
Nsigma = max(Nup - N, N - Nlo)
M = simps(mfm, massarg)
Mlo = simps(mfm_lo, massarg)
Mup = simps(mfm_up, massarg)
Msigma = max(Mup - M, M - Mlo)
ints = [N, Nsigma, Nlo, Nup, M, Msigma, Mlo, Mup]
ints = [int(round(i, -1)) for i in ints]
with open(out_file, 'w') as f:
    for x, y in zip(numbers, ints):
        print(f'{x} {y:>6}')
        f.write(f'{x} {y:>6}\n')

