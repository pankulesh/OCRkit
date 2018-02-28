#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline as spline
import argparse

# Nearest value index determination in array (for linear interpolation)
def find_nearest(array,value):
    return array[(np.abs(array-value)).argmin()]

class Magint:
    
    def __init__(self, jd, ju, n, dn):
        self.jd = jd
        self.ju = ju
        self.N = n
        self.sigmaN = dn
#Linear interpolation borders

    def evaluate_border_mass_lin(self, arr):
        idx = find_nearest(arr[0], self.jd)
        if (arr[0][idx]-self.jd)*(arr[0][idx+1]-self.jd)<0:
            j1, m1 = np.transpose(arr)[idx]
            j2, m2 = np.transpose(arr)[idx+1]
            self.md = m1 + (m2-m1)/(j2-j1)*(self.jd-j1)
        else:
            j1, m1 = np.transpose(arr)[idx-1]
            j2, m2 = np.transpose(arr)[idx]
            self.md = m1 + (m2-m1)/(j2-j1)*(self.jd-j1)
        self.sigmamd = max(abs(self.md-m1), abs(self.md-m2))
        
        idx = find_nearest(arr[0], self.ju)
        if (arr[0][idx]-self.ju)*(arr[0][idx+1]-self.ju)<0:
            j1, m1 = np.transpose(arr)[idx]
            j2, m2 = np.transpose(arr)[idx+1]
            b
            self.mu = m1 + (m2-m1)/(j2-j1)*(self.ju-j1)
        else:
            j1, m1 = np.transpose(arr)[idx-1]
            j2, m2 = np.transpose(arr)[idx]
            self.mu = m1 + (m2-m1)/(j2-j1)*(self.ju-j1)
        self.sigmamu = max(abs(self.mu-m1), abs(self.mu-m2))
    
    def evaluate_border_mass(self, massspline, arr):
        self.md = massspline(self.jd)
        if np.isnan(self.md):
            print("No appropriate mass for {}, change magnitude bounds".format(self.jd))
            quit()
        self.sigmamd = abs(find_nearest(arr, self.md)-self.md)
        self.mu = massspline(self.ju)
        if np.isnan(self.mu):
            print("No appropriate mass for {}, change magnitude bounds".format(self.ju))
            quit()
        self.sigmamu = abs(find_nearest(arr, self.mu)-self.mu)

    def evaluate_mass_kroupa(self):
            bracket1=self.mu**(-1.3)-self.md**(-1.3)
            bracket2=self.mu**(-0.3)-self.md**(-0.3)
            A = 1.3*self.N/bracket1
            self.M = A*bracket2/0.3
            dAdN=1.3/bracket1
            dAdm11=1.3*self.N/(bracket1*bracket1*self.mu**2.3)
            dAdmup=-1.3*self.N/(bracket1*bracket1*self.md**2.3)
            dMdA=bracket2/0.3
            dMdm11=-A/(0.3*self.mu**1.3)
            dMdmup=A/(0.3*self.md**1.3)
            sA = np.sqrt((dAdN*self.sigmaN)**2+(dAdm11*self.sigmamu)**2+(dAdmup*self.sigmamd)**2)
            self.sigmaM = np.sqrt((dMdA*sA)**2+(dMdm11*self.sigmamu)**2+(dMdmup*self.sigmamd)**2)
            print("M = {0} +- {1} in {2} - {3} magnitude interval".format(self.M, self.sigmaM, self.jd, self.ju))

    def evaluate_mass_mean(self):
            self.M = self.N/2*(self.mu+self.md)
            dMdN = (self.mu+self.md)/2
            self.sigmaM = np.sqrt((self.sigmaN*dMdN)**2+(self.N/2*self.sigmamu)**2+(self.N/2*self.sigmamd)**2)
            print("M = {0} +- {1} in {2} - {3} magnitude interval".format(self.M, self.sigmaM, self.jd, self.ju))

des = """Count cluster's mass from isochrone and counted number in magnitude intervals
         r, E(B-V), eE(B-V), distmod, edistmod must be in file 'photometry_NAME.txt' """
parser = argparse.ArgumentParser(description=des)
parser.add_argument('iso_in', type=str, help="Isochrone filename with two columns: mass abs_magnitude (mass(magnitude) must be one-valued)")
parser.add_argument('name', type=str, help="Cluster name (as in filenames)")
parser.add_argument('halfwidth', type=str, help="Kernel halfwidth (as in numstars filename)")
args = parser.parse_args()
name = args.name
iso_in = args.iso_in
h = args.halfwidth

param_in = "photometry_{}.txt".format(name)
with open(param_in, 'r') as f:
    f.readline()
    param = f.readline().split()
r, ebv, eebv, pr, epr, maglim, _ = (float(x) for x in param)
ebv2ag = 0.86*3.1
Ag = ebv*ebv2ag
#pr = 5*np.log10(r) - 5
fixmag = pr + Ag
sigmaJ = np.sqrt(epr**2+(ebv2ag*eebv)**2)
#param = ['1240', '0.34', '0.2']
num_in = "numstars_{0}_h={1}.txt".format(name, h)
intervals = []
num0 = np.genfromtxt(num_in, skip_header=1, usecols=(0,1,2))
intervals.append(Magint(maglim, num0[0][0], num0[0][1], num0[0][2]))

for i in range(1, len(num0)):
    intervals.append(Magint(num0[i-1][0], num0[i][0], num0[i][1]-num0[i-1][1], np.sqrt(num0[i][2]**2+num0[i-1][2]**2)))

mass = np.genfromtxt(iso_in, unpack=True)
mass_spline = spline(mass[1] + fixmag, mass[0], extrapolate=False)
plt.plot(mass[1]+fixmag, mass[0])
plt.show()
totalmass = 0
sumsq = 0

for i in intervals:
    i.evaluate_border_mass(mass_spline, mass[0])

if intervals[0].N > 20:
    intervals[0].evaluate_mass_kroupa()
    totalmass += intervals[0].M
    sumsq += intervals[0].sigmaM**2
else:
    intervals[0].evaluate_mass_mean()
    totalmass += intervals[0].M
    sumsq += intervals[0].sigmaM**2

for i in intervals[1:]:
    i.evaluate_mass_mean()
    totalmass += i.M
    sumsq += i.sigmaM**2

totalmasserror = np.sqrt(sumsq)

print("Total mass of {0} is {1} +- {2}".format(name, totalmass, totalmasserror))
