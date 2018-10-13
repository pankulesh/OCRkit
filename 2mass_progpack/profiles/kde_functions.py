#Functions for proper working of profile makers.
import numpy as np
from scipy import interpolate
from multiprocessing import Pool, cpu_count
from numpy import pi, f2py

source = """
        subroutine dens_estimator(i, step, rstar, delta, dens)
        integer i
        real*8 step, rstar, rstar2, delta, delta2, delta4, delta6
        real*8 pi, pi2, ri, ri2, bracket, bracket2
        real*8 arg, fimax, dens1, dens2, dens
cf2py   intent(out) dens
        pi = 3.14159265358979
        pi2 = pi*pi
        rstar2 = rstar*rstar
        delta2 = delta*delta
        delta4 = delta2*delta2
        delta6 = delta4*delta2
        ri = i*step
        ri2 = ri*ri
        bracket = 1.0 -(ri2+rstar2)/delta2
        bracket2 = bracket*bracket
        if (ri.lt.(delta-rstar)) then
            dens = 3.0*bracket2/pi/delta2+6.0*ri2*rstar2/pi/delta6
        else
            if (rstar.eq.0.) rstar = 0.01
            arg = (ri*ri+rstar*rstar-delta2)/(2.*ri*rstar)
            fimax=acos(arg)
            dens1 = 3.0*fimax*bracket2/pi2/delta2
            dens2 = 6.0*ri2*rstar2*fimax/pi2/delta6
            dens3 = 12.0*ri*rstar*bracket*sin(fimax)/pi2/delta4
            dens4 = 3.0*ri2*rstar2*sin(2.0*fimax)/pi2/delta6
            dens = dens1 + dens2 + dens3 + dens4
        endif
        return
        end
"""

try:
    import densest
except:
    f2py.compile(source.encode(), modulename='densest')
    import densest

dens_estimator = densest.dens_estimator

class Region:
    
    def __init__(self, ibound2):
        self.density = np.zeros(ibound2)
        self.d_mean = np.zeros(ibound2)
        self.d_disp = np.zeros(ibound2)
        self.n_tot = 0
    
    def bootstrap(self, argument, ibound, ibound2, x0, xlim, delta, step):
        de = density_estimator_1D
        nboot = 20
        density_max = np.amax(self.density)
        app_densitytck = interpolate.splrep(argument, self.density)
        dens_mags = [np.zeros(ibound) for _ in range(nboot)]
        for k in range(nboot):
        # Density estimation for "bootstrap" set of radial distances values.
            for x in (neiman_method(x0, xlim+delta, density_max, app_densitytck) for _ in range(self.n_tot)):
                rstar = x - x0
                imin, imax = make_imin_imax(rstar, delta, step, ibound)
                if imin <= ibound and imax >= 0:
                    dens_mags[k][imin:imax] += [de(i, step, rstar, delta) for i in range(imin, imax)]
        s1, s2 = (np.zeros(ibound2) for _ in range(2))

        for i in range(ibound):
            for k in range(nboot):
                s1[i] += dens_mags[k][i]
                s2[i] += dens_mags[k][i]**2

        self.d_mean = s1/nboot
        self.d_disp = np.sqrt((s2-nboot*self.d_mean**2)/(nboot-1))

def neiman_method(a, b, omega_max, tck):
    while True:
        x1 = a+(b-a)*np.random.random()
        x2 = omega_max*np.random.random()
        if x2 < interpolate.splev(x1, tck):
            return x1

def make_imin_imax(rstar, delta, step, ibound):
    if (rstar < delta):
        imin = int((rstar-delta)/step)+1
    else:
        imin = int((rstar-delta)/step)+2
    if imin < 0:
        imin = 0
    imax = int((rstar+delta)/step)+1
    if imax > ibound:
        imax = ibound
    return (imin, imax)

def density_estimator_1D(i, step, rstar, delta):
    ri = i*step
    bracket = 1 - (rstar-ri)**2/(delta**2)
    return (3*bracket/4/delta)

def edens_estimator(i, step, rstar, delta): #too slow, using fortran for this part with F2Py
    pi2 = pi**2
    rstar2 = rstar**2
    delta2 = delta**2
    delta4 = delta**4
    delta6 = delta**6
    ri = i*step
    ri2 = ri**2
    bracket = 1.0-(ri2+rstar2)/delta2
    bracket2 = bracket**2
    if ri < delta-rstar:
        return 3.0*bracket2/pi/delta2+6.0*ri2*rstar2/pi/delta6
    else:
        if rstar == 0:
            rstar = 0.01
        fimax = np.arccos((ri2+rstar2-delta2)/(2*ri*rstar))
        return 3*fimax*bracket2/pi2/delta2 + 6*ri2*rstar2*fimax/pi2/delta6 + 12*ri*rstar*bracket*np.sin(fimax)/pi2/delta4 + 3*ri2*rstar2*np.sin(2*fimax)/pi2/delta6

