# The program shows dispersion curves.
#
# Copyrignts -- GPL3 declared by habe aka Hiroshi ABE (5 Mar. 2020)

import numpy as np
from matplotlib import pyplot
import math
import cmath

N=100
omega=np.zeros(N,np.float)
refract=np.zeros(N,np.complex)
omega_p=math.sqrt(20.0)/4.0
delta=0.07

for n in range(0,N):
    omega[n]=2.0/(N-1)*n
    refract[n]=cmath.sqrt(1-omega_p*omega_p/(omega[n]*omega[n]-1+2j*delta*omega[n]))


pyplot.cla()
pyplot.grid(which="both")
pyplot.ylim([0,3])
pyplot.plot(omega,refract.real,omega,refract.imag)
pyplot.pause(0)
