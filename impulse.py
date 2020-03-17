# Impulsive Gaussian Function
#
# Copyrignts -- GPL3 declared by habe aka Hiroshi ABE (5 Mar. 2020)
#

import math
from scipy import fftpack as fft
from scipy import signal
import numpy as np
from matplotlib import pyplot
import sys

#### CONSTANTS ####

pi=3.1415926

#centre frequency of source excitation
freq=100.0e+9
omg=2.0*pi*freq
#time offset
t0=1.0e-11
#half width of Gaussian
Dt=1e-12

N=4096
t=np.zeros(N,np.float)
E=np.zeros(N,np.float)
f=np.zeros(N,np.float)
FE=np.zeros(N,np.float)
dt=10.0/freq/(N-1)
for i in range(0,N):
    t[i] = dt*(i-1)
    E[i] = math.cos(omg*(t[i]-t0))*math.exp(-(t[i]-t0)**2/(2*Dt**2))

pyplot.cla()
pyplot.grid(which="both")
#pyplot.ylim([-1.1,1.1])
pyplot.plot(t,E)
pyplot.pause(0)

df=1.0/t[N-1]
for i in range(0,N):
    f[i]=df*(N/2-i)
    
FE = fft.fftshift(fft.fft(E))

pyplot.cla()
pyplot.xlim([0,10*freq])
pyplot.yscale('log')
pyplot.plot(f,abs(FE))
pyplot.pause(0)
