# Lorentz model by TGM method (by habe 2020)
# usage:
# python tgm1d.py 20
# where 20 is the Resonant Frequency of Lorentz Media in GHz
#
# Copyrignts -- GPL3 declared by habe aka Hiroshi ABE (5 Mar. 2020)
#

import numpy as np
from matplotlib import pyplot
import math
import cmath
import sys

#### CONSTANTS ####

pi=3.1415926

#speed of light in free space
c=2.99792458e8
#permeability of free space
mu0=4.0*pi*1.0e-7
#permittivity of free space
eps0=1.0/(c*c*mu0)

#centre frequency of source excitation
freq=100.0e+9
omg=2.0*pi*freq

# Material Parameters
eps1=1.0
sig=0.0
Np=1
ff=np.array([float(sys.argv[1]),float(sys.argv[2]), float(sys.argv[3])])
omega=np.multiply(ff,2*pi*1e9)
delta=np.multiply(omega,0.1)
eps8=1.5
Np=ff.size
epss=np.array(Np*[3.0])

t0=1.0e-11
Dt=1e-12

#### NUMERICAL PARAMETERS ####
# Grids
N= 3000
L= 0.05

CFL=0.9

dx=L/(N-1)
dt=CFL*dx/c
T=L/c/CFL*2

# File to write data
F0 = open("tgm1d.data",mode='w')

x=np.zeros(N,np.float)
for i in range(0,N):
    x[i] = dx*i

# Variables    
Ey=np.array([0.0]*N)
Bz=np.zeros(N,dtype=np.float)
Fp=np.array([np.zeros(N,np.complex)]*Np)
Fm=np.array([np.zeros(N,np.complex)]*Np)
J=np.array([np.zeros(N,np.complex)]*Np)

c1 = c/math.sqrt(eps1)
c2 = c/math.sqrt(eps8)

murc1 = (c1*dt - dx)/(c1*dt + dx)
murc2 = (c2*dt - dx)/(c2*dt + dx)
#print(murc1,murc2)
# Grid number for the media boundary
N2=int(N/2)
# Grid numbers for sampling (Input, Reflection, Transparent)
n0 = int(0.1*N)
n1 = int(0.4*N)
n2 = int(0.6*N)

zp = np.zeros(Np,np.complex)
zm = np.zeros(Np,np.complex)
FpC = np.zeros(Np,np.complex)
FmC = np.zeros(Np,np.complex)
for p in range(0,Np):
    #print(p)
    omd = omega[p]*omega[p]-delta[p]*delta[p]
    zp[p] =  math.sqrt(omd) + delta[p]*1j
    zm[p] = -math.sqrt(omd) + delta[p]*1j
    FpC[p] = (cmath.exp(1j*zp[p]*dt/2) - cmath.exp(-1j*zp[p]*dt/2))/(zp[p]*(zm[p]-zp[p]))
    FmC[p] = (cmath.exp(1j*zm[p]*dt/2) - cmath.exp(-1j*zm[p]*dt/2))/(zm[p]*(zp[p]-zm[p]))
#print(FpC[0])
#print(FmC[0])

it=0
# Time Marching
for t in np.arange(0,T,dt):
    
    Bz[0:N-1] = Bz[0:N-1] - dt/dx * (Ey[1:N]-Ey[0:N-1])

    # First Mur Proc.
    FE0 = Ey[1] - murc1*Ey[0]
    FE1 = Ey[N-2] - murc2*Ey[N-1]
    
    for p in range(0,Np):
        zpdt = 1j*zp[p]*dt
        zmdt = 1j*zm[p]*dt
        zpdt2 = 1j*zp[p]*dt/2
        zmdt2 = 1j*zm[p]*dt/2
        Fp[p,:] = Fp[p,:]*cmath.exp(zpdt) + FpC[p]*Ey[:]
        Fm[p,:] = Fm[p,:]*cmath.exp(zmdt) + FmC[p]*Ey[:]
        J[p,:] = -eps0*(epss[p]-eps8)*omega[p]*omega[p]*(
            1j*zp[p]*cmath.exp(zpdt2)*Fp[p,:] +
            1j*zm[p]*cmath.exp(zmdt2)*Fm[p,:])


    # for vacuum media
    Ey[1:N2] = Ey[1:N2] - dt/(dx*eps0*mu0) * (Bz[1:N2]-Bz[0:N2-1])
    # for Lorentz Media
    Ey[N2:N-1] = Ey[N2:N-1] - dt/(dx*eps0*eps8*mu0) * (Bz[N2:N-1]-Bz[N2-1:N-2]) + J[0,N2:N-1].real/(eps0*eps8)*dt

    # Second Mur Proc.
    Ey[N-1] = FE1 + murc2*Ey[N-2]
    
    # Stimulation/Mur (switched at 4e-11 seconds)
    if(t<4e-11):
        Ey[0] = math.cos(omg*(t-t0))*math.exp(-(t-t0)**2/(2*Dt**2))
    else:
        Ey[0] = FE0 + murc1*Ey[1]
    
    s = "%e %e %e %e\n" % (t, Ey[n0], Ey[n1], Ey[n2])
    F0.write(s)

    pyplot.cla()
    pyplot.grid(which="both")
    pyplot.ylim([-1.1,1.1])
    #pyplot.plot(x,Ey,x,c*Bz)
    pyplot.plot(x,Ey)
    #pyplot.plot(x,J[0,:].real,x,J[0,:].imag)
    #pyplot.plot(x,J[0,:])
    #pyplot.plot(x,J[0,:].imag)
    pyplot.pause(0.01)
    #ims.append(im)
    #if (it%100) == 0:
    #fn = './images/' + ('%04d' % (it/100)) + '.png'
    #    fig.savefig(fn)

    print("Time = ", t)

#ani = anim.ArtistAnimation(fig, ims)
#ani.save('tgm.gif', writer='imagemagick')

F0.close()
