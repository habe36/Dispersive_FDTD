# Lorentz model by ADE method (by habe Mar. 2020)
# usage:
# python lorentz1d.py 20
# where 20 is the Resonant Frequency of Lorentz Media in GHz
#
# Copyrignts -- GPL3 declared by habe aka Hiroshi ABE (5 Mar. 2020)
#

import numpy as np
from matplotlib import pyplot
import math
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
#time offset
t0=1.0e-11
#half width of Gaussian
Dt=1e-12

# Material Parameters
eps1=1.0
sig=0.0
Np=1
omega=[2*pi*1e9*float(sys.argv[1])]
delta=np.multiply(omega,0.1)
# eps for infinity
eps8=1.5
# eps for DC
epss=[3.0]

#### NUMERICAL PARAMETERS ####
# Grids
N = 3000
L = 0.05

CFL=0.9

dx=L/(N-1)
dt=CFL*dx/c
T=L/c/CFL*2

# File to write data
F0 = open("lorentz1d.data",mode='w')

x=np.zeros(N,np.float)
for i in range(0,N):
    x[i] = dx*i

xi = np.zeros(Np,np.float)
alpha = np.zeros(Np,np.float)
gamma  = np.zeros(Np,np.float)
for p in range(0,Np):
    ot = omega[p]*dt
    alpha[p]=(2-ot**2)/(1+delta[p]*dt)
    xi[p]=(delta[p]*dt-1)/(1+delta[p]*dt)
    gamma[p]=eps0*(epss[p]-eps8)*ot**2/(1+delta[p]*dt)
print(alpha)
print(xi)
print(gamma)

#pyplot.plot(x,eps)
#pyplot.show()

Ey=np.array([0.0]*N)
Eym=np.array([0.0]*N)
Eyp=np.array([0.0]*N)
Bz=np.zeros(N,dtype=np.float)
J=np.array([np.zeros(N,np.float)]*Np)
Jm=np.array([np.zeros(N,np.float)]*Np)
Jp=np.array([np.zeros(N,np.float)]*Np)
G=0
for p in range(0,Np):
    G=G+gamma[p]
G=0.5*G
D=2*eps0*eps8 + G + sig*dt
C1=G/D
C2=(2*eps0*eps8-sig*dt)/D
C3=2*dt/D
print(C1,C2,C3)

c1 = c/math.sqrt(eps1)
c2 = c/math.sqrt(eps8)

murc1 = (c1*dt - dx)/(c1*dt + dx)
murc2 = (c2*dt - dx)/(c2*dt + dx)
#print(murc1,murc2)
# Grid number for the media boundary
N2 = int(N/2.0)
# Grid numbers for sampling (Input, Reflection, Transparent)
n0 = int(0.1*N)
n1 = int(0.4*N)
n2 = int(0.6*N)

# Time Marching
for t in np.arange(0,T,dt):
    
    Bz[0:N-1] = Bz[0:N-1] - dt/dx * (Ey[1:N]-Ey[0:N-1])

    # First Mur Proc.
    FE0 = Ey[1] - murc1*Ey[0]
    FE1 = Ey[N-2] - murc2*Ey[N-1]
    
    JJ = np.zeros(N,np.float)
    for p in range(0,Np):
        JJ[:] = JJ[:] + (1+alpha[p])*J[p,:] + xi[p]*Jm[p,:]
    # for vacuum media
    Eyp[1:N2] = Ey[1:N2] - dt/(dx*eps0*mu0) * (Bz[1:N2]-Bz[0:N2-1])
    # for Lorentz media
    Eyp[N2:N-1] = C1*Eym[N2:N-1] + C2*Ey[N2:N-1] + C3*( -1.0/(mu0*dx) *(Bz[N2:N-1]-Bz[N2-1:N-2]) -0.5*JJ[N2:N-1] )
    
    # Stimulation/Mur (switched at 4e-11 seconds)
    if(t<4e-11):
        Eyp[0] = math.cos(omg*(t-t0))*math.exp(-(t-t0)**2/(2*Dt**2))
    else:
        Eyp[0] = FE0 + murc1*Eyp[1]

    # Second Mur Proc.
    Eyp[N-1] = FE1 + murc2*Eyp[N-2]

    for p in range(0,Np):
        Jp[p,1:N-1] = alpha[p]*J[p,1:N-1] + xi[p]*Jm[p,1:N-1] + gamma[p]*(Eyp[1:N-1]-Eym[1:N-1])/(2*dt)
    #print(Jp)
    
    Jm[:,1:N-1] = J[:,1:N-1]
    J[:,1:N-1] = Jp[:,1:N-1]

    Eym[:] = Ey[:]
    Ey[:] = Eyp[:]
    #print("Time = ", t)
    s = "%e %e %e %e\n" % (t, Eyp[n0], Eyp[n1], Eyp[n2])
    F0.write(s)
    
    pyplot.cla()
    pyplot.grid(which="both")
    pyplot.ylim([-1.1,1.1])
    #pyplot.plot(x,Ey,x,c*Bz)
    pyplot.plot(x,Ey)
    #pyplot.plot(x,J[0,:])
    pyplot.pause(0.01)

F0.close()
