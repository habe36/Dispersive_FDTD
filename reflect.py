# Calculate Reflection Coefficient in Fourier space from simulated data
# Usage:
# python reflect.py tgm1d.data
# where tgm1d.data can be lorentz1d.data
# This program generate tgm1d.data.fft data which contains FFTed data.
# To display use 'reflection.gp'.
# if you want more resolution, change NN values to bigger value, saying 16*N.
# this will padd zero data to input data and reflected data.
#
# Copyrignts -- GPL3 declared by habe aka Hiroshi ABE (5 Mar. 2020)
#
from scipy import fftpack as fft
from scipy import signal
from matplotlib import pyplot as plt
import numpy as np
import csv
import sys

fin=sys.argv[1]
print(fin)

d0=[]
d1=[]
d2=[]
d3=[]
with open(fin,mode='r') as f:
    csv_reader = csv.reader(f, delimiter=' ')
    for row in csv_reader:
        d0.append(row[0])
        d1.append(row[1])
        d2.append(row[2])
        d3.append(row[3])

N=len(d0)
print(N)
NN=N
x = np.zeros(NN,np.float)
e0 = np.zeros(NN,np.float)
e1 = np.zeros(NN,np.float)
e2 = np.zeros(NN,np.float)

han = np.hamming(NN)

dx=float(d0[2])-float(d0[1])

for n in range(0,NN):
    x[n] = dx*(n-1)
    
for n in range(0,NN):
    if(x[n] < 1e-10):
        e0[n] = d1[n]
        e1[n] = 0.0
    elif(x[n] < 3e-10):
        e0[n] = 0.0
        e1[n] = d2[n]
    else:
        e0[n] = 0.0
        e1[n] = 0.0

#e0 = han*e0
#e1 = han*e1
        
fe0 = fft.fftshift(fft.fft(e0))
fe1 = fft.fftshift(fft.fft(e1))
ff=np.zeros(NN,np.float)

fout=fin + '.fft'
FO = open(fout,mode='w')

df=1.0/x[NN-1]
f=np.zeros(NN,np.float)

for n in range(0,NN):
    if( fe0[n]!=0.0 ):
        ff[n] = abs(fe1[n]/fe0[n])
    else:
        ff[n] = 0.0

    f[n] = df*(NN/2-n)
    s="%e %e %e %e\n" % (f[n],abs(fe0[n]),abs(fe1[n]),ff[n])
    FO.write(s)

FO.close()

plt.cla()
plt.grid(which="both")
#plt.ylim([-1.1,1.1])
plt.plot(x,e0,x,e1)
plt.show()
plt.pause(0)
plt.plot(f,abs(fe0),f,abs(fe1))
plt.show()
plt.pause(0)
plt.ylim([0.0,1.1])
plt.xlim([0.0,2.0e11])
plt.plot(f,ff)
plt.show()
