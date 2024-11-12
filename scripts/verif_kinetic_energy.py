
import pandas as pd
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

import struct
import sys
import os

species=["Electron", "Proton", "Positon", "Neutron", "Photon", "ion3", "ion4", "ion5", "ion6", "ion7"]

if len( sys.argv ) != 3:
    print( "Usage: python3 verif_kinetic_energy.py path/to/species/data/ path/to/field/data/" )
    exit()

filename = sys.argv[1]
field_filename = sys.argv[2]

#####


file = open(filename, 'rb')

content = file.read()

k = 0

n_it = struct.unpack('i',content[k:k+4])[0] ; k+= 4
n_sp = struct.unpack('i',content[k:k+4])[0] ; k+= 4

kin_energy=np.empty((n_sp,n_it))

for s in range(n_sp):
    for i in range(n_it):
        kin_energy[s][i] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8

x = np.arange(n_it)

for s in range(n_sp):
    plt.plot(x, kin_energy[s], label=species[s], marker='+')

#####

file = open(field_filename, 'rb')

content = file.read()

k = 0

n_it = struct.unpack('i',content[k:k+4])[0] ; k+= 4
x = np.arange(n_it)

Ex=np.empty(n_it)
Ey=np.empty(n_it)
Ez=np.empty(n_it)
Bx=np.empty(n_it)
By=np.empty(n_it)
Bz=np.empty(n_it)
for i in range(n_it):
    Ex[i] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
for i in range(n_it):
    Ey[i] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
for i in range(n_it):
    Ez[i] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
for i in range(n_it):
    Bx[i] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
for i in range(n_it):
    By[i] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
for i in range(n_it):
    Bz[i] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
field_sum=Ex+Ey+Ez+Bx+By+Bz

print(Ex[n_it-1])
print(Ey[n_it-1])
print(Ez[n_it-1])
print(Bx[n_it-1])
print(By[n_it-1])
print(Bz[n_it-1])



plt.plot(x, Ex, label='Ex', marker='+')
plt.plot(x, Ey, label='Ey', marker='+')
plt.plot(x, Ez, label='Ez', marker='+')
plt.plot(x, Bx, label='Bx', marker='+')
plt.plot(x, By, label='By', marker='+')
plt.plot(x, Bz, label='Bz', marker='+')
plt.plot(x, field_sum, label='sum', marker='+')

######

kin_energy_total=np.zeros((n_it))
for s in range(n_sp):
    kin_energy_total[:] += kin_energy[s][:]

plt.plot(x,kin_energy_total, label="Total", marker='+')

######

plt.legend()
plt.show()
