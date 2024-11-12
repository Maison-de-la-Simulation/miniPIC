
import pandas as pd
import numpy as np
from numpy import linalg as LA

import matplotlib.pyplot as plt
import math


import sys
import os
import glob
import struct


# ______________________________________________________________________________
# Read command line arguments

file_path = ""

if len(sys.argv) > 1:
    file_path = sys.argv[1]
else:
    sys.exit("Please, specify a valid path to a binary file as a first command line argument.")


# ______________________________________________________________________________
# Read all the binary file
binary_files = sorted(glob.glob(file_path + '/cloud_s00_*'))

number_of_iterations=len(binary_files)

tab_id = np.zeros(number_of_iterations)
tab_w  = np.zeros(number_of_iterations)
tab_x  = np.zeros(number_of_iterations)
tab_y  = np.zeros(number_of_iterations)
tab_z  = np.zeros(number_of_iterations)
tab_px = np.zeros(number_of_iterations)
tab_py = np.zeros(number_of_iterations)
tab_pz = np.zeros(number_of_iterations)


it=0
for binary_file in binary_files:

    print(binary_file)

    file = open(binary_file, 'rb')

    content = file.read()

    k = 0

    particle_number = struct.unpack('I',content[k:k+4])[0] ; k+= 4

#    id = np.zeros(particle_number)
#    w  = np.zeros(particle_number)
#    x  = np.zeros(particle_number)
#    y  = np.zeros(particle_number)
#    z  = np.zeros(particle_number)
#    px = np.zeros(particle_number)
#    py = np.zeros(particle_number)
#    pz = np.zeros(particle_number)

    #for ip in range(particle_number):
    
    #tab_id[ip] = ip;
    tab_w[it]  = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_x[it]  = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_y[it]  = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_z[it]  = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_px[it] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_py[it] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_pz[it] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8

    it=it+1

####################
#     THEORIC
####################

# ______________________________________________________________________________
# Params for the simulation

dx=1/32
dy=1/32
dz=1/32

q=-1.
m=1.

dt= 0.9 * np.sqrt(1 / (1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz)));
print("dt= {}".format(dt))
iterations=np.arange(0,number_of_iterations)

# ______________________________________________________________________________
# Compute theoric

X0=np.array([tab_x[0], tab_y[0], tab_z[0]])
normX0=LA.norm(X0)
print("X0: {}".format(X0))
print("Norme de X0: {}".format(normX0))

p0=np.array([tab_px[0], tab_py[0], tab_pz[0]])
normp0=LA.norm(p0)
print("p0: {}".format(p0))
print("Norme de p0: {}".format(normp0))

B0=np.array([0.0, 0.0, -9.0])
normB0=LA.norm(B0)
print("B0: {}".format(B0))
print("Norme de B0: {}".format(normB0))

Rth = m*normp0/(abs(q)*normB0)
print("R theoric: {}".format(Rth))

gamma=np.sqrt(1+np.sum((p0)**2))
print("Gamma, Lorentz: {}".format(gamma))
w=(abs(q)/m)*(normB0/gamma)
print("Omega, angular frequency: {}".format(w))
period=(2*math.pi)/w
print("Periode: {}".format(period))



normp=np.zeros(number_of_iterations)
R=np.zeros(number_of_iterations)
xth=np.zeros(number_of_iterations)
yth=np.zeros(number_of_iterations)

xth0 = tab_x[0] - Rth
yth0 = tab_y[0]

xlim = [xth0 - 1.5*Rth, xth0 + 1.5*Rth]
ylim = [yth0 - 1.5*Rth, yth0 + 1.5*Rth]

for i,it in enumerate(iterations):
    t = dt * it 

    normp[i] = LA.norm(np.array([tab_px[i],tab_py[i],tab_pz[i]]))

    xth[i] = xth0 + Rth * math.cos(w*t)
    yth[i] = yth0 + Rth * math.sin(w*t)

#error_abs = np.sqrt((tab_x - xth)**2 + (tab_y - yth)**2)

Rs =  np.sqrt((tab_x - xth0)**2 + (tab_y - yth0)**2)
print(Rs)

error_abs = np.abs(Rs - Rth)
error_rel = error_abs / Rth

print("tab_x[0]: {}".format(tab_x))
print("tab_y[0]: {}".format(tab_y))
print("xth: {}".format(xth))
print("yth: {}".format(yth))


###########
fig1,ax=plt.subplots(1,1)
ax.plot(iterations, normp)
ax.set_xlabel("t")
ax.set_ylabel("||p||")
ax.set_ylim([normp0-1,normp0+1])



fig2, ax = plt.subplots(1,1)
circleRth=plt.Circle((xth0, yth0), Rth, color='green', fill=False)
ax.add_patch(circleRth)
scat = ax.scatter(x=tab_x,y=tab_y,marker='x', label='Compute')
scat = ax.scatter(x=xth,y=yth,marker='+', label='Theoric')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.legend()
ax.grid()



fig3, axs = plt.subplots(1,2)
axs[0].plot(iterations, tab_x, marker='+', label='Compute')
axs[0].plot(iterations, xth, label='Theoric')
axs[0].set_xlabel("t")
axs[0].set_ylabel("x")
axs[0].set_xlim([0,42])
axs[0].set_ylim(xlim)
axs[0].legend()
axs[1].plot(iterations, tab_y, marker='+', label='Compute')
axs[1].plot(iterations, yth, label='Theoric')
axs[1].set_xlabel("t")
axs[1].set_ylabel("y")
axs[1].set_xlim([0,42])
axs[1].set_ylim(ylim)
axs[1].legend()



fig4,ax=plt.subplots(1,1)
ax.plot(iterations, error_abs)
ax.set_xlabel("t")
ax.set_ylabel("error_abs")
fig5,ax=plt.subplots(1,1)
ax.plot(iterations, error_rel)
ax.set_xlabel("t")
ax.set_ylabel("error_rel")



plt.show()
