
import pandas as pd
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

import sys
import os
import glob
import struct


# ______________________________________________________________________________
# Read command line arguments

file_path = "/Users/mathieu/Codes/minipic/kokkos/tests/Ecst/diags"

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


x_round = 0
y_round = 0
z_round = 0

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
#
#    for ip in range(particle_number):

        #id[ip] = ip;
    tab_w[it]  = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_x[it]  = struct.unpack('d',content[k:k + 8])[0]  ; k+= 8
    tab_y[it]  = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_z[it]  = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_px[it] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_py[it] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
    tab_pz[it] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8

    tab_x[it] += x_round
    tab_y[it] += y_round
    tab_z[it] += z_round

    if (it > 0):
        # Periodic conditions
        if (tab_x[it-1] >= tab_x[it]):
            x_round += 1
            tab_x[it] += 1
        if (tab_y[it-1] >= tab_y[it]):
            y_round += 1
            tab_y[it] += 1
        if (tab_z[it-1] >= tab_z[it]):
            z_round += 1
            tab_z[it] += 1

    print(" {} {} {}".format(tab_x[it], tab_y[it], tab_z[it]))

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

dt= 0.9 * np.sqrt(1 / (1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz)))
iterations=np.arange(0,number_of_iterations)

print("dt: {}".format(dt))

# ______________________________________________________________________________
# Compute theoric

X0=np.array([tab_x[0], tab_y[0], tab_z[0]])
print("X0: {}".format(X0))
p0_correct=np.array([tab_px[0], tab_py[0], tab_pz[0]])
p0=np.array([0., 0., 0.])
normp0=LA.norm(p0)
print("p0: {}".format(p0))
print("normp0: {}".format(normp0))
E0=np.array([-0.05, -0.01, -0.08])
normE0=LA.norm(E0)
print("E0: {}".format(E0))
print("normE0: {}".format(normE0))

direction = E0 / normE0

Fl0=(q/m)*E0
print("Lorentz force: {}".format(Fl0))

u0=p0/m
print("u0: {}".format(u0))
C = X0 - (m/(q*normE0)) * np.sqrt(1+np.sum(u0**2)) * direction
print("C: {}".format(C))
print(direction)
print(u0)

print( (m/(q*normE0)) * np.sqrt(1+np.sum(u0**2)) * direction + C )


pth=np.zeros((3,number_of_iterations))
xth=np.zeros((3,number_of_iterations))

for i,it in enumerate(iterations):

    #print(i, it)

    t = dt * it * 10
    
    pth[:,i] = Fl0*t+p0_correct
    #print(pth[:,i])

    u=(p0+q*E0*t)/m
    xth[:,i] = (m/(q*normE0)) * np.sqrt(1+np.sum(u**2)) * direction + C
    

print("tab_x: {}".format(tab_x))
print("xth[0,:]: {}".format(xth[0,:]))

print("tab_px: {}".format(tab_px))
print("pth[0,:]: {}".format(pth[0,:]))


error_abs_x = abs(tab_x - xth[0])
error_abs_y = abs(tab_y - xth[1])
error_abs_z = abs(tab_z - xth[2])
error_rel_x = abs(tab_x - xth[0]) / xth[0]
error_rel_y = abs(tab_y - xth[1]) / xth[1]
error_rel_z = abs(tab_z - xth[2]) / xth[2]

error_abs_px = abs(tab_px - pth[0])
error_abs_py = abs(tab_py - pth[1])
error_abs_pz = abs(tab_pz - pth[2])
error_rel_px = abs(tab_px - pth[0]) / pth[0]
error_rel_py = abs(tab_py - pth[1]) / pth[1]
error_rel_pz = abs(tab_pz - pth[2]) / pth[2]


# ______________________________________________________________________________
# Plot graph

fig, axs = plt.subplots(1,3)
axs[0].plot(iterations, tab_px, marker='+', label='Compute')
axs[0].plot(iterations, pth[0], label='Theoric')
axs[0].set_xlabel("t")
axs[0].set_ylabel("px")
axs[1].plot(iterations, tab_py, marker='+', label='Compute')
axs[1].plot(iterations, pth[1], label='Theoric')
axs[1].set_xlabel("t")
axs[1].set_ylabel("py")
axs[2].plot(iterations, tab_pz, marker='+', label='Compute')
axs[2].plot(iterations, pth[2], label='Theoric')
axs[2].set_xlabel("t")
axs[2].set_ylabel("pz")

axs[0].legend()
axs[1].legend()
axs[2].legend()


fig = plt.figure(figsize=(12, 8))
gs = plt.GridSpec(1, 3)
ax0 = plt.subplot(gs[:,0:1])
ax1 = plt.subplot(gs[:,1:2])
ax2 = plt.subplot(gs[:,2:3])

ax0.plot(iterations, tab_x, marker='+', label='Compute')
ax0.plot(iterations, xth[0], label='Theoric')
ax0.set_xlabel("t")
ax0.set_ylabel("x")

ax1.plot(iterations, tab_y, marker='+', label='Compute')
ax1.plot(iterations, xth[1], label='Theoric')
ax1.set_xlabel("t")
ax1.set_ylabel("y")

ax2.plot(iterations, tab_z, marker='+', label='Compute')
ax2.plot(iterations, xth[2], label='Theoric')
ax2.set_xlabel("t")
ax2.set_ylabel("z")

ax0.legend()
ax1.legend()
ax2.legend()

if True:
    fig2, axs = plt.subplots(1,3)
    axs[0].plot(iterations, error_abs_x, marker='+')
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("error_abs_x")
    axs[1].plot(iterations, error_abs_y, marker='+')
    axs[1].set_xlabel("t")
    axs[1].set_ylabel("error_abs_y")
    axs[2].plot(iterations, error_abs_z, marker='+')
    axs[2].set_xlabel("t")
    axs[2].set_ylabel("error_abs_z")

if True:
    fig3, axs = plt.subplots(1,3)
    axs[0].plot(iterations, error_rel_x, marker='+')
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("error_rel_x")
    axs[1].plot(iterations, error_rel_y, marker='+')
    axs[1].set_xlabel("t")
    axs[1].set_ylabel("error_rel_y")
    axs[2].plot(iterations, error_rel_z, marker='+')
    axs[2].set_xlabel("t")
    axs[2].set_ylabel("error_rel_z")

if True:
    fig4, axs = plt.subplots(1,3)
    axs[0].plot(iterations, error_abs_px, marker='+')
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("error_abs_px")
    axs[1].plot(iterations, error_abs_py, marker='+')
    axs[1].set_xlabel("t")
    axs[1].set_ylabel("error_abs_py")
    axs[2].plot(iterations, error_abs_pz, marker='+')
    axs[2].set_xlabel("t")
    axs[2].set_ylabel("error_abs_pz")

if True:
    fig5, axs = plt.subplots(1,3)
    axs[0].plot(iterations, error_rel_px, marker='+')
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("error_rel_px")
    axs[1].plot(iterations, error_rel_py, marker='+')
    axs[1].set_xlabel("t")
    axs[1].set_ylabel("error_rel_py")
    axs[2].plot(iterations, error_rel_pz, marker='+')
    axs[2].set_xlabel("t")
    axs[2].set_ylabel("error_rel_pz")


plt.show()


