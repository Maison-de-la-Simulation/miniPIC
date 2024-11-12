
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import struct
import sys
import os

path=""

if len( sys.argv ) == 1:
    print( "Usage: python3 print_fields.py path/to/data/directory it" )
    exit()

path = sys.argv[1]
it = int(sys.argv[2])


def read_file(path):
    file = open(path, "rb")

    print(path)
    
    content = file.read()

    k = 0


    # n nodes
    nx = struct.unpack('i',content[k:k+4])[0] ; k+= 4
    ny = struct.unpack('i',content[k:k+4])[0] ; k+= 4
    nz = struct.unpack('i',content[k:k+4])[0] ; k+= 4

    print("nx{}".format(nx))
    print("ny{}".format(ny))
    print("nz{}".format(nz))
    
    size = nx*ny*nz

    Field = np.array(struct.unpack('{}d'.format(size),content[k:k + 8*size])) ; k+= 8*size

    Field_map = np.reshape(Field, (nx, ny, nz))
    print("Map sum = {}".format(Field_map.sum()))
    
    x = np.arange(nx)
    y = np.arange(ny)
    z = np.arange(nz)

    return Field_map, x, y

Ex_map, Ex_x, Ex_y = read_file("{}/Ex_field.{:05d}.data".format(path,it))
Ey_map, Ey_x, Ey_y = read_file("{}/Ey_field.{:05d}.data".format(path,it))
Ez_map, Ez_x, Ez_y = read_file("{}/Ez_field.{:05d}.data".format(path,it))
Bx_map, Bx_x, Bx_y = read_file("{}/Bx_field.{:05d}.data".format(path,it))
By_map, By_x, By_y = read_file("{}/By_field.{:05d}.data".format(path,it))
Bz_map, Bz_x, Bz_y = read_file("{}/Bz_field.{:05d}.data".format(path,it))



def symetric(im):
    clim = np.array(im.get_clim())
    max = 0.8*np.max(np.absolute(clim))
    im.set_clim([-max,max])


fig, axs = plt.subplots(2,3)

im00 = axs[0,0].pcolormesh(Ex_x,Ex_y, Ex_map[:,:,4].T, cmap="RdBu")
im01 = axs[0,1].pcolormesh(Ey_x,Ey_y, Ey_map[:,:,4].T, cmap="RdBu")
im02 = axs[0,2].pcolormesh(Ez_x,Ez_y, Ez_map[:,:,4].T, cmap="RdBu")
im10 = axs[1,0].pcolormesh(Bx_x,Bx_y, Bx_map[:,:,4].T, cmap="RdBu")
im11 = axs[1,1].pcolormesh(By_x,By_y, By_map[:,:,4].T, cmap="RdBu")
im12 = axs[1,2].pcolormesh(Bz_x,Bz_y, Bz_map[:,:,4].T, cmap="RdBu")

#symetric(im00)
#symetric(im01)
#symetric(im02)
#symetric(im10)
#symetric(im11)
#symetric(im12)

cb00 = plt.colorbar(im00, ax=axs[0,0])
cb01 = plt.colorbar(im01, ax=axs[0,1])
cb02 = plt.colorbar(im02, ax=axs[0,2])
cb10 = plt.colorbar(im10, ax=axs[1,0])
cb11 = plt.colorbar(im11, ax=axs[1,1])
cb12 = plt.colorbar(im12, ax=axs[1,2])

axs[0,0].set_xlabel("x")
axs[0,0].set_ylabel("y")
axs[0,0].set_title("Ex")
axs[0,1].set_xlabel("x")
axs[0,1].set_ylabel("y")
axs[0,1].set_title("Ey")
axs[0,2].set_xlabel("x")
axs[0,2].set_ylabel("y")
axs[0,2].set_title("Ez")

axs[1,0].set_xlabel("x")
axs[1,0].set_ylabel("y")
axs[1,0].set_title("Bx")
axs[1,1].set_xlabel("x")
axs[1,1].set_ylabel("y")
axs[1,1].set_title("By")
axs[1,2].set_xlabel("x")
axs[1,2].set_ylabel("y")
axs[1,2].set_title("Bz")


plt.show()
