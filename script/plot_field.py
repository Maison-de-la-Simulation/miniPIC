# __________________________________________________________
#
# This script is used to plot the particle binning diag
# __________________________________________________________


#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation

import struct
import sys
import os

import argparse

import lib.minipic_diag

path=""

# Use argparse to parse arguments
# -h for help
# -f for file name
# -c for colormap

parser = argparse.ArgumentParser(description='Plot particle binning diag')

parser.add_argument('-f', '--file', type=str, help="Path toward the diag binning file")
parser.add_argument('-c', '--colormap', type=str, help="colormap to use", default="seismic")

args = parser.parse_args()

# Check if the file exists
if not os.path.isfile(args.file):
    print("File does not exist")
    sys.exit()

# Get the binning dimension
dim = lib.minipic_diag.get_diag_dimension(args.file)

# Get the iteration number from the file name (diags_yyyyy.bin)
iteration = int(args.file.split('_')[-1].split('.')[0])

# Utiliser les arguments
print(" File path: ", args.file)
print(" Colormap: ", args.colormap)
print(" Iteration: ", iteration)

x_axis, x, y_axis, y, z_axis, z, data_name, data_map = lib.minipic_diag.read_3d_diag(args.file)

xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]
zmin = z[0]
zmax = z[-1]

print(" data name : {}".format(data_name))
print(" x_axis : {}".format(x_axis))
print(" x_min : {}".format(xmin))
print(" x_max : {}".format(xmax))
x_ncells = len(x)
print(" x_ncells : {}".format(x_ncells))
# print(" x {}".format(x))

print(" y_axis : {}".format(y_axis))
print(" y_min : {}".format(ymin))
print(" y_max : {}".format(ymax))
y_ncells = len(y)
print(" y_ncells : {}".format(y_ncells))
# print(" y {}".format(y))

print(" z_axis : {}".format(z_axis))
print(" z_min : {}".format(zmin))
print(" z_max : {}".format(zmax))
z_ncells = len(z)
print(" z_ncells : {}".format(z_ncells))
# print(" z {}".format(z))

fig = plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(4, 4)
axs = fig.add_subplot(gs[0:4, 0:4])

img=axs.pcolormesh(x, y, data_map[:,:,int(z_ncells/2)].T, cmap=args.colormap)

axs.set_title("Field {} at iteration {}".format(data_name, iteration))

axs.set_xlabel(x_axis)
axs.set_ylabel(y_axis)

fig.colorbar(img, ax=axs)

#axs.set_xscale("log")
#axs.set_yscale("log")
#axs.set_zscale("log")

fig.tight_layout()

plt.show()
