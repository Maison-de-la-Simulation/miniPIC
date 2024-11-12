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
parser.add_argument('-c', '--colormap', type=str, help="colormap to use", default="plasma")

args = parser.parse_args()

# Utiliser les arguments
print(" File path: ", args.file)
print(" Colormap: ", args.colormap)

# Check if the file exists
if not os.path.isfile(args.file):
    print("File does not exist")
    sys.exit()

# Get the binning dimension

dim = lib.minipic_diag.get_diag_dimension(args.file);

print(" diag dimension : {}".format(dim))

if (dim == 1):

    # Read file 1D
    x_axis, x_min, x_max, x, data_name, data = lib.minipic_diag.read_1d_diag(args.file)

    print(" data name : {}".format(data_name))
    print(" x_axis : {}".format(x_axis))
    print(" x_min : {}".format(x_min))
    print(" x_max : {}".format(x_max))
    x_ncells = len(x)
    print(" x_ncells : {}".format(x_ncells))
    print(" x {}".format(x))

    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(4, 4)
    axs = fig.add_subplot(gs[0:4, 0:4])

    axs.plot(x, data)

    axs.set_xlabel(x_axis)
    axs.set_ylabel(data_name)

elif (dim == 2):

    x_axis, x, y_axis, y, data_name, data_map = lib.minipic_diag.read_3d_diag(args.file)

    print(" data name : {}".format(data_name))
    print(" x_axis : {}".format(x_axis))
    # print(" x_min : {}".format(x_min))
    # print(" x_max : {}".format(x_max))
    x_ncells = len(x)
    print(" x_ncells : {}".format(x_ncells))
    print(" x {}".format(x))

    print(" y_axis : {}".format(y_axis))
    # print(" y_min : {}".format(y_min))
    # print(" y_max : {}".format(y_max))
    y_ncells = len(y)
    print(" y_ncells : {}".format(y_ncells))
    print(" y {}".format(y))

    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(4, 4)
    axs = fig.add_subplot(gs[0:4, 0:4])

    img=axs.pcolormesh(x, y, data_map[:,:].T, cmap=args.colormap)

    axs.set_xlabel(x_axis)
    axs.set_ylabel(y_axis)

    fig.colorbar(img, ax=axs)

    #axs.set_xscale("log")
    #axs.set_yscale("log")
    #axs.set_zscale("log")

    fig.tight_layout()

elif (dim == 3):

    x_axis, x, y_axis, y, z_axis, z, data_name, data_map = lib.minipic_diag.read_3d_diag(args.file)

    print(" data name : {}".format(data_name))
    print(" x_axis : {}".format(x_axis))
    # print(" x_min : {}".format(x_min))
    # print(" x_max : {}".format(x_max))
    x_ncells = len(x)
    print(" x_ncells : {}".format(x_ncells))
    print(" x {}".format(x))

    print(" y_axis : {}".format(y_axis))
    # print(" y_min : {}".format(y_min))
    # print(" y_max : {}".format(y_max))
    y_ncells = len(y)
    print(" y_ncells : {}".format(y_ncells))
    print(" y {}".format(y))

    print(" z_axis : {}".format(z_axis))
    # print(" z_min : {}".format(z_min))
    # print(" z_max : {}".format(z_max))
    z_ncells = len(z)
    print(" z_ncells : {}".format(z_ncells))
    print(" z {}".format(z))

    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(4, 4)
    axs = fig.add_subplot(gs[0:4, 0:4])

    img=axs.pcolormesh(x, y, data_map[:,:,int(z_ncells/2)].T, cmap=args.colormap)

    axs.set_xlabel(x_axis)
    axs.set_ylabel(y_axis)

    fig.colorbar(img, ax=axs)

    #axs.set_xscale("log")
    #axs.set_yscale("log")
    #axs.set_zscale("log")

    fig.tight_layout()

plt.show()
