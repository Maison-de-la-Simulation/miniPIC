# ______________________________________________________________________________
#
# Read and plot the particles at a specific iteration
#
# ______________________________________________________________________________

import sys
import struct
import numpy as np
from matplotlib import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
import lib.minipic_diag as minipic_diag

# from mayavi import mlab

# ______________________________________________________________________________
# RCParams - personalize the figure output

rcParams['figure.facecolor'] = 'w'
rcParams['font.size'] = 15
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['axes.labelsize'] = 15

rcParams['xtick.major.size'] = 10
rcParams['ytick.major.size'] = 10

rcParams['xtick.minor.size'] = 5
rcParams['ytick.minor.size'] = 5

rcParams['axes.linewidth'] = 1.5

rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5

# ______________________________________________________________________________
# Read command line arguments

parser = argparse.ArgumentParser(description='Plot particle cloud')

parser.add_argument('-f', '--file', type=str, help="Path toward the diag cloud file")
# parser.add_argument('-c', '--colormap', type=str, help="colormap to use", default="plasma")

args = parser.parse_args()

file_path = args.file

# ______________________________________________________________________________
# Read the binary file

particle_number, id, w, x, y, z, px, py, pz = minipic_diag.read_particle_cloud(file_path)

radius = 0.01

print(" Number of particles: {}".format(particle_number))


# vx = np.array(struct.unpack('{}d'.format(particle_number),content[k:k + 8*particle_number])) ; k+= 8*particle_number
# vy = np.array(struct.unpack('{}d'.format(particle_number),content[k:k + 8*particle_number])) ; k+= 8*particle_number
# vz = np.array(struct.unpack('{}d'.format(particle_number),content[k:k + 8*particle_number])) ; k+= 8*particle_number

# mass = np.array(struct.unpack('{}d'.format(particle_number),content[k:k + 8*particle_number])) ; k+= 8*particle_number

# energy = 0.5 * mass * ( vx*vx + vy*vy + vz*vz)

# ______________________________________________________________________________
# Figure and plot using Matplotlib

fig = figure(figsize=(12, 8))
ax = fig.add_subplot(projection='3d')

im = ax.scatter(x, y, z, c=id, marker='o', alpha = 1) # s = np.sqrt(radius)
cb = colorbar(im)
cb.set_label("Energy")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax.set_title("Particles from {}".format(file_path))

fig.tight_layout()

# ______________________________________________________________________________
# Figure and plot using Mayavi

# fig = mlab.figure(size=(600,600))

# glyph = mlab.points3d(x, y, z, energy, resolution=16) #scale_factor=radius

# mlab.axes()
# mlab.colorbar()

# mlab.show()

show()
