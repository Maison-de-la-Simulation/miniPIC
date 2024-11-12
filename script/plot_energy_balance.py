# ______________________________________________________________________________
#
# Plot the energy balance
#
# ______________________________________________________________________________

import sys
import numpy as np
import glob
from matplotlib import *
from matplotlib.pyplot import *
import argparse

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

parser = argparse.ArgumentParser(description='Print 3D diag')

parser.add_argument('-f', '--folder', type=str, help="Path toward diags folder", default="diags")

args = parser.parse_args()

# Utiliser les arguments
print("Diags folder path: ", args.folder)

file_path = args.folder

scalars = {}

# ______________________________________________________________________________
# Read the species scalar files

# Look for all species scalar files
species_files = glob.glob(file_path + '/species_*')

for species_file in species_files:

    print(species_file)

    name = species_file[2:].split(".")[0]

    scalars[name] =  {}

    file = open( species_file, 'r')

    lines = file.readlines()

    number_of_lines = len(lines)

    scalars[name]["iterations"] = np.zeros(number_of_lines-1)
    scalars[name]["number_of_particles"] = np.zeros(number_of_lines-1)
    scalars[name]["kinetic_energy"] = np.zeros(number_of_lines-1)

    for iline, line in enumerate(lines[1:]):
        content = line.split(" ")
        scalars[name]["iterations"][iline] = int(content[0])
        scalars[name]["number_of_particles"][iline] = int(content[1])
        scalars[name]["kinetic_energy"][iline] = float(content[2])

# ______________________________________________________________________________
# Read the field scalars

file = open(file_path + '/fields.txt', 'r')

lines = file.readlines()

number_of_lines = len(lines)
number_of_iterations = number_of_lines-2

dt = float((lines[0].split(":"))[1])

print(dt)

scalars["iterations"] = np.zeros(number_of_iterations)

scalars["Ex"] = np.zeros(number_of_iterations)
scalars["Ey"] = np.zeros(number_of_iterations)
scalars["Ez"] = np.zeros(number_of_iterations)

scalars["Bx"] = np.zeros(number_of_iterations)
scalars["By"] = np.zeros(number_of_iterations)
scalars["Bz"] = np.zeros(number_of_iterations)

for iline, line in enumerate(lines[2:]):
    content = line.split(" ")
    scalars["iterations"][iline] = int(content[0])
    scalars["Ex"][iline] = float(content[1])
    scalars["Ey"][iline] = float(content[2])
    scalars["Ez"][iline] = float(content[3])
    scalars["Bx"][iline] = float(content[4])
    scalars["By"][iline] = float(content[5])
    scalars["Bz"][iline] = float(content[6])

# ______________________________________________________________________________
# Total energy 

scalars["total_energy"] = scalars["Ex"] + scalars["Ey"] + scalars["Ez"] + scalars["Bx"] + scalars["By"] + scalars["Bz"]

for species_file in species_files:

    name = species_file[2:].split(".")[0]

    scalars["total_energy"]  += scalars[name]["kinetic_energy"]

# ______________________________________________________________________________
# Figure and plot

fig = figure(figsize=(12, 6))

gs = GridSpec(3, 3)
ax0 = subplot(gs[:,:])

# for species_file in species_files:

#     name = species_file[2:].split(".")[0]

#     im0 = ax0.plot(scalars[name]["iterations"],scalars[name]["number_of_particles"])

for species_file in species_files:

    name = species_file[2:].split(".")[0]

    im0 = ax0.plot(scalars[name]["iterations"]*dt,scalars[name]["kinetic_energy"],label=name)

im0 = ax0.plot(scalars["iterations"]*dt,scalars["Ex"],label="Ex")
im0 = ax0.plot(scalars["iterations"]*dt,scalars["Ey"],label="Ey")
im0 = ax0.plot(scalars["iterations"]*dt,scalars["Ez"],label="Ez")

im0 = ax0.plot(scalars["iterations"]*dt,scalars["Bx"],label="Bx")
im0 = ax0.plot(scalars["iterations"]*dt,scalars["By"],label="By")
im0 = ax0.plot(scalars["iterations"]*dt,scalars["Bz"],label="Bz")

im0 = ax0.plot(scalars["iterations"]*dt,scalars["total_energy"] ,label="total_energy", color='k')

ax0.set_xlabel("Time")
ax0.set_ylabel("Energy")

ax0.set_yscale('log')

#ax0.set_title("Plants {}".format(iteration))

ax0.legend(loc="best")

fig.tight_layout()

show()