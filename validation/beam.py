# ______________________________________________________________________
#
# Validation script for the `beam` benchmark
#
# ______________________________________________________________________

import os
import math
import lib.minipic_diag
import lib.minipic_ci as minipic_ci
import numpy as np

def validate(evaluate=True, threshold=1e-10):

    # ______________________________________________________________________
    # Check output files are created

    # list of output files
    output_file_list = []

    # Add *.vtk files
    for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz", 
                  "diag_x_y_z_d_s00", "diag_x_y_z_d_s01"]:
        for it in range(0,550,50):
            file = "{}_{:03d}.vtk".format(field, it)
            output_file_list.append(file)

    # Add *.bin files
    for field in ["diag_w_gamma_s00", "diag_w_gamma_s01"]:
        for it in range(0,550,50):
            file = "{}_{:03d}.bin".format(field, it)
            output_file_list.append(file)

    # Add scalars
    output_file_list.append("fields.txt")
    output_file_list.append("species_00.txt")
    output_file_list.append("species_01.txt")

    # Check that all output files exist
    for file in output_file_list:
        if (not(os.path.exists("diags/"+file))):
            minipic_ci.error('File {} not generated'.format(file))
        
    # ______________________________________________________________________
    # Check scalars

    print(" > Check scalars")



    # Check final scalar values for fields

    reference_data = [500, 4.124448552038559e-11, 7.163173980414326e-11, 7.844505241044072e-11, 1.0753209235313e-13, 8.019575880396739e-11, 7.667103357828862e-11]

    with open("diags/fields.txt", "r") as f:
        lines = f.readlines()

        last_line = lines[-1].split(" ")

        iteration = int(last_line[0])
        Ex = float(last_line[1])
        Ey = float(last_line[2])
        Ez = float(last_line[3])
        Bx = float(last_line[4])
        By = float(last_line[5])
        Bz = float(last_line[6])

    # print all field values
    print("   - Final field scalar ({}): Ex = {}, Ey = {}, Ez = {}, Bx = {}, By = {}, Bz = {}".format(iteration, Ex, Ey, Ez, Bx, By, Bz))

    if (evaluate):

        minipic_ci.evaluate(iteration, reference_data[0], reference_data[0], '==', 'Last iteration in fields.txt is not correct')

        minipic_ci.evaluate(Ex, reference_data[1], threshold, 'relative', 'Ex value at final iteration in fields.txt is not correct')  
        minipic_ci.evaluate(Ey, reference_data[2], threshold, 'relative', 'Ey value at final iteration in fields.txt is not correct')  
        minipic_ci.evaluate(Ez, reference_data[3], threshold, 'relative', 'Ez value at final iteration in fields.txt is not correct')  
        minipic_ci.evaluate(Bx, reference_data[4], threshold, 'relative', 'Bx value at final iteration in fields.txt is not correct')  
        minipic_ci.evaluate(By, reference_data[5], threshold, 'relative', 'By value at final iteration in fields.txt is not correct')  
        minipic_ci.evaluate(Bz, reference_data[6], threshold, 'relative', 'Bz value at final iteration in fields.txt is not correct')

    else:

        print("    - reference_data = [{}, {}, {}, {}, {}, {}, {}]".format(iteration, Ex, Ey, Ez, Bx, By, Bz)) 

    # Check initial scalar for species

    reference_data = [[0, 17171, 1.486505893401623e-05], 
                      [0, 17171, 0.02727053664199493]]

    for ispecies in range(2):

        with open("diags/species_{:02d}.txt".format(ispecies), "r") as f:
            lines = f.readlines()

            last_line = lines[1].split(" ")

            iteration = int(last_line[0])
            particles = float(last_line[1])
            energy = float(last_line[2])

        print("    - Initial scalar for species {}: {}, {}, {}".format(ispecies, iteration, particles, energy))

        if (evaluate):

            if (reference_data[ispecies][0] != iteration):
                minipic_ci.error('First iteration in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(particles, reference_data[ispecies][1], reference_data[ispecies][1], '==', 'Number of particles in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(energy, reference_data[ispecies][2], threshold, 'relative', 'Kinetic energy in species_{:02d}.txt is not correct'.format(ispecies)) 

    # Check final scalar for species

    reference_data = [[500, 17171, 1.315416203546459e-05], [500, 17171, 0.02727052464789659]]

    for ispecies in range(2):

        with open("diags/species_{:02d}.txt".format(ispecies), "r") as f:
            lines = f.readlines()

            last_line = lines[-1].split(" ")

            iteration = int(last_line[0])
            particles = float(last_line[1])
            energy = float(last_line[2])

        print("    - Final scalar for species {}: {}, {}, {}".format(ispecies, iteration, particles, energy))

        if (evaluate):

            if (reference_data[ispecies][0] != iteration):
                minipic_ci.error('Last iteration in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(particles, reference_data[ispecies][1], reference_data[ispecies][1], '==', 'Number of particles in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(energy, reference_data[ispecies][2], threshold, 'relative', 'Kinetic energy in species_{:02d}.txt is not correct'.format(ispecies)) 

    # ______________________________________________________________________
    # Check gamma spectrums

    reference_sum_data = [
        [2.6341788307225324e-05, 2.6340193071521304e-05, 2.630707755958522e-05, 2.620567052764382e-05, 2.5977739794666403e-05, 2.5745016811351395e-05, 2.544705187191869e-05, 2.516323961620942e-05, 2.4961479614754423e-05, 2.4719498510823074e-05],
        [2.6328534922736266e-05, 2.632853714034016e-05, 2.6328533189048132e-05, 2.6328544284653243e-05, 2.6328542661491538e-05, 2.6328563090011446e-05, 2.6328573517162e-05, 2.632857026984209e-05, 2.6328586809892927e-05, 2.6328556641940606e-05]
    ]

    print(" > Checking gamma spectrums")

    for ispecies in range(2):

        new_data = []

        for i,it in enumerate(range(0,500,50)):

            file = "diag_w_gamma_s{:02d}_{:03d}.bin".format(ispecies,it)

            x_axis_name, x_min, x_max, x_data, data_name, data = lib.minipic_diag.read_1d_diag("diags/" + file)

            new_data.append(np.sum(data * x_data))

        print("    - For species {}: {}".format(ispecies, new_data))

        if (evaluate):

            for i,it in enumerate(range(0,500,50)):
                minipic_ci.evaluate(new_data[i], reference_sum_data[ispecies][i], 1e-9, 'relative', 'Gamma spectrum not similar at iteration {} for species {}'.format(it, ispecies))







