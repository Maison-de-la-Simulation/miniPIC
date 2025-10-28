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
    # Reference data

    reference_dict = {'scalar': {'fields': [500, 4.176915389526767e-11, 6.84235342733258e-11, 7.384372445331327e-11, 1.490359461146354e-13, 7.824479096513594e-11, 7.512299160674653e-11], 'species': [[[0, 17171.0, 1.486505893401618e-05], [0, 17171.0, 0.02727053664199492]], [[500, 17171.0, 1.315434546122977e-05], [500, 17171.0, 0.02727052355889318]]]}, 'gamma_spectrum': [[2.6341788307225324e-05, 2.6340187152554154e-05, 2.6306898149346666e-05, 2.6206688009942488e-05, 2.5975705446921505e-05, 2.5749079793833593e-05, 2.5452293402771052e-05, 2.517019636457644e-05, 2.4977679280317923e-05, 2.4735294342510844e-05], [2.6328534922736266e-05, 2.632853644089038e-05, 2.6328533755778843e-05, 2.6328544155090046e-05, 2.6328544972441783e-05, 2.6328567899607564e-05, 2.632856920474773e-05, 2.6328571264628347e-05, 2.632857386819752e-05, 2.6328558916765818e-05]]}

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

    reference_data = reference_dict["scalar"]["fields"]

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

        reference_dict["scalar"]["fields"] = [iteration, Ex, Ey, Ez, Bx, By, Bz]

    # Check initial scalar for species

    reference_data = reference_dict["scalar"]["species"][0]

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

        else:

            reference_dict["scalar"]["species"][0][ispecies] = [iteration, particles, energy]

    # Check final scalar for species

    reference_data = reference_dict["scalar"]["species"][1]

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

        else:

            reference_dict["scalar"]["species"][1][ispecies] = [iteration, particles, energy]

    # ______________________________________________________________________
    # Check gamma spectrums

    reference_sum_data = reference_dict["gamma_spectrum"]

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

        else:
                
                reference_dict["gamma_spectrum"][ispecies] = new_data

    if (not(evaluate)):

        print("Reference data:")

        # import json
        # print(json.dumps(reference_dict, indent=4))
        print(reference_dict)







