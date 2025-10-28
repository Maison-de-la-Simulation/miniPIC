# ______________________________________________________________________
#
# Validation script for the `default` benchmark
#
# ______________________________________________________________________

import os
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
                  "diag_x_y_z_d_s00", "diag_x_y_z_d_s01", 
                  "diag_px_py_pz_d_s00", "diag_px_py_pz_d_s01"]:
        for it in range(0,300,50):
            file = "{}_{:03d}.vtk".format(field, it)
            output_file_list.append(file)

    # Add *.bin files
    for field in ["cloud_s00", "cloud_s01", 
                  "diag_w_gamma_s00", "diag_w_gamma_s01"]:
        for it in range(0,300,50):
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

    print(" > Checking scalars")

    # Check initial scalar for species

    reference_data = [ [0, 262144, 1.499802569489315e-07], [0, 262144, 1.500982120183934e-07]]

    for ispecies in range(2):

        with open("diags/species_{:02d}.txt".format(ispecies), "r") as f:
            lines = f.readlines()

            last_line = lines[1].split(" ")

            iteration = int(last_line[0])
            particles = float(last_line[1])
            energy = float(last_line[2])

        print("    - Initial scalar for species {}: {}, {}, {}".format(ispecies, iteration, particles, energy))

        if (evaluate):

            minipic_ci.evaluate(iteration, reference_data[ispecies][0], reference_data[ispecies][0], '==', 'First iteration in species_{}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(particles, reference_data[ispecies][1], reference_data[ispecies][1], '==', 'Number of particles in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(energy, reference_data[ispecies][2], threshold, 'relative', 'Kinetic energy in species_{:02d}.txt is not correct'.format(ispecies)) 

    # Check final scalar for species

    reference_data = [ [300, 262144, 1.49980167722712e-07], 
                       [300, 262144, 1.500981415474678e-07]]

    for ispecies in range(2):

        with open("diags/species_{:02d}.txt".format(ispecies), "r") as f:
            lines = f.readlines()

            last_line = lines[-1].split(" ")

            iteration = int(last_line[0])
            particles = float(last_line[1])
            energy = float(last_line[2])

        print("    - Final scalar for species {}: {}, {}, {}".format(ispecies, iteration, particles, energy))

        if (evaluate):

            minipic_ci.evaluate(iteration, reference_data[ispecies][0], reference_data[ispecies][0], '==', 'Last iteration in species_{}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(particles, reference_data[ispecies][1], reference_data[ispecies][1], '==', 'Number of particles in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(energy, reference_data[ispecies][2], threshold, 'relative', 'Kinetic energy in species_{:02d}.txt is not correct'.format(ispecies)) 
        

    # Check final scalar values for fields

    reference_data = [300, 8.707252968548611e-16, 8.187202771374761e-16, 8.133778516563632e-16, 1.014322708204923e-15, 2.738517219101327e-15, 2.932533327560403e-15]

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

    print("    - Field values at final iteration {}: \n {}, {}, {}, {}, {}, {}".format(iteration, Ex, Ey, Ez, Bx, By, Bz))

    if (evaluate):

        minipic_ci.evaluate(iteration, reference_data[0], reference_data[0], '==',  'Last iteration in fields.txt is not correct')

        # We check that the fields do not explode (numerical instability)
        minipic_ci.evaluate(Ex, reference_data[1], threshold, 'relative', 'Ex value at it {} in fields.txt is not correct'.format(iteration))   
        minipic_ci.evaluate(Ey, reference_data[2], threshold, 'relative', 'Ey value at it {} in fields.txt is not correct'.format(iteration))
        minipic_ci.evaluate(Ez, reference_data[3], threshold, 'relative', 'Ez value at it {} in fields.txt is not correct'.format(iteration))
        minipic_ci.evaluate(Bx, reference_data[4], threshold, 'relative', 'Bx value at it {} in fields.txt is not correct'.format(iteration))
        minipic_ci.evaluate(By, reference_data[5], threshold, 'relative', 'By value at it {} in fields.txt is not correct'.format(iteration))
        minipic_ci.evaluate(Bz, reference_data[6], threshold, 'relative', 'Bz value at it {} in fields.txt is not correct'.format(iteration))
        


    # ______________________________________________________________________
    # Check gamma spectrums

    reference_sum_data = [[1.0150306498218023e-05, 1.0150305070953287e-05, 1.015031026575343e-05, 1.0150307734880672e-05, 1.0150311192601487e-05, 1.015030843474158e-05], 
                          [1.000157589162598e-05, 1.0001575891523015e-05, 1.000157589136796e-05, 1.000157589118226e-05, 1.0001575890954698e-05, 1.0001575890766535e-05]]

    print(" > Checking gamma spectrums")

    for ispecies in range(2):

        new_data = []

        for i,it in enumerate(range(0,300,50)):

            file = "diag_w_gamma_s{:02d}_{:03d}.bin".format(ispecies,it)

            x_axis_name, x_min, x_max, x_data, data_name, data = lib.minipic_diag.read_1d_diag("diags/" + file)

            new_data.append(np.sum(np.abs(data * x_data)))

        print("    - For species {}: {}".format(ispecies, new_data))

        if (evaluate):
            for i,it in enumerate(range(0,300,50)):
                minipic_ci.evaluate(new_data[i], reference_sum_data[ispecies][i], 1e-13, 'relative', 'Gamma spectrum for species {} at iteration {} not similar'.format(ispecies,it))

    # ______________________________________________________________________
    # Check initial cloud file (particle initialization)

    reference_data = [[131076.1602721234, 131072.11534091408, 131069.52464953822, 21014.376991466845, 21028.564131718493, 21014.480007866434],
                      [131076.1602721234, 131072.11534091408, 131069.52464953822, 488.33937456089905, 488.375137441201, 488.2761036710393]]

    print(" > Checking initial cloud file")

    for ispecies in range(2):

        file = "cloud_s{:02}_000.bin".format(ispecies)

        particle_number, id, w, x, y, z, px, py, pz = lib.minipic_diag.read_particle_cloud("diags/"+file)

        x_sum = np.sum(np.abs(x))
        y_sum = np.sum(np.abs(y))
        z_sum = np.sum(np.abs(z))

        px_sum = np.sum(np.abs(px))
        py_sum = np.sum(np.abs(py))
        pz_sum = np.sum(np.abs(pz))

        print("   - For Species {}: {}, {}, {}, {}, {}, {}".format(ispecies, x_sum, y_sum, z_sum, px_sum, py_sum, pz_sum))

        if (evaluate):
            minipic_ci.evaluate(x_sum, reference_data[ispecies][0], threshold, 'relative', 'Sum over x positions not similar')
            minipic_ci.evaluate(y_sum, reference_data[ispecies][1], threshold, 'relative', 'Sum over y positions not similar')
            minipic_ci.evaluate(z_sum, reference_data[ispecies][2], threshold, 'relative', 'Sum over z positions not similar')

            minipic_ci.evaluate(px_sum, reference_data[ispecies][3], threshold, 'relative', 'Sum over px not similar')
            minipic_ci.evaluate(py_sum, reference_data[ispecies][4], threshold, 'relative', 'Sum over py not similar')
            minipic_ci.evaluate(pz_sum, reference_data[ispecies][5], threshold, 'relative', 'Sum over pz not similar')

    # ______________________________________________________________________
    # Check final cloud file (particle initialization)

    reference_data = [[131059.99739212358, 131023.15182992323, 131142.13555391456, 21014.42145328789, 21028.616998630845, 21014.458067973428],
                      [131100.95742997795, 130988.17625004012, 131067.20585065625, 488.3403714801243, 488.3735558866029, 488.2763737032677]]

    print(" > Checking final cloud file")

    for ispecies in range(2):

        file = "cloud_s{:02}_300.bin".format(ispecies)

        particle_number, id, w, x, y, z, px, py, pz = lib.minipic_diag.read_particle_cloud("diags/"+file)

        x_sum = np.sum(np.abs(x))
        y_sum = np.sum(np.abs(y))
        z_sum = np.sum(np.abs(z))

        px_sum = np.sum(np.abs(px))
        py_sum = np.sum(np.abs(py))
        pz_sum = np.sum(np.abs(pz))

        print("    - For Species {}: {}, {}, {}, {}, {}, {}".format(ispecies, x_sum, y_sum, z_sum, px_sum, py_sum, pz_sum))

        if (evaluate):
            minipic_ci.evaluate(x_sum, reference_data[ispecies][0], threshold, 'relative', 'Sum over x positions not similar')
            minipic_ci.evaluate(y_sum, reference_data[ispecies][1], threshold, 'relative', 'Sum over y positions not similar')
            minipic_ci.evaluate(z_sum, reference_data[ispecies][2], threshold, 'relative', 'Sum over z positions not similar')

            minipic_ci.evaluate(px_sum, reference_data[ispecies][3], threshold, 'relative', 'Sum over px not similar')
            minipic_ci.evaluate(py_sum, reference_data[ispecies][4], threshold, 'relative', 'Sum over py not similar')
            minipic_ci.evaluate(pz_sum, reference_data[ispecies][5], threshold, 'relative', 'Sum over pz not similar')



