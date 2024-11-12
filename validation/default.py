# ______________________________________________________________________
#
# Validation script for the `default` benchmark
#
# ______________________________________________________________________

import os
import lib.minipic_diag
import lib.minipic_ci as minipic_ci
import numpy as np

def validate(evaluate=True, threshold = 1e-10):

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

    # Check final scalar values for fields

    reference_data = [300, 8.676516835993083e-16, 8.308212089099558e-16, 8.176656652584353e-16, 1.06349695228065e-15, 3.156277623217432e-15, 2.740691983778254e-15]


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

    if (evaluate):
        minipic_ci.evaluate(iteration, reference_data[0], reference_data[0], '==',  'Last iteration in fields.txt is not correct')

    print("Field values: \n {}, {}, {}, {}, {}, {}".format(Ex, Ey, Ez, Bx, By, Bz))

    # We check that the fields do not explode (numerical instability)
    if (evaluate):
        minipic_ci.evaluate(Ex, reference_data[1], threshold, 'relative', 'Ex value in fields.txt is not correct')   
        minipic_ci.evaluate(Ey, reference_data[2], threshold, 'relative', 'Ey value in fields.txt is not correct')  
        minipic_ci.evaluate(Ez, reference_data[3], threshold, 'relative', 'Ez value in fields.txt is not correct')  
        minipic_ci.evaluate(Bx, reference_data[4], threshold, 'relative', 'Bx value in fields.txt is not correct')  
        minipic_ci.evaluate(By, reference_data[5], threshold, 'relative', 'By value in fields.txt is not correct')  
        minipic_ci.evaluate(Bz, reference_data[6], threshold, 'relative', 'Bz value in fields.txt is not correct')  

    # Check initial scalar for species

    reference_data = [ [0, 262144, 1.499802569489304e-07], [0, 262144, 1.500982120183948e-07]]

    for ispecies in range(2):

        with open("diags/species_{:02d}.txt".format(ispecies), "r") as f:
            lines = f.readlines()

            last_line = lines[1].split(" ")

            iteration = int(last_line[0])
            particles = float(last_line[1])
            energy = float(last_line[2])

        print("Species {}: {}, {}, {}".format(ispecies, iteration, particles, energy))

        if (evaluate):
            minipic_ci.evaluate(iteration, reference_data[ispecies][0], reference_data[ispecies][0], '==', 'First iteration in species_{}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(particles, reference_data[ispecies][1], reference_data[ispecies][1], '==', 'Number of particles in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(energy, reference_data[ispecies][2], threshold, 'relative', 'Kinetic energy in species_{:02d}.txt is not correct'.format(ispecies)) 

    # Check final scalar for species

    reference_data = [ [300, 262144, 1.499801491562789e-07], [300, 262144, 1.500981343380884e-07]]

    for ispecies in range(2):

        with open("diags/species_{:02d}.txt".format(ispecies), "r") as f:
            lines = f.readlines()

            last_line = lines[-1].split(" ")

            iteration = int(last_line[0])
            particles = float(last_line[1])
            energy = float(last_line[2])

        print("Species {}: {}, {}, {}".format(ispecies, iteration, particles, energy))

        if (evaluate):
            minipic_ci.evaluate(iteration, reference_data[ispecies][0], reference_data[ispecies][0], '==', 'Last iteration in species_{}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(particles, reference_data[ispecies][1], reference_data[ispecies][1], '==', 'Number of particles in species_{:02d}.txt is not correct'.format(ispecies))

            minipic_ci.evaluate(energy, reference_data[ispecies][2], threshold, 'relative', 'Kinetic energy in species_{:02d}.txt is not correct'.format(ispecies)) 

    # ______________________________________________________________________
    # Check gamma spectrums

    reference_sum_data = [[1.0150306498218023e-05, 1.0150305460964983e-05, 1.0150309455244676e-05, 1.0150307973949819e-05, 1.0150311383418524e-05, 1.015030659799281e-05], 
                          [1.000157589162598e-05, 1.0001575891509972e-05, 1.0001575891335875e-05, 1.0001575891120025e-05, 1.0001575890858374e-05, 1.0001575890629277e-05]]

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

    reference_data = [[131076.16027212347, 131072.11534091394, 131069.52464953832, 21014.37699146685, 21028.564131718496, 21014.48000786643],
                      [131076.16027212347, 131072.11534091394, 131069.52464953832, 488.33937456089893, 488.37513744120105, 488.27610367103927]]

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

        print("    - For Species {}: {}, {}, {}, {}, {}, {}".format(ispecies, x_sum, y_sum, z_sum, px_sum, py_sum, pz_sum))

        if (evaluate):
            minipic_ci.evaluate(x_sum, reference_data[ispecies][0], 1e-9, 'relative', 'Sum over x positions not similar')
            minipic_ci.evaluate(y_sum, reference_data[ispecies][1], 1e-9, 'relative', 'Sum over y positions not similar')
            minipic_ci.evaluate(z_sum, reference_data[ispecies][2], 1e-9, 'relative', 'Sum over z positions not similar')

            minipic_ci.evaluate(px_sum, reference_data[ispecies][3], 1e-9, 'relative', 'Sum over px not similar')
            minipic_ci.evaluate(py_sum, reference_data[ispecies][4], 1e-9, 'relative', 'Sum over py not similar')
            minipic_ci.evaluate(pz_sum, reference_data[ispecies][5], 1e-9, 'relative', 'Sum over pz not similar')

    # ______________________________________________________________________
    # Check final cloud file (particle initialization)

    reference_data = [[131052.0194511577, 131022.16752569123, 131145.1132868646, 21014.411780135964, 21028.626451922173, 21014.459008521262],
                      [131100.95740321468, 130988.17619900177, 131067.20581941119, 488.34032316893814, 488.37354869420835, 488.2763679441135]]

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
            minipic_ci.evaluate(x_sum, reference_data[ispecies][0], 1e-9, 'relative', 'Sum over x positions not similar')
            minipic_ci.evaluate(y_sum, reference_data[ispecies][1], 1e-9, 'relative', 'Sum over y positions not similar')
            minipic_ci.evaluate(z_sum, reference_data[ispecies][2], 1e-9, 'relative', 'Sum over z positions not similar')

            minipic_ci.evaluate(px_sum, reference_data[ispecies][3], 1e-9, 'relative', 'Sum over px not similar')
            minipic_ci.evaluate(py_sum, reference_data[ispecies][4], 1e-9, 'relative', 'Sum over py not similar')
            minipic_ci.evaluate(pz_sum, reference_data[ispecies][5], 1e-9, 'relative', 'Sum over pz not similar')


