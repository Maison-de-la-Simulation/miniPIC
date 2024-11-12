# ______________________________________________________________________
#
# Validation script for the `Bcst` benchmark
#
# ______________________________________________________________________

import os
import numpy as np
import lib.minipic_ci as minipic_ci
import lib.minipic_diag as minipic_diag

def validate(evaluate=True, threshold=1e-10):

    number_of_iterations = 50000

    # ______________________________________________________________________
    # Check output files are created

    # list of output files
    output_file_list = []

    # Create output file list
    for field in ["cloud_s00"]:
        for it in range(0,number_of_iterations,1000):

            file = "{}_{:05d}.bin".format(field, it)
            output_file_list.append(file)

    # Add scalars
    output_file_list.append("fields.txt")
    output_file_list.append("species_00.txt")

    # Check that all output files exist
    for file in output_file_list:
        if (not(os.path.exists("diags/"+file))):
            raise ValueError('File {} not generated'.format(file))

    # ______________________________________________________________________
    # Check final scalar for species

    reference_data = [50000, 1.0, 0.003123105625616999]

    with open("diags/species_00.txt", "r") as f:
        lines = f.readlines()

        last_line = lines[-1].split(" ")

        iteration = int(last_line[0])
        particles = float(last_line[1])
        energy = float(last_line[2])

    if (evaluate):

        if (reference_data[0] != iteration):
            minipic_ci.error('Last iteration in species_00.txt is not correct'.format(iteration, reference_data[0]))

        minipic_ci.evaluate(particles, reference_data[1], reference_data[1], '==', 'Number of particles in species_00.txt is not correct'.format(particles, reference_data[1]))

        minipic_ci.evaluate(energy, reference_data[2], threshold, 'relative', 'Kinetic energy in species_00.txt is not correct'.format(energy, reference_data[2])) 

    else:

        print("reference_data = [{}, {}, {}]".format(iteration, particles, energy))

    # ______________________________________________________________________
    # Check cloud files

    x_sum_ref = 27.258704826095418
    y_sum_ref = 25.45718322703324
    z_sum_ref = 25.499999999999996
    px_sum_ref = 126.92720347227889
    py_sum_ref = 127.28848117864518
    pz_sum_ref = 0.0

    nb_files = int(number_of_iterations/1000)

    x_array = np.zeros(nb_files)
    y_array = np.zeros(nb_files)
    z_array = np.zeros(nb_files)

    px_array = np.zeros(nb_files)
    py_array = np.zeros(nb_files)
    pz_array = np.zeros(nb_files)

    for i,it in enumerate(range(0, number_of_iterations, 1000)):

        file = "cloud_s00_{:05d}.bin".format(it)

        particle_number, id, w, x, y, z, px, py, pz = minipic_diag.read_particle_cloud("diags/"+file)

        x_array[i] = x[0]
        y_array[i] = y[0]
        z_array[i] = z[0]

        px_array[i] = px[0]
        py_array[i] = py[0]
        pz_array[i] = pz[0]

    x_sum = np.sum(np.abs(x_array))
    y_sum = np.sum(np.abs(y_array))
    z_sum = np.sum(np.abs(z_array))

    px_sum = np.sum(np.abs(px_array))
    py_sum = np.sum(np.abs(py_array))
    pz_sum = np.sum(np.abs(pz_array))

    if (evaluate):

        minipic_ci.evaluate(x_sum, x_sum_ref, threshold, 'relative', 'Sum over x positions not similar')
        minipic_ci.evaluate(y_sum, y_sum_ref, threshold, 'relative', 'Sum over y positions not similar')
        minipic_ci.evaluate(z_sum, z_sum_ref, threshold, 'relative', 'Sum over z positions not similar')

        minipic_ci.evaluate(px_sum, px_sum_ref, threshold, 'relative', 'Sum over px not similar')
        minipic_ci.evaluate(py_sum, py_sum_ref, threshold, 'relative', 'Sum over py not similar')
        minipic_ci.evaluate(pz_sum, pz_sum_ref, threshold, 'absolute', 'Sum over pz not similar')

    else:

        print("x_sum_ref = {}".format(x_sum))
        print("y_sum_ref = {}".format(y_sum))
        print("z_sum_ref = {}".format(z_sum))

        print("px_sum_ref = {}".format(px_sum))
        print("py_sum_ref = {}".format(py_sum))
        print("pz_sum_ref = {}".format(pz_sum))
    