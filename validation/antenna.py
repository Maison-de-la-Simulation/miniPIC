# ______________________________________________________________________
#
# Validation script for the `Bcst` benchmark
#
# ______________________________________________________________________

import os
import numpy as np
import lib.minipic_ci as minipic_ci
import lib.minipic_diag as minipic_diag

def validate(evaluate=True, threshold = 1e-10):

    # ______________________________________________________________________
    # Check output files are created

    # list of output files
    output_file_list = []

    # Create output file list
    for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
        for it in range(0,700,100):

            file = "{}_{:03d}.bin".format(field, it)
            output_file_list.append(file)

    # Add scalars
    output_file_list.append("fields.txt")

    # Check that all output files exist
    for file in output_file_list:
        if (not(os.path.exists("diags/"+file))):
            raise ValueError('File {} not generated'.format(file))
        
    # ______________________________________________________________________
    # Check scalars

    print(" > Check scalars")

    # Check final scalar values for fields

    reference_data = [600, 8.392314275494247e-06, 7.685688704219648e-08, 9.162999693074787e-04, 8.581201889410053e-06, 9.523643303672555e-04, 1.993340067202490e-34, 0.000000000000000e+00, 0.000000000000000e+00, 1.877328812863773e-37]

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

        minipic_ci.evaluate(iteration, reference_data[0], threshold, 'relative', 'Last iteration in fields.txt is not correct')

        minipic_ci.evaluate(Ex, reference_data[1], threshold, 'relative', 'Ex value in fields.txt is not correct')   
        minipic_ci.evaluate(Ey, reference_data[2], threshold, 'relative', 'Ey value in fields.txt is not correct')  
        minipic_ci.evaluate(Ez, reference_data[3], threshold, 'relative', 'Ez value in fields.txt is not correct')  
        minipic_ci.evaluate(Bx, reference_data[4], threshold, 'relative', 'Bx value in fields.txt is not correct')  
        minipic_ci.evaluate(By, reference_data[5], threshold, 'relative', 'By value in fields.txt is not correct')  
        # minipic_ci.evaluate(Bz, reference_data[6], 1e-2, 'relative', 'Bz value in fields.txt is not correct') ===> too small to be relevant

    # ______________________________________________________________________
    # Check field 

    print(" > Check fields")

    reference_data = {
        'Ex': [0.0, 82.16528908415819, 344.34858930166956, 590.8495077993598, 557.1785762453475, 835.1262141764287],
        'Ey': [0.0, 9.971267188545653, 38.4454341611614, 65.35670493378103, 59.67453062284768, 76.50264530492844],
        'Ez' : [0.0, 774.7241409632107, 3237.01821064664, 5666.167748229646, 6216.341002673892, 9765.210440891475],
        'Bx' : [0.0, 86.91713443298572, 355.1215926140225, 614.3956604220939, 634.2442079543871, 881.1441186493583],
        'By' : [0.0, 784.3391276892934, 3246.3762073832377, 5627.157839639507, 6428.848152848923, 9793.622315394427],
        # 'Bz' : [0.0, 8.68335830083499e-14, 8.109979554106818e-13, 2.2881775959027207e-12, 3.4172712318143877e-12, 4.161059825829489e-12]
    }

    for field  in reference_data.keys():

        new_data = []

        for i,it in enumerate(range(0,600,100)):

            file = "{}_{:03d}.bin".format(field,it)

            x_axis_name, x_data, y_axis_name, y_data, z_axis_name, z_data, data_name, data_map = minipic_diag.read_3d_diag("diags/" + file)

            #print(len(x_data), len(y_data), len(z_data), np.shape(data_map))

            xv, yv, zv = np.meshgrid(x_data, y_data, z_data, indexing='ij')

            new_data.append(np.sum(np.abs(data_map) * xv * yv * zv))

        print("    - For field {}: {}".format(field, new_data))

        if (evaluate):
            for i,it in enumerate(range(100,600,100)):
                minipic_ci.evaluate(new_data[i+1], reference_data[field][i+1], threshold, 'relative', '{} field not similar'.format(field))