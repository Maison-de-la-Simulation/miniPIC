# ______________________________________________________________________
#
# Functions to read minipic diags
#
# ______________________________________________________________________

import numpy as np
import struct


# ________________________________
#
# Get file dimension
# ________________________________
def get_diag_dimension(path):

    with open(path, 'rb') as file:
        content = file.read()
        k = 0
        # Get the diag dimension
        dim = struct.unpack('i',content[k:k+4])[0] ; k+= 4

    return dim


# ________________________________
#
# Read file 1D
# ________________________________
def read_1d_diag(path):

    file = open(path, 'rb')

    content = file.read()

    k = 0

    # Get the diag dimension
    dim = struct.unpack('i',content[k:k+4])[0] ; k+= 4

    if (dim != 1):
        raise Exception('diag dimension should be 1 not {}.'.format(dim))

    data_name = np.array(struct.unpack('16s',content[k:k+16]))[0].decode("utf-8").strip() ; k+=16

    x_axis_name = np.array(struct.unpack('8s',content[k:k+8]))[0].decode("utf-8").strip() ; k+= 8
    x_min = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    x_max = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    x_ncells = struct.unpack('i',content[k:k+4])[0] ; k+= 4
    dx = (x_max-x_min)/x_ncells
    x_data = np.linspace(x_min+0.5*dx, x_max-0.5*dx, x_ncells)

    data = np.array(struct.unpack('{}d'.format(x_ncells),content[k:k + 8*x_ncells])) ; k+= 8*x_ncells

    return x_axis_name, x_min, x_max, x_data, data_name, data

# Read 2D structured grid 

# ________________________________
#
# Read 3D structured grid
# ________________________________

def read_3d_diag(path):

    # Read file 3D
    file = open(path, 'rb')

    content = file.read()

    k = 0

    # Get the diag dimension
    dim = struct.unpack('i',content[k:k+4])[0] ; k+= 4

    if (dim != 3):
        raise Exception('diag dimension should be 3 not {}.'.format(dim))

    data_name = np.array(struct.unpack('16s',content[k:k+16]))[0].decode("utf-8").strip() ; k+=16

    x_axis_name = np.array(struct.unpack('8s',content[k:k+8]))[0].decode("utf-8").strip() ; k+= 8
    x_min = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    x_max = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    x_ncells = struct.unpack('i',content[k:k+4])[0] ; k+= 4

    y_axis_name = np.array(struct.unpack('8s',content[k:k+8]))[0].decode("utf-8").strip() ; k+= 8
    y_min = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    y_max = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    y_ncells = struct.unpack('i',content[k:k+4])[0] ; k+= 4

    z_axis_name = np.array(struct.unpack('8s',content[k:k+8]))[0].decode("utf-8").strip() ; k+= 8
    z_min = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    z_max = struct.unpack('d',content[k:k+8])[0] ; k+= 8
    z_ncells = struct.unpack('i',content[k:k+4])[0] ; k+= 4

    dx = (x_max-x_min)/x_ncells
    dy = (y_max-y_min)/y_ncells
    dz = (z_max-z_min)/z_ncells
    x_data = np.linspace(x_min+0.5*dx, x_max-0.5*dx, x_ncells)
    y_data = np.linspace(y_min+0.5*dy, y_max-0.5*dy, y_ncells)
    z_data = np.linspace(z_min+0.5*dz, z_max-0.5*dz, z_ncells)

    size = x_ncells * y_ncells * z_ncells
    raw_data = np.array(struct.unpack('{}d'.format(size),content[k:k + 8*size])) ; k+= 8*size

    data_map = np.reshape(raw_data, (x_ncells, y_ncells, z_ncells))

    return x_axis_name, x_data, y_axis_name, y_data, z_axis_name, z_data, data_name, data_map

# ________________________________
#
# Read particle cloud
# ________________________________
def read_particle_cloud(path):

    file = open(path, 'rb')

    content = file.read()

    k = 0

    particle_number = struct.unpack('I',content[k:k+4])[0] ; k+= 4

    id = np.zeros(particle_number)
    w = np.zeros(particle_number)
    x = np.zeros(particle_number)
    y = np.zeros(particle_number)
    z = np.zeros(particle_number)
    px = np.zeros(particle_number)
    py = np.zeros(particle_number)
    pz = np.zeros(particle_number)

    for ip in range(particle_number):

        id[ip] = ip
        w[ip] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
        x[ip] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
        y[ip] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
        z[ip] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
        px[ip] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
        py[ip] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8
        pz[ip] = struct.unpack('d',content[k:k + 8])[0] ; k+= 8

    return particle_number, id, w, x, y, z, px, py, pz

# ________________________________
#
# Read timers
#
# Format: json
# ________________________________

def read_timers(path):
    
    import json

    # check if the file exists
    if not os.path.isfile(path):
        raise Exception('File {} does not exist'.format(path))

    with open(path) as f:
        data = json.load(f)

    return data

