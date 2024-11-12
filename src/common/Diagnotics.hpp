
/* _____________________________________________________________________ */
//! \file Diagnostics.hpp

//! \brief Diagnostics functions to output data on disks
//!
/* _____________________________________________________________________ */

#pragma once

#include <csignal>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Params.hpp"
#include "Patch.hpp"

namespace Diags {

// _____________________________________________________________________
//! \brief Initialization of the diagnostics
//! -> Create the output directory
// _____________________________________________________________________
void initialize(Params &params) {
  // Create the output directory
  std::string command = "mkdir -p " + params.output_directory;
  system(command.c_str());
}

// _____________________________________________________________________
//! \brief output a structured array of any dimension in a binary file

//! WARNING: This function is not optimized
//! This function should be put in subdomain or in a specific diag class

//! \param[in] file_name - name of the file on disk
//! \param[in] projected_parameter - name of the data parameter
//! \param[in] axis - axis names
//! \param[in] n_cells number of cells for each axis
//! \param[in] min - min value for each axis
//! \param[in] max - max value for each axis
//! \param[in] data pointer to the raw data
// _____________________________________________________________________
void output_binary_structured_grid(std::string file_name,
                                   std::string projected_parameter,
                                   std::vector<std::string> axis,
                                   std::vector<int> n_cells,
                                   std::vector<double> min,
                                   std::vector<double> max,
                                   double *__restrict__ data) {

  const unsigned int dim = axis.size();
  unsigned int data_size = 1;

  std::ofstream binary_file(file_name, std::ios::out | std::ios::binary);
  // char * data_ptr = (char *) data.data();
  binary_file.write((char *)&dim, sizeof(int));

  char projected_data_short_name[16] = "";
  snprintf(projected_data_short_name, 16, "%s", projected_parameter.c_str());
  // sprintf(projected_data_short_name, "%16c", projected_parameter.c_str());

  binary_file.write((char *)projected_data_short_name, sizeof(char) * 16);

  for (int idim = 0; idim < dim; ++idim) {

    char axis_short_name[8] = "";
    snprintf(axis_short_name, 8, "%s", axis[idim].c_str());

    binary_file.write((char *)&axis_short_name, sizeof(char) * 8);
    binary_file.write((char *)&min[idim], sizeof(double));
    binary_file.write((char *)&max[idim], sizeof(double));
    binary_file.write((char *)&n_cells[idim], sizeof(int));

    data_size *= n_cells[idim];
  }

  binary_file.write((char *)data, sizeof(double) * data_size);

  binary_file.close();
}

// _____________________________________________________________________
//! \brief output a structured array of any dimension in a binary file

//! WARNING: This function is not optimized
//! This function should be put in subdomain or in a specific diag class

//! \param[in] file_name - name of the file on disk
//! \param[in] projected_parameter - name of the data parameter
//! \param[in] axis - axis names
//! \param[in] min - min value for each axis
//! \param[in] max - max value for each axis
// _____________________________________________________________________
void output_vtk_structured_grid(std::string file_name,
                                std::string projected_parameter,
                                std::vector<std::string> axis,
                                std::vector<int> n_cells,
                                std::vector<double> min,
                                std::vector<double> delta,
                                double *__restrict__ data) {

  const unsigned int dim = axis.size();
  unsigned int data_size = 1;
  for (int idim = 0; idim < dim; ++idim) {
    data_size *= n_cells[idim];
  }

  std::ofstream binary_file(file_name, std::ios::out);

  binary_file << "# vtk DataFile Version 2.0\n"
              << projected_parameter.c_str() << "\n"
              << "ASCII\n";
  binary_file << "DATASET STRUCTURED_POINTS\n";

  // Add the dimensions
  binary_file << "DIMENSIONS ";
  for (int idim = dim - 1; idim >= 0; --idim) {
    binary_file << n_cells[idim] << " ";
  }
  binary_file << "\n";

  // Add the aspect ratio
  binary_file << "SPACING ";
  for (int idim = dim - 1; idim >= 0; --idim) {
    binary_file << delta[idim] << " ";
  }
  binary_file << "\n";

  // Origin
  binary_file << "ORIGIN ";
  for (int idim = dim - 1; idim >= 0; --idim) {
    binary_file << min[idim] << " ";
  }
  binary_file << "\n";

  // data
  binary_file << "POINT_DATA " << data_size << "\n"
              << "SCALARS " << projected_parameter.c_str() << " double 1\n"
              << "LOOKUP_TABLE default\n";

  for (unsigned int i = 0; i < data_size; ++i) {
    binary_file << data[i] << " ";
  }
  binary_file.close();
}

// _____________________________________________________________________
//! \brief Particle binning diagnostic for any dimension (1d, 2d, 3d, and more)

//! NOTE: if for a given axis, max <= min, then the min and max values are automatically computed

//! WARNING: This function is not optimized
//! This function should be put in subdomain or in a specific diag class

//! \param[in] diag_name - string used to initiate the diag name
//! \param[in] patches - vector of patches
//! \param[in] projected_parameter - property projected on the grid, can be
//! `weight`
//! \param[in] axis - axis to use for the grid. The axis vector size
//! determines the dimension of the diag. Axis can be `gamma`, `weight`, `x`,
//! `y`, `z`, `px`, `py`, `pz`
//! \param[in] n_cells - number of cells for each axis
//! \param[in] min - min value for each axis
//! \param[in] max - max value for each axis
//! \param[in] is - species
//! \param[in] it - iteration
//! \param[in] format - (optional argument) - determine the output format, can
//! be `binary` (default) or `vtk`
// _____________________________________________________________________
void particle_binning(std::string diag_name,
                      Params &params,
                      std::vector<Patch> &patches,
                      std::string projected_parameter,
                      std::vector<std::string> axis,
                      std::vector<int> n_cells,
                      std::vector<double> min_in,
                      std::vector<double> max_in,
                      int is,
                      int it,
                      std::string format = "binary",
                      bool verbose       = false) {

  if (verbose) {
    std::cout << " create particle binning: " << diag_name << std::endl;
    std::cout << " - it: " << it << std::endl;
  }

  // Get diag dimension
  const unsigned int dim = axis.size();

  // compute the total size of the diag grid
  unsigned int diag_data_size = n_cells[0];
  for (int idim = 1; idim < dim; ++idim) {
    diag_data_size *= n_cells[idim];
  }

  // Allocate diag grid
  std::vector<double> diag_data(diag_data_size, 0.);

  // Use an int code for the axis instead of string
  int axis_code[3];
  for (int idim = 0; idim < dim; ++idim) {

    if (axis[idim] == "gamma") {
      axis_code[idim] = 0;
    } else if (axis[idim] == "x") {
      axis_code[idim] = 1;
    } else if (axis[idim] == "y") {
      axis_code[idim] = 2;
    } else if (axis[idim] == "z") {
      axis_code[idim] = 3;
    } else if (axis[idim] == "px") {
      axis_code[idim] = 4;
    } else if (axis[idim] == "py") {
      axis_code[idim] = 5;
    } else if (axis[idim] == "pz") {
      axis_code[idim] = 6;
    }
  }

  // Use a code for the projected parameter
  int data_code;
  if (projected_parameter == "weight") {
    data_code = 0;
  } else if (projected_parameter == "density") {
    data_code = 1;
  } else if (projected_parameter == "particles") {
    data_code = 2;
  }

  // Check the min max values
  double min_value[3];
  double max_value[3];

  for (int idim = 0; idim < dim; ++idim) {

    min_value[idim] = min_in[idim];
    max_value[idim] = max_in[idim];

    // We compute the boundaries in this case
    // if min_in <= max_in
    if (min_value[idim] >= max_value[idim]) {

      if (axis_code[idim] == 0) {
        min_value[idim] = 1e9;
        max_value[idim] = -1e9;
      }

      for (int i_patch = 0; i_patch < patches.size(); ++i_patch) {
        const unsigned int n_particles = patches[i_patch].particles_m[is].size();
        for (unsigned int ip = 0; ip < n_particles; ip++) {

          mini_float value;

          if (axis_code[idim] == 0) {
            value = sqrt(1 +
                         patches[i_patch].particles_m[is].mx_h(ip) *
                           patches[i_patch].particles_m[is].mx_h(ip) +
                         patches[i_patch].particles_m[is].my_h(ip) *
                           patches[i_patch].particles_m[is].my_h(ip) +
                         patches[i_patch].particles_m[is].mz_h(ip) *
                           patches[i_patch].particles_m[is].mz_h(ip));
          } else if (axis_code[idim] == 1) {
            value = patches[i_patch].particles_m[is].x_h(ip);
          } else if (axis_code[idim] == 2) {
            value = patches[i_patch].particles_m[is].y_h(ip);
          } else if (axis_code[idim] == 3) {
            value = patches[i_patch].particles_m[is].z_h(ip);
          } else if (axis_code[idim] == 4) {
            value = patches[i_patch].particles_m[is].mx_h(ip);
          } else if (axis_code[idim] == 5) {
            value = patches[i_patch].particles_m[is].my_h(ip);
          } else if (axis_code[idim] == 6) {
            value = patches[i_patch].particles_m[is].mz_h(ip);
          }

          min_value[idim] = min(min_value[idim], static_cast<double>(value));
          max_value[idim] = max(max_value[idim], static_cast<double>(value));
        }
      } // loop on patches

      max_value[idim] *= 1.01;

      if (verbose) {
        std::cout << " - Auto min/max for axis " << axis[idim] << " " << min_value[idim] << " "
                  << max_value[idim] << std::endl;
      }
    }
  }

  // Delta
  std::vector<double> delta(dim);
  for (int idim = 0; idim < dim; ++idim) {
    delta[idim] = (max_value[idim] - min_value[idim]) / n_cells[idim];
  }

  // Project the particle properties of each patch in diag_data
  for (int i_patch = 0; i_patch < patches.size(); ++i_patch) {

    // get number of particles to project in the current patch
    const unsigned int n_particles = patches[i_patch].particles_m[is].size();

    // Compute data
    for (int ip = 0; ip < n_particles; ip++) {

      mini_float value[3];
      bool inside_diag_data = true;

      for (int idim = 0; idim < dim; ++idim) {

        if (axis_code[idim] == 0) {
          value[idim] = sqrt(
            1 +
            patches[i_patch].particles_m[is].mx_h(ip) * patches[i_patch].particles_m[is].mx_h(ip) +
            patches[i_patch].particles_m[is].my_h(ip) * patches[i_patch].particles_m[is].my_h(ip) +
            patches[i_patch].particles_m[is].mz_h(ip) * patches[i_patch].particles_m[is].mz_h(ip));
        } else if (axis_code[idim] == 1) {
          value[idim] = patches[i_patch].particles_m[is].x_h(ip);
        } else if (axis_code[idim] == 2) {
          value[idim] = patches[i_patch].particles_m[is].y_h(ip);
        } else if (axis_code[idim] == 3) {
          value[idim] = patches[i_patch].particles_m[is].z_h(ip);
        } else if (axis_code[idim] == 4) {
          value[idim] = patches[i_patch].particles_m[is].mx_h(ip);
        } else if (axis_code[idim] == 5) {
          value[idim] = patches[i_patch].particles_m[is].my_h(ip);
        } else if (axis_code[idim] == 6) {
          value[idim] = patches[i_patch].particles_m[is].mz_h(ip);
        }

        if (value[idim] < min_value[idim] || value[idim] >= max_value[idim]) {
          // Compute global cell index
          inside_diag_data = false;
        }
      }

      if (inside_diag_data) {

        int index[3];
        unsigned int global_index = 0;

        // For each axis, we compute the value to be projected and the grid
        // index
        for (int idim = 0; idim < dim; ++idim) {

          // Compute global cell index
          index[idim] = static_cast<int>(floor((value[idim] - min_value[idim]) / delta[idim]));

          // Compute the global 1d index in diag_data
          global_index = global_index * n_cells[idim] + index[idim];
        }

        if (data_code == 0) {
          diag_data[global_index] += static_cast<double>(patches[i_patch].particles_m[is].w_h(ip));
        } else if (data_code == 1) {
          diag_data[global_index] +=
            static_cast<double>(patches[i_patch].particles_m[is].w_h(ip) * params.inv_cell_volume);
        } else if (data_code == 2) {
          diag_data[global_index] += 1;
        }
      }
    } // end for each particles

  } // end for each patch

  // _______________________________________________________________
  // Save the diag_data

  // Output name template
  std::string template_string = "%s_s%02d_%0" + std::to_string(params.max_it_digits) + "d";

  if (format == "binary" or format == "bin") {
    template_string += ".bin";
  } else if (format == "vtk") {
    template_string += ".vtk";
  }

  char buffer[64];
  snprintf(buffer, 64, template_string.c_str(), diag_name.c_str(), is, it);

  std::string file_name(params.output_directory + "/" + buffer);

  // Binary format
  if (format == "binary" or format == "bin") {

    Diags::output_binary_structured_grid(file_name,
                                         projected_parameter,
                                         axis,
                                         n_cells,
                                         {min_value[0], min_value[1], min_value[2]},
                                         {max_value[0], max_value[1], max_value[2]},
                                         diag_data.data());

    // Old style VTK format
  } else if (format == "vtk") {

    output_vtk_structured_grid(file_name,
                               projected_parameter,
                               axis,
                               n_cells,
                               {min_value[0], min_value[1], min_value[2]},
                               {delta[0], delta[1], delta[2]},
                               diag_data.data());
  }
}

// _____________________________________________________________________
//! \brief output the field grids

//! NOTE: Create a binary file per field component

//! \param[in] Params - global constant parameters
//! \param[in] em - species
//! \param[in] it - iteration
//! \param[in] format - (optional argument) - determine the output format, can
//! be `binary` (default) or `vtk`
// _____________________________________________________________________
void fields(Params &params, ElectroMagn &em, unsigned int it, std::string format = "binary") {

  // Output name template
  std::string template_string = "%s_%0" + std::to_string(params.max_it_digits) + "d";

  if (format == "binary" or format == "bin") {
    template_string += ".bin";
  } else if (format == "vtk") {
    template_string += ".vtk";
  }

  char buffer[64];

  std::vector<Field<mini_float> *> field_list =
    {&em.Ex_m, &em.Ey_m, &em.Ez_m, &em.Bx_m, &em.By_m, &em.Bz_m};

  for (auto ifield = 0; ifield < field_list.size(); ++ifield) {

#if defined(__SHAMAN__)
    const int nx = field_list[ifield]->nx();
    const int ny = field_list[ifield]->ny();
    const int nz = field_list[ifield]->nz();

    double *tmp_grid = new double[field_list[ifield]->size()];

    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        for (int iz = 0; iz < nz; ++iz) {
          int i       = ix * ny * nz + iy * nz + iz;
          tmp_grid[i] = static_cast<double>(field_list[ifield]->operator()(ix, iy, iz));
        }
      }
    }
#else
    double *tmp_grid = field_list[ifield]->get_raw_pointer(minipic::host);
#endif

    // Binary format
    if (format == "binary" or format == "bin") {
      snprintf(buffer, 64, template_string.c_str(), field_list[ifield]->name_m.c_str(), it);
      std::string file_name(params.output_directory + "/" + buffer);

      output_binary_structured_grid(
        file_name,
        field_list[ifield]->name_m,
        {"x", "y", "z"},
        {field_list[ifield]->nx_m, field_list[ifield]->ny_m, field_list[ifield]->nz_m},
        {params.inf_x, params.inf_y, params.inf_z},
        {params.sup_x, params.sup_y, params.sup_z},
        tmp_grid);

      // Old style vtk format
    } else if (format == "vtk") {

      snprintf(buffer, 64, template_string.c_str(), field_list[ifield]->name_m.c_str(), it);
      std::string file_name(params.output_directory + "/" + buffer);

      output_vtk_structured_grid(
        file_name,
        field_list[ifield]->name_m.c_str(),
        {"x", "y", "z"},
        {field_list[ifield]->nx_m, field_list[ifield]->ny_m, field_list[ifield]->nz_m},
        {params.inf_x, params.inf_y, params.inf_z},
        {params.dx, params.dy, params.dz},
        tmp_grid);
    }

#if defined(__SHAMAN__)

    delete[] tmp_grid;

#endif
  }
}

// _____________________________________________________________________
//! \brief Particle cloud diagnostic
//!        Dump all particle properties in a binary or vtk file.

//! NOTE: This function should be put in subdomain or in a specific diag class

//! \param[in] diag_name - string used to initiate the diag name
//! \param[in] patches - vector of patches
//! \param[in] is - species
//! \param[in] it - iteration
//! \param[in] format - output format, can be "binary" (default) or "vtk". This
//! defines the file extension, respectively `.bin` or `.vtk`
// _____________________________________________________________________
void particle_cloud(std::string diag_name,
                    Params &params,
                    std::vector<Patch> &patches,
                    int is,
                    int it,
                    std::string format = "binary") {

  unsigned int number_of_particles = 0;
  for (int i_patch = 0; i_patch < patches.size(); i_patch++) {
    number_of_particles += patches[i_patch].particles_m[is].size();
  }

  // Output name template
  // binary: %s_s%02d_%0<max_it_digits>d.bin
  // vtk:    %s_s%02d_%0<max_it_digits>d.vtk
  // where %s is the diag_name, %02d is the species index
  // and %0<max_it_digits>d the iteration number
  std::string template_string = "%s_s%02d_%0" + std::to_string(params.max_it_digits) + "d";

  if (format == "binary" or format == "bin") {

    template_string += ".bin";

    char buffer[64];
    snprintf(buffer, 64, template_string.c_str(), diag_name.c_str(), is, it);

    std::string file_name(params.output_directory + "/" + buffer);

    std::ofstream binary_file(file_name, std::ios::out | std::ios::binary);

    binary_file.write((char *)&number_of_particles, sizeof(number_of_particles));

    if (!binary_file) {
      ERROR(" Error while creating the file :" << file_name)
      std::raise(SIGABRT);
    }

    for (int i_patch = 0; i_patch < patches.size(); i_patch++) {
      const unsigned int n_particles = patches[i_patch].particles_m[is].size();
      for (int ip = 0; ip < n_particles; ++ip) {

        binary_file.write((char *)(&patches[i_patch].particles_m[is].w_h(ip)), sizeof(double));

        binary_file.write((char *)(&patches[i_patch].particles_m[is].x_h(ip)), sizeof(double));
        binary_file.write((char *)(&patches[i_patch].particles_m[is].y_h(ip)), sizeof(double));
        binary_file.write((char *)(&patches[i_patch].particles_m[is].z_h(ip)), sizeof(double));

        binary_file.write((char *)(&patches[i_patch].particles_m[is].mx_h(ip)), sizeof(double));
        binary_file.write((char *)(&patches[i_patch].particles_m[is].my_h(ip)), sizeof(double));
        binary_file.write((char *)(&patches[i_patch].particles_m[is].mz_h(ip)), sizeof(double));

      } // End particle loop
    }   // End patch loop

    // Cleaning
    binary_file.close();

  } else if (format == "vtk") {

    template_string += ".vtk";

    char buffer[64];
    snprintf(buffer, 64, "%s_s%02d_%05d.vtk", diag_name.c_str(), is, it);

    std::string file_name(params.output_directory + "/" + buffer);

    std::ofstream vtk_file(file_name, std::ios::out | std::ios::trunc);

    if (!vtk_file) {
      std::cerr << " Error while creating the file :" << file_name << std::endl;
      std::raise(SIGABRT);
    }

    vtk_file << "# vtk DataFile Version 3.0" << std::endl;
    vtk_file << "vtk output" << std::endl;
    vtk_file << "ASCII" << std::endl;
    vtk_file << "DATASET POLYDATA" << std::endl;

    // Particle positions
    vtk_file << std::endl;
    vtk_file << "POINTS " << number_of_particles << " float" << std::endl;

    for (auto i_patch = 0; i_patch < patches.size(); ++i_patch) {
      const unsigned int n_particles = patches[i_patch].particles_m[is].size();
      for (auto ip = 0; ip < patches[i_patch].particles_m[is].size(); ++ip) {

        vtk_file << patches[i_patch].particles_m[is].z_h(ip) << " "
                 << patches[i_patch].particles_m[is].y_h(ip) << " "
                 << patches[i_patch].particles_m[is].x_h(ip) << std::endl;

      } // End particle loop
    }   // End patch loop

    // Construction of the weight
    vtk_file << std::endl;
    vtk_file << "POINT_DATA " << number_of_particles << std::endl;
    vtk_file << "SCALARS weight float" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;

    for (int i_patch = 0; i_patch < patches.size(); i_patch++) {
      for (int ip = 0; ip < patches[i_patch].particles_m[is].size(); ++ip) {

        vtk_file << patches[i_patch].particles_m[is].w_h(ip) << " ";

      } // End particle loop
    }   // End patch loop

    vtk_file << std::endl;

    // Construction of the energy
    vtk_file << std::endl;
    vtk_file << "SCALARS gamma float" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (int i_patch = 0; i_patch < patches.size(); i_patch++) {

      for (int ip = 0; ip < patches[i_patch].particles_m[is].size(); ++ip) {

        const mini_float gamma =
          1 /
          sqrt(
            1. +
            patches[i_patch].particles_m[is].mx_h(ip) * patches[i_patch].particles_m[is].mx_h(ip) +
            patches[i_patch].particles_m[is].my_h(ip) * patches[i_patch].particles_m[is].my_h(ip) +
            patches[i_patch].particles_m[is].mz_h(ip) * patches[i_patch].particles_m[is].mz_h(ip));

        vtk_file << gamma << " ";
      }
    }
    vtk_file << std::endl;

    // Construction of the momentum vector
    vtk_file << std::endl;
    vtk_file << "VECTORS momentum float" << std::endl;
    for (int i_patch = 0; i_patch < patches.size(); i_patch++) {
      for (int ip = 0; ip < patches[i_patch].particles_m[is].size(); ++ip) {
        vtk_file << patches[i_patch].particles_m[is].mx_h(ip) << " ";
        vtk_file << patches[i_patch].particles_m[is].my_h(ip) << " ";
        vtk_file << patches[i_patch].particles_m[is].mz_h(ip) << " ";
      }
    }

    vtk_file.close();
  }
}

// _____________________________________________________________________
//! \brief Particle scalars
//!        Compute reduced 0D scalars for particles

//! NOTE: Create an ascii file called `species_*.txt` where `*` is the species
//! number `is`.

//! NOTE: If the parameter `it` is equal to 0, the file is reinitialized
//! (ios::trunc). Else, the new data is added at the end of the file.

//! NOTE: The first line of the file is a header describing all columns.

//! \param[in] Params - global constant parameters
//! \param[in] patches - vector of patches
//! \param[in] is - species
//! \param[in] it - iteration
// _____________________________________________________________________
void scalars(Params &params, std::vector<Patch> &patches, unsigned int is, unsigned int it) {

  // Particle scalars _________________________________________

  unsigned int number_of_particles  = 0;
  mini_float species_kinetic_energy = 0;

  for (unsigned int i_patch = 0; i_patch < patches.size(); ++i_patch) {

    species_kinetic_energy += patches[i_patch].particles_m[is].get_kinetic_energy(minipic::device);

    number_of_particles += patches[i_patch].particles_m[is].size();

  } // end loop patch

  // output scalars file _____________________________________

  std::ofstream scalars_file;

  char buffer[64];
  snprintf(buffer, 64, "species_%02d.txt", is);

  std::string file_name(params.output_directory + "/" + buffer);

  if (it == 0) {
    scalars_file.open(file_name, std::ios::out | std::ios::trunc);
  } else {
    scalars_file.open(file_name, std::ios::out | std::ios::app);
  }

  if (!scalars_file.is_open()) {
    std::cerr << " Error while accessing the scalars file." << std::endl;
    std::raise(SIGABRT);
  }

  // Header
  if (it == 0) {
    scalars_file << "iter number_of_particles E_kin \n";
  }

  // Data
  scalars_file << it;
  scalars_file << " " << number_of_particles;
  scalars_file << " " << std::scientific << std::setprecision(15) << species_kinetic_energy;
  scalars_file << "\n";

  scalars_file.close();
}

// _____________________________________________________________________
//! \brief Field scalars
//!        Compute reduced 0D scalars for Fields

//! NOTE: Create an ascii file called `fields.txt`.

//! NOTE: If the parameter `it` is equal to 0, the file is reinitialized
//! (ios::trunc). Else, the new data is added at the end of the file.

//! NOTE: The first line of the file is a header describing all columns.

//! \param[in] Params - global constant parameters
//! \param[in] em - electromagnetic object containing all fields
//! \param[in] it - iteration
// _____________________________________________________________________
void scalars(Params &params, ElectroMagn &em, unsigned int it) {

  // compute Field scalars _________________________________________

  mini_float Ex_energy = 0.5 * em.Ex_m.sum(2, minipic::device) * params.cell_volume;
  mini_float Ey_energy = 0.5 * em.Ey_m.sum(2, minipic::device) * params.cell_volume;
  mini_float Ez_energy = 0.5 * em.Ez_m.sum(2, minipic::device) * params.cell_volume;

  mini_float Bx_energy = 0.5 * em.Bx_m.sum(2, minipic::device) * params.cell_volume;
  mini_float By_energy = 0.5 * em.By_m.sum(2, minipic::device) * params.cell_volume;
  mini_float Bz_energy = 0.5 * em.Bz_m.sum(2, minipic::device) * params.cell_volume;

  mini_float Jx_energy = 0.5 * em.Jx_m.sum(2, minipic::device) * params.cell_volume;
  mini_float Jy_energy = 0.5 * em.Jy_m.sum(2, minipic::device) * params.cell_volume;
  mini_float Jz_energy = 0.5 * em.Jz_m.sum(2, minipic::device) * params.cell_volume;

  // #endif

  // output Field scalars _________________________________________

  std::ofstream scalars_file;

  std::string file_name(params.output_directory + "/fields.txt");

  if (it == 0) {
    scalars_file.open(file_name, std::ios::out | std::ios::trunc);
  } else {
    scalars_file.open(file_name, std::ios::out | std::ios::app);
  }

  if (!scalars_file.is_open()) {
    std::cerr << " Error while accessing the file `fields.txt`." << std::endl;
    std::raise(SIGABRT);
  }

  // Header
  if (it == 0) {
    scalars_file << "dt: " << params.dt << "\n";
    scalars_file << "iter Ex Ey Ez Bx By Bz Jx Jy Jz\n";
  }

  // Data
  scalars_file << it;
  scalars_file << " " << std::scientific << std::setprecision(15) << Ex_energy;
  scalars_file << " " << std::scientific << std::setprecision(15) << Ey_energy;
  scalars_file << " " << std::scientific << std::setprecision(15) << Ez_energy;

  scalars_file << " " << std::scientific << std::setprecision(15) << Bx_energy;
  scalars_file << " " << std::scientific << std::setprecision(15) << By_energy;
  scalars_file << " " << std::scientific << std::setprecision(15) << Bz_energy;

  scalars_file << " " << std::scientific << std::setprecision(15) << Jx_energy;
  scalars_file << " " << std::scientific << std::setprecision(15) << Jy_energy;
  scalars_file << " " << std::scientific << std::setprecision(15) << Jz_energy;
  scalars_file << "\n";

  scalars_file.close();
}

} // namespace Diags