/* _____________________________________________________________________ */
//! \file Shaman_analysis.hpp

//! \brief Analysis of the simulation numerical precision with Shaman

/* _____________________________________________________________________ */

#include "Backend.hpp"
#include "SubDomain.hpp"

#ifndef SHAMAN_ANALYSIS_H
#define SHAMAN_ANALYSIS_H

auto shaman_summary(SubDomain &subdomain) -> void {

  // ____________________________________________________
  // For the full domain, get the particle statistics

  // x, y, z, mx, my, mz, gamma_inv
  double error[7][3];
  double rel_error[7][3];
  double digits[7][3];
  for (int i = 0; i < 7; i++) {
    error[i][0] = 0;
    error[i][1] = 0;
    error[i][2] = 0;

    rel_error[i][0] = 0;
    rel_error[i][1] = 0;
    rel_error[i][2] = 0;

    digits[i][0] = 15;
    digits[i][1] = 0;
    digits[i][2] = 0;
  }

  for (int idx_patch = 0; idx_patch < subdomain.patches_.size(); idx_patch++) {
    // for each species
    for (int idx_species = 0; idx_species < subdomain.patches_[idx_patch].n_species_m;
         idx_species++) {

      unsigned int n_particles = subdomain.patches_[idx_patch].particles_m[idx_species].size();

      if (n_particles > 0) {

        // For each particle in this species
        for (int ip = 1; ip < n_particles; ip++) {

          // Get the particle error
          const double particle_error[7] = {
            (subdomain.patches_[idx_patch].particles_m[idx_species].x_[ip].error),
            (subdomain.patches_[idx_patch].particles_m[idx_species].y_[ip].error),
            (subdomain.patches_[idx_patch].particles_m[idx_species].z_[ip].error),
            (subdomain.patches_[idx_patch].particles_m[idx_species].mx_[ip].error),
            (subdomain.patches_[idx_patch].particles_m[idx_species].my_[ip].error),
            (subdomain.patches_[idx_patch].particles_m[idx_species].mz_[ip].error),
            (subdomain.patches_[idx_patch].particles_m[idx_species].gamma_inv_[ip].error)};

          // Get the relative error
          const double particle_rel_error[7] = {
            (subdomain.patches_[idx_patch].particles_m[idx_species].x_[ip].error /
             static_cast<double>(subdomain.patches_[idx_patch].particles_m[idx_species].x_[ip])),
            (subdomain.patches_[idx_patch].particles_m[idx_species].y_[ip].error /
             static_cast<double>(subdomain.patches_[idx_patch].particles_m[idx_species].y_[ip])),
            (subdomain.patches_[idx_patch].particles_m[idx_species].z_[ip].error /
             static_cast<double>(subdomain.patches_[idx_patch].particles_m[idx_species].z_[ip])),
            (subdomain.patches_[idx_patch].particles_m[idx_species].mx_[ip].error /
             static_cast<double>(subdomain.patches_[idx_patch].particles_m[idx_species].mx_[ip])),
            (subdomain.patches_[idx_patch].particles_m[idx_species].my_[ip].error /
             static_cast<double>(subdomain.patches_[idx_patch].particles_m[idx_species].my_[ip])),
            (subdomain.patches_[idx_patch].particles_m[idx_species].mz_[ip].error /
             static_cast<double>(subdomain.patches_[idx_patch].particles_m[idx_species].mz_[ip])),
            (subdomain.patches_[idx_patch].particles_m[idx_species].gamma_inv_[ip].error /
             static_cast<double>(
               subdomain.patches_[idx_patch].particles_m[idx_species].gamma_inv_[ip]))};

          // Get the particle digits
          const double particle_digits[7] = {
            (subdomain.patches_[idx_patch].particles_m[idx_species].x_[ip].digits()),
            (subdomain.patches_[idx_patch].particles_m[idx_species].y_[ip].digits()),
            (subdomain.patches_[idx_patch].particles_m[idx_species].z_[ip].digits()),
            (subdomain.patches_[idx_patch].particles_m[idx_species].mx_[ip].digits()),
            (subdomain.patches_[idx_patch].particles_m[idx_species].my_[ip].digits()),
            (subdomain.patches_[idx_patch].particles_m[idx_species].mz_[ip].digits()),
            (subdomain.patches_[idx_patch].particles_m[idx_species].gamma_inv_[ip].digits())};

          // For each component
          for (int i = 0; i < 7; i++) {

            // Update the min max and mean error
            error[i][0] = std::min(error[i][0], particle_error[i]);
            error[i][1] = std::max(error[i][1], particle_error[i]);
            error[i][2] += abs(particle_error[i]);

            rel_error[i][0] = std::min(rel_error[i][0], particle_rel_error[i]);
            rel_error[i][1] = std::max(rel_error[i][1], particle_rel_error[i]);
            rel_error[i][2] += abs(particle_rel_error[i]);

            // Update the min max and mean digits
            digits[i][0] = std::min(digits[i][0], particle_digits[i]);
            digits[i][1] = std::max(digits[i][1], particle_digits[i]);
            digits[i][2] += abs(particle_digits[i]);
          }
        }
      }
    }
  }

  // Compute the mean error
  for (int i = 0; i < 7; i++) {
    error[i][2] /= subdomain.get_total_number_of_particles();
    rel_error[i][2] /= subdomain.get_total_number_of_particles();
    digits[i][2] /= subdomain.get_total_number_of_particles();
  }

  // Print the error

  std::cout << endl;
  std::cout << (" ------------------------------------ \n");
  std::cout << (" SHAMAN\n");
  std::cout << (" ------------------------------------ \n");

  std::cout << endl;
  std::cout << " Particles absolute error: " << std::endl;
  std::cout << " Error on x: " << error[0][0] << " " << error[0][1] << " " << error[0][2]
            << std::endl;
  std::cout << " Error on y: " << error[1][0] << " " << error[1][1] << " " << error[1][2]
            << std::endl;
  std::cout << " Error on z: " << error[2][0] << " " << error[2][1] << " " << error[2][2]
            << std::endl;
  std::cout << " Error on mx: " << error[3][0] << " " << error[3][1] << " " << error[3][2]
            << std::endl;
  std::cout << " Error on my: " << error[4][0] << " " << error[4][1] << " " << error[4][2]
            << std::endl;
  std::cout << " Error on mz: " << error[5][0] << " " << error[5][1] << " " << error[5][2]
            << std::endl;
  std::cout << " Error on gamma_inv: " << error[6][0] << " " << error[6][1] << " " << error[6][2]
            << std::endl;

  std::cout << endl;
  std::cout << " Particles relative error: " << std::endl;
  std::cout << " Error on x: " << rel_error[0][0] << " " << rel_error[0][1] << " "
            << rel_error[0][2] << std::endl;
  std::cout << " Error on y: " << rel_error[1][0] << " " << rel_error[1][1] << " "
            << rel_error[1][2] << std::endl;
  std::cout << " Error on z: " << rel_error[2][0] << " " << rel_error[2][1] << " "
            << rel_error[2][2] << std::endl;
  std::cout << " Error on mx: " << rel_error[3][0] << " " << rel_error[3][1] << " "
            << rel_error[3][2] << std::endl;
  std::cout << " Error on my: " << rel_error[4][0] << " " << rel_error[4][1] << " "
            << rel_error[4][2] << std::endl;
  std::cout << " Error on mz: " << rel_error[5][0] << " " << rel_error[5][1] << " "
            << rel_error[5][2] << std::endl;
  std::cout << " Error on gamma_inv: " << rel_error[6][0] << " " << rel_error[6][1] << " "
            << rel_error[6][2] << std::endl;

  std::cout << endl;
  std::cout << " Particles digits: " << std::endl;
  std::cout << " Digits on x: " << setprecision(5) << digits[0][0] << " " << digits[0][1] << " "
            << digits[0][2] << std::endl;
  std::cout << " Digits on y: " << setprecision(5) << digits[1][0] << " " << digits[1][1] << " "
            << digits[1][2] << std::endl;
  std::cout << " Digits on z: " << setprecision(5) << digits[2][0] << " " << digits[2][1] << " "
            << digits[2][2] << std::endl;
  std::cout << " Digits on mx: " << setprecision(5) << digits[3][0] << " " << digits[3][1] << " "
            << digits[3][2] << std::endl;
  std::cout << " Digits on my: " << setprecision(5) << digits[4][0] << " " << digits[4][1] << " "
            << digits[4][2] << std::endl;
  std::cout << " Digits on mz: " << setprecision(5) << digits[5][0] << " " << digits[5][1] << " "
            << digits[5][2] << std::endl;
  std::cout << " Digits on gamma_inv: " << setprecision(5) << digits[6][0] << " " << digits[6][1]
            << " " << digits[6][2] << std::endl;

  // ____________________________________________________
  // For the full domain, get the field statistics

  // Field reference list Ex, Ey, Ez, Bx, By, Bz
  std::vector<Field<mini_float> *> field_list = {&subdomain.em_.Ex_m,
                                                 &subdomain.em_.Ey_m,
                                                 &subdomain.em_.Ez_m,
                                                 &subdomain.em_.Bx_m,
                                                 &subdomain.em_.By_m,
                                                 &subdomain.em_.Bz_m};

  // Field error
  double field_error[6][3];
  double field_rel_error[6][3];
  double field_digits[6][3];

  for (int i = 0; i < 6; i++) {
    field_error[i][0] = 0;
    field_error[i][1] = 0;
    field_error[i][2] = 0;

    field_rel_error[i][0] = 0;
    field_rel_error[i][1] = 0;
    field_rel_error[i][2] = 0;

    field_digits[i][0] = 15;
    field_digits[i][1] = 0;
    field_digits[i][2] = 0;
  }

  // for each field component
  for (int ifield = 0; ifield < 6; ifield++) {
    for (unsigned int ix = 0; ix < field_list[ifield]->nx(); ix++) {
      for (unsigned int iy = 0; iy < field_list[ifield]->ny(); iy++) {
        for (unsigned int iz = 0; iz < field_list[ifield]->nz(); iz++) {

          // Get the field error
          const double local_field_error = (field_list[ifield]->operator()(ix, iy, iz).error);

          // Get the relative error
          const double local_field_rel_error =
            (field_list[ifield]->operator()(ix, iy, iz).error /
             static_cast<double>(field_list[ifield]->operator()(ix, iy, iz)));

          // Get the field digits
          const double local_field_digits = (field_list[ifield]->operator()(ix, iy, iz).digits());

          // Update the min max and mean error
          field_error[ifield][0] = std::min(field_error[ifield][0], local_field_error);
          field_error[ifield][1] = std::max(field_error[ifield][1], local_field_error);
          field_error[ifield][2] += abs(local_field_error);

          field_rel_error[ifield][0] = std::min(field_rel_error[ifield][0], local_field_rel_error);
          field_rel_error[ifield][1] = std::max(field_rel_error[ifield][1], local_field_rel_error);
          field_rel_error[ifield][2] += abs(local_field_rel_error);

          // Update the min max and mean digits
          field_digits[ifield][0] = std::min(field_digits[ifield][0], local_field_digits);
          field_digits[ifield][1] = std::max(field_digits[ifield][1], local_field_digits);
          field_digits[ifield][2] += abs(local_field_digits);
        }
      }
    }

    // Compute the mean error

    unsigned int n_cells =
      field_list[ifield]->nx() * field_list[ifield]->ny() * field_list[ifield]->nz();

    field_error[ifield][2] /= n_cells;
    field_rel_error[ifield][2] /= n_cells;
    field_digits[ifield][2] /= n_cells;

  } // end for each field component

  // print the error

  std::cout << endl;
  std::cout << " Fields analysis: " << std::endl << std::endl;

  std::cout << " Absolute error: " << std::endl;
  std::cout << (" Component |     min      |     max     |     mean    |\n");
  std::cout << (" ----------|--------------|-------------|------------ |\n");

  // print the min, max, mean respecting the array order adding "|" between columns
  for (int i = 0; i < 6; i++) {
    std::cout << " " << std::setw(9) << field_list[i]->name_m << " | " << std::setw(6)
              << field_error[i][0] << " | " << std::setw(6) << field_error[i][1] << " | "
              << std::setw(6) << field_error[i][2] << " |" << std::endl;
  }

  std::cout << " Relative error: " << std::endl;
  std::cout << (" Component |     min      |     max     |     mean    |\n");
  std::cout << (" ----------|--------------|-------------|------------ |\n");

  for (int i = 0; i < 6; i++) {
    std::cout << " " << std::setw(9) << field_list[i]->name_m << " | " << std::setw(6)
              << field_rel_error[i][0] << " | " << std::setw(6) << field_rel_error[i][1] << " | "
              << std::setw(6) << field_rel_error[i][2] << " |" << std::endl;
  }

  std::cout << " digits: " << std::endl;
  std::cout << (" Component |     min      |     max     |     mean    |\n");
  std::cout << (" ----------|--------------|-------------|------------ |\n");

  for (int i = 0; i < 6; i++) {
    std::cout << " " << std::setw(9) << field_list[i]->name_m << " | " << std::setw(6)
              << field_digits[i][0] << " | " << std::setw(6) << field_digits[i][1] << " | "
              << std::setw(6) << field_digits[i][2] << " |" << std::endl;
  }
}

#endif