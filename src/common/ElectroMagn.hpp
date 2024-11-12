/* _____________________________________________________________________ */
//! \file ElectroMagn.hpp

//! \brief Structure to store global current and EM fields

/* _____________________________________________________________________ */

#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>

#include "Field.hpp"
#include "Params.hpp"

//! Structure to store global current and EM fields
class ElectroMagn {
public:
  ElectroMagn() {}

  //! n nodes on global primal grid (N Cells+1)
  int nx_p_m, ny_p_m, nz_p_m;
  //! n nodes on global dual grid (N Primal+1 | N Cells+2)
  int nx_d_m, ny_d_m, nz_d_m;
  //! Delta in each dimenion | cell size
  mini_float dx_m, dy_m, dz_m, cell_volume_m;
  //! Invert certain value to avoid div
  mini_float inv_dx_m, inv_dy_m, inv_dz_m, inv_cell_volume_m;

  /*
    FIELDS on Yee lattice | Staggered grid
    Ex -> (ix+1/2)*dx, iy*dy      , iz*dz       same for Jx
    Ey -> ix*dx      , (iy+1/2)*dy, iz*dz       same for Jy
    Ez -> ix*dx      , iy*dy      , (iz+1/2)*dz same for Jz

    Bx -> ix*dx      , (iy+1/2)*dy, (iz+1/2)*dz
    By -> (ix+1/2)*dx, iy*dy      , (iz+1/2)*dz
    Bz -> (ix+1/2)*dx, (iy+1/2)*dy, iz*dz
  */

  //! Global electric field
  Field<mini_float> Ex_m, Ey_m, Ez_m;

  //! Global courant total field
  Field<mini_float> Jx_m, Jy_m, Jz_m;

  //! Global magnetic field
  Field<mini_float> Bx_m, By_m, Bz_m;

  // ____________________________________________________________________________
  //
  //! \brief Memory allocation for global fields and init parameters
  //! \param params global parameters
  //! \param backend backend properties
  // ____________________________________________________________________________
  void allocate(const Params &params, Backend &backend) {
    nx_p_m = params.nx_p;
    ny_p_m = params.ny_p;
    nz_p_m = params.nz_p;
    nx_d_m = params.nx_d;
    ny_d_m = params.ny_d;
    nz_d_m = params.nz_d;

    dx_m          = params.dx;
    dy_m          = params.dy;
    dz_m          = params.dz;
    cell_volume_m = params.cell_volume;

    inv_dx_m          = params.inv_dx;
    inv_dy_m          = params.inv_dy;
    inv_dz_m          = params.inv_dz;
    inv_cell_volume_m = params.inv_cell_volume;

    Jx_m.allocate(nx_d_m, ny_p_m + 2, nz_p_m + 2, backend, 0.0, 1, 0, 0, "Jx");
    Jy_m.allocate(nx_p_m + 2, ny_d_m, nz_p_m + 2, backend, 0.0, 0, 1, 0, "Jy");
    Jz_m.allocate(nx_p_m + 2, ny_p_m + 2, nz_d_m, backend, 0.0, 0, 0, 1, "Jz");

    // Jx_m.allocate(nx_d_m, ny_p_m, nz_p_m, 0.0, 1, 0, 0, "Jx");
    // Jy_m.allocate(nx_p_m, ny_d_m, nz_p_m, 0.0, 0, 1, 0, "Jy");
    // Jz_m.allocate(nx_p_m, ny_p_m, nz_d_m, 0.0, 0, 0, 1, "Jz");

    Ex_m.allocate(nx_d_m, ny_p_m, nz_p_m, backend, params.E0_[0], 1, 0, 0, "Ex");
    Ey_m.allocate(nx_p_m, ny_d_m, nz_p_m, backend, params.E0_[1], 0, 1, 0, "Ey");
    Ez_m.allocate(nx_p_m, ny_p_m, nz_d_m, backend, params.E0_[2], 0, 0, 1, "Ez");

    Bx_m.allocate(nx_p_m, ny_d_m, nz_d_m, backend, params.B0_[0], 0, 1, 1, "Bx");
    By_m.allocate(nx_d_m, ny_p_m, nz_d_m, backend, params.B0_[1], 1, 0, 1, "By");
    Bz_m.allocate(nx_d_m, ny_d_m, nz_p_m, backend, params.B0_[2], 1, 1, 0, "Bz");

    // Load all field to the device
    sync(minipic::host, minipic::device);
  }

  // ____________________________________________________________________________
  //
  //! Reset all currents grid
  // ____________________________________________________________________________
  template <class T_space> void reset_currents(const T_space space) {

    // Jx_m.reset();
    // Jy_m.reset();
    // Jz_m.reset();

    Jx_m.reset(space);
    Jy_m.reset(space);
    Jz_m.reset(space);
  }

  // ____________________________________________________________________________
  //! \brief sync host <-> device
  // ____________________________________________________________________________
  template <class T_from, class T_to> void sync(const T_from from, const T_to to) {

    Ex_m.sync(from, to);
    Ey_m.sync(from, to);
    Ez_m.sync(from, to);

    Bx_m.sync(from, to);
    By_m.sync(from, to);
    Bz_m.sync(from, to);

    Jx_m.sync(from, to);
    Jy_m.sync(from, to);
    Jz_m.sync(from, to);
  }

  // ____________________________________________________________________________
  //! \brief Print a slice of the global current grid, debug function
  // ____________________________________________________________________________
  void print_current_slice(int slice_idx) {
    std::cout << "Jx:\n";
    for (unsigned int ix = 0; ix < nx_d_m; ix++) {
      for (unsigned int iy = 0; iy < ny_p_m; iy++) {
        // printf("\t%3.4e", Jx_m(ix, iy, slice_idx));
        std::cout << Jx_m(ix, iy, slice_idx);
      }
      std::cout << ("\n");
    }
    std::cout << ("Jy:\n");
    for (unsigned int ix = 0; ix < nx_p_m; ix++) {
      for (unsigned int iy = 0; iy < ny_d_m; iy++) {
        // printf ("\t%3.4e", Jy_m(ix, iy, slice_idx));
        std::cout << Jy_m(ix, iy, slice_idx);
      }
      std::cout << ("\n");
    }
    std::cout << ("Jz:\n");
    for (unsigned int ix = 0; ix < nx_p_m; ix++) {
      for (unsigned int iy = 0; iy < ny_p_m; iy++) {
        // printf ("\t%3.4e", Jz_m(ix, iy, slice_idx));
        std::cout << Jz_m(ix, iy, slice_idx);
      }
      std::cout << ("\n");
    }
  }
};

