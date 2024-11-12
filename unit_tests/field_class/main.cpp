// _________________________________________________________________________
//
// Unit test - Field class
//
//! \brief test the Field class
//
// _________________________________________________________________________

#include "Headers.hpp"
#include "Params.hpp"
#include "Backend.hpp"
#include "Field.hpp"

#include <iomanip>
#include <cassert>


int main(int argc, char *argv[]) {


    std::cout << " ___________________________________________________________ " << std::endl;
    std::cout << "|                                                           |" << std::endl;
    std::cout << "|                      Field class                          |" << std::endl;
    std::cout << "|___________________________________________________________|" << std::endl;

    Params params;

    Backend backend;

    backend.init(argc, argv, params);
    {

        // backend.info();

        int nx = 256;
        int ny = 234;
        int nz = 1024;
        int nynz=ny*nz;
        int size = nx*ny*nz;

        std::string name = "my_field";

        std::cout << " > Test constructor" << std::endl;

        Field<mini_float> field;
        
        field.allocate(nx, ny, nz, backend, 0, 0, 0, 0, name);

        // _________________________________________________________________________
        // Get size

        std::cout << " > Test size method" << std::endl;
        std::cout << "   - nx: " << field.nx() << " " << nx << std::endl;
        std::cout << "   - ny: " << field.ny() << " " << ny << std::endl;
        std::cout << "   - nz: " << field.nz() << " " << nz << std::endl;
        std::cout << "   - size: " << field.size() << " " << nx*ny*nz << std::endl;

        assert(field.size() == nx*ny*nz);

        // _________________________________________________________________________
        // Test fill method

        std::cout << " > Test fill method" << std::endl;

        field.fill(2.0, minipic::device);
        field.fill(2.0, minipic::host);

        std::cout << "   - sum on device: " << field.sum(2,minipic::device) 
                << "   - sum on host: " << field.sum(2,minipic::host)
                << std::endl;

        assert(field.sum(2,minipic::device) == nx*ny*nz*4.0);
        assert(field.sum(2,minipic::host) == nx*ny*nz*4.0);

        // _________________________________________________________________________
        // Kernel on host

        std::cout << " > Test kernel on host" << std::endl;

#if defined (__MINIPIC_KOKKOS__)

        field.reset(minipic::host);

        field_t host_view = field.data_m.h_view;

        typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>> host_mdrange_policy;
        Kokkos::parallel_for(
        host_mdrange_policy({0, 0, 0},
                        {field.nx(), field.ny(), field.nz()}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
            host_view(ix, iy, iz) = ix - iy + iz;
        });

#elif defined (__MINIPIC_KOKKOS_UNIFIED__)

        field.reset(minipic::host);

        field_t host_view = field.data_m;

        typedef Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<3>> host_mdrange_policy;

        Kokkos::parallel_for(
        host_mdrange_policy({0, 0, 0},
                        {field.nx(), field.ny(), field.nz()}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
            host_view(ix, iy, iz) = ix - iy + iz;
        });

#elif defined(__MINIPIC_SYCL__)

        field.reset(minipic::host);

        double * host_data = field.host_data_;

        for (auto ix = 0; ix < nx; ++ix) {
            for (auto iy = 0; iy < ny; ++iy) {
                for (auto iz = 0; iz < nz; ++iz) {
                    auto index = ix*ny*nz + iy*nz + iz;
                    host_data[index] = ix - iy + iz;
                }
            }
        }

        // backend.sycl_queue_->submit([&](sycl::handler& cgh) {
        // cgh.parallel_for(sycl::range<3>(nx,ny,nz), [=](sycl::id<3> idx) {
        //     int ix = idx[0];
        //     int iy = idx[1];
        //     int iz = idx[2];
        //     int index = ix*ny*nz + iy*nz + iz;
        //     device_data[index] = index*index;
        // });
        // });
        // backend.sycl_queue_->wait();

#elif defined(__MINIPIC_OPENACC__)

        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    field(ix,iy,iz) = ix - iy + iz;
                }
            }
        }

#elif defined(__MINIPIC_STDPAR__)

        mini_float *const ptrv1 = field.get_raw_pointer(minipic::host);
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    auto index = ix*ny*nz + iy*nz + iz;
                    ptrv1[index] = ix - iy + iz;
                }
            }
        }

#endif

        mini_float expected_sum = 0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    expected_sum += pow(ix - iy + iz,2);
                }
            }
        }

        auto field_host_sum = field.sum(2,minipic::host);
        auto error_host = std::abs(field_host_sum - expected_sum) / expected_sum;

        std::cout << std::setprecision(15) << "   - expected sum: " << expected_sum << std::endl;
        std::cout << std::setprecision(15) << "   - sum on host: " << field_host_sum << " with error: " << error_host
                << std::endl;

        assert( error_host < 1e-12);

        // _________________________________________________________________________
        // test host->device transfers

        std::cout << " > Test host->device transfers" << std::endl;

        field.sync(minipic::host, minipic::device);

        auto field_device_sum = field.sum(2,minipic::device);

        auto error_device = std::abs(field_device_sum - expected_sum) / expected_sum;

        std::cout << std::setprecision(15) << "   - sum on device: " << field_device_sum << " with error: " << error_device << std::endl;

        assert( error_device < 1e-12);

        // _________________________________________________________________________
        // Kernel on device

        std::cout << " > Test kernel on device" << std::endl;

#if defined (__MINIPIC_KOKKOS__)

        field.reset(minipic::device);

        device_field_t view = field.data_m.d_view;

        typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
        Kokkos::parallel_for(
        mdrange_policy({0, 0, 0},
                        {field.nx(), field.ny(), field.nz()}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
            view(ix, iy, iz) = ix + iy + iz;
        });

#elif defined (__MINIPIC_KOKKOS_UNIFIED__)

        field.reset(minipic::device);

        device_field_t view = field.data_m;

        typedef Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<3>> mdrange_policy;
        Kokkos::parallel_for(
        mdrange_policy({0, 0, 0},
                        {field.nx(), field.ny(), field.nz()}),
        KOKKOS_LAMBDA(const int ix, const int iy, const int iz) {
            view(ix, iy, iz) = ix + iy + iz;
        });

#elif defined(__MINIPIC_SYCL__)

        field.reset(minipic::device);

        double * device_data = field.device_data_;

        backend.sycl_queue_->submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::range<3>(nx,ny,nz), [=](sycl::id<3> idx) {
            int ix = idx[0];
            int iy = idx[1];
            int iz = idx[2];
            int index = ix*ny*nz + iy*nz + iz;
            device_data[index] = ix + iy + iz;
        });
        });
        backend.sycl_queue_->wait();

#elif defined(__MINIPIC_OPENACC__)

        #pragma acc parallel present(field)
        #pragma acc loop gang worker vector collapse(3)
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    field(ix,iy,iz) = ix + iy + iz;
                }
            }
        }
#elif defined(__MINIPIC_STDPAR__)

        mini_float *const ptrv1 = field.get_raw_pointer(minipic::device);
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), size, [=](int i) 
        {
            int ix       = i/ nynz;
            int iy       = (i - ix * nynz) / nz;
            const int iz = i - ix * nynz - iy * nz;
            ptrv1[i]=ix+iy+iz; });


#endif

        // get the right sum
        expected_sum = 0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    expected_sum += pow(ix + iy + iz,2);
                }
            }
        }

        mini_float device_sum = field.sum(2,minipic::device);

        std::cout << "   - sum on device: " << device_sum << "  - expected: " << expected_sum
                << std::endl;

        assert(device_sum == expected_sum);

        // _________________________________________________________________________
        // Device to host transfer

        std::cout << " > Test device->host transfers" << std::endl;        

#if defined(__MINIPIC_STDPAR__) || defined(__MINIPIC_KOKKOS_UNIFIED__)

#else
        field.reset(minipic::host);
#endif
        field.sync(minipic::device, minipic::host);

        mini_float host_sum = field.sum(2,minipic::host);

        std::cout << "   - sum on host: " << host_sum << "  - expected: " << expected_sum
                << std::endl;

        assert(host_sum == expected_sum);

        // delete [] data;
        // delete &field;
        
    }

    // _________________________________________________________________________
    // Destructor

    std::cout << " > Test destructor" << std::endl;

    backend.finalize();

}
