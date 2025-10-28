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
        std::cout << std::endl;

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
        std::cout << std::endl;
        
        // _________________________________________________________________________
        // Kernel on host

        std::cout << " > Test kernel on host" << std::endl;

        mini_float expected_sum = 0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    expected_sum += pow(ix - iy + iz,2);
                    field(ix,iy,iz) = ix - iy + iz;
                }
            }
        }

        auto field_host_sum = field.sum(2,minipic::host);
        auto error_host = std::abs(field_host_sum - expected_sum) / expected_sum;

        std::cout << std::setprecision(15) << "   - expected sum: " << expected_sum << std::endl;
        std::cout << std::setprecision(15) << "   - sum on host: " << field_host_sum << " with error: " << error_host
                << std::endl;

        assert( error_host < 1e-12);
        std::cout << std::endl;

        auto field_device_sum=0.0;

        // _________________________________________________________________________
        // test host->device transfers
        std::cout << " > Test host->device transfers " << std::endl;
        #if defined (__MINIPIC_CUDA__)
        std::cout << " > Test reset on device " << std::endl;
        field.reset(minipic::device);
        field_device_sum = field.sum(1,minipic::device);
        std::cout << std::setprecision(15) << "   - sum on device before sync: " << field_device_sum << " ( expected : 0 )" << std::endl;

        #endif

        field.sync(minipic::host, minipic::device);
        field_device_sum = field.sum(2,minipic::device);
        auto error_device = std::abs(field_device_sum - expected_sum) / expected_sum;
        std::cout << std::setprecision(15) << "   - sum on device after sync : " << field_device_sum << " with error: " << error_device << std::endl;

        assert( error_device < 1e-12);
        std::cout << std::endl;


        // _________________________________________________________________________
        // Kernel on device
       std::cout << " > Test kernel on device" << std::endl;

        // get the right sum
        expected_sum = 0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    expected_sum += pow(ix - iy + iz,2);
                }
            }
        }

        mini_float device_sum = field.sum(2,minipic::device);

        std::cout << "   - sum on device: " << device_sum << "  - expected: " << expected_sum
                << std::endl;

        assert(device_sum == expected_sum);
        std::cout << std::endl;

        // _________________________________________________________________________
        // Device to host transfer
        std::cout << " > Test device->host transfers " << std::endl;  

        #if defined (__MINIPIC_CUDA__)
        std::cout << " > Test reset on host " << std::endl;
        field.reset(minipic::host);
        field_host_sum = field.sum(1,minipic::host);
        std::cout << std::setprecision(15) << "   - sum on host before sync : " << field_host_sum << " ( expected : 0 )" << std::endl;
        #endif

        field.sync(minipic::device, minipic::host);
        mini_float host_sum = field.sum(2,minipic::host);
        std::cout << "   - sum on host: " << host_sum << "  - expected: " << expected_sum
                << std::endl;
        assert(host_sum == expected_sum);
        std::cout << std::endl;
        
        #if !defined(__MINIPIC_CUDA__)
        std::cout << " > Test reset on host and device " << std::endl;
        field.reset(minipic::host);
        field.reset(minipic::device);
        
        field_host_sum = field.sum(1,minipic::host);
        field_device_sum = field.sum(1,minipic::device);

        std::cout << std::setprecision(15) << "   - sum on host after reset : " << field_host_sum << " ( expected : 0 )" << std::endl;
        std::cout << std::setprecision(15) << "   - sum on device after reset : " << field_device_sum << " ( expected : 0 )" << std::endl;
        #endif


        // delete [] data;
        // delete &field;
        
    }

    // _________________________________________________________________________
    // Destructor
    std::cout << std::endl;
    std::cout << " > Test destructor" << std::endl;

    backend.finalize();

}
