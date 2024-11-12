// _________________________________________________________________________
//
// Unit test - Particles class
//
//! \brief test the particles class
//
// _________________________________________________________________________

#include "Headers.hpp"
#include "Particles.hpp"
#include "Vector.hpp"

#include <cassert>
#include <cmath>


int main(int argc, char *argv[]) {

    std::cout << " ___________________________________________________________ " << std::endl;
    std::cout << "|                                                           |" << std::endl;
    std::cout << "|                      Particles class                         |" << std::endl;
    std::cout << "|___________________________________________________________|" << std::endl;

    Params params;

    Backend backend;

    backend.init(argc, argv, params);
    {

        // _________________________________________________________________________
        // Parameters
        int nbparticle = 20;

        // _________________________________________________________________________
        // Test constructor

        std::cout << " > Test constructor" << std::endl;
        Particles<mini_float> part;

        //Test allocate 
        part.allocate(1, 1, 1, nbparticle, 10, backend);

        // _________________________________________________________________________
        // Test size method

        std::cout << " > Test size method" << std::endl;
        std::cout << "   - size: " << part.size() << " (expected: " << nbparticle << ")" << std::endl;
        assert(part.size() == nbparticle);

        // _________________________________________________________________________
        // Test size method

        std::cout << " > Test attributs" << std::endl;
        std::cout << "   - charge q: " << part.charge_m << " (expected: " << 1 << ")" << std::endl;
        std::cout << "   - mass m: " << part.mass_m << " (expected: " << 1 << ")" << std::endl;
        std::cout << "   - temperature t: " << part.temperature_m << " (expected: " << 1 << ")" << std::endl;
        for (int ip = 0; ip < nbparticle; ++ip) 
        {
            part.set(ip, 2.0, ip, ip+2, ip+3, ip, ip*2, ip*3);
        }
        std::cout << " > Test set method on host" << std::endl;
        std::cout << "particule 0 : x = " <<part.x_[0]<<"(0), y = "<< part.y_[0] <<"(2), z = "<<part.z_[0]<< "(3), mx= "<<part.mx_[0]<<"(0), my = "<< part.my_[0] <<"(0), mz = "<<part.mz_[0]<< "(0), weight= "<<part.weight_[0]<<"(2)" << std::endl;
        std::cout << "particule 5 : x = " <<part.x_[5]<<"(5), y = "<< part.y_[5] <<"(7), z = "<<part.z_[5]<< "(8), mx= "<<part.mx_[5]<<"(5), my = "<< part.my_[5] <<"(10), mz = "<<part.mz_[5]<< "(15), weight= "<<part.weight_[5]<<"(2)" << std::endl;
        std::cout << "particule 10 : x = " <<part.x_[10]<<"(10), y = "<< part.y_[10] <<"(12), z = "<<part.z_[10]<< "(13), mx= "<<part.mx_[10]<<"(10), my = "<< part.my_[10] <<"(20), mz = "<<part.mz_[10]<< "(30), weight= "<<part.weight_[10]<<"(2)" << std::endl;
        std::cout << "particule 19 : x = " <<part.x_[19]<<"(19), y = "<< part.y_[19] <<"(21), z = "<<part.z_[19]<< "(22), mx= "<<part.mx_[19]<<"(19), my = "<< part.my_[19] <<"(38), mz = "<<part.mz_[19]<< "(57), weight= "<<part.weight_[19]<<"(2)" << std::endl;
    
        // _________________________________________________________________________
        // Test get_kinetic_energy method

        std::cout << " > Test get_kinetic_energy method" << std::endl;

         mini_float kinetic_energy = 0;

        for (auto ip = 0; ip < nbparticle; ++ip) {
        const mini_float gamma = sqrt(1. + part.mx_(ip) * part.mx_(ip) + part.my_(ip) * part.my_(ip) + part.mz_(ip) * part.mz_(ip));
        kinetic_energy += part.weight_(ip) * (gamma - 1.);
        }

        auto kinetic_energy_host = part.get_kinetic_energy(minipic::host);
        auto kinetic_energy_device = part.get_kinetic_energy(minipic::device);

        std::cout << "   - kinetic_energy on host: " << kinetic_energy_host << " (expected: " << kinetic_energy << ")" << std::endl;
        std::cout << "   - kinetic_energy on device: " << kinetic_energy_device << " (expected: " << kinetic_energy<< ")" << std::endl;

        assert(std::fabs(kinetic_energy_host - kinetic_energy) < 1e-10);
        assert(std::fabs(kinetic_energy_device - kinetic_energy) < 1e-10);

        std::cout << " > Test copy from host to device" << std::endl;

        // init values on host
        for (int ip = 0; ip < nbparticle; ++ip) 
        {
            part.set(ip, 2.0, ip, ip+6, ip+9, ip, ip*10, ip*5);
        }

        mini_float kinetic_energy_ref = 0; 
        for (auto ip = 0; ip < nbparticle; ++ip) {

            const mini_float mxl = ip;
            const mini_float myl = ip*2;
            const mini_float mzl = ip*3;

            const mini_float gamma = sqrt(1. + mxl*mxl + myl*myl + mzl*mzl);
            kinetic_energy_ref += 2.0 * (gamma - 1.);
        }

        // copy to device
        part.sync(minipic::host, minipic::device);

        kinetic_energy_host = part.get_kinetic_energy(minipic::host);
        kinetic_energy_device = part.get_kinetic_energy(minipic::device);


        
        std::cout << "   - kinetic_energy on host: " << kinetic_energy_host << " (expected: " << kinetic_energy_ref << ")" << std::endl;
        std::cout << "   - kinetic_energy on device: " << kinetic_energy_device << " (expected: " << kinetic_energy_ref<< ")" << std::endl;
        
        assert(std::fabs(kinetic_energy_device - kinetic_energy_host) < 1e-10);

        // _________________________________________________________________________
        std::cout << " > Test copy from device to host" << std::endl;

#if defined(__MINIPIC_KOKKOS_COMMON__)

    #if defined(__MINIPIC_KOKKOS__)

        auto x = part.x_.data_.d_view;
        auto y = part.y_.data_.d_view;
        auto z = part.z_.data_.d_view;

        auto mx = part.mx_.data_.d_view;
        auto my = part.my_.data_.d_view;
        auto mz = part.mz_.data_.d_view;

    #elif defined(__MINIPIC_KOKKOS_UNIFIED__)

        auto x = part.x_.data_;
        auto y = part.y_.data_;
        auto z = part.z_.data_;

        auto mx = part.mx_.data_;
        auto my = part.my_.data_;
        auto mz = part.mz_.data_;

    #endif

        Kokkos::parallel_for("init_on_device",Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, nbparticle),KOKKOS_LAMBDA(const int& i) {

            x(i) = i;
            y(i) = i+2;
            z(i) = i+3;

            mx(i) = i;
            my(i) = i*2;
            mz(i) = i*3;
        });

        Kokkos::fence();

#elif defined(__MINIPIC_STDPAR__)

        // init values on device
        mini_float * mx = part.mx_.get_raw_pointer(minipic::device);
        mini_float * my = part.my_.get_raw_pointer(minipic::device);
        mini_float * mz = part.mz_.get_raw_pointer(minipic::device);
        mini_float * x = part.x_.get_raw_pointer(minipic::device);
        mini_float * y = part.y_.get_raw_pointer(minipic::device);
        mini_float * z = part.z_.get_raw_pointer(minipic::device);

        part.weight_.fill(2.0);
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), nbparticle, [=](int i) {mx[i]=i; });
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), nbparticle, [=](int i) {my[i]=i*2; });
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), nbparticle, [=](int i) {mz[i]=i*3; });
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), nbparticle, [=](int i) {x[i]=i; });
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), nbparticle, [=](int i) {y[i]=i+2; });
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), nbparticle, [=](int i) {z[i]=i+3; });

#endif

        part.sync(minipic::device, minipic::host);

        kinetic_energy_host = part.get_kinetic_energy(minipic::host);
        kinetic_energy_device = part.get_kinetic_energy(minipic::device);

        kinetic_energy_ref = 0; 
        for (auto ip = 0; ip < nbparticle; ++ip) {

            const mini_float mxl = ip;
            const mini_float myl = ip*2;
            const mini_float mzl = ip*3;

            const mini_float gamma = sqrt(1. + mxl*mxl + myl*myl + mzl*mzl);
            kinetic_energy_ref += 2 * (gamma - 1.);
        }
        
        std::cout << "   - kinetic_energy on host: "  <<kinetic_energy_host << " (expected: " << kinetic_energy_ref << ")" << std::endl;
        std::cout << "   - kinetic_energy on device: " << kinetic_energy_device << " (expected: " << kinetic_energy_ref<< ")" << std::endl;

        assert(std::fabs(kinetic_energy_device - kinetic_energy_host) < 1e-10);

        // _________________________________________________________________________
        // Test clear method
        std::cout << " > Test clear method" << std::endl;
        part.clear();

        std::cout << "   - size: " << part.size() << " (expected: 0)" << std::endl;
        
        kinetic_energy_host = part.get_kinetic_energy(minipic::host);
        kinetic_energy_device = part.get_kinetic_energy(minipic::device);
        std::cout << "   - kinetic_energy on host: " << kinetic_energy_host << " (expected: " << 0 << ")" << std::endl;
        std::cout << "   - kinetic_energy on device: " << kinetic_energy_device << " (expected: " << 0<< ")" << std::endl;

        std::cout << " > Test destructor" << std::endl;
    }
    backend.finalize();

}
