// _________________________________________________________________________
//
// Unit test - Vector class
//
//! \brief test the vector class
//
// _________________________________________________________________________

#include "Headers.hpp"
#include "Field.hpp"
#include "Vector.hpp"

#include <cassert>


int main(int argc, char *argv[]) {

    std::cout << " ___________________________________________________________ " << std::endl;
    std::cout << "|                                                           |" << std::endl;
    std::cout << "|                      Vector class                         |" << std::endl;
    std::cout << "|___________________________________________________________|" << std::endl;

    Params params;

    Backend backend;

    backend.init(argc, argv, params);
    {

        // _________________________________________________________________________
        // Parameters

        int size = 100;
        #if defined(__MINIPIC_SYCL__)
        auto sycl_queue_ptr = backend.sycl_queue_;
        #endif

        // _________________________________________________________________________
        // Test constructor

        std::cout << " > Test constructor" << std::endl;

        Vector <mini_float> v1(size, 1.0, backend);
        //Vector <int> v1(size, 1.0, backend);

        // _________________________________________________________________________
        // Test size method

        std::cout << " > Test size method" << std::endl;
        std::cout << "   - size: " << v1.size() << " (expected: " << size << ")" << std::endl;

        assert(v1.size() == size);

        // _________________________________________________________________________
        // Test fill method

        std::cout << " > Test fill method" << std::endl;

        v1.fill(2.0);
        //v1.fill(19);

        auto sum_host = v1.sum(1, minipic::host);

        std::cout << " - sum on host: " << sum_host << " (expected: " << size*2.0 << ")" << std::endl;

        auto sum_device = v1.sum(1, minipic::device);

        std::cout << " - sum on device: " << sum_device << " (expected: " << size*2.0 << ")" << std::endl;

        assert(sum_host == size*2.0);
        assert(sum_device == size*2.0);

        // _________________________________________________________________________
        // test kernel on host

#if defined(__MINIPIC_KOKKOS_UNIFIED__)

        auto view = v1.data_;

        Kokkos::parallel_for("init_on_host",Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, size),KOKKOS_LAMBDA(const int& i) {
            view(i) = i;
        });

#else

        for (auto i = 0; i < v1.size(); ++i) {
            v1[i] = i;
        }


#endif

        // _________________________________________________________________________
        // test sum on host


        std::cout << " > Test sum on host" << std::endl;

        sum_host = v1.sum(1, minipic::host);

        std::cout << " - sum on host: " << sum_host << std::endl;

        // _________________________________________________________________________
        // test copy from host to device

        std::cout << " > Test copy from host to device" << std::endl;

        // copy to device
        v1.sync(minipic::host, minipic::device);

        sum_device = v1.sum(1, minipic::device);

        std::cout << " - sum on device: " << sum_device << std::endl;

        assert(sum_host == sum_device);

        // _________________________________________________________________________
        // Test kernel on device

        std::cout << " > Test kernel on device" << std::endl;
        
        // init values on device
#if defined(__MINIPIC_SYCL__)

        mini_float *const __restrict__ device_data = v1.get_raw_pointer(minipic::device);

        sycl::range<1> n_particles {static_cast<size_t>(size)};

        // fill on device
        sycl_queue_ptr->parallel_for(n_particles, [=](sycl::id<1> i) {
        device_data[i] = i*i;
        });

        sycl_queue_ptr->wait();

#elif defined(__MINIPIC_KOKKOS__)

        device_vector_t device_view = v1.data_.d_view;

        Kokkos::parallel_for(size, KOKKOS_LAMBDA(const int i) {
            device_view(i) = i*i;
        });
        Kokkos::fence();

#elif defined(__MINIPIC_KOKKOS_UNIFIED__)

        device_vector_t device_view = v1.data_;

        Kokkos::parallel_for(size, KOKKOS_LAMBDA(const int i) {
            device_view(i) = i*i;
        });
        Kokkos::fence();

#elif defined(__MINIPIC_STDPAR__)
        mini_float *const ptrv1 = v1.get_raw_pointer(minipic::device);
        std::for_each_n(std::execution::par_unseq, counting_iterator(0), size, [=](int i) {ptrv1[i]=i*i; });
#else
        for (int i = 0; i < size; i++) {
            v1[i] = i*i;
        }
#endif

        // Compute the reference on host
        mini_float reference_sum = 0.0;
        for (int i = 0; i < size; i++) {
            reference_sum += i*i;
        }


        // _________________________________________________________________________
        // Test sum on device

        sum_device = v1.sum(1, minipic::device);

        auto error_device = std::abs((sum_device - reference_sum) / reference_sum);

        std::cout << " - sum on device: " << sum_device << " with error: " << error_device << std::endl;


        // _________________________________________________________________________
        // Test copy from device to host

        std::cout << " > Test copy from device to host" << std::endl;

        v1.sync(minipic::device, minipic::host);

        // _________________________________________________________________________
        // Test sum on host

        sum_host = v1.sum(1, minipic::host);

        auto error_host = std::abs((sum_host - reference_sum) / reference_sum);

        std::cout << " - sum on host: " << sum_host << " with error: " << error_host << std::endl;

        assert(sum_host == sum_device);

    }

    std::cout << " > Test destructor" << std::endl;

    backend.finalize();

}
