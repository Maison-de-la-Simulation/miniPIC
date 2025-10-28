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
#include <iomanip>

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

        std::cout << " > Test constructor on Host" << std::endl;
        auto sum_host = v1.sum(1, minipic::host);
        std::cout << " - sum on host: " << sum_host << " (expected: " << size*1.0 << ")" << std::endl;

        std::cout << " > Test constructor on Device" << std::endl;
        auto sum_device = v1.sum(1, minipic::device);
        std::cout << " - sum on device: " << sum_device << " (expected: " << size*1.0 << ")" << std::endl;

        assert(sum_host == size*1.0);
        assert(sum_device == size*1.0);

        std::cout << std::endl;


        // _________________________________________________________________________
        // Test size method

        std::cout << " > Test size method" << std::endl;
        std::cout << "   - size: " << v1.size() << " (expected: " << size << ")" << std::endl;

        assert(v1.size() == size);
        std::cout<< std::endl;

        // _________________________________________________________________________
        // Test fill method
        std::cout << " > Test fill method on Host" << std::endl;
        v1.fill(2.0, minipic::host);
        sum_host = v1.sum(1, minipic::host);
        std::cout << " - sum on host: " << sum_host << " (expected: " << size*2.0 << ")" << std::endl;

        std::cout << " > Test fill method on Device" << std::endl;
        v1.fill(2.0, minipic::device);
        sum_device = v1.sum(1, minipic::device);
        std::cout << " - sum on device: " << std::setprecision(15) << sum_device << " (expected: " << size*2.0 << ")" << std::endl;

        assert(sum_host == size*2.0);
        assert(sum_device == size*2.0);
        //assert(std::abs(sum_device - size * 3.0) < 1e-5);

        std::cout << std::endl;


        // _________________________________________________________________________
        // test kernel on host


        for (auto i = 0; i < v1.size(); ++i) {
            v1[i] = i;
        }

        // _________________________________________________________________________
        // test sum on host
        std::cout << " > Test sum on host" << std::endl;
        sum_host = v1.sum(1, minipic::host);
        std::cout << " - sum on host: " << sum_host << " (expected: \"4950\")" << std::endl;
#if defined (__MINIPIC_CUDA__)
        sum_device = v1.sum(1, minipic::device);
        std::cout << " - sum on device: " << sum_device << " (expected: " << size*3.0 << ")" << std::endl;
#else
#endif        
        std::cout << std::endl;

        // _________________________________________________________________________
        // test copy from host to device

        std::cout << " > Test copy from host to device" << std::endl;

        // copy to device
        v1.sync(minipic::host, minipic::device);
        sum_device = v1.sum(1, minipic::device);
        std::cout << " - sum on device: " << sum_device << " (expected: \"4950\")" << std::endl;
        assert(sum_host == sum_device);
        std::cout << std::endl;


        // _________________________________________________________________________
        // Test kernel on device

        std::cout << " > Test kernel on device" << std::endl;
        #if (__MINIPIC_STDPAR__)
        mini_float* v1_d= v1.get_raw_pointer(minipic::device);
        std::for_each(std::execution::par_unseq, counting_iterator(0),
                              counting_iterator(size),
                  [=](int idx) {
                      v1_d[idx] = idx * idx;      
                  });
        #else
        // init values on device
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
        std::cout << " - sum on host: " << sum_host << " (expected: \"4950\")" << std::endl;
        std::cout << std::endl;

        // _________________________________________________________________________
        // Test copy from device to host
        std::cout << " > Test copy from device to host" << std::endl;

        v1.sync(minipic::device, minipic::host);
        sum_host = v1.sum(1, minipic::host);
        error_device = std::abs((sum_host - reference_sum) / reference_sum);
        std::cout << " - sum on host: " << sum_host << " (expected: " << reference_sum << ")" << " with error: " << error_device << std::endl;
        assert(sum_host == sum_device);
        std::cout << std::endl;

    }

    std::cout << " > Test destructor" << std::endl;

    backend.finalize();

}

