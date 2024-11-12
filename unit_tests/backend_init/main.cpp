// _________________________________________________________________________
//
// Unit test - Backend initialization
//
//! \brief test the Backend class initialization for different backends
//
// _________________________________________________________________________

#include "Backend.hpp"
#include "Params.hpp"

int main(int argc, char *argv[]) {

    std::cout << " ___________________________________________________________ " << std::endl;
    std::cout << "|                                                           |" << std::endl;
    std::cout << "|                      Backend init                         |" << std::endl;
    std::cout << "|___________________________________________________________|" << std::endl;

    Params params;

    Backend backend;

    backend.init(argc, argv, params);
    {
        backend.info();
    }

    backend.finalize();

    return 0;
}
