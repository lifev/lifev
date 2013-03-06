#include <lifev/core/LifeV.hpp>

using namespace LifeV;

int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    std::cout << "End Result: TEST PASSED" << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
