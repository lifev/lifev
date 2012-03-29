===================================
==  OUT OF DATE as of 9.03.2012  ==
==  NEEDS UPDATING !             ==
===================================

For configuring and building LifeV with CMake the following variables are of importance
(provide them as arguments to the cmake command -D VARIABLE_NAME:VARIABLE_TYPE=VALUE):

* CMAKE_CXX_COMPILER (STRING) - set this to the desired C++ compiler (g++, c++, mpicxx, mpixlcxx etc.)
        By default, this is set to "mpicxx"
* CMAKE_CXX_FLAGS (STRING) - provide a string containing additional compiler flags (this variable
        can be used to override the default settings, like optimization level, of the build type
* CMAKE_BUILD_TYPE (STRING) - set this to Debug ( -Wall -g) or Release (-O2).
        By default, this is set to "Debug"

* LIFE_EXTRA_LINK_LINE (STRING) - use this to specify flags for the linker or append additional libraries
        to the link line

* CMAKE_INSTALL_PREFIX (PATH) - set this to the desired install location for LifeV

Third party library dependencies:

* ParMETIS (required) - if it is installed in a non-standard location, you will need to set
        PARMETIS_ROOT (PATH) to the appropriate location

* MPI - if you plan on using the MPI wrappers (mpicxx, mpixlcxx etc.), then you need to set
        LIFE_USE_MPI_WRAPPER (BOOL), and no additional information is required about the MPI installation.
        Otherwise, the scripts will try to determine the location of MPI headers and libraries.
        In case of a non standard installation, please provide the location of MPI through the
        MPI_INCLUDE_PATH (PATH) and MPI_LIBRARIES (STRING) variables.

* Trilinos (required) - if the Trilinos libraries aren't installed in a standard location,
        please set Trilinos_ROOT (PATH) to point to the installation

* UMFPACK (optional) - if planning to compile LifeV with UMFPACK support, please set LIFE_USE_UMFPACK
        (BOOL) to TRUE. In the event of a non standard UMFPACK installation, please provide a hint
        through UMFPACK_ROOT (PATH). Please note that the scripts expect that both UMFPACK, UFConfig and
        AMD are installed in the "include" and "lib" subdirectories of UMFPACK_ROOT.

* LAPACK and BLAS (required) - in the case of non-standard location and name for the LAPACK and BLAS
        libraries, they can be specified through the LAPACK_LIBRARIES (STRING) and BLAS_LIBRARIES (STRING)
        variables.

* Boost (required) - in the case of a non-standard Boost installation, please provide it's location
        through the BOOST_ROOT (PATH) variable

* HDF5 (optional) - the HDF5 support in LifeV can be activated with LIFE_USE_HDF5 (BOOL).
        Non standard paths can be provided through HDF5_ROOT (PATH)

* QHULL (optional) - the QHULL support in LifeV is activated with LIFE_USE_QHULL (BOOL).
        Non standard paths can be provided with QHULL_ROOT (PATH)

All the optional third party libraries are disabled by default.

After the configuration phase, compile the libraries with "make", tests can be compiled individually with
"make <NAME_OF_TEST>" or collectively with "make all_tests". Running the testsuite is done by issuing
"ctest".



The configure scripts respond to other default CMake variables. These aren't needed in most cases, but
additional information about them can be found on the official CMake website.


The collection of CMake scripts for LifeV is work in progress.

----------------------------------------------------------------------------
Radu Popescu, EPFL - CMCS
<radu.popescu@epfl.ch>
