# CMake options
include (CheckIncludeFiles)
cmake_minimum_required (VERSION 2.6)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Default compilers
set (CMAKE_C_COMPILER "mpicc" CACHE STRING "Default C compiler")
set (CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "Default C++ compiler")

# Default build type
set (CMAKE_BUILD_TYPE "Debug" CACHE STRING "Default build type")

# Default value for TPL flags
set (LIFE_USE_HDF5 FALSE CACHE STRING "By default, don't use HDF5")
set (LIFE_USE_QHULL FALSE CACHE STRING "By default, don't use QHULL")
set (LIFE_USE_SuperLU FALSE CACHE STRING "By default, don't use SuperLU")
set (LIFE_USE_SuperLUDist FALSE CACHE STRING "By default, don't use SuperLUDist")
set (LIFE_USE_UMFPACK FALSE CACHE STRING "By default, don't use UMFPACK")
set (LIFE_USE_MPI_WRAPPER TRUE CACHE STRING "By default, use the MPI wrappers")
set (LIFE_USE_SPIRIT_PARSER TRUE CACHE STRING "By default, use the Boost Spirit parser")

# Put some other variables in the cache so they can be set with ccmake
set (PARMETIS_ROOT "PARMETIS_ROOT-NOTFOUND" CACHE STRING "Location of ParMETIS")
set (UMFPACK_ROOT "UMFPACK_ROOT-NOTFOUND" CACHE STRING "Location of UMFPACK")
set (HDF5_ROOT "HDF5_ROOT-NOTFOUND" CACHE STRING "Location of HDF5")
set (QHULL_ROOT "QHULL_ROOT-NOTFOUND" CACHE STRING "Location of QHULL")
set (SuperLU_ROOT "SuperLU_ROOT-NOTFOUND" CACHE STRING "Location of SuperLU")
set (SuperLUDist_ROOT "SuperLUDist_ROOT-NOTFOUND" CACHE STRING "Location of SuperLUDist")
set (BOOST_ROOT "BOOST_ROOT-NOTFOUND" CACHE STRING "Location of Boost")

# Project options
project (LifeV)
set (PROJECT_SHORT_NAME "lifev")
set (LifeV_VERSION_MAJOR 2)
set (LifeV_VERSION_MINOR 0)
set (LifeV_VERSION_PATCH 0)
set (PROJECT_BUGREPORT "lifev-user@googlegroups.com")
set (PROJECT_HOMEPAGE "www.lifev.org")

# Languages
enable_language (C CXX Fortran)

# Here we add CXXFLAGS, LDFLAGS that are set by the user
# It's a peculiar hack: needed to reset the CMAKE_LANG_FLAGS
# variables in order to place user specified flags at the end
# (this is useful for overriding an optimization level and
# possibly other things)
set (CMAKE_C_FLAGS_SAVED ${CMAKE_C_FLAGS})
set (CMAKE_C_FLAGS "")
set (CMAKE_CXX_FLAGS_SAVED ${CMAKE_CXX_FLAGS})
set (CMAKE_CXX_FLAGS "")

set (CMAKE_C_FLAGS_DEBUG "-Wall -g")
set (CMAKE_C_FLAGS_RELEASE "-O2")
set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_SAVED}")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS_SAVED}")
set (CMAKE_CXX_FLAGS_DEBUG "-Wall -g")
set (CMAKE_CXX_FLAGS_RELEASE "-O2")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_SAVED}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_SAVED}")

add_definitions (-DTHREEDIM)

# Dependencies

# The list of dependencies will end up being quite long, as it will contain all
# the dependencies of LifeV. Current list:
#   * MPI (v)
#   * Trilinos (v)
#   * ParMETIS (v)
#   * UMFPACK (should be optional) (v)
#   * HDF5 (should be optional) (v)
#   * BLAS, LAPACK (v)
#   * Boost (v)
#   * Qhull (should be optional) (v)

# MPI
if (NOT LIFE_USE_MPI_WRAPPER)
    find_package (MPI REQUIRED)
    if (MPI_FOUND)
        include_directories (${MPI_CXX_INCLUDE_PATH})
        set (LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
    endif()
endif()

# Trilinos
find_package (Trilinos COMPONENTS anasazi belos ml ifpack aztecoo amesos zoltan epetraext teuchos epetra REQUIRED)
if (Trilinos_FOUND)
    if (Trilinos_amesos_FOUND)
        set (HAVE_TRILINOS_AMESOS TRUE)
    endif()
    if (Trilinos_anasazi_FOUND)
        set (HAVE_TRILINOS_ANASAZI TRUE)
    endif()
    if (Trilinos_aztecoo_FOUND)
        set (HAVE_TRILINOS_AZTECOO TRUE)
    endif()
    if (Trilinos_belos_FOUND)
        set (HAVE_TRILINOS_BELOS TRUE)
    endif()
    if (Trilinos_epetra_FOUND)
        set (HAVE_TRILINOS_EPETRA TRUE)
    endif()
    if (Trilinos_epetraext_FOUND)
        set (HAVE_TRILINOS_EPETRAEXT TRUE)
    endif()
    if (Trilinos_ifpack_FOUND)
        set (HAVE_TRILINOS_IFPACK TRUE)
    endif()
    if (Trilinos_ml_FOUND)
        set (HAVE_TRILINOS_ML TRUE)
    endif()
    if (Trilinos_teuchos_FOUND)
        set (HAVE_TRILINOS_TEUCHOS TRUE)
    endif()
    if (Trilinos_zoltan_FOUND)
        set (HAVE_TRILINOS_ZOLTAN TRUE)
    endif()
    include_directories (${Trilinos_INCLUDE_DIR})
    set (LIBS ${LIBS} ${Trilinos_LIBRARIES})
endif()

# ParMETIS
find_package (PARMETIS REQUIRED)
if (PARMETIS_FOUND)
    include_directories (${PARMETIS_INCLUDE_DIR})
    set (LIBS ${LIBS} ${PARMETIS_LIBRARIES})
endif()

# UMFPACK (optional)
if (LIFE_USE_UMFPACK)
    find_package (UMFPACK)
    if (UMFPACK_FOUND)
        set (HAVE_UMFPACK TRUE)
        set (HAVE_UMFPACK_H TRUE)
        set (HAVE_AMD_H TRUE)
        include_directories (${UMFPACK_INCLUDE_DIR})
        set (LIBS ${LIBS} ${UMFPACK_LIBRARIES})
    endif()
endif(LIFE_USE_UMFPACK)

# SuperLU (optional)
if (LIFE_USE_SuperLU)
    find_package (SuperLU)
    if (SuperLU_FOUND)
        set (HAVE_SuperLU TRUE)
        include_directories (${SuperLU_INCLUDE_DIRS})
        set (LIBS ${LIBS} ${SuperLU_LIBRARIES})
    endif()
endif(LIFE_USE_SuperLU)

# SuperLUDist (optional)
if (LIFE_USE_SuperLUDist)
    find_package (SuperLUDist)
    if (SuperLUDist_FOUND)
        set (HAVE_SuperLUDist TRUE)
        include_directories (${SuperLUDist_INCLUDE_DIRS})
        set (LIBS ${LIBS} ${SuperLUDist_LIBRARIES})
    endif()
endif(LIFE_USE_SuperLUDist)

# LAPACK
find_package (LAPACK REQUIRED)
if (LAPACK_FOUND)
    include_directories (${LAPACK_INCLUDE_DIR})
    set (LIBS ${LIBS} ${LAPACK_LIBRARIES})
endif()

# BLAS
find_package (BLAS REQUIRED)
if (BLAS_FOUND)
    include_directories (${BLAS_INCLUDE_DIR})
    set (LIBS ${LIBS} ${BLAS_LIBRARIES})
endif()

# HDF5 (optional)
if (LIFE_USE_HDF5)
    if (HDF5_ROOT)
        set (HDF5_INCLUDE_DIRS ${HDF5_ROOT}/include)
        set (HDF5_INCLUDE_DIR ${HDF5_ROOT}/include)
        set (HDF5_LIBRARIES "-L${HDF5_ROOT}/lib -lhdf5 -lz -lm")
    endif()
    find_package (HDF5)
    if (HDF5_FOUND)
        set (HAVE_HDF5 TRUE)
        include_directories (${HDF5_INCLUDE_DIRS})
        link_directories (${HDF5_LIBRARY_DIRS})
        set (LIBS ${LIBS} ${HDF5_LIBRARIES})
    endif()
endif(LIFE_USE_HDF5)

# Boost
find_package (Boost REQUIRED)
if (Boost_FOUND)
    set (HAVE_BOOST TRUE)
    set (HAVE_BOOST_VERSION ${Boost_VERSION})
    if (Boost_MINOR_VERSION GREATER 39)
        set (HAVE_BOOST_GT_1_39 TRUE)
    endif(Boost_MINOR_VERSION GREATER 39)
    if (Boost_MINOR_VERSION GREATER 40)
        set (HAVE_BOOST_SPIRIT_QI TRUE)
    endif(Boost_MINOR_VERSION GREATER 40)
    include_directories (${Boost_INCLUDE_DIR})
    link_directories (${Boost_LIBRARY_DIRS})
    set (LIBS ${LIBS} ${Boost_LIBRARIES})
    if (NOT LIFE_USE_SPIRIT_PARSER)
        set (DISABLE_SPIRIT_PARSER TRUE)
    endif()
endif()

# QHULL (optional)
if (LIFE_USE_QHULL)
    find_package (QHULL)
    if (QHULL_FOUND)
        set (HAVE_QHULL TRUE)
        include_directories (${QHULL_INCLUDE_DIRS})
        set (LIBS ${LIBS} ${QHULL_LIBRARIES})
    endif()
endif(LIFE_USE_QHULL)

# Make sure that all headers are in the include path
include_directories (${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
set (LIFEV_PACKAGES core array mesh alg fem filters functions solver)
foreach (component ${LIFEV_PACKAGES})
    include_directories (${PROJECT_SOURCE_DIR}/life/life${component})
    set (LIBS life${component} ${LIBS})
endforeach(component)

# Here we add the user-define value of LIBS
set (LIBS ${LIBS} ${LIFE_EXTRA_LINK_LINE})

# lifeconfig.h
string (COMPARE EQUAL ${CMAKE_BUILD_TYPE} "Debug" HAVE_LIFEV_DEBUG)
string (COMPARE EQUAL ${CMAKE_BUILD_TYPE} "Release" HAVE_LIFEV_OPTO3)

check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(memory.h HAVE_MEMORY_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(unistd.h HAVE_UNISTD_H)

configure_file(
  "${PROJECT_SOURCE_DIR}/lifeconfig_cmake.h.in"
  "${PROJECT_BINARY_DIR}/lifeconfig.h"
)
install(FILES "${PROJECT_BINARY_DIR}/lifeconfig.h"
        DESTINATION include)

# The following subdirectories are part of the project
add_subdirectory(life)
add_subdirectory(examples)

include (CTest)
enable_testing()
add_subdirectory(testsuite)

# Documentation (not yet working, don't activate)
if (LIFE_BUILD_DOCS)
    add_subdirectory(doc)
endif()
