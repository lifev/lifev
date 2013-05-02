cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Default compilers
set(CMAKE_C_COMPILER "mpicc" CACHE STRING "Default C compiler")
set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "Default C++ compiler")

# Default build type
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default build type")

set(APPLIFEV_USE_MPI_WRAPPER ON CACHE BOOL "By default, use the MPI wrappers")

# Project options
project(AppLifeV)
set(PROJECT_SHORT_NAME "app-lifev")

# Languages
enable_language(C CXX)

# MPI
if(NOT APPLIFEV_USE_MPI_WRAPPER)
    find_package (MPI REQUIRED)
    if(MPI_FOUND)
        include_directories(${MPI_CXX_INCLUDE_PATH})
        set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
    endif()
endif()

# LIFEV
find_package(LIFEV REQUIRED)
message(STATUS "LifeV include dirs: " "${LIFEV_INCLUDE_DIRS}")
include_directories(${LIFEV_INCLUDE_DIRS})
message(STATUS "LifeV libraries: " "${LIFEV_LIBRARIES}")
set(LIBS ${LIBS} ${LIFEV_LIBRARIES})

# dummy test application
add_subdirectory(dummy)

