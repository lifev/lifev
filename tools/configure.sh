#!/bin/bash

MY_LIB_DIR=/u/cmcs/rpopescu/lib/opt

     

cmake \
    -D LifeV_ENABLE_TESTS:BOOL=ON \
    -D LifeV_ENABLE_STRONG_CXX_COMPILE_WARNINGS:BOOL=OFF \
    -D LifeV_ENABLE_ALL_PACKAGES:BOOL=ON \
    -D LifeV_ENABLE_OneDFSI:BOOL=ON \
    -D LifeV_ENABLE_ZeroDimensional:BOOL=ON \
    -D LifeV_ENABLE_FSI:BOOL=ON \
    -D HDF5_INCLUDE_DIRS:PATH=$MY_LIB_DIR/hdf5/include \
    -D HDF5_LIBRARY_DIRS:PATH=$MY_LIB_DIR/hdf5/lib \
    -D ParMETIS_INCLUDE_DIRS:PATH=$MY_LIB_DIR/parmetis/include \
    -D ParMETIS_LIBRARY_DIRS:PATH=$MY_LIB_DIR/parmetis/lib \
    -D Trilinos_INCLUDE_DIRS:PATH=$MY_LIB_DIR/trilinos_test/include \
    -D Trilinos_LIBRARY_DIRS:PATH=$MY_LIB_DIR/trilinos_test/lib \
    -D CMAKE_INSTALL_PREFIX:PATH=$MY_LIB_DIR/new_life \
    $* \
    ../../
