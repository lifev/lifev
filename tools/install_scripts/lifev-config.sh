#! /bin/bash

LIFEV_SRC_DIR=${HOME}/lifev-env/lifev
LIFEV_BUILD_DIR=${HOME}/lifev-env/lifev-build
LIFEV_INSTALL_DIR=${HOME}/lifev-env/lifev-install

TRILINOS_HOME=${HOME}/lifev-env/install/lib/trilinos

mkdir -p ${LIFEV_BUILD_DIR}
cd ${LIFEV_BUILD_DIR}

rm -f CMakeCache.txt
cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=${LIFEV_INSTALL_DIR} \
  -D CMAKE_BUILD_TYPE:STRING=Release \
 \
  -D Trilinos_INCLUDE_DIRS:PATH=${TRILINOS_HOME}/include \
  -D Trilinos_LIBRARY_DIRS:PATH=${TRILINOS_HOME}/lib \
  \
  -D LifeV_ENABLE_ALL_PACKAGES:BOOL=ON \
  -D LifeV_ENABLE_TESTS:BOOL=ON \
  \
  $* \
  ${LIFEV_SRC_DIR} || return 1
