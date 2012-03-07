#! /bin/bash
# this script generates the list of include files and source files 
# to be used with TriBITS build system

for dir in $(find . -type d); do
  echo "$dir"
  cd $dir
  echo "INCLUDE_DIRECTORIES(\${CMAKE_CURRENT_SOURCE_DIR})" >> CMakeLists.txt
  echo >> CMakeLists.txt
  echo "SET(HEADERS \${HEADERS}" >> CMakeLists.txt
  ls -1 *.hpp >> CMakeLists.txt
  echo ")" >> CMakeLists.txt
  echo >> CMakeLists.txt
  echo "SET(SOURCES \${SOURCES}"  >> CMakeLists.txt
  ls -1 *.cpp >> CMakeLists.txt
  echo ")" >> CMakeLists.txt
  cd ..
done
