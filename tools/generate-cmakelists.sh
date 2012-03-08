#! /bin/bash
# this script generates the list of include files and source files 
# to be used with TriBITS build system

for find_dir in $(find ./* -type d -not -path '*/.*/*' -not -path "./cmake" -not -path "./testsuite" -not -path "./example"); do
  dir=$(echo $find_dir | sed -e "s:\./::g")
  echo "$dir"
  echo "SET(${dir}_HEADERS" > $dir/CMakeLists.txt
  ls -1 $dir/*.hpp >> $dir/CMakeLists.txt
  echo "CACHE INTERNAL \"\")" >> $dir/CMakeLists.txt
  echo >> $dir/CMakeLists.txt
  echo "SET(${dir}_SOURCES"  >> $dir/CMakeLists.txt
  ls -1 $dir/*.cpp >> $dir/CMakeLists.txt
  echo "CACHE INTERNAL \"\")" >> $dir/CMakeLists.txt
done
