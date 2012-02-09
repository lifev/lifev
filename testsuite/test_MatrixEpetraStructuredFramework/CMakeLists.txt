# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_MatrixEpetraStructuredFramework EXCLUDE_FROM_ALL main.cpp )
target_link_libraries (test_MatrixEpetraStructuredFramework ${LIBS})

add_test (NAME LifeArray::MatrixEpetraStructuredFramework COMMAND mpirun -np 2 test_MatrixEpetraStructuredFramework )
