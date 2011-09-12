# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_fsi EXCLUDE_FROM_ALL main.cpp ud_functions.cpp)
target_link_libraries (test_fsi ${LIBS})

add_custom_target (test_fsi_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_fsi DEPENDS test_fsi_links)

add_test (NAME LifeSolver::FSI COMMAND mpirun -np 4 test_fsi)
