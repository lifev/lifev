# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_darcy EXCLUDE_FROM_ALL main.cpp darcy.cpp user_fun.cpp)
target_link_libraries (test_darcy ${LIBS})

add_custom_target (test_darcy_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_darcy DEPENDS test_darcy_links)

add_test (NAME LifeSolver::Darcy COMMAND mpirun -np 2 test_darcy)
