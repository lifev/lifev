# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_hyperbolic EXCLUDE_FROM_ALL main.cpp hyperbolic.cpp)
target_link_libraries (test_hyperbolic ${LIBS})

add_custom_target (test_hyperbolic_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_hyperbolic DEPENDS test_hyperbolic_links)

add_test (NAME LifeSolver::Hyperbolic COMMAND mpirun -np 1 test_hyperbolic)
