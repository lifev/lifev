# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_timeAdvance1 EXCLUDE_FROM_ALL main.cpp timeAdvance.cpp)
target_link_libraries (test_timeAdvance1 ${LIBS})

add_custom_target (test_timeAdvance1_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_timeAdvance1 DEPENDS test_timeAdvance1_links)

add_test (NAME LifeSolver::TimeAdvanceOrder1 COMMAND mpirun -np 2 test_timeAdvance1)
