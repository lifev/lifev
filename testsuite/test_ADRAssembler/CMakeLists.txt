# Add the executable
add_executable (test_ADRAssembler EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_ADRAssembler ${LIBS})

add_custom_target (test_ADRAssembler_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_ADRAssembler DEPENDS test_ADRAssembler_links)

# Add the executable
add_executable (test_ADRAssembler_2d EXCLUDE_FROM_ALL main_2d.cpp)
target_link_libraries (test_ADRAssembler_2d ${LIBS})

add_custom_target (test_ADRAssembler_2d_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data_2d;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_ADRAssembler_2d DEPENDS test_ADRAssembler_2d_links)

add_test (NAME LifeSolver::ADRAssembler COMMAND mpirun -np 2 test_ADRAssembler)
add_test (NAME LifeSolver::ADRAssembler_2D COMMAND mpirun -np 2 test_ADRAssembler_2d)
