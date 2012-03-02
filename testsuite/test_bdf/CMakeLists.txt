# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_bdf EXCLUDE_FROM_ALL main.cpp test_bdf.cpp)
target_link_libraries (test_bdf ${LIBS})

add_custom_target (test_bdf_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data1;data1
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data2;data2
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data3;data3
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_bdf DEPENDS test_bdf_links)

# Add the executable
add_executable (test_bdf_2d EXCLUDE_FROM_ALL main.cpp test_bdf_2d.cpp)
target_link_libraries (test_bdf_2d ${LIBS})

add_custom_target (test_bdf_2d_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data_2d;data_2d
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_bdf_2d DEPENDS test_bdf_2d_links)

add_test (NAME LifeSolver::BDF_Order1 COMMAND mpirun -np 2 test_bdf -f data1)
add_test (NAME LifeSolver::BDF_Order2 COMMAND mpirun -np 2 test_bdf -f data2)
add_test (NAME LifeSolver::BDF_Order3 COMMAND mpirun -np 2 test_bdf -f data3)
add_test (NAME LifeSolver::BDF_2D COMMAND mpirun -np 2 test_bdf_2d -f data_2d)
