# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_partition EXCLUDE_FROM_ALL main_partition.cpp)
target_link_libraries (test_partition ${LIBS})

add_executable (test_cylinder EXCLUDE_FROM_ALL main_cylinder.cpp cylinder.cpp)
target_link_libraries (test_cylinder ${LIBS})

add_custom_target (test_partition_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data_partition;data_partition
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data_cylinder;data_cylinder
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh;run_test.sh
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_partition DEPENDS test_partition_links)
add_dependencies (test_cylinder DEPENDS test_partition_links)

add_test (NAME LifeSolver::Partition_Offline COMMAND mpirun -np 1 test_partition -f data_partition)
add_test (NAME LifeSolver::Partition_Online COMMAND mpirun -np 4 test_cylinder -f data_cylinder)
