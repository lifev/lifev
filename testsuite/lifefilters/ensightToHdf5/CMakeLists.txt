# Add the executable
add_executable (test_filters EXCLUDE_FROM_ALL main.cpp ensightToHdf5.cpp)
target_link_libraries (test_filters ${LIBS})

add_custom_target (test_filters_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria;Mesh
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR};importDir
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_filters DEPENDS test_filters_links)

add_test (NAME LifeFilters::EnsightToHDF5 COMMAND mpirun -np 4 test_filters)
