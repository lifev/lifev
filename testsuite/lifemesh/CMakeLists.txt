# Add the executable
add_executable (test_mesh EXCLUDE_FROM_ALL test_mesh.cpp)
target_link_libraries (test_mesh ${LIBS})

add_custom_target (test_mesh_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/mesh++/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_mesh DEPENDS test_mesh_links)

add_test (NAME LifeMesh::Mesh COMMAND test_mesh)
