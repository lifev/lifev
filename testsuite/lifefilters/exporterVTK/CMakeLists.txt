# Add the executable
add_executable (testExport EXCLUDE_FROM_ALL mainExport.cpp)
target_link_libraries (testExport ${LIBS})

add_executable (testImport EXCLUDE_FROM_ALL mainImport.cpp)
target_link_libraries (testImport ${LIBS})

add_custom_target (exporter_vtk_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria;Mesh
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (testExport DEPENDS exporter_vtk_links)
add_dependencies (testImport DEPENDS exporter_vtk_links)

# add_test (NAME LifeFilters::TestExport COMMAND testExport)
# add_test (NAME LifeFilters::TestImport COMMAND testImport)
