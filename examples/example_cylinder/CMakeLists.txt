# Add the executable
add_executable (example_cylinder EXCLUDE_FROM_ALL main.cpp cylinder.cpp)
target_link_libraries (example_cylinder ${LIBS})

add_custom_target (example_cylinder_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria;Mesh
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   VERBATIM)

add_dependencies (example_cylinder DEPENDS example_cylinder_links)

