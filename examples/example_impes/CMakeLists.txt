# Add the executable
add_executable (example_impes EXCLUDE_FROM_ALL main.cpp impes.cpp user_fun.cpp)
target_link_libraries (example_impes ${LIBS})

add_custom_target (example_impes_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria;Mesh
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   VERBATIM)

add_dependencies (example_impes DEPENDS example_impes_links)

