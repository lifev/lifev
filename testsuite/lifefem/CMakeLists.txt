# Add the executable
add_executable (test_currentFE EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_currentFE ${LIBS})

add_executable (test_quadRule EXCLUDE_FROM_ALL test_quadrule.cpp SetOfFun.cpp)
target_link_libraries (test_quadRule ${LIBS})

add_executable (test_interpolate EXCLUDE_FROM_ALL test_interpolate.cpp)
target_link_libraries (test_interpolate ${LIBS})

add_custom_target (lifefem_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_currentFE DEPENDS lifefem_links)
add_dependencies (test_quadRule DEPENDS lifefem_links)
add_dependencies (test_interpolate DEPENDS lifefem_links)

add_test (NAME LifeFem::CurrentFE COMMAND test_currentFE)
add_test (NAME LifeFem::QuadRule COMMAND test_quadRule)
add_test (NAME LifeFem::Interpolate COMMAND test_interpolate)
