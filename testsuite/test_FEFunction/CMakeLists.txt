# Add the executable
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_FEFunction EXCLUDE_FROM_ALL main.cpp fefct.cpp user_fun.cpp)
target_link_libraries (test_FEFunction ${LIBS})

add_custom_target (test_FEFunction_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_FEFunction DEPENDS test_FEFunction_links)

add_test (NAME LifeFem::FEFunction COMMAND test_FEFunction)
