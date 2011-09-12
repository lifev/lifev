# Add the executable
add_executable (test_ADRAssembler EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_ADRAssembler ${LIBS})

add_custom_target (test_ADRAssembler_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_ADRAssembler DEPENDS test_ADRAssembler_links)

add_test (NAME LifeSolver::ADRAssembler COMMAND test_ADRAssembler)
