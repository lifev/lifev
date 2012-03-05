# Add the executable
add_executable (test_getPot EXCLUDE_FROM_ALL main_getPot.cpp)
target_link_libraries (test_getPot ${LIBS})

add_custom_target (test_getPot_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data_getPot;data_getPot
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_getPot DEPENDS test_getPot_links)

add_test (NAME LifeCore::GetPot COMMAND test_getPot)
