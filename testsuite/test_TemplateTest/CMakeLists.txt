# Add the executable
add_executable (test_TemplateTest EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_TemplateTest ${LIBS})

add_custom_target (test_TemplateTest_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data.txt;data.txt
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_TemplateTest DEPENDS test_TemplateTest_links)

add_test (NAME LifeCore::TemplateTest COMMAND test_TemplateTest)
