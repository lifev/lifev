# Add the executable
add_executable (test_ZeroDimensional EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_ZeroDimensional ${LIBS})

add_custom_target (test_ZeroDimensional_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/RythmosAztecOOParams.xml;RythmosAztecOOParams.xml
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/RythmosAztecOOParamsLowsf.xml;RythmosAztecOOParamsLowsf.xml
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/RythmosBelosParams.xml;RythmosBelosParams.xml
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data.dat;data.dat
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/circuitFile.dat;circuitFile.dat
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/inputFile.dat;inputFile.dat
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_ZeroDimensional DEPENDS test_ZeroDimensional_links)

add_test (NAME LifeSolver::ZeroDimensional COMMAND test_ZeroDimensional)
