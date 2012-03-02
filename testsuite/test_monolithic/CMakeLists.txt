# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_monolithic EXCLUDE_FROM_ALL main.cpp ud_functions.cpp flowConditions.cpp lumpedHeart.cpp)
target_link_libraries (test_monolithic ${LIBS})

add_custom_target (test_monolithic_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/dataCE;dataCE
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/dataGCE;dataGCE
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria/;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_monolithic DEPENDS test_monolithic_links)

add_test (NAME LifeSolver::MonolithicCE COMMAND mpirun -np 2 test_monolithic -f dataCE)
add_test (NAME LifeSolver::MonolithicGCE COMMAND mpirun -np 2 test_monolithic -f dataGCE)
