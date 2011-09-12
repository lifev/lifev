# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_LevelSetSolver EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_LevelSetSolver ${LIBS})

add_custom_target (test_LevelSetSolver_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data;data
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_LevelSetSolver DEPENDS test_LevelSetSolver_links)

add_test (NAME LifeSolver::LevelSetSolver COMMAND mpirun -np 3 test_LevelSetSolver)
