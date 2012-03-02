# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_navierStokes EXCLUDE_FROM_ALL mainRossEthierSteinman.cpp ../../life/lifefunctions/RossEthierSteinmanDec.cpp)
target_link_libraries (test_navierStokes ${LIBS})

add_custom_target (test_navierStokes_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/dataRossEthierSteinman;dataRossEthierSteinman
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_navierStokes DEPENDS test_navierStokes_links)

# Add the executable
add_executable (test_navierStokes_2d EXCLUDE_FROM_ALL mainKimMoin.cpp ../../life/lifefunctions/KimMoin.cpp)
target_link_libraries (test_navierStokes_2d ${LIBS})

add_custom_target (test_navierStokes_2d_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/dataKimMoin;dataKimMoin
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh;Mesh
                   VERBATIM)

# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_navierStokes_2d DEPENDS test_navierStokes_2d_links)

# add_test (NAME LifeSolver::navierStokes COMMAND mpirun -np 2 test_navierStokes -f dataRossEthierSteinman)
add_test (NAME LifeSolver::navierStokes COMMAND mpirun -np 2 test_navierStokes
-f dataRossEthierSteinman)
add_test (NAME LifeSolver::navierStokes_2D COMMAND mpirun -np 2 test_navierStokes_2d -f dataKimMoin)
