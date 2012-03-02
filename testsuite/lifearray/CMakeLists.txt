# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_meshentitycontainer EXCLUDE_FROM_ALL test_meshentitycontainer.cpp)
target_link_libraries (test_meshentitycontainer ${LIBS})

add_executable (test_vectorsmall EXCLUDE_FROM_ALL test_vectorsmall.cpp)
target_link_libraries (test_vectorsmall ${LIBS})

add_executable (test_ghosthandler EXCLUDE_FROM_ALL test_ghosthandler.cpp)
target_link_libraries (test_ghosthandler ${LIBS})

add_custom_target (test_ghosthandler_links
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data_ghosthandler;data_ghosthandler
                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/testsuite/data/mesh/inria;Mesh
                   VERBATIM)
                   
# Fake dependency to ensure that links are made prior to running the test
add_dependencies (test_ghosthandler DEPENDS test_ghosthandler_links)

add_test (NAME LifeArray::MeshEntityContainer COMMAND mpirun -np 1 test_meshentitycontainer)
add_test (NAME LifeArray::VectorSmall COMMAND mpirun -np 1 test_vectorsmall)
add_test (NAME LifeArray::GhostHandler COMMAND mpirun -np 2 test_ghosthandler)
