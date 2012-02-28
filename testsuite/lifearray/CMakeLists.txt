# Add the executable
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_executable (test_meshentitycontainer EXCLUDE_FROM_ALL test_meshentitycontainer.cpp)
target_link_libraries (test_meshentitycontainer ${LIBS})

add_executable (test_vectorsmall EXCLUDE_FROM_ALL test_vectorsmall.cpp)
target_link_libraries (test_vectorsmall ${LIBS})

add_test (NAME LifeArray::MeshEntityContainer COMMAND mpirun -np 1 test_meshentitycontainer)
add_test (NAME LifeArray::VEctorSmall COMMAND mpirun -np 1 test_vectorsmall)
