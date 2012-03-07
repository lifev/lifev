# Add the executable
add_executable (test_VectorContainer EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_VectorContainer ${LIBS})

add_test (NAME LifeArray::VectorContainer COMMAND test_VectorContainer)
