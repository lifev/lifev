# Add the executable
add_executable (test_Parser EXCLUDE_FROM_ALL main.cpp)
target_link_libraries (test_Parser ${LIBS})

add_test (NAME LifeCore::Parser COMMAND test_Parser)
