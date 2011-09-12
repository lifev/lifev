set (EXAMPLE_LIST example_cylinder)

foreach (component ${EXAMPLE_LIST})
    add_subdirectory(${component})
endforeach()

add_custom_target (examples
                   DEPENDS ${EXAMPLE_LIST})
