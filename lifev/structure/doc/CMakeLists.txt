
SET(TEX_FILES
  documentation
)

CONFIGURE_FILE(ltxcompile.cmake
               ltxcompile @ONLY)

FOREACH(TEX_FILE ${TEX_FILES})
  ADD_CUSTOM_TARGET(pdf_${TEX_FILE} echo
    DEPENDS ${TEX_FILE}.pdf
  )
  ADD_DEPENDENCIES(pdf pdf_${TEX_FILE})

  ADD_CUSTOM_COMMAND(
    OUTPUT    ${TEX_FILE}.pdf
    COMMAND   /bin/sh ltxcompile
    ARGS      ${TEX_FILE}.tex
    DEPENDS   ${TEX_FILE}.tex
  )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${TEX_FILE}.pdf
    DESTINATION share/doc/lifev/structure
    COMPONENT doc
    OPTIONAL
  )
ENDFOREACH()

