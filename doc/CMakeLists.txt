
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  # copy LifeV logo
  FILE(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/common/LifeV.png
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/api/html)

  # due to doxygen bug, these files must be copied
  FILE(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tex/lifev.bib
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  FILE(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/tex/biblio.bib
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # generate Doxyfile
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.cmake
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

  # add doc target
  add_custom_target(doc
                    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen" VERBATIM
                   )
  install(
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/api
    DESTINATION share/doc/lifev
    COMPONENT doc
    OPTIONAL
  )

endif(DOXYGEN_FOUND)

FIND_PACKAGE(LATEX)

IF(PDFLATEX_COMPILER)

  ADD_CUSTOM_TARGET(pdf echo)

  ADD_SUBDIRECTORY(CodingGuidelines)

  ADD_SUBDIRECTORY(manual)

ENDIF()
