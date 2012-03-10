
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  # copy LifeV logo
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/common/LifeV.png
                 ${CMAKE_CURRENT_BINARY_DIR}/api/html/LifeV.png COPYONLY)
  # this can be substituted only with cmake >= 2.8
  # file(COPY yourDir DESTINATION yourDestination)
 
  # generate Doxyfile
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.cmake
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

  # add doc target
  add_custom_target(doc
                    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen" VERBATIM
                   )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/api DESTINATION share/doc/lifev)

endif(DOXYGEN_FOUND)

FIND_PACKAGE(LATEX)

IF(PDFLATEX_COMPILER)

  ADD_CUSTOM_TARGET(pdf echo)

  ADD_SUBDIRECTORY(CodingGuidelines)

ENDIF()
