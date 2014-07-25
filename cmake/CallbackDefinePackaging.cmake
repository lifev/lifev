INCLUDE(TribitsGlobalMacros)

MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING)

  #MESSAGE("TRIBITS_REPOSITORY_DEFINE_PACKAGING() called for LifeV!")

  GET_FILENAME_COMPONENT(LifeV_SOURCE_PATH ${LifeV_SOURCE_DIR} PATH)

  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    /.git/
    ".gitignore"
  )
  
  #removing any packages not enabled from the tarball
  set(ENABLED_FLAG OFF)
  set(INCLUDE_EMPTY TRUE)
  TRIBITS_GET_ENABLED_LIST(${PROJECT_NAME}_PACKAGES ${PROJECT_NAME} ${ENABLED_FLAG} ${INCLUDE_EMPTY} 
    NON_ENABLED_PACKAGES NUM_NON_ENABLED)
  STRING(REPLACE " " ";" NON_ENABLED_PACKAGES "${NON_ENABLED_PACKAGES}")

  FOREACH(TRIBITS_PACKAGE ${NON_ENABLED_PACKAGES})
    LIST(FIND ${PROJECT_NAME}_PACKAGES ${TRIBITS_PACKAGE} PACKAGE_IDX)
    LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    
    #checking if we have a relative path to the package's files. Since the exclude is a
    #regular expression any "../" will be interpretted as <any char><any char>/ which
    #would never match the package's actual directory. There isn't a direct way in cmake
    #to convert a relative path into an absolute path with string operations so as a way
    #of making sure that we get the correct path of the package we use a find_path for the
    #CMakeLists.txt file for the package. Since the package has to have this file to work
    #correctly it should be guaranteed to be there.
    STRING(REGEX MATCH "[.][.]/" IS_RELATIVE_PATH ${PACKAGE_DIR})
    IF("${IS_RELATIVE_PATH}" STREQUAL "")
      SET(CPACK_SOURCE_IGNORE_FILES ${LifeV_SOURCE_PATH}/${PACKAGE_DIR} ${CPACK_SOURCE_IGNORE_FILES})
    ELSE()
      FIND_PATH(ABSOLUTE_PATH CMakeLists.txt PATHS ${LifeV_SOURCE_PATH}/${PACKAGE_DIR} NO_DEFAULT_PATH)
      IF("${ABSOLUTE_PATH}" STREQUAL "ABSOLUTE_PATH-NOTFOUND")
        MESSAGE(AUTHOR_WARNING "Relative path found for disabled package ${TRIBITS_PACKAGE} but package was missing a CMakeLists.txt file. This disabled package will likely not be excluded from a source release")
      ENDIF()
      SET(CPACK_SOURCE_IGNORE_FILES ${ABSOLUTE_PATH} ${CPACK_SOURCE_IGNORE_FILES})
    ENDIF()
  ENDFOREACH()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Exclude files when building source packages")
    FOREACH(item ${CPACK_SOURCE_IGNORE_FILES})
      MESSAGE(${item})
    ENDFOREACH()
  ENDIF()

  # The CPACK_RESOURCE_FILE_[LICENSE|README] files must end in one of
  # .txt .rtf .html. Copy the pertinant file to the binary directory with
  # a .txt extension. This is only the case with the PackageMaker 
  # generator, but it doesn't hurt to do it for other generators as
  # well.
  MACRO(COPY_INSTALLER_RESOURCE _varname _source _destination)
    SET("${_varname}" "${_destination}")
    IF (EXISTS "${_destination}")
      FILE(REMOVE_RECURSE "${_destination}")
    ENDIF ()
    CONFIGURE_FILE(
      "${_source}" 
      "${_destination}" 
      COPYONLY)
  ENDMACRO()
  COPY_INSTALLER_RESOURCE(LifeV_README
    "${LifeV_SOURCE_DIR}/README.md"
    "${LifeV_BINARY_DIR}/README.md")
  COPY_INSTALLER_RESOURCE(LifeV_LICENSE
    "${LifeV_SOURCE_DIR}/LICENSE"
    "${LifeV_BINARY_DIR}/LICENSE.txt")

  SET(CPACK_PACKAGE_DESCRIPTION "LifeV provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "lifev-setup-${LifeV_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "LifeV ${LifeV_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "LifeV ${LifeV_VERSION}")
  SET(CPACK_PACKAGE_NAME "lifev")
  SET(CPACK_PACKAGE_VENDOR "EPFL - CMCS, MOX Polimi, INRIA REO, ESTIME, Sc. Comp. Emory Univ.")
  SET(CPACK_PACKAGE_VERSION "${LifeV_VERSION}")
  SET(CPACK_RESOURCE_FILE_README "${LifeV_README}")
  SET(CPACK_RESOURCE_FILE_LICENSE "${LifeV_LICENSE}")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "lifev-source-${LifeV_VERSION}")
  SET(CPACK_COMPONENTS_ALL ${LifeV_PACKAGES} Unspecified)
  
  set(ENABLED_FLAG ON)
  set(INCLUDE_EMPTY FALSE)
  TRIBITS_GET_ENABLED_LIST( LifeV_PACKAGES LifeV ${ENABLED_FLAG}
    ${INCLUDE_EMPTY} ENABLED_PACKAGES NUM_ENABLED)
  string(REPLACE " " ";" ENABLED_PACKAGES "${ENABLED_PACKAGES}")
  
  #message("ENABLED PACKAGES: ${ENABLED_PACKAGES} ${NUM_ENABLED}")
  FOREACH(PKG ${ENABLED_PACKAGES})
    IF(NOT "${${PKG}_LIB_REQUIRED_DEP_PACKAGES}" STREQUAL "")
        string(TOUPPER ${PKG} UPPER_PKG)
        #message("${UPPER_PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
        SET(CPACK_COMPONENT_${UPPER_PKG}_DEPENDS ${${PKG}_LIB_REQUIRED_DEP_PACKAGES})
    ENDIF()
    #message("${PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
  ENDFOREACH()

  INCLUDE(CPack)

ENDMACRO()
