
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  multirate
  hyperbolic 
 )

IF(TPL_ENABLE_HDF5)
    ADD_SUBDIRECTORIES(
        partition_io
        )
ENDIF()

IF(LifeV_${PACKAGE_NAME}_ENABLE_SPIRIT_PARSER)
  ADD_SUBDIRECTORY(parser)
ENDIF()
