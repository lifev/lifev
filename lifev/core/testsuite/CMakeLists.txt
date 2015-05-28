
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  adr_assembler
  array
  bdf
  fe_function
  fem
  filter
  hyperbolic
  interpolation
  linear_solver
  linear_solver_preconditioner
  matrix_epetra_structured_framework
  mesh
  repeated_mesh
  region_marker_id
  template_test
  verify_solution
  vector_container
 )

IF(TPL_ENABLE_HDF5)
    ADD_SUBDIRECTORIES(
        offline_partition_io
        )
ENDIF()

IF(LifeV_${PACKAGE_NAME}_ENABLE_SPIRIT_PARSER)
  ADD_SUBDIRECTORY(parser)
ENDIF()
