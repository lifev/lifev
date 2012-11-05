
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
  linear_solver
  matrix_epetra_structured_framework
  mesh
  repeated_mesh
  region_marker_id
  partition_io
  template_test
  vector_container
)

IF(LifeV_${PACKAGE_NAME}_ENABLE_SPIRIT_PARSER)
  ADD_SUBDIRECTORY(parser)
ENDIF()
