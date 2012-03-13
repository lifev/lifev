
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  adr_assembler
  array
  bdf
  exporter_vtk
  exporter_ensight_to_hdf5
  fe_function
  filter
  hyperbolic
  matrix_epetra_structured_framework
  mesh
  parser
  template_test
  vector_container
)
