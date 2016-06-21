
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  static_graph
  mt_assembly
  ADR_1D
  ADR_2D
  vectorial_ADR_2D
  interpolate_gradient_2D
  blocks_2D
  boundary_integrals
  #VMS
  #VMS_Cylinder
  repeated_mesh_2D
)
