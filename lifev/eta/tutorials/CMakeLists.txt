
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  1_ETA_laplacian
  2_ETA_ADR
  3_ETA_rhs_and_value
  4_ETA_vectorial_laplacian
  5_ETA_debug_expressions
  6_ETA_functor
  7_ETA_blocks
  8_ETA_block_manip
  9_ETA_QRProvider
  10_ETA_QR_Advanced
  9_ETA_gradient_interpolation
  11_ETA_LaplacianPhiI
  12_ETA_LaplacianPhiJ
  13_ETA_LaplacianVector
)
