
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  example_partitionMesh
  example_computationNorm
  example_checkingFibersDirection
  example_computePrincipalTensions
  example_bodyForces
  example_anisotropicTraction
  example_evaluateActivationMultimechanism
  example_evaluatingScalarVectorialTensorialQuantitiesUsingETA
)

IF(EXPORTVECTORS_DEBUG)
  ADD_SUBDIRECTORIES(
  example_checkigForceTerm
)
ENDIF()

IF(COLORING_MESH)
  ADD_SUBDIRECTORIES(
  example_creatingDamagedZone
)
ENDIF()

IF(COMPUTATION_JACOBIAN)
  ADD_SUBDIRECTORIES(
    example_computingJacobian
)
ENDIF()
