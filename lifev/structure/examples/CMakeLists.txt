
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  example_InflationExtension
  example_CookMembrane
  example_tractionWithSymmetry
  example_computePrincipalTensions
  example_partitionMesh
  example_computationNorm
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


IF(ENABLE_ANISOTROPIC_LAW)
  ADD_SUBDIRECTORIES(
    example_bodyForces
    example_anisotropicTraction
    example_multimechanismTraction
    example_evaluateActivationMultimechanism
    example_evaluatingScalarVectorialTensorialQuantitiesUsingETA
)
ENDIF()
