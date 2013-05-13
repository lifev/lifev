
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  example_InflationExtension
  example_CookMembrane
  example_creatingDamagedZone
  example_tractionWithSymmetry
  example_principalTensionsInflationExtensions
  example_partitionMesh
  example_computationNorm
)

IF(EXPORTVECTORS_DEBUG)
  ADD_SUBDIRECTORIES(
  example_checkigForceTerm
)
ENDIF()

IF(ENABLE_ANISOTROPIC_LAW)
  ADD_SUBDIRECTORIES(
    example_bodyForces
)
ENDIF()
