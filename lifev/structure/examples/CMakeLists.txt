
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
  example_bodyForces
)

IF(EXPORTVECTORS_DEBUG)
  ADD_SUBDIRECTORIES(
  example_checkigForceTerm
)
ENDIF()
