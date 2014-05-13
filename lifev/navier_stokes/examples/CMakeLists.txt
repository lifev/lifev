
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  resistanceBCs
  oseen_assembler
  cavity_ns
  cavity_stokes
  cylinder
  TestCases
)

IF(LifeV_ENABLE_ETA)
ADD_SUBDIRECTORIES(
  ETpseudo_fsi
)
ENDIF()
