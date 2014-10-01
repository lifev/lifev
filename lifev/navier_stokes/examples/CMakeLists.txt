
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  cavity_ns
  cavity_stokes
  cylinder
  TestCases
)

IF(LifeV_ENABLE_ETA)
ADD_SUBDIRECTORIES(
  #ETpseudo_fsi
  externalFluidFlow
  #oseen_assembler
  #resistanceBCs
  #computeWSS
  externalFluidFlow_FORTESTING
  #computeAeroCoeff
  Turek_benchmark
)
ENDIF()
