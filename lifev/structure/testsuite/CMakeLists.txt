
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  time_advance_i
  time_advance_ii
  principalTensions
  structuralsolver
  )


IF(ENABLE_ANISOTROPIC_LAW)
  ADD_SUBDIRECTORIES(
    anisotropicLaw
    evaluateNodalETA
    )
ENDIF()
