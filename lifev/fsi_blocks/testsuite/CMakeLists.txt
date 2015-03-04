INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
    fsi_tube
    fsi_challenge
    #fsi_challenge_steady
    fsi_couplingMaps
    fsi_couplingMapsVectorial
)