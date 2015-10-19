INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
    fsi_tube
    fsi_restart
    Bypass
    #fsi_bypass
    #fsi_challenge
    #fsi_challenge_steady
    #fsi_couplingMapsVectorial
)