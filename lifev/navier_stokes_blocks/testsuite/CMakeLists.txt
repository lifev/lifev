
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
	semi_implicit_time_dependent_navier_stokes
	nonlinear_steady_navier_stokes
	nonlinear_time_dependent_navier_stokes
	Turek_Reynolds_20_semi_implicit
        Turek_Reynolds_20
)
