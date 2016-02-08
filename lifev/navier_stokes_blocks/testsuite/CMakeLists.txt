
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
	semi_implicit_time_dependent_navier_stokes
	nonlinear_steady_navier_stokes
	nonlinear_time_dependent_navier_stokes
	Turek_Reynolds_20_semi_implicit
	Turek_Reynolds_100_semi_implicit
    Turek_Reynolds_20
    non_conforming_fluids
        evaluate_at_quadrature_points
)
