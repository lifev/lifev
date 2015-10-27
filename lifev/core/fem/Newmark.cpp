#include <lifev/core/fem/Newmark.hpp>

namespace LifeV
{

Newmark::Newmark():
	M_beta(0.0),
    M_gamma(0.0),
    M_timeStep(0.0)
{}

Newmark::~Newmark()
{}

void
Newmark::initialize(const vectorPtr_Type& state, const vectorPtr_Type& first_derivative, const vectorPtr_Type& second_derivative)
{

	// Previous timestep - given
	M_old_state.reset( new vector_Type ( *state, Unique ) );
	M_old_first_derivative.reset( new vector_Type ( *first_derivative, Unique ) );
	M_old_second_derivative.reset( new vector_Type ( *second_derivative, Unique ) );

	M_old_state->zero();
	M_old_first_derivative->zero();
	M_old_second_derivative->zero();

	// Current timestep - initialized to zero
	M_current_state.reset( new vector_Type ( *state, Unique ) );
	M_current_first_derivative.reset( new vector_Type ( *first_derivative, Unique ) );
	M_current_second_derivative.reset( new vector_Type ( *second_derivative, Unique ) );

	M_current_state->zero();
	M_current_first_derivative->zero();
	M_current_second_derivative->zero();

	// Instantiate csi
	M_csi.reset( new vector_Type ( *state, Unique ) );
}

void
Newmark::compute_csi( )
{
	ASSERT( M_beta != 0.0, "Beta coefficient has not been set in Newmark");

	ASSERT( M_timeStep != 0.0, "Beta coefficient has not been set in Newmark");

	M_csi->zero();

	*M_csi = ( ( ( 1.0/(M_timeStep*M_timeStep*M_beta) ) * ( *M_old_state ) ) +
			   ( ( 1.0/(M_timeStep*M_beta) ) * ( *M_old_first_derivative ) ) +
			   ( ( (1.0 - 2.0 * M_beta)/(2.0 * M_beta) ) * (*M_old_second_derivative) ) );
 }


void
Newmark::shift( const vectorPtr_Type& state )
{
	M_current_state->zero();

	*M_current_state = *state;

	ASSERT( M_beta != 0.0, "Beta coefficient has not been set in Newmark");

	ASSERT( M_gamma != 0.0, "Gamma coefficient has not been set in Newmark");

	ASSERT( M_timeStep != 0.0, "Timestep has not been set in Newmark");

	M_current_first_derivative->zero();

	M_current_second_derivative->zero();

	*M_current_second_derivative = ( ( 1.0/(M_timeStep*M_timeStep*M_beta) * (*M_current_state) ) - *M_csi );

	*M_current_first_derivative = ( (*M_old_first_derivative) +
			                        ( M_timeStep * M_gamma * (*M_current_second_derivative) ) +
			                        ( M_timeStep * ( 1.0 - M_gamma ) * (*M_old_second_derivative) ) );

	// Shift

	*M_old_state = *M_current_state;

	*M_old_first_derivative = *M_current_first_derivative;

	*M_old_second_derivative = *M_current_second_derivative;
}

}
