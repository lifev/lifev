#include <lifev/core/fem/BDFSecondOrderDerivative.hpp>

namespace LifeV
{

BDFSecondOrderDerivative::BDFSecondOrderDerivative():
    M_BDForder(0),
    M_states(),
    M_timeStep(0.0)
{}

BDFSecondOrderDerivative::BDFSecondOrderDerivative(const UInt orderBDF):
    M_BDForder(orderBDF),
    M_states(),
    M_timeStep(0.0)
{}

BDFSecondOrderDerivative::~BDFSecondOrderDerivative()
{}

void
BDFSecondOrderDerivative::setBDForder(const UInt order)
{
	M_BDForder  = order;
}

void
BDFSecondOrderDerivative::setTimeStep(const Real dt)
{
    M_timeStep = dt;
}

void
BDFSecondOrderDerivative::initialize(const std::vector<vector_Type> InitialData)
{
	ASSERT( M_BDForder != 0, "Order of the BDF scheme has not been set, please use BDFSecondOrderDerivative::setBDForder(const UInt order)");

	M_sizeStencil = M_BDForder + 1;

    for ( int i = 0; i < M_sizeStencil; ++i )
        M_states.push_back(InitialData[i]);
}

void
BDFSecondOrderDerivative::shift(const vector_Type newVector)
{
    switch (M_BDForder) {
        case 1: // size 1: we just replace the last entry by the new vector
        	M_states[M_sizeStencil-2] = M_states[M_sizeStencil-1];
            M_states[M_sizeStencil-1] = newVector;
            break;
        case 2: // size 2
        	M_states[M_sizeStencil-3] = M_states[M_sizeStencil-2];
            M_states[M_sizeStencil-2] = M_states[M_sizeStencil-1];
            M_states[M_sizeStencil-1] = newVector;
            break;
        default:
            break;
    }
}

Real
BDFSecondOrderDerivative::massCoefficient() // Note: there is no division by dt^2!!!!
{
    switch (M_BDForder) {
        case 1:
            return 1.0;
            break;
        case 2:
            return 2.0;
            break;
        default:
            break;
    }
    RETURN_UNDEFINED;
}

Real
BDFSecondOrderDerivative::coefficientFirstDerivative() // Note: there is no division by dt!!!!
{
    switch (M_BDForder) {
        case 1:
            return 1.0;
            break;
        case 2:
            return 1.5;
            break;
        default:
            break;
    }
    RETURN_UNDEFINED;
}

void
BDFSecondOrderDerivative::first_der_old_dts(vector_Type& vec_old_timesteps) // Note: there is no division by dt!!!!
{
	ASSERT( M_timeStep != 0, "Timestep has not been set, please use BDFSecondOrderDerivative::setTimeStep(const Real dt) ");

    switch (M_BDForder) {
        case 1:
        	vec_old_timesteps = -1.0 * M_states[M_sizeStencil-1]; // vec = -1.0 * d_n
            break;
        case 2:
        	vec_old_timesteps = -2.0 * M_states[M_sizeStencil-1] + 0.5 * M_states[M_sizeStencil-2]; // vec = -2.0 * d_n + 0.5 * d_{n-1}
            break;
        default:
            break;
    }
}

void
BDFSecondOrderDerivative::second_der_old_dts(vector_Type& vec_old_timesteps) // Note: there is no division by dt^2!!!!
{
	ASSERT( M_timeStep != 0, "Timestep has not been set, please use BDFSecondOrderDerivative::setTimeStep(const Real dt) ");

    switch (M_BDForder) {
        case 1:
        	vec_old_timesteps = -2.0 * M_states[M_sizeStencil-1] + 1.0 * M_states[M_sizeStencil-2]; // vec = -2 * d_n + 1.0 * d_{n-1}
            break;
        case 2:
        	vec_old_timesteps = -5.0 * M_states[M_sizeStencil-1] + 4.0 * M_states[M_sizeStencil-2] - 1.0 * M_states[M_sizeStencil-3]; // vec = -5 * d_n + 4.0 * d_{n-1} - 1.0 * d_{n-2}
            break;
        default:
            break;
    }
}

}
