#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>

namespace LifeV
{

TimeAndExtrapolationHandler::TimeAndExtrapolationHandler():
    M_BDForder(0),
    M_maximumExtrapolationOrder(0),
    M_states(),
    M_timeStep(0.0)
{}

TimeAndExtrapolationHandler::TimeAndExtrapolationHandler(const UInt orderBDF, const UInt maximumOrderExtrapolation):
    M_BDForder(orderBDF),
    M_maximumExtrapolationOrder(maximumOrderExtrapolation),
    M_states(),
    M_timeStep()
{}

TimeAndExtrapolationHandler::~TimeAndExtrapolationHandler()
{}

void
TimeAndExtrapolationHandler::setBDForder(const UInt order)
{
	M_BDForder  = order;
}

void
TimeAndExtrapolationHandler::setMaximumExtrapolationOrder(const UInt order)
{
	M_maximumExtrapolationOrder  = order;
}

void
TimeAndExtrapolationHandler::setTimeStep(const Real dt)
{
    M_timeStep = dt;
}

void
TimeAndExtrapolationHandler::initialize(const std::vector<vector_Type> InitialData)
{
	ASSERT( M_BDForder != 0, "Order of the BDF scheme has not been set, please use TimeAndExtrapolationHandler::setBDForder(const UInt order)");

	ASSERT( M_maximumExtrapolationOrder != 0, "Maximum order for extrapolation has not been set, please use TimeAndExtrapolationHandler::setMaximumExtrapolationOrder(const UInt order)");

	M_sizeStencil = ( M_maximumExtrapolationOrder > M_BDForder ) ? M_maximumExtrapolationOrder : M_BDForder;

    ASSERT( InitialData.size() == M_sizeStencil, "Wrong initial data dimension, it has to be of size equal max(M_BDForder, M_maximumExtrapolationOrder)");

    for ( int i = 0; i < M_sizeStencil; ++i )
        M_states.push_back(InitialData[i]);
}

void
TimeAndExtrapolationHandler::shift(const vector_Type newVector)
{
    switch (M_sizeStencil) {
        case 1: // size 1: we just replace the last entry by the new vector
            M_states[M_sizeStencil-1] = newVector;
            break;
        case 2: // size 2
            M_states[M_sizeStencil-2] = M_states[M_sizeStencil-1];
            M_states[M_sizeStencil-1] = newVector;
            break;
        case 3: // size 3
        	M_states[M_sizeStencil-3] = M_states[M_sizeStencil-2];
        	M_states[M_sizeStencil-2] = M_states[M_sizeStencil-1];
        	M_states[M_sizeStencil-1] = newVector;
        	break;
        default:
            break;
    }
}

void
TimeAndExtrapolationHandler::extrapolate(const UInt order, vector_Type& extrapolation)
{
    ASSERT( order <= M_maximumExtrapolationOrder, "Order of extrapolation is higher than the maximum order previously set");
    switch (order) {
        case 1:
            extrapolation = M_states[M_sizeStencil-1]; // u_star = u_n
            break;
        case 2:
            extrapolation = 2*M_states[M_sizeStencil-1] - M_states[M_sizeStencil-2]; // u_star = 2*u_n - u_{n-1}
            break;
        case 3:
        	extrapolation = 3*M_states[M_sizeStencil-1] - 3*M_states[M_sizeStencil-2] + M_states[M_sizeStencil-3]; // u_star = 3*u_n - 3*u_{n-1} + u_{n-2}
        	break;
        default:
            break;
    }
}

void
TimeAndExtrapolationHandler::rhsContribution(vector_Type& rhs_bdf)
{
	ASSERT( M_timeStep != 0, "Timestep has not been set, please use TimeAndExtrapolationHandler::setTimeStep(const Real dt) ");

    switch (M_BDForder) {
        case 1:
            rhs_bdf = 1/M_timeStep*M_states[M_sizeStencil-1]; // u_rhs = 1/dt*u_n
            break;
        case 2:
            rhs_bdf = 1/M_timeStep*(2*M_states[M_sizeStencil-1] - (1.0/2.0)*M_states[M_sizeStencil-2]); // u_rhs = 1/dt*(2*u_n - 0.5*u_{n-1})
            break;
        case 3:
        	rhs_bdf = 1/M_timeStep*(3*M_states[M_sizeStencil-1] - (3.0/2.0)*M_states[M_sizeStencil-2] + (1.0/3.0)*M_states[M_sizeStencil-3]); // u_rhs = 1/dt*(3*u_n - 3/2*u_{n-1} + 1/3*u_{n-2})
        	break;
        default:
            break;
    }
}

Real
TimeAndExtrapolationHandler::alpha()
{
    switch (M_BDForder) {
        case 1:
            return 1.0;
            break;
        case 2:
            return 1.5;
            break;
        case 3:
        	return 11.0/6.0;
        	break;
        default:
            break;
    }
    RETURN_UNDEFINED;
}

}
