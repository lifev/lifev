#ifndef TIMEHANDLERQUADPTS_H
#define TIMEHANDLERQUADPTS_H 1

/*
 *  author: DAVIDE FORTI, davide.forti@epfl.ch
 *  Lightweighted class to Handle the time advancing scheme (based on BDF approximation of the time derivative and for the extrapolation).
 *  It can be used for first derivatives.
 *
 *  The formulas implemented in the methods extrapolate and rhsContribution are taken from
 *  "Algebraic fractional-step schemes with spectral methods for the incompressible Navier–Stokes equations"
 *  by Paola Gervasio, Fausto Saleri and Alessandro Veneziani. (see pag. 3 for details)
 */

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

template <UInt DIM>
class TimeAndExtrapolationHandlerQuadPts
{

    typedef std::vector<std::vector<VectorSmall<DIM>>> vector_Type;

    typedef std::shared_ptr<vector_Type> vectorPtr_Type;

public:

    // empty constructor
    TimeAndExtrapolationHandlerQuadPts();

    // constructor taking as first input the order of the Time Advancing BDF scheme, as second the maximum order of extrapolation desired (maximum 3)
    TimeAndExtrapolationHandlerQuadPts(const UInt orderBDF);

    // empty destructor
    ~TimeAndExtrapolationHandlerQuadPts();

    // setup of the class taking the order of the method
    void setBDForder(const UInt order);

    // initialize the time handler class with initial data, need a vector of M_order vectorEpetra
    void initialize(const std::vector<vector_Type> InitialData);

    // shift - to be used when a timestep is solved
    void shift(const vector_Type newVector);

    // getter for the state
    std::vector<vector_Type> state()
    {
        return M_states;
    }

    // set the timestep
    void setTimeStep(const Real dt);

    // get the part from the discretization of the time derivative that goes to the right hand side
    void rhsContribution(vector_Type& rhs_bdf);

    // get the value of alpha that should go in front of the term u_(n+1) (see the paper cited at the beginning of the doc)
    Real alpha();
    
    // extrapolation
    void extrapolate(vector_Type& extrapolation);

private:

    // order of the BDF scheme used
    UInt M_BDForder;

    // vector with the state variable at time n, n-1, .., n-(p-1)
    std::vector<vector_Type> M_states;

    // variable that contains the bigger value between M_BDForder and M_maximumExtrapolationOrder
    UInt M_sizeStencil;

    // timestep
    Real M_timeStep;

};

template <UInt DIM>
TimeAndExtrapolationHandlerQuadPts<DIM>::TimeAndExtrapolationHandlerQuadPts():
M_BDForder(0),
M_states(),
M_timeStep(0.0)
{}

template <UInt DIM>
TimeAndExtrapolationHandlerQuadPts<DIM>::TimeAndExtrapolationHandlerQuadPts(const UInt orderBDF):
M_BDForder(orderBDF),
M_states(),
M_timeStep()
{}

template <UInt DIM>
TimeAndExtrapolationHandlerQuadPts<DIM>::~TimeAndExtrapolationHandlerQuadPts()
{}

template <UInt DIM>
void
TimeAndExtrapolationHandlerQuadPts<DIM>::setBDForder(const UInt order)
{
    M_BDForder  = order;
}

template <UInt DIM>
void
TimeAndExtrapolationHandlerQuadPts<DIM>::setTimeStep(const Real dt)
{
    M_timeStep = dt;
}

template <UInt DIM>
void
TimeAndExtrapolationHandlerQuadPts<DIM>::initialize(const std::vector<vector_Type> InitialData)
{
    ASSERT( M_BDForder != 0, "Order of the BDF scheme has not been set, please use TimeAndExtrapolationHandler::setBDForder(const UInt order)");
    
    M_sizeStencil = M_BDForder;
    
    ASSERT( InitialData.size() == M_sizeStencil, "Wrong initial data dimension, it has to be of size equal max(M_BDForder, M_maximumExtrapolationOrder)");
    
    for ( int i = 0; i < M_sizeStencil; ++i )
        M_states.push_back(InitialData[i]);
}

template <UInt DIM>
void
TimeAndExtrapolationHandlerQuadPts<DIM>::shift(const vector_Type newVector)
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

template <UInt DIM>
void
TimeAndExtrapolationHandlerQuadPts<DIM>::rhsContribution(vector_Type& rhs_bdf)
{
    ASSERT( M_timeStep != 0, "Timestep has not been set, please use TimeAndExtrapolationHandler::setTimeStep(const Real dt) ");
    
    // nota: M_states.size() -> orderBDF
    //       M_states[0].size() -> numElements per process
    //       M_states[0][0].size() -> numQuadPts per element
    //       DIM is the size of the vector small
    
    switch (M_BDForder) {
            case 1:
            for (int i = 0 ; i < M_states[0].size() ; ++i ) // loop elements
            {
                for (int j = 0 ; j < M_states[0][0].size(); ++j ) // loop quadrature points
                {
                    for (int k = 0 ; k < DIM; ++k ) // loop quadrature points
                    {
                        rhs_bdf[i][j](k) = 1/M_timeStep*M_states[M_sizeStencil-1][i][j](k); // u_rhs = 1/dt*u_n
                    }
                }
            }
            break;
            case 2:
            for (int i = 0 ; i < M_states[0].size() ; ++i ) // loop elements
            {
                for (int j = 0 ; j < M_states[0][0].size(); ++j ) // loop quadrature points
                {
                    for (int k = 0 ; k < DIM; ++k ) // loop quadrature points
                    {
                        rhs_bdf[i][j](k) = 1/M_timeStep*(2*M_states[M_sizeStencil-1][i][j](k) - (1.0/2.0)*M_states[M_sizeStencil-2][i][j](k));
                        // u_rhs = 1/dt*(2*u_n - 0.5*u_{n-1})
                    }
                }
            }
            break;
            case 3:
            for (int i = 0 ; i < M_states[0].size() ; ++i ) // loop elements
            {
                for (int j = 0 ; j < M_states[0][0].size(); ++j ) // loop quadrature points
                {
                    for (int k = 0 ; k < DIM; ++k ) // loop quadrature points
                    {
                        rhs_bdf[i][j](k) = 1/M_timeStep*(3*M_states[M_sizeStencil-1][i][j](k) - (3.0/2.0)*M_states[M_sizeStencil-2][i][j](k) +
                                                         (1.0/3.0)*M_states[M_sizeStencil-3][i][j](k) );
                        // u_rhs = 1/dt*(3*u_n - 3/2*u_{n-1} + 1/3*u_{n-2})
                    }
                }
            }
            break;
        default:
            break;
    }
}

template <UInt DIM>
void
TimeAndExtrapolationHandlerQuadPts<DIM>::extrapolate(vector_Type& extrapolation)
{
    ASSERT( M_timeStep != 0, "Timestep has not been set, please use TimeAndExtrapolationHandler::setTimeStep(const Real dt) ");
    
    switch (M_BDForder) {
        case 1:
            for (int i = 0 ; i < M_states[0].size() ; ++i ) // loop elements
            {
                for (int j = 0 ; j < M_states[0][0].size(); ++j ) // loop quadrature points
                {
                    for (int k = 0 ; k < DIM; ++k ) // loop quadrature points
                    {
                        extrapolation[i][j](k) = M_states[M_sizeStencil-1][i][j](k); // u_star = u_n
                    }
                }
            }
            break;
        case 2:
            for (int i = 0 ; i < M_states[0].size() ; ++i ) // loop elements
            {
                for (int j = 0 ; j < M_states[0][0].size(); ++j ) // loop quadrature points
                {
                    for (int k = 0 ; k < DIM; ++k ) // loop quadrature points
                    {
                        extrapolation[i][j](k) = 2*M_states[M_sizeStencil-1][i][j](k) - M_states[M_sizeStencil-2][i][j](k);
                        // u_star = 2*u_n - u_{n-1}
                    }
                }
            }
            break;
        case 3:
            for (int i = 0 ; i < M_states[0].size() ; ++i ) // loop elements
            {
                for (int j = 0 ; j < M_states[0][0].size(); ++j ) // loop quadrature points
                {
                    for (int k = 0 ; k < DIM; ++k ) // loop quadrature points
                    {
                        extrapolation[i][j](k) = 3.0*M_states[M_sizeStencil-1][i][j](k) - 3.0*M_states[M_sizeStencil-2][i][j](k) +
                                                     M_states[M_sizeStencil-3][i][j](k);
                        // u_star = 3*u_n - 3*u_{n-1} + u_{n-2}
                    }
                }
            }
            break;
        default:
            break;
    }
}
    
template <UInt DIM>
Real
TimeAndExtrapolationHandlerQuadPts<DIM>::alpha()
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
}
    
} // end namespace LifeV

#endif
