#ifndef TIMEHANDLER_H
#define TIMEHANDLER_H 1

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

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{

class TimeAndExtrapolationHandler
{

    typedef VectorEpetra vector_Type;

    typedef std::shared_ptr<vector_Type> vectorPtr_Type;

public:

    // empty constructor
    TimeAndExtrapolationHandler();

    // constructor taking as first input the order of the Time Advancing BDF scheme, as second the maximum order of extrapolation desired (maximum 3)
    TimeAndExtrapolationHandler(const UInt orderBDF, const UInt maximumOrderExtrapolation);

    // empty destructor
    ~TimeAndExtrapolationHandler();

    // setup of the class taking the order of the method
    void setBDForder(const UInt order);

    // set the maximum extrapolation order
    void setMaximumExtrapolationOrder(const UInt order);

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

    // extrapolation
    void extrapolate(const UInt order, vector_Type& extrapolation);

    // get the part from the discretization of the time derivative that goes to the right hand side
    void rhsContribution(vector_Type& rhs_bdf);

    // get the value of alpha that should go in front of the term u_(n+1) (see the paper cited at the beginning of the doc)
    Real alpha();

private:

    // order of the BDF scheme used
    UInt M_BDForder;

    // maximum order used for extrapolation
    UInt M_maximumExtrapolationOrder;

    // vector with the state variable at time n, n-1, .., n-(p-1)
    std::vector<vector_Type> M_states;

    // variable that contains the bigger value between M_BDForder and M_maximumExtrapolationOrder
    UInt M_sizeStencil;

    // timestep
    Real M_timeStep;

};

} // end namespace LifeV

#endif
