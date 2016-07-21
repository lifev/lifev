#ifndef BDFSECONDORDERDERIVATIVE_H
#define BDFSECONDORDERDERIVATIVE_H 1

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

class BDFSecondOrderDerivative
{

    typedef VectorEpetra vector_Type;

    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

public:

    // empty constructor
    BDFSecondOrderDerivative();

    // constructor taking as input the order of the Time Advancing BDF scheme
    BDFSecondOrderDerivative(const UInt orderBDF);

    // empty destructor
    ~BDFSecondOrderDerivative();

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

    // used in the coupling of the velocities for FSI
    void first_der_old_dts( vector_Type& vec_old_timesteps);

    // used in the computation of the residual
    void second_der_old_dts( vector_Type& vec_old_timesteps);

    // coefficient in front of d_{n+1}
    Real massCoefficient();

    // coefficient first derivative in front of d_{n+1}
    Real coefficientFirstDerivative();

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

} // end namespace LifeV

#endif
