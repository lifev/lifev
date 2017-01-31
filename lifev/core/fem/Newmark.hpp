#ifndef TIMEHANDLERNEWMARK_H
#define TIMEHANDLERNEWMARK_H 1

/*
 *  author: DAVIDE FORTI, davide.forti@epfl.ch
 *  Lightweighted class to Handle the time advancing scheme (based on Newmark).
 *
 */

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{

class Newmark
{

    typedef VectorEpetra vector_Type;

    typedef std::shared_ptr<vector_Type> vectorPtr_Type;

public:

    // empty constructor
    Newmark();

    // empty destructor
    ~Newmark();

    void initialize ( const vectorPtr_Type& state, const vectorPtr_Type& first_derivative, const vectorPtr_Type& second_derivative);

    void compute_csi( );

    void shift( const vectorPtr_Type& state );

    void restart( const vectorPtr_Type& state, const vectorPtr_Type& first_derivative, const vectorPtr_Type& second_derivative );

    // set the value of beta
    void set_beta( const Real beta ) { M_beta = beta; };

    // set the value of beta
    void set_gamma( const Real gamma ) { M_gamma = gamma; };

    // set the value of beta
    void set_timestep( const Real timestep ) { M_timeStep = timestep; };

    // set the value of beta
    Real get_beta( ) { return M_beta; };

    // set the value of beta
    Real get_gamma( ) { return M_gamma; };

    // set the value of beta
    Real get_timestep( ) { return M_timeStep; };

    vectorPtr_Type const& get_csi() const { return M_csi; };

    vectorPtr_Type const& old_first_derivative() const { return M_old_first_derivative; };

    vectorPtr_Type const& old_second_derivative() const { return M_old_second_derivative; };

private:

    // timestep
    Real M_timeStep;

    // Coefficients. For second order derivatives the method is unconditionally stable for M_beta >= 0.25 and M_gamma = 0.5
    Real M_beta;
    Real M_gamma;

    // Vector needed by Newmark
    vectorPtr_Type M_csi;

    // vectors for the current approximation
    vectorPtr_Type M_current_state;
    vectorPtr_Type M_current_first_derivative;
    vectorPtr_Type M_current_second_derivative;

    // vectors for the approximation at the previous time step
    vectorPtr_Type M_old_state;
    vectorPtr_Type M_old_first_derivative;
    vectorPtr_Type M_old_second_derivative;

};

} // end namespace LifeV

#endif
