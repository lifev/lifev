//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief   File containing a class for an easy handling of different order time
 *           discretizations/extrapolations BDF based
 *
 *  @author Alessandro Veneziani <ale@mathcs.emory.edu>
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *  @author C. Winkelmann
 *  @author Mauro Perego
 *  @author Umberto Villa <uvilla@emory.edu>
 *
 *  @contributor Claudia Colciago
 *  @mantainer Claudia Colciago
 */


#ifndef _BDF_VARIABLE_TIMESTEP_H
#define _BDF_VARIABLE_TIMESTEP_H

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
const UInt bdfMaxOrder = 6;
typedef boost::numeric::ublas::vector<Real> ScalarVector;


//!TimeAdvanceBDFVariableStep - Backward differencing formula time discretization with non-uniform time step
/*!

    A differential equation of the form

    \f$ M u' = A u + f \f$

    is discretized in time as

    \f$ M p'(t_{k+1}) = A u_{k+1} + f_{k+1} \f$

    where p denotes the polynomial of order n in t that interpolates
    \f$ (t_i,u_i) \f$ for i = k-n+1,...,k+1.

    The approximative time derivative \f$ p'(t_{k+1}) \f$ is a linear
    combination of state vectors \f$ u_i \f$:

    \f$ p'(t_{k+1}) = \frac{1}{\Delta t} (\alpha_0 u_{k+1} - \sum_{i=0}^n \alpha_i u_{k+1-i} )\f$

    Thus we have

    \f$ \frac{\alpha_0}{\Delta t} M u_{k+1} = A u_{k+1} + f + M \bar{p} \f$

    with

    \f$ \bar{p} = \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i u_{k+1-i} \f$

    This class stores the n last state vectors in order to be able to
    calculate \f$ \bar{p} \f$. It also provides \f$ \alpha_i \f$
    and can extrapolate the the new state from the n last states with a
    polynomial of order n-1:

    \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$

    The class allows to change the time-steps, recomputing the coefficients.
    to change \f$ \Delta t \f$ at each time-step,
    use shift_right(uCurrent, timeStepNew) instead of shift_right(uCurrent).

    Other methods that could be helpful while adapting in time are:
    set_deltat, store_unk and restore_unk
*/
template< typename FEVectorType = VectorEpetra >
class TimeAdvanceBDFVariableStep
{
public:

    //! @name Public Types
    //@{
    typedef ScalarVector                         container_Type;
    typedef FEVectorType                         feVector_Type;
    typedef boost::shared_ptr< feVector_Type >   feVectorPtr_Type;
    typedef std::vector< feVectorPtr_Type >      feVectorPtrContainer_Type;
    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    TimeAdvanceBDFVariableStep();

    //!Destructor
    ~TimeAdvanceBDFVariableStep();
    //@}

    //! @name Methods
    //@{

    //!Setup the TimeAdvanceBDFVariableStep with order n
    /*!
        @param order order of the BDF
     */
    void setup ( const UInt order );

    //! Initialize all the entries of the unknown vector to be derived with the
    //! vector u0 (duplicated if startup=0)
    /*!
        If startup = 0, it initializes all the entries of the unknown vectors with a given function
        The array of initial conditions needed by the selected BDF is
        initialized as follows: _unk=[ u0(t0), u0(t0-dt), u0(t0-2*dt), ...]
        When startUp = true, only M_unknown[0] is initialized (with u0).
        At the first step, class Bfd computes the coefficients of BDF1,
        at the second step, the coefficient of BDF2, and so on,
        until it reaches the wanted method.
        Second order of accuracy is obtained using this procedure to startup BDF2.
        Using the startup procedure with BDF3, the accuracy order is limited to the second order;
        anyway it is better to use the startup procedure instead of taking u_{-2} = u_{-1} = u_0.

        @param u0 initial condition vector
        @param timeStep time step
        @param startup
     */
    void setInitialCondition ( feVector_Type u0, Real const timeStep, bool startup = 0 );

    //! Initialize all the entries of the unknown vector to be derived with a
    //! set of vectors u0s
    /*!
        @param u0s initial condition vectors
        @param timeStep time step
     */
    void setInitialCondition ( std::vector<feVector_Type> u0s, Real const timeStep );

    //!Initialize all the entries of the unknonwn vectors with a given function
    /*!
        The array of initial conditions needed by the selected BDF is
        initialized as follows: M_unknown=[ u0Function(t0),
        u0Function(t0-dt), u0Function(t0-2*dt), ...]
        For the space dependence of the initial conditions we need informations
        on:
        @param u0Function the function we want to interpolate
        @param u0Vector corrected mapped vector, it contains the u0Function(t0) in output
        @param feSpace finite element space
        @param t0 initial time (t0) and the
        @param timeStep time step (for solutions before the initial instant)

        Based on the NavierStokesHandler::initialize by M. Fernandez
     */
    template<typename FunctionType, typename FESpaceType>
    void setInitialCondition ( const FunctionType& u0Function, feVector_Type& u0Vector, FESpaceType& feSpace,
                               Real t0, Real timeStep );

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula divided by dt
    /*!
        @return a feVector_Type, by default feVector_Type = VectorEpetra
     */

    feVector_Type rhsContribution() const;

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula divided by timeStep (backward compatibility version, will be discontinued)
    /*!
        @param timeStep time step, 1 by default
        @return a feVector_Type, by default feVector_Type = VectorEpetra
     */
    feVector_Type rhsContributionTimeStep ( Real timeStep = 1 ) const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    feVector_Type extrapolateSolution() const;

    //! Update the vectors of the previous time steps by shifting on the right
    //! the old values.
    /*!
        @param uCurrent current (new) value of the state vector
     */
    void shiftRight ( feVector_Type const& uCurrent );

    //! Update the vectors of the previous time steps by shifting on the right
    //! the old values.
    /*!
        Set the the new time step and shift right the old time steps.
        @param uCurrent current (new) value of the state vector
        @param timeStepNew new value of the time step
     */
    void shiftRight ( feVector_Type const& uCurrent, Real timeStepNew );

    //! Save the current vector M_unknowns and the current vector M_timeStep
    void storeSolution();

    //! Restore the vector M_unknowns and the vector M_timeStep with the ones saved with store_unk()
    void restoreSolution();

    //! It is equivalent to do : storeSolution() + setTimeStep()
    /*!
        @param timeStep new time step
     */
    void restoreSolution ( Real timeStep );

    //! Show informations about the BDF
    void showMe() const;

    //@}

    //! @name Set Methods
    //@{
    //! Replace the current time-step with timeStep and computes the coefficients a_i and beta_i as functions of M_timeStep.
    /*!
        @param timeStep time step
     */
    void setTimeStep ( Real timeStep );

    //@}


    //! @name Get Methods
    //@{
    //! Return the i-th coefficient of the time derivative alpha_i
    /*!
        @param i index of the coefficient
     */
    const Real& coefficientDerivative ( UInt i ) const;

    //! Return the i-th coefficient of the time derivative alpha_i divided by dt
    /*!
        @param i index of the coefficient
     */
    Real coefficientDerivativeOverTimeStep ( UInt i ) const;

    //! Return the i-th coefficient of the time extrapolation beta_i
    /*!
        @param i index of the coefficient
     */
    const Real& coefficientExtrapolation ( UInt i ) const;

    //! Return the vector of the time steps, ordered starting from the most recent one.
    const container_Type& timeStepVector() const
    {
        return M_timeStep;
    }

    //! Return a vector with the last n state vectors
    const feVectorPtrContainer_Type& stateVector() const
    {
        return M_unknowns;
    }

    const feVectorPtrContainer_Type& stateVectorBack() const
    {
        return M_unknownsBack;
    }

    //@}

private:

    //! Private methods
    //@{
    //! Computes the coefficients a_i and beta_i as functions of _M_delta_t.
    /*!
       Arbitrary order, variable time step BDF coefficients.
       Reference: Comincioli pag. 618-619
       Matlab implementation:
          function [alpha0, alpha, beta] = compute_coef(dts, order)
          assert( (length(dts) == order) );
          rho = dts(1)./cumsum(dts);

          alpha0 = sum(rho);
          alpha = zeros(order,1);
          beta = zeros(order,1);

          for j = 1:order
               ind = setdiff([1:order],j);
               beta(j) = 1/prod(1-rho(ind)/rho(j));
               alpha(j) = rho(j)*beta(j);
          end

          end
     */
    void computeCoefficient();

    //@}

    //! Order of the BDF derivative/extrapolation: the time-derivative
    //! coefficients vector has size n+1, the extrapolation vector has size n
    UInt M_order;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    container_Type M_alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    container_Type M_beta;

    //! container_Type \f$ \delta_t \f$ of time intervals
    container_Type M_timeStep;
    container_Type M_timeStepBack;

    //! Last n state vectors
    feVectorPtrContainer_Type M_unknowns;
    feVectorPtrContainer_Type M_unknownsBack;

};

// ===================================================
// Constructors & Destructor
// ===================================================
template<typename FEVectorType>
TimeAdvanceBDFVariableStep<FEVectorType>::TimeAdvanceBDFVariableStep()
    :
    M_order ( 0 ),
    M_alpha ( M_order + 1 ),
    M_beta ( M_order ),
    M_timeStep ( ScalarVector ( M_order, 1. ) )
{}


template<typename FEVectorType>
TimeAdvanceBDFVariableStep<FEVectorType>::~TimeAdvanceBDFVariableStep()
{}

// ===================================================
// Methods
// ===================================================
template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::setup ( const UInt order )
{
    M_order = order;
    M_alpha.resize ( M_order + 1 );
    M_beta.resize ( M_order );
    M_timeStep = ScalarVector ( M_order, 1. );

    if ( order <= 0 || order > bdfMaxOrder )
    {
        std::ostringstream errorMessage;
        errorMessage << "Error: wrong BDF order\n"
                     << " you want to use BDF order " << order << "\n"
                     << " we support BDF order from 1 to " << bdfMaxOrder << "\n";
        throw std::invalid_argument ( errorMessage.str() );
    }
    computeCoefficient();
    M_unknowns.reserve ( order );
}


template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::setInitialCondition ( feVector_Type u0, Real const timeStep, bool startup )
{
    M_unknowns.resize ( 0 );
    for ( UInt i = 0; i < M_order; i++ )
    {
        feVectorPtr_Type feVectorPtr ( new feVector_Type ( u0 ) );
        M_unknowns.push_back ( feVectorPtr );
        M_timeStep[ i ] = timeStep;
    }

    if ( startup )
    {
        for ( UInt  i = 1; i < M_order; i++ )
        {
            M_timeStep[ i ] = 1e20;
        }
    }

    computeCoefficient();

    return ;
}


template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::setInitialCondition ( std::vector<feVector_Type> u0s, Real const timeStep )
{
    UInt n0 ( u0s.size() );
    // Check if u0s has the right dimensions
    ASSERT ( n0 >= M_order, "Initial data are not enough for the selected BDF" );

    M_unknowns.resize (0);
    for ( UInt i = 0; i < M_order; i++ )
    {
        feVectorPtr_Type tmp ( new feVector_Type ( u0s[ i ] ) );
        M_unknowns.push_back ( tmp );
        M_timeStep[ i ] = timeStep;
    }

    computeCoefficient();

    // if n0>n, only the first n inital data will be considered
    if ( n0 > M_order )
    {
        std::cout << "The initial data set is larger than needed by the BDF."
                  << std::endl;
        std::cout << "Only the first " << M_order << " data will be considered. "
                  << std::endl;
    }

    return ;
}


template<typename FEVectorType>
template<typename FunctionType, typename FESpaceType>
void TimeAdvanceBDFVariableStep<FEVectorType>::setInitialCondition ( const FunctionType& u0Function, feVector_Type& u0Vector,
                                                                     FESpaceType& feSpace, Real t0, Real timeStep )
{
    M_unknowns.resize ( 0 );

    for ( UInt i = 0 ; i < M_order; ++i )
    {
        feVectorPtr_Type tmp ( new feVector_Type ( u0Vector ) );
        M_unknowns.push_back ( tmp );
        feSpace.interpolate ( static_cast<typename FESpaceType::function_Type> ( u0Function ), *M_unknowns[ i ], t0 - i * timeStep );
        M_timeStep[ i ] = timeStep;
    }
    u0Vector = *M_unknowns[ 0 ];

    computeCoefficient();

    return ;
}


template<typename FEVectorType>
typename TimeAdvanceBDFVariableStep<FEVectorType>::feVector_Type
TimeAdvanceBDFVariableStep<FEVectorType>::rhsContribution() const
{
    feVector_Type uTimeDerivative ( *M_unknowns[ 0 ] );
    uTimeDerivative *= M_alpha[ 1 ] / M_timeStep[ 0 ];

    for ( UInt i = 1; i < M_order; ++i )
    {
        uTimeDerivative += ( M_alpha[ i + 1 ] / M_timeStep[ 0 ] ) * *M_unknowns[ i ];
    }

    return uTimeDerivative;
}


template<typename FEVectorType>
typename TimeAdvanceBDFVariableStep<FEVectorType>::feVector_Type
TimeAdvanceBDFVariableStep<FEVectorType>::rhsContributionTimeStep ( Real timeStep ) const
{
    feVector_Type uTimeDerivative ( *M_unknowns[ 0 ] );
    uTimeDerivative *= M_alpha[ 1 ] / timeStep;

    for ( UInt i = 1; i < M_order; ++i )
    {
        uTimeDerivative += ( M_alpha[ i + 1 ] / timeStep ) * *M_unknowns[ i ];
    }

    return uTimeDerivative;
}


template<typename FEVectorType>
typename TimeAdvanceBDFVariableStep<FEVectorType>::feVector_Type
TimeAdvanceBDFVariableStep<FEVectorType>::extrapolateSolution() const
{
    feVector_Type uExtrapolated ( *M_unknowns[ 0 ] );
    uExtrapolated *= M_beta[ 0 ];

    for ( UInt i = 1; i < M_order; ++i )
    {
        uExtrapolated += M_beta[ i ] * *M_unknowns[ i ];
    }

    return uExtrapolated;
}


template<typename FEVectorType>
void
TimeAdvanceBDFVariableStep<FEVectorType>::shiftRight ( feVector_Type const& uCurrent )
{
    typedef typename feVectorPtrContainer_Type::iterator bdfContainerIterator_Type;
    bdfContainerIterator_Type it = M_unknowns.end() - 1;
    bdfContainerIterator_Type itMinusOne = M_unknowns.end() - 1;
    bdfContainerIterator_Type itBegin = M_unknowns.begin();

    for ( ; it != itBegin; --it )
    {
        itMinusOne--;
        *it = *itMinusOne;
    }
    feVectorPtr_Type tmp ( new feVector_Type ( uCurrent ) );
    *itBegin = tmp;

    for ( UInt i = M_order - 1; i > 0; i-- )
    {
        M_timeStep[ i ] = M_timeStep[ i - 1 ];
    }

    computeCoefficient();
}

template<typename FEVectorType>
void
TimeAdvanceBDFVariableStep<FEVectorType>::shiftRight ( feVector_Type const& uCurrent, Real timeStepNew )
{

    typedef typename feVectorPtrContainer_Type::iterator bdfContainerIterator_Type;
    bdfContainerIterator_Type it = M_unknowns.end() - 1;
    bdfContainerIterator_Type itMinusOne = M_unknowns.end() - 1;
    bdfContainerIterator_Type itBegin = M_unknowns.begin();

    for ( ; it != itBegin; --it )
    {
        itMinusOne--;
        *it = *itMinusOne;
    }

    feVectorPtr_Type tmp ( new feVector_Type ( uCurrent ) );
    *itBegin = tmp;

    for ( UInt i = M_order - 1; i > 0; i-- )
    {
        M_timeStep[ i ] = M_timeStep[ i - 1 ];
    }
    M_timeStep[ 0 ] = timeStepNew;

    computeCoefficient();
}

template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::storeSolution()
{
    typedef typename feVectorPtrContainer_Type::iterator bdfContainerIterator_Type;

    M_unknownsBack.resize ( 0 );
    for (bdfContainerIterator_Type it = M_unknowns.begin(); it < M_unknowns.end(); ++it)
    {
        M_unknownsBack.push_back ( new feVectorPtr_Type ( **it ) );
    }

    M_timeStepBack = M_timeStep;
}


template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::restoreSolution()
{
    typedef typename feVectorPtrContainer_Type::iterator bdfContainerIterator_Type;

    M_unknowns.resize ( 0 );
    for ( bdfContainerIterator_Type it = M_unknownsBack.begin(); it < M_unknownsBack.end(); ++it )
    {
        M_unknowns.push_back ( new feVectorPtr_Type ( **it ) );
    }


    M_timeStep = M_timeStepBack;

    computeCoefficient();
}

template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::restoreSolution ( Real timeStep )
{
    typedef typename feVectorPtrContainer_Type::iterator bdfContainerIterator_Type;

    M_unknowns.resize ( 0 );
    for ( bdfContainerIterator_Type it = M_unknownsBack.begin(); it < M_unknownsBack.end(); ++it )
    {
        M_unknowns.push_back ( new feVectorPtr_Type ( **it ) );
    }

    M_timeStep = M_timeStepBack;
    M_timeStep[ 0 ] = timeStep;

    computeCoefficient();
}

template<typename FEVectorType>
void
TimeAdvanceBDFVariableStep<FEVectorType>::showMe() const
{
    std::cout << "*** BDF Time discretization of order " << M_order << " ***"
              << std::endl;
    std::cout << "    Coefficients: " << std::endl;
    for ( UInt i = 0; i < M_order + 1; ++i )
        std::cout << "       alpha(" << i << ") = " << M_alpha[ i ]
                  << std::endl;
    for ( UInt i = 0; i < M_order; ++i )
        std::cout << "       beta (" << i << ") = " << M_beta[ i ]
                  << std::endl;

    return ;
}

// ===================================================
//Set Methods
// ===================================================
template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::setTimeStep ( Real timeStep )
{
    M_timeStep[ 0 ] = timeStep;
    computeCoefficient();
}

// ===================================================
// Get Methods
// ===================================================
template<typename FEVectorType>
const Real&
TimeAdvanceBDFVariableStep<FEVectorType>::coefficientDerivative ( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT ( i < M_order + 1,
             "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return M_alpha[ i ];
}


template<typename FEVectorType>
Real
TimeAdvanceBDFVariableStep<FEVectorType>::coefficientDerivativeOverTimeStep ( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT ( i < M_order + 1,
             "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return M_alpha[ i ] / M_timeStep[ 0 ];
}


template<typename FEVectorType>
const Real&
TimeAdvanceBDFVariableStep<FEVectorType>::coefficientExtrapolation ( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT (  i < M_order,
              "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return M_beta[ i ];
}

// ===================================================
// Private Methods
// ===================================================
template<typename FEVectorType>
void TimeAdvanceBDFVariableStep<FEVectorType>::computeCoefficient()
{
    container_Type rho ( ScalarVector ( M_order, M_timeStep[ 0 ] ) );

    container_Type cumulativeSumTimeStep ( M_order );
    std::partial_sum ( M_timeStep.begin(), M_timeStep.end(), cumulativeSumTimeStep.begin() );

    std::transform ( rho.begin(), rho.end(), cumulativeSumTimeStep.begin(),
                     rho.begin(), std::divides<double>() );

    M_alpha[ 0 ] = boost::numeric::ublas::sum ( rho );

    for ( UInt j = 0; j < M_order; ++j )
    {
        Real tmp ( 1.0 );

        if ( rho[j] != 0 ) //rho[j] may be 0 with the start-up procedure
        {
            for ( UInt kk = 0; kk < M_order; ++kk )
            {
                if ( kk != j )
                {
                    tmp *= ( 1 - rho[ kk ] / rho[ j ] );
                }
            }
            M_alpha[ j + 1 ] = rho[ j ] / tmp;
            M_beta[ j ]    = 1. / tmp;
        }
        else //we don't have enough initial data to start a higher order BDF.
        {
            M_alpha[ j + 1 ] = 0;
            M_beta[ j ] = 0;
        }
    }
}

} //Namespace LifeV

#endif  /*_BDF_VARIABLE_TIMESTEP_H*/
