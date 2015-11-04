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
    @file
    @brief File containing a class to  deal the time advancing scheme.
    A class for an easy handling of different order time
    discretizations/extrapolations BDF based for first and second order problem

    @date 09-2010

    @author Simone Deparis  <simone.deparis@epfl.ch>
    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */


#ifndef TIMEADVANCEBDF_H
#define TIMEADVANCEBDF_H 1




#include <lifev/core/fem/TimeAdvance.hpp>

namespace LifeV
{
const UInt BDF_MAX_ORDER = 5;

//!class TimeAdvanceBDF - Backward differencing formula time discretization for the first and the second order problem in time.
/*!
 @author Simone Deparis  <simone.deparis@epfl.ch>
 @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

<ol>
<li> First order problem

A differential equation of the form

\f[ M \dot{u} + A u = f \f]

is discretizated in time as

\f[ M V_{k+1}  + A U_{k+1} = f_{k+1} \f]

where V denotes the polynomial of order n in t that interpolates
\f$(t_i,u_i)\f$ for \f$i = k-n+1,...,k+1\f$.

The approximative time derivative \f$ V \f$ is a linear
combination of state vectors \f$u_i\f$:

\f[ V_{k+1} = \frac{1}{\Delta t} (\alpha_0 U_{k+1} - \sum_{i=0}^n \alpha_i U_{k+1-i} )\f]

Thus we have

\f[ \frac{\alpha_0}{\Delta t} M U_{k+1} = A U_{k+1} + f + M f_V \f]

with

\f[ f_V= \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i U_{k+1-i} \f]

This class stores the n last state vectors in order to be able to
calculate \f$ f_V \f$. It also provides \f$\alpha_i\f$
and can extrapolate the the new state from the \f$n\f$ last states with a
polynomial of order \f$n-1\f$:

\f[ U_{k+1} \approx \sum_{i=0}^{n-1} \beta_i U_{k-i} \f]
</li>
<li>  Second order problem

A differential equation of the form

\f[ M \ddot{u}= D( u, \dot{u}) + A(u) + f \f]

is discretized in time as

\f[ M W_{k+1} = A(U^*) U_{k+1} +D(U^*, V^*) V_{k+1} + f_{k+1} \f]

where \f$W\f$ and \f$V\f$ denotes the polynomial of order \f$n+1\f$ and order \f$ n\f$,
while \f$U^*\f$ and \f$V^*\f$ are suitable extrapolations.

The velocity vector, as for first order  problem, is:

\f[ V_{k+1} = \frac{1}{\Delta t} (\alpha_0 U_{k+1} - \sum_{i=0}^n \alpha_i U_{k+1-i} )\f]

while the acceleration vector \f$W^{n+1}\f$ is:

\f[ W_{k+1} = \frac{1}{\Delta t^2} (\xi_0 U_{k+1} - \sum_{i=0}^{n+1} \xi_i U_{k+1-i} )\f]

Thus we have

\f[ \frac{\xi_0}{\Delta t^2} M U_{k+1} + \frac{\alpha_0}{\Delta t } D U_{k+1} + A U_{k+1} + f + M f_W + D f_V \f]

with

\f[ f_V= \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i U_{k+1-i} \f]

and

\f[  f_W =\frac{1}{\Delta t^2} \sum_{i=1}^{n+1} \xi_i U_{ k+1 - i}  \f]

It also provides \f$\alpha_i\f$  and can extrapolate the the new state from the \f$n\f$ last states with a
polynomial of order \f$n-1\f$:

\f[ U_{k+1} \approx  U^*= \sum_{i=0}^{n-1} \beta_i U_{k-i} \f]

and  \f$V^*\f$ in following way:

\f[  V_{k+1} \approx V^*=\sum_{i=0}^p \frac{\beta_i^V}{\Delta t} U^{n-i} = W^{n+1}+O(\Delta t^p),  \f]
</li>
</ol>
*/

template<typename feVectorType = VectorEpetra >
class TimeAdvanceBDF:
    public  TimeAdvance < feVectorType >
{
public:

    //!@name Public Types
    //@{

    //! class super
    typedef TimeAdvance< feVectorType >                    super;
    //! type of template
    typedef typename super::feVector_Type                  feVector_Type;

    //! container of feVector
    typedef typename super::feVectorContainer_Type         feVectorContainer_Type;

    //! container of pointer of feVector;
    typedef typename super::feVectorContainerPtr_Type      feVectorContainerPtr_Type;

    //! iterator;
    typedef typename feVectorContainerPtr_Type::iterator   feVectorContainerPtrIterate_Type;

    //! container of pointer of feVector;
    typedef typename super::feVectorSharedPtrContainer_Type        feVectorSharedPtrContainer_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty  Constructor

    TimeAdvanceBDF();

    //! Destructor
    virtual ~TimeAdvanceBDF() {}

    //@}

    //! @name Methods
    //@{

    //!Update the state vector
    /*! Update the vectors of the previous time steps by shifting on the right  the old values.
      @param solution current (new) value of the state vector
    */
    void shiftRight (const feVector_Type&  solution );


    //! Update the right hand side \f$ f_V \f$ of the time derivative formula
    /*!
      Return the right hand side \f$ f_V \f$ of the time derivative formula
      @param timeStep defined the  time step need to compute the
      @returns rhsV
    */
    void RHSFirstDerivative (const Real& timeStep, feVectorType& rhsContribution ) const;

    //! Update the right hand side \f$ f_W \f$ of the time derivative formula
    /*!
      Sets and Returns the right hand side \f$ f_W \f$ of the time derivative formula
      @param timeStep defined the  time step need to compute the \f$ f_W \f$
      @returns rhsW
    */
    void updateRHSSecondDerivative (const Real& timeStep = 1 );

    //!Show the properties  of temporal scheme
    void showMe (std::ostream& output = std::cout ) const;
    //@}

    //!@name Set Methods
    //@{

    //! Initialize the parameters of time advance scheme
    /*!
      Initialize parameters of time advance scheme;
      @param  order define the order of BDF;
      @param  orderDerivative  define the order of derivate;
    */
    void setup ( const UInt& order, const UInt& orderDerivative = 1 );

    //! Initialize the parameters of time advance scheme used in TimeAdvanceNewmark
    /*!
      @note: this setup does not run for BDF class;
    */
    void setup ( const  std::vector<Real>&  /*coefficients*/, const  UInt& /*orderDerivative*/)
    {
        ERROR_MSG ("use setup for TimeAdvanceNewmark but the time advance scheme is BDF");
    }

    //! Initialize the StateVector
    /*!
      Initialize all the entries of the unknown vector to be derived with the vector x0 (duplicated).
      this class is virtual because used in BDF;
      @param x0 is the initial unk;
    */
    void setInitialCondition ( const feVector_Type& x0);

    //! Initialize the StateVector used in TimeAdvanceNewmark
    void setInitialCondition (const feVector_Type& /* x0*/, const feVector_Type& /*v0*/ )
    {
        ERROR_MSG ( "this method  is not yet implemented" );
    }

    //! Initialize the StateVector used in TimeAdvanceNewmark
    void setInitialCondition (const feVector_Type& /*x0*/, const feVector_Type& /*v0*/, const feVector_Type& /* w0*/ )
    {
        ERROR_MSG ( "this method  is not yet implemented" );
    }

    //! Initialize all the entries of the unknown vector to be derived with a
    /*! set of vectors x0
      @note: this is taken as a copy (not a reference), since x0 is resized inside the method.
    */
    void setInitialCondition (const feVectorSharedPtrContainer_Type& x0 );

    //@}

    //!@name Get Methods
    //@{

    //!Return the \f$i\f$-th coefficient of the unk's extrapolation
    /*!
      @param \f$i\f$ index of  extrapolation coefficient
      @returns beta
    */
    inline Real coefficientExtrapolation (const UInt& i ) const;

    //! Return the \f$i\f$-th coefficient of the velocity's extrapolation
    /*!
      @param \f$i\f$ index of velocity's extrapolation  coefficient
      @returns betaFirstDerivative
    */
    inline Real coefficientExtrapolationFirstDerivative (const UInt& i ) const;

    //! Compute the polynomial extrapolation of solution
    /*!
      Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
      \f$u^{n+1}\f$ defined by the n stored state vectors
    */

    void extrapolation (feVector_Type& extrapolation) const;

    //! Compute the polynomial extrapolation of velocity
    /*!
      Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
      \f$u^{n+1}\f$ defined by the n stored state vectors
    */
    void extrapolationFirstDerivative (feVector_Type& extrapolation) const;

    //! Return the current velocity
    feVectorType firstDerivative()  const;

    //!Return the current acceleration
    feVectorType secondDerivative() const;

    //@}

};

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename feVectorType>
TimeAdvanceBDF <feVectorType> :: TimeAdvanceBDF() :
    super()
{
}

// ===================================================
// Methods
// ===================================================
template<typename feVectorType>
void
TimeAdvanceBDF<feVectorType>::shiftRight (feVector_Type const&  solution )
{
    ASSERT ( this->M_unknowns.size() == this->M_size,
             "M_unknowns.size() and  M_size must be equal" );

    feVectorContainerPtrIterate_Type it   = this->M_unknowns.end() - 1;
    feVectorContainerPtrIterate_Type itm1 = this->M_unknowns.end() - 1;
    feVectorContainerPtrIterate_Type itb  = this->M_unknowns.begin();

    delete *itm1;

    for ( ; it != itb; --it )
    {
        itm1--;
        *it = *itm1;
    }

    *itb = new feVector_Type (solution);

}

template<typename feVectorType>
void
TimeAdvanceBDF<feVectorType>::RHSFirstDerivative (const Real& timeStep, feVectorType& rhsContribution ) const
{

    rhsContribution *= (this->M_alpha[ 1 ] / timeStep);

    for ( UInt i = 1; i < this->M_order; ++i )
    {
        rhsContribution += (this->M_alpha[ i + 1 ] / timeStep) * *this->M_unknowns[ i ];
    }
}


template<typename feVectorType>
void
TimeAdvanceBDF<feVectorType>::updateRHSSecondDerivative (const Real& timeStep )
{
    ASSERT ( this->M_orderDerivative == 2 ,
             " M_orderDerivative must be equal two" );

    feVectorContainerPtrIterate_Type it  = this->M_rhsContribution.end() - 1;

    *it = new feVector_Type (*this->M_unknowns[ 0 ]);

    ** it *= this->M_xi[ 1 ] / (timeStep * timeStep);

    for ( UInt i = 1; i < this->M_order + 1; ++i )
    {
        **it += ( this->M_xi[ i + 1 ] / (timeStep * timeStep) ) * *this->M_unknowns[ i ];
    }
}

template<typename feVectorType>
void
TimeAdvanceBDF<feVectorType>::showMe (std::ostream& output) const
{
    output << "*** BDF Time discretization of order " << this->M_order << " maximum order of derivate " << this->M_orderDerivative << " ***"
           << std::endl;
    output << "    Coefficients: " << std::endl;
    for ( UInt i = 0; i < this->M_order + 1; ++i )
        output << "       alpha(" << i << ") = " << this->M_alpha[ i ]
               << std::endl;
    for ( UInt i = 0; i < this->M_order; ++i )
        output << "       beta (" << i << ") = " << this->M_beta[ i ]
               << std::endl;
    if (this->M_orderDerivative == 2)
    {
        for ( UInt i = 0; i < this->M_order + this->M_orderDerivative; ++i )
        {
            output << "     xi(" << i << ") = " << this->M_xi[ i ]  << std::endl;
        }
        for ( UInt i = 0;  i < this->M_order + 1; ++i  )
            output << "       beta of the extrapolation of the first derivative ("
                   << i << ") = " << this->M_betaFirstDerivative[ i ]
                   << std::endl;
    }
    output << "Delta Time : " << this->M_timeStep << "\n";
    output << "*************************************\n";
    return ;
}

// ===================================================
// Set Methods
// ==================================================

template<typename feVectorType>
void
TimeAdvanceBDF<feVectorType>::setup ( const UInt& order, const UInt& orderDerivative)
{
    if ( order <= 0 || order > BDF_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
             << " you want to use BDF order " << order << "\n"
             << " we support BDF order from 1 to " << BDF_MAX_ORDER << "\n";
        throw std::invalid_argument ( __ex.str() );
    }

    this->M_order = order ;
    this->M_orderDerivative = orderDerivative ;
    this->M_size = order ;
    this->M_alpha.resize ( order + 1 );
    this->M_xi.resize ( order + 2 );
    this->M_beta.resize ( order );
    this->M_betaFirstDerivative.resize ( order + 1 );
    this->M_coefficientsSize = order + orderDerivative;

    switch ( order )
    {
        case 1:
            this->M_alpha[ 0 ] = 1.; // Backward Euler
            this->M_alpha[ 1 ] = 1.;
            this->M_beta[ 0 ] = 1.; // u^{n+1} \approx u^n
            this->M_xi[ 0 ] = 1.;
            this->M_xi[ 1 ] = 2.;
            this->M_xi[ 2 ] = -1.;
            this->M_betaFirstDerivative[0] = 2.;
            this->M_betaFirstDerivative[0] = -1.;
            break;
        case 2:
            this->M_alpha[ 0 ] = 3. / 2.;
            this->M_alpha[ 1 ] = 2.;
            this->M_alpha[ 2 ] = -1. / 2.;
            this->M_beta[ 0 ] = 2.;
            this->M_beta[ 1 ] = -1.;
            this->M_xi[ 0 ] =  2.;
            this->M_xi[ 1 ] = 5.;
            this->M_xi[ 2 ] =  -4.;
            this->M_xi[ 3 ] = 1.;
            this->M_betaFirstDerivative[0] = 3.;
            this->M_betaFirstDerivative[1] = -3.;
            this->M_betaFirstDerivative[2] = 1.;
            break;
        case 3:
            this->M_alpha[ 0 ] = 11. / 6.;
            this->M_alpha[ 1 ] = 3.;
            this->M_alpha[ 2 ] = -3. / 2.;
            this->M_alpha[ 3 ] = 1. / 3.;
            this->M_beta[ 0 ] = 3.;
            this->M_beta[ 1 ] = -3.;
            this->M_beta[ 2 ] = 1.;
            this->M_xi[ 0 ] =  35. / 12.;
            this->M_xi[ 1 ] =  26. / 3.;
            this->M_xi[ 2 ] =  -19. / 2;
            this->M_xi[ 3 ] =  14. / 3.;
            this->M_xi[ 4 ] =  -11. / 12.;
            this->M_betaFirstDerivative[0] = 4.;
            this->M_betaFirstDerivative[ 1 ] = -6.;
            this->M_betaFirstDerivative[ 2 ] = 4.;
            this->M_betaFirstDerivative[ 3 ] = -1.;
            break;
        case 4:
            this->M_alpha[ 0 ] = 25. / 12.;
            this->M_alpha[ 1 ] = 4.;
            this->M_alpha[ 2 ] = -3. ;
            this->M_alpha[ 3 ] = 4. / 3.;
            this->M_alpha[ 4 ] = - 1 / 4.;
            this->M_beta[ 0 ] = 4.;
            this->M_beta[ 1 ] = -6.;
            this->M_beta[ 2 ] = 4.;
            this->M_beta[ 3 ] = -1;
            this->M_xi[ 0 ] =  15. / 4.;
            this->M_xi[ 1 ] =  77. / 6.;
            this->M_xi[ 2 ] =  -107. / 6.;
            this->M_xi[ 3 ] =  13.;
            this->M_xi[ 4 ] =  -61. / 12.;
            this->M_xi[ 5 ] =  5. / 6.;
            this->M_betaFirstDerivative[ 0 ] = 5.;
            this->M_betaFirstDerivative[ 1 ] = -10.;
            this->M_betaFirstDerivative[ 2 ] = 10.;
            this->M_betaFirstDerivative[ 3 ] = -5.;
            this->M_betaFirstDerivative[ 4 ] = 1.;
            break;
        case 5:
            this->M_alpha[ 0 ] = 137. / 60.;
            this->M_alpha[ 1 ] = 5.;
            this->M_alpha[ 2 ] = -5. ;
            this->M_alpha[ 3 ] = 10. / 3.;
            this->M_alpha[ 4 ] = -5. / 4.;
            this->M_alpha[ 5 ] =  1. / 5.;
            this->M_beta[ 0 ] = 5.;
            this->M_beta[ 1 ] = -10.;
            this->M_beta[ 2 ] = 10.;
            this->M_beta[ 3 ] = -5.;
            this->M_beta[ 4 ] = 1.;
            this->M_xi[ 0 ] =  203. / 45.;
            this->M_xi[ 1 ] =  87. / 5.;
            this->M_xi[ 2 ] =  -117. / 4.;
            this->M_xi[ 3 ] =  254. / 9.;
            this->M_xi[ 4 ] =  -33. / 2.;
            this->M_xi[ 5 ] =  27. / 5.;
            this->M_xi[ 6 ] =  -137. / 180.;
            this->M_betaFirstDerivative[ 0 ] = 6.;
            this->M_betaFirstDerivative[ 1 ] = -15. / 2.;
            this->M_betaFirstDerivative[ 2 ] = 20. / 3.;
            this->M_betaFirstDerivative[ 3 ] = -15. / 4.;
            this->M_betaFirstDerivative[ 4 ] = 6. / 5.;
            this->M_betaFirstDerivative[ 5 ] = -1.;
            break;
    }

    if ( this->M_orderDerivative == 2 )
    {
        this->M_size++;
    }
    this->M_unknowns.reserve ( this->M_size);
    this->M_rhsContribution.reserve (2);
}

template<typename feVectorType>
void TimeAdvanceBDF<feVectorType>::setInitialCondition ( const feVector_Type& x0)
{
    feVectorContainerPtrIterate_Type iter     = this->M_unknowns.begin();
    feVectorContainerPtrIterate_Type iter_end = this->M_unknowns.end();

    for ( ; iter != iter_end; iter++ )
    {
        delete *iter;
        *iter = new feVector_Type (x0);
    }

    for ( UInt i (this->M_unknowns.size() ) ; i < this->M_order; i++ )
    {
        this->M_unknowns.push_back (new feVector_Type (x0) );
    }

    feVector_Type zero (x0);
    zero *= 0;
    this->setInitialRHS (zero);

    ASSERT ( this->M_unknowns.size() == this->M_order,
             "M_unknowns.size() and  M_order must be equal" );
}

template<typename feVectorType>
void TimeAdvanceBDF<feVectorType>::setInitialCondition ( const feVectorSharedPtrContainer_Type& x0)
{
    UInt n0 = x0.size();

    ASSERT ( n0 != 0, "Initial data are not enough for the selected BDF" );

    feVectorContainerPtrIterate_Type iter     = this->M_unknowns.begin();
    feVectorContainerPtrIterate_Type iter_end = this->M_unknowns.end();

    UInt i (0);

    for ( ; iter != iter_end && i < n0 ; ++iter, ++i )
    {
        delete *iter;
        *iter = new feVector_Type (*x0[i]);
    }

    for ( i = this->M_unknowns.size() ; i < this->M_order && i < n0; ++i )
    {
        this->M_unknowns.push_back (new feVector_Type (*x0[i]) );
    }

    if (this->M_orderDerivative == 1)
    {
        for ( i = this->M_unknowns.size() ; i < this->M_order; ++i )
        {
            this->M_unknowns.push_back (new feVector_Type (*x0[n0 - 1]) );
        }

        ASSERT ( this->M_unknowns.size() == this->M_order,
                 "M_unknowns.size() and  M_order must be equal" );
    }
    if (this->M_orderDerivative == 2)
    {
        for ( i = this->M_unknowns.size() ; i < this->M_order + 1; ++i )
        {
            this->M_unknowns.push_back (new feVector_Type (*x0[n0 - 1]) );
        }

        ASSERT ( this->M_unknowns.size() == this->M_order + 1,
                 "M_unknowns.size() and  M_order must be equal" );
    }

    //!initialize zero
    feVector_Type zero (*x0[0]);
    zero *= 0;
    this->setInitialRHS (zero);
}

// ===================================================
// Get Methods
// ===================================================

template<typename feVectorType>
inline Real
TimeAdvanceBDF<feVectorType>::coefficientExtrapolation (const UInt& i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT ( i < this->M_order,
             "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return this->M_beta[ i ];
}

template<typename feVectorType>
inline double
TimeAdvanceBDF<feVectorType>::coefficientExtrapolationFirstDerivative (const UInt& i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT ( i < this->M_order,
             "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return this->M_betaFirstDerivative[ i ];
}

template<typename feVectorType>
void
TimeAdvanceBDF<feVectorType>::extrapolation (feVector_Type& extrapolation) const
{
    extrapolation = this->M_beta[ 0 ] * (*this->M_unknowns[ 0 ]);

    for ( UInt i = 1; i < this->M_order; ++i )
    {
        extrapolation += this->M_beta[ i ] * (*this->M_unknowns[ i ]);
    }
}

template<typename feVectorType>
void
TimeAdvanceBDF<feVectorType>::extrapolationFirstDerivative (feVector_Type& extrapolation) const
{
    ASSERT ( this->M_orderDerivative == 2,
             "extrapolationFirstDerivative: this method must be used with the second order problem." )

    extrapolation = this->M_betaFirstDerivative[ 0 ] * (*this->M_unknowns[ 0 ]);

    for ( UInt i = 1; i < this->M_order; ++i )
    {
        extrapolation += this->M_betaFirstDerivative[ i ] * (*this->M_unknowns[ i ]);
    }
}

template<typename feVectorType>
feVectorType
TimeAdvanceBDF<feVectorType>::firstDerivative() const
{
    // Note by Cristiano Malossi - 13/02/2013
    //
    // This method is valid only if
    //
    // 1) shiftRight has been called
    // 2) updateRHSFirstDerivative has NOT been called
    //
    // TODO In case both shiftRight and updateRHSFirstDerivative have been called
    // a possible implementation is the one commented below.
    //
    //    feVector_Type derivative( *this->M_unknowns[ 0 ] * this->M_alpha[ 0 ] );
    //    for ( UInt i = 1; i <= this->M_order; ++i )
    //    {
    //        derivative -= *this->M_unknowns[ i ] * this->M_alpha[ i ];
    //    }
    //
    //    return derivative / this->M_timeStep;
    //
    // Before going in this direction the design of the TimeAdvance needs to be discussed.
    // The same consideration is valid for the second derivative.
    return *this->M_unknowns[0] * this->M_alpha[ 0 ] / this->M_timeStep
           - ( *this->M_rhsContribution[0] );
}

template<typename feVectorType>
feVectorType
TimeAdvanceBDF<feVectorType>::secondDerivative() const
{
    return (*this->M_unknowns[0]) * this->M_xi[ 0 ] / (this->M_timeStep * this->M_timeStep) -  ( *this->M_rhsContribution[1] );
}


// ===================================================
// Macros
// ===================================================

//! define the BDF factory; this class runs only the default template parameter.
inline
TimeAdvance<VectorEpetra>* createBDF()
{
    return new TimeAdvanceBDF<VectorEpetra>();
}

namespace
{
static bool registerBDF = TimeAdvanceFactory::instance().registerProduct ( "BDF",  &createBDF );
}

}  //Namespace LifeV

#endif  /*TIMEADVANCEBDF */
