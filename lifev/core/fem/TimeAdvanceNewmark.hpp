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
    This class consider \f$\theta\f$-method for first order problems and
    TimeAdvanceNewmark scheme for the second order problems.

    @date 09-2010

    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
*/

#ifndef TimeAdvanceNewmark_H
#define TimeAdvanceNewmark_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
//! TimeAdvanceNewmark  - Class to deal the \f$theta\f$-method and TimeAdvanceNewmark scheme
/*!
  @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

  This class can be used to approximate problems of the first order and the second order in time.
  In the first case the temporal scheme is a theta-method, while in the second case is a TimeAdvanceNewmark scheme.

  This class defines the state vector \f$X^{n+1}\f$, a suitable  extrapolation of vector \f$X^*\f$ and   opportune coefficients used to determinate \f$X^{n+1}\f$ and \f$X^*\f$.
<ol>
<li> First order problem:

   \f[ M u' + A(u) = f \f]

   the state vector is \f[X^{n+1} = (U^{n+1},V^{n+1}, U^{n}, V^{n})\f].
   where \f$U\f$ is an approximation of \f$u\f$ and \f$V\f$  of \f$\dot{u}\f$.

   We consider a parameter \f$\theta\f$, and we apply the following theta method:

   \f[ U^{n+1} = U^{n} + \Delta t  \theta V^{n+1} + (1 − \theta) V^n,  \f]

    so the approximation of velocity at timestep \f$n+1\f$ is:

    \f[  V^{ n+1} = 1/ (\theta * \Delta t) (U^{n+1}-U^n)+( 1 -  1 / \theta) V^n;\f]

    We can  linearize non-linear term \f$A(u)\f$  as \f$A(U^*) U^{n+1}\f$,  where the \f$U^*\f$ becomes:

    \f[ U^∗ = U^n + \Delta t V^n, \f]

    The coefficients \f$\alpha\f$, \f$beta\f$ depend on \f$theta\f$,
     and we can use a explicit \f$theta\f$-method.
     </li>
    <li> Second order Method:

     \f[ M \ddot{u} + D(u, \dot{u})+ A(u) = f \f]

     the state vector is \f[X^{n+1} = (U^{n+1},V^{n+1}, W^{n+1}, U^{n}, V^{n}, W^{n}).\f]
     where U is an approximation of \f$u\f$ and \f$V\f$  of \f$\dot{u}\f$ and \f$W\f$ of \f$\ddot{u}\f$ .

     We introduce the parameters (\f$\theta\f$, \f$\gamma\f$)  and we apply the following TimeAdvanceNewmark method:

     \f[ U^{n+1} = U^{n} + \Delta t  \theta V^{n+1} + (1 − \theta) V^n,  \f]

     so the approximation of velocity at time step \f$n+1\f$ is

     \f[ V^{ n+1} = \beta / (\theta * \Delta t) (U^{n+1}-U^n)+( 1 - \gamma/ \theta) V^n+ (1- \gamma /(2\theta)) \Delta t W^n;\f]

     and we determine \f$W^{n+1}\f$ by

     \f[ W^{n+1} =1/(\theta \Delta t^2) (U^{n+1}-U^{n}) - 1 / \Delta t V^n + (1-1/(2\theta)) W^n \f].

     We can  linearize non-linear term \f$A(u)\f$  as \f$A(U^*) U^{n+1}\f$ and  \f$D(u, \dot{u})\f$ as \f$D(U^*, V^*) V^{n+1}\f$,
     where \f$U^*\f$ is given by

     \f[ U^∗ = U^n + \Delta t V^n \f]

      and  \f$V^*\f$ is

      \f[  V^* =V^n+ \Delta t W^n  \f]

      The coefficients \f$\xi\f$, \f$\alpha\f$, \f$\beta\f$, \f$\beta_V\f$ depend on \f$\theta\f$ and \f$\gamma\f$.
      </li>
      </ol>
*/

template<typename feVectorType = VectorEpetra >
class TimeAdvanceNewmark:
        public TimeAdvance < feVectorType >
{
public:

  //! @name Public Types
  //@{

  // class super;
  typedef TimeAdvance< feVectorType >                    super;
  // type of template
  typedef typename super::feVector_Type                  feVector_Type;
  // container of feVector
  typedef typename super::feVectorContainer_Type         feVectorContainer_Type;

  // container of pointer of feVector;
  typedef typename super::feVectorContainerPtr_Type      feVectorContainerPtr_Type;

  // iterator;
  typedef typename feVectorContainerPtr_Type::iterator   feVectorContainerPtrIterate_Type;

  // container of pointer of feVector;
  typedef typename super::feVectorSharedPtrContainer_Type        feVectorSharedPtrContainer_Type;

  //@}

  //! @name Constructor & Destructor
  //@{

  //! Empty  Constructor
  TimeAdvanceNewmark();

  //! Destructor
  virtual ~TimeAdvanceNewmark() {}

  //@}

  //! @name Methods
  //@{

  //!Update the state vector
  /*! Update the vectors of the previous time steps by shifting on the right  the old values.
    @param solution current (new) value of the state vector
  */
  void shiftRight(const feVector_Type& solution);

  //! Update the right hand side \f$ f_V \f$ of the time derivative formula
  /*!
    Return the right hand side \f$ f_V \f$ of the time derivative formula
    @param timeStep defined the  time step need to compute the
    @returns rhsV
  */
  void RHSFirstDerivative(const Real& timeStep, feVectorType& rhsContribution, int const shift = 0 ) const;

  //! Update the right hand side \f$ f_W \f$ of the time derivative formula
  /*!
    Set and Return the right hand side \f$ f_W \f$ of the time derivative formula
    @param timeStep defined the  time step need to compute the \f$ f_W \f$
  */
  void updateRHSSecondDerivative(const Real& timeStep = 1 );

  //!Show the properties  of temporal scheme
  void showMe(std::ostream& output = std::cout) const;

  //@}

  //!@name Set Methods
  //@{

  //! Initialize the parameters of time advance scheme
  /*!
    @param  order define the order of BDF;
    @param  orderDerivatve  define the order of derivative;
  */
  void setup ( const UInt& order,  const  UInt& orderDerivative )
  {
    ERROR_MSG("use setup for TimeAdvanceBDF but the time advance scheme is Newmark");
  }

  //! Initialize the parameters of time advance scheme
  /*!
    @param  coefficients define the TimeAdvanceNewmark's coefficients (\theta, \gamma);
    @param  orderDerivative  define the order of derivative;
  */
  void setup(const std::vector<Real>& coefficients, const  UInt& orderDerivative);

  //! Initialize the StateVector
  /*!
    Initialize all the entries of the unknown vector to be derived with the vector x0 (duplicated).
    this class is virtual because used in bdf;
    @param x0 is the initial solution;
  */
  void setInitialCondition( const feVector_Type& x0);

  //! initialize the state vector
  /*!
    Initialize all the entries of the unknown vector to be derived with the vector x0, v0 (duplicated).
    this class is virtual because used in \f$\theta\f$-method scheme;
    @param x0 is the initial unk;
    @param v0 is the initial velocity
  */
  void setInitialCondition( const feVector_Type& x0, const feVector_Type& v0 );

  //! initialize the state vector
  /*!
    Initialize all the entries of the unknown vector to be derived with the vector x0, v0,w0 (duplicated).
    this class is virtual because used in Newamrk scheme;
    @param x0 is the initial solution;
    @param v0 is the initial velocity
    @param w0 is the initial acceleration
  */
  void setInitialCondition(const feVector_Type& x0, const feVector_Type& v0, const feVector_Type& w0 );

  //! Initialize the state vector
  /*! Initialize all the entries of the unknown vector to be derived with a
    set of vectors x0
    note: this is taken as a copy (not a reference), since x0 is resized inside the method.
  */
  void setInitialCondition( const feVectorSharedPtrContainer_Type& x0);

  //@}

  //!@name Get Methods
  //@{

  //!Return the \f$i\f$-th coefficient of the unk's extrapolation
  /*!
    @param i index of  extrapolation coefficient
    @returns beta
  */
   Real coefficientExtrapolation(const  UInt& i ) const;

  //! Return the \f$i\f$-th coefficient of the velocity's extrapolation
  /*!
    @param \f$i\f$ index of the coefficient of the first derivative
    @returns beta
  */
  Real coefficientExtrapolationFirstDerivative(const UInt& i ) const;

  //! Compute the polynomial extrapolation of solution
  /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
  */
  //feVectorType   extrapolation() const;
  void extrapolation(feVector_Type& extrapolation) const;

  //! Compute the polynomial extrapolation of velocity
  /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
  */
  void extrapolationFirstDerivative(feVector_Type& extrapolation) const;

  //! Return the current velocity
  feVector_Type firstDerivative()  const
  {
    return( *this->M_unknowns[1]);
  }

  //!Return the current acceleration
  feVector_Type secondDerivative() const
  {
    return  *this->M_unknowns[2];
  }

  //@}

private:

  //! Coefficient of TimeAdvanceNewmark: \f$theta\f$
  Real M_theta;

  //! Coefficient of TimeAdvanceNewmark: \f$\gamma\f$
  Real M_gamma;

};

// ==================================================
// Constructors & Destructor
// ==================================================
template<typename feVectorType>
TimeAdvanceNewmark <feVectorType> ::TimeAdvanceNewmark():super()
{
}

// ===================================================
// Methods
// ===================================================
template<typename feVectorType>
void TimeAdvanceNewmark <feVectorType>::shiftRight(const feVector_Type& solution)
{
  ASSERT (  this->M_timeStep != 0 ,  "M_timeStep must be different to 0");

  feVectorContainerPtrIterate_Type it   =  this->M_unknowns.end();
  feVectorContainerPtrIterate_Type itb1 =  this->M_unknowns.begin() +  this->M_size/2;
  feVectorContainerPtrIterate_Type itb  =  this->M_unknowns.begin();

  for ( ; itb1 != it; ++itb1, ++itb)
    *itb1 = *itb;

  itb  =  this->M_unknowns.begin();

  // insert unk in unknowns[0];
  *itb = new feVector_Type( solution );

  itb++;

  //update velocity
  feVector_Type velocityTemp(solution);
  velocityTemp *=  this->M_alpha[0] / this->M_timeStep;
  velocityTemp -= *this->M_rhsContribution[0];

  // update unknows[1] with velocityTemp is current velocity
  *itb = new feVector_Type(velocityTemp);

  if ( this->M_orderDerivative == 2 )
    {
      itb++;

      //update acceleration
      feVector_Type accelerationTemp( solution );
      accelerationTemp *= this->M_xi[ 0 ] / ( this->M_timeStep * this->M_timeStep);
      accelerationTemp -= * this->M_rhsContribution[ 1 ];

      *itb = new feVector_Type( accelerationTemp );
    }
  return;
}

template<typename feVectorType>
void
TimeAdvanceNewmark<feVectorType>::RHSFirstDerivative(const Real& timeStep, feVectorType& rhsContribution, int const shift ) const
{
    rhsContribution *=  (this->M_alpha[ 1 ] / timeStep) ;

    Real timeStepPower(1.); // was: std::pow( timeStep, static_cast<Real>(i - 1 ) )

    for (UInt i= 1; i  < this->M_firstOrderDerivativeSize; ++i )
    {
        rhsContribution += ( this->M_alpha[ i+1 ] * timeStepPower ) *  (* this->M_unknowns[ i-shift ]);
        timeStepPower *= timeStep;
    }
}

template<typename feVectorType>
void
TimeAdvanceNewmark<feVectorType>::updateRHSSecondDerivative(const Real& timeStep )
{
  feVectorContainerPtrIterate_Type it =  this->M_rhsContribution.end()-1;

  *it = new feVector_Type(*this->M_unknowns[0]);

  **it *=  this->M_xi[ 1 ] /(timeStep * timeStep) ;

  for ( UInt i = 1;  i < this->M_secondOrderDerivativeSize; ++i )
    **it += ( this->M_xi[ i+1 ] * std::pow(timeStep, static_cast<Real>(i - 2) ) ) * ( *this->M_unknowns[ i ]);
}

template<typename feVectorType>
void
TimeAdvanceNewmark<feVectorType>::showMe(std::ostream& output ) const
{
    output << "*** TimeAdvanceNewmark discretization maximum order of derivate "
	   << this->M_orderDerivative<< " ***"<< std::endl;
    output <<" Coefficients : "      <<  std::endl;
    output <<" theta :        "      << M_theta<<"\n"
	   <<" gamma :  "            <<  M_gamma<<"\n"
	   <<" size unknowns :"      << this->M_size<<"\n";

    for ( UInt i = 0; i <  this->M_alpha.size(); ++i )
      output << "  alpha(" << i << ") = " <<  this->M_alpha[ i ]<< std::endl;

    if (this->M_orderDerivative == 2)
      {
        for ( UInt i = 0; i <  this->M_xi.size(); ++i )
	  output << "       xi(" << i << ") = " <<  this->M_xi[ i ] << std::endl;
      }

    output << "Delta Time : "<< this->M_timeStep<<"\n";
    output << "*************************************\n";
}

// ===================================================
// Set Methods
// ===================================================

template<typename feVectorType>
void
TimeAdvanceNewmark<feVectorType>::setup(const std::vector<Real>& coefficients, const  UInt& orderDerivative)
{
  //initialize theta
  M_theta = coefficients[0];

  //initilialize gamma
  M_gamma = coefficients[1];

  //initialize Order Derivative
  this->M_orderDerivative= orderDerivative;

  // If theta equal 0, explicit meta method
  if (M_theta == 0)
    {
      ASSERT (this->M_orderDerivative == 2,  "theta is 0 must be different from 0 in TimeAdvanceNewmark");
      this->M_size = 4;
      this->M_alpha[ 0 ] =  1;
      this->M_alpha[ 1 ] =  1;
      this->M_alpha[ 2 ] =  1;
      this->M_order  = 1;
    }
  else
    {
      if (this->M_orderDerivative == 1 )  // Theta method
        {
	  this->M_gamma = 1;
	  //  unknown vector's  dimension;
	  this->M_size = 4;
	  this->M_alpha.resize(3);
	  this->M_xi.resize(3);
	  this->M_beta.resize(3);
	  this->M_alpha[ 0 ] =  M_gamma / M_theta;
	  this->M_alpha[ 1 ] =  M_gamma / M_theta;
	  this->M_alpha[ 2 ] =  M_gamma / M_theta - 1.0;
	  this->M_beta[0] = 1;
	  this->M_beta[1] = 1;
	  this->M_beta[2] = 0.5;
	  this->M_xi[0]   = 0;
	  this->M_xi[1]   = 0;
	  this->M_xi[2]   = 0;
	  this->M_coefficientsSize = 3;
        }
      else     //TimeAdvanceNewmarkMethod
        {
	  //unknown vector's dimension
	  this->M_size = 6 ;
	  this->M_alpha.resize(4);
	  this->M_xi.resize(4);
	  this->M_beta.resize(3);
	  this->M_betaFirstDerivative.resize(3);
	  //initialitation alpha coefficients
	  this->M_alpha[ 0 ] =  M_gamma / M_theta;
	  this->M_alpha[ 1 ] =  M_gamma / M_theta;
	  this->M_alpha[ 2 ] =  M_gamma / M_theta - 1.0;
	  this->M_alpha[ 3 ] = M_gamma / (2.0 * M_theta) -1.0;

	  //initialitation xi coefficients
	  this->M_xi[ 0 ] =  1. / M_theta;
	  this->M_xi[ 1 ] =  1. / M_theta;
	  this->M_xi[ 2 ] =  1. / M_theta;
	  this->M_xi[ 3 ] =  1. / ( 2.0 * M_theta )-1.0;


	  //initialitation extrap coefficients
	  this->M_beta[ 0 ] = 1;
	  this->M_beta[ 1 ] = 1;
	  this->M_beta[ 2 ] = 0.5;
	  this->M_betaFirstDerivative[ 0 ] = 0;
	  this->M_betaFirstDerivative[ 1 ] = 1;
	  this->M_betaFirstDerivative[ 2 ] = 1;

	  this->M_coefficientsSize  = 4;
        }
      this->M_unknowns.reserve(this->M_size);
      this-> M_rhsContribution.reserve(2);
      // check order  scheme
      if (this->M_alpha[0] == 0.5)
	this->M_order = 2;
      else
	this->M_order = 1;

      this->M_firstOrderDerivativeSize  =  static_cast<Real>(this->M_size) / 2.0;
      this->M_secondOrderDerivativeSize = static_cast<Real>(this->M_size) / 2.0;
    }
}

template<typename feVectorType>
void TimeAdvanceNewmark<feVectorType>::setInitialCondition( const feVector_Type& x0)
{
  feVectorContainerPtrIterate_Type iter     = this->M_unknowns.begin();
  feVectorContainerPtrIterate_Type iter_end = this->M_unknowns.end();

  feVector_Type zero(x0);
  zero *=0;

  for ( ; iter != iter_end; ++iter )
    {
      delete *iter;
      *iter = new feVector_Type(zero, Unique);
    }

  for ( UInt i(this->M_unknowns.size()) ; i < this->M_size; ++i )
    this->M_unknowns.push_back(new feVector_Type(x0));

  feVectorContainerPtrIterate_Type iterRhs     = this->M_rhsContribution.begin();
  feVectorContainerPtrIterate_Type iterRhs_end = this->M_rhsContribution.end();

  for ( ; iterRhs !=iterRhs_end ; ++iterRhs)
    {
      delete *iterRhs;
      *iterRhs = new feVector_Type(zero, Unique);
    }
}

template<typename feVectorType>
void TimeAdvanceNewmark<feVectorType>::setInitialCondition( const feVector_Type& x0, const feVector_Type& v0)
{
  feVectorContainerPtrIterate_Type iter       = this->M_unknowns.begin();
  feVectorContainerPtrIterate_Type iter_end   = this->M_unknowns.end();

  //!initialize zero
  feVector_Type zero(x0);
  zero *=0;

  for ( ; iter != iter_end; ++iter )
    {
      delete *iter;
      *iter = new feVector_Type(zero, Unique);
    }

  this->M_unknowns.push_back(new feVector_Type(x0));
  this->M_unknowns.push_back(new feVector_Type(v0));

  for (UInt i=0; i<4; ++i )
    {
      this->M_unknowns.push_back(new feVector_Type(zero, Unique));
      *this->M_unknowns[i+2]  *=0;
    }
  this->setInitialRHS(zero);
}

template<typename feVectorType>
void TimeAdvanceNewmark<feVectorType>::setInitialCondition( const feVector_Type& x0, const feVector_Type& v0, const feVector_Type& w0)
{
  feVectorContainerPtrIterate_Type iter       = this->M_unknowns.begin();
  feVectorContainerPtrIterate_Type iter_end   = this->M_unknowns.end();

  //!initialize zero
  feVector_Type zero(x0);
  zero *=0;

  for ( ; iter != iter_end; ++iter )
    {
      delete *iter;
      *iter = new feVector_Type(zero, Unique);
    }

  this->M_unknowns.push_back(new feVector_Type(x0));
  this->M_unknowns.push_back(new feVector_Type(v0));
  this->M_unknowns.push_back(new feVector_Type(w0));

  for (UInt i=0; i<3; ++i )
    {
      this->M_unknowns.push_back(new feVector_Type(zero, Unique));
      *this->M_unknowns[i+3] *=0;
    }
  this->setInitialRHS(zero);
}

template<typename feVectorType>
void TimeAdvanceNewmark<feVectorType>::setInitialCondition( const feVectorSharedPtrContainer_Type& x0)
{
  const UInt n0 = x0.size();

  ASSERT( n0 != 0, "vector null " );

  feVectorContainerPtrIterate_Type iter     = this->M_unknowns.begin();
  feVectorContainerPtrIterate_Type iter_end = this->M_unknowns.end();

  UInt i(0);

  //!initialize zero
  feVector_Type zero(*x0[0]);
  zero *=0;

  for ( ; iter != iter_end && i< n0 ; ++iter, ++i )
    {
      delete *iter;
      *iter = new feVector_Type(*x0[i]);
    }

  for ( i = this->M_unknowns.size() ; i < this->M_size && i< n0; ++i )
    this->M_unknowns.push_back(new feVector_Type(*x0[i]));

  for ( i = this->M_unknowns.size() ; i < this->M_size; i++ )
    this->M_unknowns.push_back(new feVector_Type( *x0[ n0-1 ] ) );

  this->setInitialRHS(zero);
}

// ===================================================
// Get Methods
// ===================================================

template<typename feVectorType>
Real
TimeAdvanceNewmark<feVectorType>::coefficientExtrapolation(const UInt& i) const
{
  ASSERT ( i <  3 ,  "coeff_der i must equal 0 or 1 because U^*= U^n + timeStep*V^n + timeStep^2 / 2 W^n");

  switch (i)
  {
  case 0:
	  return this->M_beta(i);
  case 1:
	  return this->M_beta(i)*this->M_timeStep;
  case 2:
	  return this->M_beta(i)*this->M_timeStep*this->M_timeStep;
  default:
	  ERROR_MSG ("coeff_der i must equal 0 or 1 because U^*= U^n + timeStep*V^n + timeStep^2 / 2 W^n");
  }
  return 1;
}

template<typename feVectorType>
Real
TimeAdvanceNewmark<feVectorType>::coefficientExtrapolationFirstDerivative(const UInt& i ) const
{
    switch (i)
    {
    case 0:
        return this->M_betaFirstDerivative(i);
    case 1:
        return this->M_betaFirstDerivative(i)*this->M_timeStep;
    case 2:
        return this->M_betaFirstDerivative(i)*this->M_timeStep*this->M_timeStep;
    default:
        ERROR_MSG ("coeff_der i must equal 0 or 1 because U^*= U^n + timeStep*V^n + timeStep^2 / 2 W^n");
    }
	  return 1;
}

template<typename feVectorType>
void
TimeAdvanceNewmark<feVectorType>::extrapolation(feVector_Type& extrapolation) const
{
  extrapolation += this->M_timeStep * ( *this->M_unknowns[ 1 ]);

  if ( this->M_orderDerivative == 2 )
    extrapolation += ( this->M_timeStep * this->M_timeStep ) / 2.0 * ( *this->M_unknowns[2]);
}


template<typename feVectorType>
void
TimeAdvanceNewmark<feVectorType>::extrapolationFirstDerivative(feVector_Type& extrapolation) const
{
  ASSERT ( this->M_orderDerivative == 2,
	   "extrapolationFirstDerivative: this method must be used with the second order problem." )

  extrapolation = *this->M_unknowns[1];
  extrapolation += this->M_timeStep * ( *this->M_unknowns[ 2 ]);
}

// ===================================================
// Macros
// ===================================================

//! define the TimeAdvanceNewmark;  this class runs only the default template parameter.
inline
TimeAdvance< VectorEpetra >* createTimeAdvanceNewmark() { return new TimeAdvanceNewmark< VectorEpetra >(); }

namespace
{
  static bool registerTimeAdvanceNewmark = TimeAdvanceFactory::instance().registerProduct( "Newmark",  &createTimeAdvanceNewmark);
}

}
#endif  /*TimeAdvanceNewmark*/
