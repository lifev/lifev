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
    @date
 
    @author Simone Deparis  <simone.deparis@epfl.ch>
    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 
    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */


#ifndef TIMEADVANCEBDF_H
#define TIMEADVANCEBDF_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <life/lifearray/EpetraVector.hpp>
#include <life/lifefem/timeAdvance_template.hpp>

namespace LifeV
{
const UInt BdfT_MAX_ORDER = 4;

  //!class BdfT - Backward differencing formula time discretization for the first and the second order problem in time.
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

template<typename feVectorType = EpetraVector >
class BdfT:
        public  TimeAdvance < feVectorType >
{
public:

    //!@name Public Types
    //@{

     typedef TimeAdvance< feVectorType >                   super;
     typedef typename super::feVectorContainer_Type        feVectorContainer_Type ;
     typedef typename super::feVectorContainerPtr_Type     feVectorContainerPtr_Type;
     typedef typename feVectorContainerPtr_Type::iterator  feVectorContainerPtrIterate_Type;

    //@}

     //! @name Constructor & Destructor
    //@{

    //! Empty  Constructor
    
    BdfT();
     //! Constructor
     /*! 
     @param  order of the BDF
     */
    //  BdfT( const UInt& order);

    /*! Constructor
     @param order: is accurancy's order of the BDF,
     @param orderDerivate: is the maximum order of derivate
     */
    //BdfT( const UInt& order, const  UInt& orderDerivate );

     //! Destructor
     ~BdfT() {}

   //@}

    //! @name Methods
    //@{
  
     //!Update the state vector
     /*! Update the vectors of the previous time steps by shifting on the right  the old values.
     @param solution current (new) value of the state vector
     */
     void shiftRight(const feVectorType&  solution );

     //! Update the right hand side \f$ f_V \f$ of the time derivative formula
     /*!
     Set and Return the right hand side \f$ f_V \f$ of the time derivative formula
     @param timeStep defined the  time step need to compute the
     @returns rhsV
     */
     feVectorType updateRHSFirstDerivative(const Real& timeStep = 1 );
   
     //! Update the right hand side \f$ f_W \f$ of the time derivative formula
     /*
     Sets and Returns the right hand side \f$ f_W \f$ of the time derivative formula
     @param timeStep defined the  time step need to compute the \f$ f_W \f$
     @returns rhsW
     */
     feVectorType updateRHSSecondDerivative(const Real& timeStep = 1 );

     //!Show the properties  of temporal scheme
     void showMe() const;
    //@}

    //!@name Set Methods
    //@{

     //! Initialize the parameters of time advance scheme
     /*
     Initialize parameters of time advance scheme;
     @param  order define the order of BDF;
     @param  orderDerivate  define the order of derivate;
     */
    void setup ( const UInt& order, const UInt& orderDerivate = 1 );

    //! Initialize the parameters of time advance scheme used in Newmark

    void setup ( const  std::vector<Real>&  coefficients, const  UInt& orderDerivate);

     //! Initialize the StateVector
     /*!
     Initialize all the entries of the unknown vector to be derived with the vector x0 (duplicated).
     this class is virtual because used in BDF;
     @param x0 is the initial unk;
     */
    void setInitialCondition( const feVectorType& x0);

    //! Initialize the StateVector used in Newmark
    void setInitialCondition(const feVectorType& x0, const feVectorType& v0 );

    //! Initialize the StateVector used in Newmark
    void setInitialCondition(const feVectorType& x0, const feVectorType& v0, const feVectorType&  w0 );
   
    //! Initialize all the entries of the unknown vector to be derived with a
    //! set of vectors x0
    //! note: this is taken as a copy (not a reference), since x0 is resized inside the method.
    void setInitialCondition(const feVectorContainer_Type& x0 );
  
    //@}

    //!@name Get Methods
    //@{

    //!Return the \f$i\f$-th coefficient of the unk's extrapolation
    /*!
    @param \f$i\f$ index of  extrapolation coefficient
    @returns beta
    */
    Real coefficientExtrapolation(const UInt& i ) const;
   
    //! Return the \f$i\f$-th coefficient of the velocity's extrapolation
    /*!
    @param \f$i\f$ index of velocity's extrapolation  coefficient
    @returns betaVelocity
    */
    Real coefficientExtrapolationVelocity(const UInt& i ) const;
   
    //! Compute the polynomial extrapolation of solution
    /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
    */ 
    feVectorType extrapolation() const;

    //! Compute the polynomial extrapolation of velocity
    /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
    */
    feVectorType extrapolationVelocity() const;

    //! Return the current velocity
    feVectorType velocity()  const;

    //!Return the current accelerate
   feVectorType accelerate() const;

   //@}

};

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename feVectorType>
BdfT <feVectorType> :: BdfT() :
        super()
{
}

// ===================================================
// Methods
// ===================================================
template<typename feVectorType>
void
BdfT<feVectorType>::shiftRight(feVectorType const&  solution )
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

    *itb = new feVectorType(solution);

}

template<typename feVectorType>
feVectorType
BdfT<feVectorType>::updateRHSFirstDerivative(const Real& timeStep )
{
    feVectorContainerPtrIterate_Type it  = this->M_rhsContribution.begin();

    feVectorType fv ( *this->M_unknowns[ 0 ] );

    fv *= this->M_alpha[ 1 ] / timeStep;

    for ( UInt i = 1; i < this->M_order; ++i )
    {
        fv += (this->M_alpha[ i + 1 ] / timeStep) * *this->M_unknowns[ i ];
    }

    *it = new feVectorType( fv );
    return fv;
}

template<typename feVectorType>
feVectorType
BdfT<feVectorType>::updateRHSSecondDerivative(const Real& timeStep )
{
    ASSERT ( this->M_orderDerivate== 2 ,
             " M_orderDerivatemust be equal two" );

    feVectorContainerPtrIterate_Type it  = this->M_rhsContribution.end()-1;

    feVectorType fw(*this->M_unknowns[ 0 ]);

    fw *= this->M_xi[ 1 ] / (timeStep*timeStep);

    for ( UInt i = 1; i < this->M_order + 1; ++i )
        fw += ( this->M_xi[ i + 1 ] / (timeStep*timeStep) ) * *this->M_unknowns[ i ];

    *it = new feVectorType( fw );

    return fw;
}

template<typename feVectorType>
void
BdfT<feVectorType>::showMe() const
{
    std::cout << "*** BDF Time discretization of order " << this->M_order << " maximum order of derivate "<< this->M_orderDerivate<< " ***"
              << std::endl;
    std::cout << "    Coefficients: " << std::endl;
    for ( UInt i = 0; i < this->M_order + 1; ++i )
        std::cout << "       alpha(" << i << ") = " << this->M_alpha[ i ]
                  << std::endl;
    for ( UInt i = 0; i < this->M_order; ++i )
        std::cout << "       beta (" << i << ") = " << this->M_beta[ i ]
                  << std::endl;
    if (this->M_orderDerivate==2)
    {
        for ( UInt i = 0; i < this->M_order + this->M_orderDerivate; ++i )
            std::cout << "     xi(" << i << ") = " << this->M_xi[ i ]  << std::endl;
        for ( UInt i = 0;  i < this->M_order+1; ++i  )
            std::cout << "       beta Velocity (" << i << ") = " << this->M_betaVelocity[ i ]
                      << std::endl;
    }
    std::cout <<"Delta Time : "<<this->M_timeStep<<"\n";
    std::cout <<"*************************************\n";
    return ;
}

// ===================================================
// Set Methods
// ==================================================

template<typename feVectorType>
void
BdfT<feVectorType>::setup( const UInt& order, const UInt& orderDerivate)
{
    if ( order <= 0 || order > BdfT_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
        << " you want to use BDF order " << order << "\n"
        << " we support BDF order from 1 to " << BdfT_MAX_ORDER << "\n";
        throw std::invalid_argument( __ex.str() );
    }

    this->M_order = order ;
    this->M_orderDerivate= orderDerivate ;
    this->M_size = order ;
    this->M_alpha.resize( order + 1 );
    this->M_xi.resize( order + 2 );
    this->M_beta.resize( order );
    this->M_betaVelocity.resize( order + 1 );
    this->M_coefficientsSize = order + orderDerivate ;

    switch ( order )
    {
    case 1:
        this->M_alpha[ 0 ] = 1.; // Backward Euler
        this->M_alpha[ 1 ] = 1.;
        this->M_beta[ 0 ] = 1.; // u^{n+1} \approx u^n
        this->M_xi[ 0 ] = 1.;
        this->M_xi[ 1 ] = 2;
        this->M_xi[ 2 ] = -1;
        this->M_betaVelocity[0]=2.;
        this->M_betaVelocity[1]=-1;
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
        this->M_betaVelocity[0]=3.;
        this->M_betaVelocity[ 1 ] = -3.;
        this->M_betaVelocity[ 2 ] = 1.;
        break;
    case 3:
        this->M_alpha[ 0 ] = 11. / 6.;
        this->M_alpha[ 1 ] = 3.;
        this->M_alpha[ 2 ] = -3. / 2.;
        this->M_alpha[ 3 ] = 1. / 3.;
        this->M_beta[ 0 ] = 3.;
        this->M_beta[ 1 ] = -3.;
        this->M_beta[ 2 ] = 1.;
        this->M_xi[ 0 ] =  35./12.;
        this->M_xi[ 1 ] =  26./3.;
        this->M_xi[ 2 ] =  -19./2;
        this->M_xi[ 3 ] =  14./3.;
        this->M_xi[ 4 ] =  -11./12.;
        this->M_betaVelocity[ 0 ] = 4.;
        this->M_betaVelocity[ 1 ] = -6.;
        this->M_betaVelocity[ 2 ] = 4.;
        this->M_betaVelocity[ 3 ] = -1;
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
        break;
    }

    if ( this->M_orderDerivate== 2 )
        this->M_size++;
    this->M_unknowns.reserve( this->M_size);
    this->M_rhsContribution.reserve(2);
}

template<typename feVectorType>
void
BdfT<feVectorType>::setup ( const  std::vector<Real>&  /*coefficients*/,  const  UInt& /*orderDerivate*/)
{
    ERROR_MSG("use setup for Newmark but the time advance scheme is BDF");
}

template<typename feVectorType>
void BdfT<feVectorType>::setInitialCondition( const feVectorType& x0)
{
    feVectorContainerPtrIterate_Type iter     = this->M_unknowns.begin();
    feVectorContainerPtrIterate_Type iter_end = this->M_unknowns.end();

    for ( ; iter != iter_end; iter++ )
    {
        delete *iter;
        *iter = new feVectorType(x0);
    }

    for ( UInt i(this->M_unknowns.size()) ; i < this->M_order; i++ )
        this->M_unknowns.push_back(new feVectorType(x0));

    feVectorType zero(x0);
    zero *=0;
    this->setInitialRHS(zero);
   
    ASSERT ( this->M_unknowns.size() == this->M_order,
             "M_unknowns.size() and  M_order must be equal" );
}

template<typename feVectorType>
void  BdfT<feVectorType>::setInitialCondition( const feVectorType& /*x0*/, const feVectorType& /*v0*/)
{
    ERROR_MSG( "this method  is not yet implemented" );
}

template<typename feVectorType>
void BdfT<feVectorType>::setInitialCondition( const feVectorType& /*x0*/, const feVectorType& /*v0*/, const feVectorType& /*w0*/)
{
    ERROR_MSG( "this method  is not yet implemented" );
}

template<typename feVectorType>
void BdfT<feVectorType>::setInitialCondition( const feVectorContainer_Type& x0)
{
    UInt n0 = x0.size();

    ASSERT( n0 != 0, "Initial data are not enough for the selected BDF" );

    feVectorContainerPtrIterate_Type iter     = this->M_unknowns.begin();
    feVectorContainerPtrIterate_Type iter_end = this->M_unknowns.end();

    UInt i(0);

    for ( ; iter != iter_end && i< n0 ; ++iter, ++i )
    {
        delete *iter;
        *iter = new feVectorType(x0[i]);
    }

    for ( i = this->M_unknowns.size() ; i < this->M_order && i< n0; ++i )
        this->M_unknowns.push_back(new feVectorType(x0[i]));

    if (this->M_orderDerivate == 1)
    {
        for ( i = this->M_unknowns.size() ; i < this->M_order; ++i )
            this->M_unknowns.push_back(new feVectorType(x0[n0-1]));

        ASSERT ( this->M_unknowns.size() == this->M_order,
                 "M_unknowns.size() and  M_order must be equal" );
    }
    if (this->M_orderDerivate == 2)
    {
        for ( i = this->M_unknowns.size() ; i < this->M_order + 1; ++i )
            this->M_unknowns.push_back(new feVectorType(x0[n0-1]));

        ASSERT ( this->M_unknowns.size() == this->M_order + 1,
                 "M_unknowns.size() and  M_order must be equal" );
    }

    //!initialize zero 
    feVectorType zero(x0[0]);
    zero *=0;
    this->setInitialRHS(zero);
}

// ===================================================
// Get Methods
// ===================================================  

template<typename feVectorType>
Real
BdfT<feVectorType>::coefficientExtrapolation(const UInt& i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i < this->M_order,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return this->M_beta[ i ];
}

template<typename feVectorType>
double
BdfT<feVectorType>::coefficientExtrapolationVelocity (const UInt& i ) const
{
      // Pay attention: i is c-based indexed
    ASSERT( i < this->M_order+1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return this->M_betaVelocity[ i +1];
}

template<typename feVectorType>
feVectorType
BdfT<feVectorType>::extrapolation() const
{
    feVectorType ue(*this->M_unknowns[ 0 ]);
    ue *= this->M_beta[ 0 ];

    for ( UInt i = 1; i < this->M_order; ++i )
    {
        ue += this->M_beta[ i ] * *this->M_unknowns[ i ];
    }

    return ue;
}

template<typename feVectorType>
feVectorType
BdfT<feVectorType>::extrapolationVelocity() const
{
    feVectorType velocity(*this->M_unknowns[ 0 ]);
    velocity *= this->M_betaVelocity[ 0 ];

    for ( UInt i = 1; i < this->M_order+1; ++i )
        velocity += this->M_betaVelocity[ i ] * *this->M_unknowns[ i ];

    return velocity;
}

template<typename feVectorType>
feVectorType
BdfT<feVectorType>::velocity()  const
{
    feVectorType velocity(*this->M_unknowns[0]);
    velocity  *= this->M_alpha[ 0 ] / this->M_timeStep;
    velocity  -= (*this->M_rhsContribution[0]);
    return velocity;
}

template<typename feVectorType>
feVectorType
BdfT<feVectorType>::accelerate() const
{
    feVectorType accelerate(*this->M_unknowns[0]);
    accelerate  *= this->M_xi[ 0 ]  /  (this->M_timeStep*this->M_timeStep);
    accelerate  -=  ( *this->M_rhsContribution[1] );
    return accelerate;
}

// ===================================================
// Macros
// ===================================================

//! define the BDF factory
inline
TimeAdvance<EpetraVector>* createBDF() { return new BdfT<EpetraVector>(); }

namespace
{
static bool registerBDF = TimeAdvanceFactory::instance().registerProduct( "BDF",  &createBDF );
}

}  //Namespace LifeV

#endif  /*TIMEADVANCEBDF */
