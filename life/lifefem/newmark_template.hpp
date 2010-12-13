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
    Newmark scheme for the second order problems.
 
    @date
 
    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
*/

#ifndef TIMEADVANCENEWMARK_H
#define TIMEADVANCENEWMARK_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <math.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <life/lifefem/timeAdvance_template.hpp>
#include <life/lifearray/EpetraVector.hpp>

namespace LifeV
{

//! Newmark  - Class to deal the \f$theta\f$-method and Newmark scheme
/*!
  @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

  This class can be used to approximate problems of the first order and the second order in time.
  In the first case the temporal scheme is a theta-method, while in the second case is a Newmark scheme.

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

     We introduce the parameters (\f$\theta\f$, \f$\gamma\f$)  and we apply the following Newmark method:

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

template<typename VectorType = EpetraVector >
class Newmark :
        public TimeAdvance < VectorType >
{
public:

    //! @name Public Types
    //@{

    typedef TimeAdvance< VectorType >                 super;
    typedef typename super::vectorVector_Type         vectorVector_Type ;
    typedef typename super::vectorVectorPtr_Type      vectorVectorPtr_Type;
    typedef typename vectorVectorPtr_Type::iterator   vectorVectorPtrIterator_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty  Constructor
    Newmark();
    
     //! Destructor
    ~Newmark() {}
  
    //@}

    //! @name Methods
    //@{

    //!Update the state vector
    /*! Update the vectors of the previous time steps by shifting on the right  the old values.
    @param solution current (new) value of the state vector
    */
    void shift_right(const VectorType& solution);
    
    //! Update the right hand side \f$ f_V \f$ of the time derivative formula
    /*!
    Set and Return the right hand side \f$ f_V \f$ of the time derivative formula
    @param dt defined the  time step need to compute the
    */
    VectorType time_der(const Real& dt = 1 );

    //! Update the right hand side \f$ f_W \f$ of the time derivative formula
    /*
    Set and Return the right hand side \f$ f_W \f$ of the time derivative formula
    @param dt defined the  time step need to compute the \f$ f_W \f$
    */
    VectorType time_derOrder2(const Real& dt = 1 );

    //!Show the properties  of temporal scheme
    void showMe() const;

    //@}

    //!@name Set Methods
    //@{

     //! Initialize the parameters of time advance scheme
     /*
     Initialize parameters of time advance scheme used in Newmark scheme
     @param  coefficients define the Newmark's coefficients
     @param  orderDerivate  define the order of derivate;
     */
     void setup (const  std::vector<Real>&  coefficients, const  UInt& orderDerivate);

    //! initialize parameters of time advance scheme;
     /*
     Initialize parameters of time advance scheme used in BDF;
     @param  order define the order of BDF;
     @param  orderDerivate  define the order of derivate;
     */
    void setup ( const UInt& order, const  UInt& orderDerivate);
    
     //! Initialize the StateVector
     /*!
     Initialize all the entries of the unknown vector to be derived with the vector x0 (duplicated).
     this class is virtual because used in bdf;
     @param x0 is the initial solution;
     */
     void initialize_unk( const VectorType& x0);

     //! initialize the state vector
     /*!
     Initialize all the entries of the unknown vector to be derived with the vector x0, v0 (duplicated).
     this class is virtual because used in \f$\theta\f$-method scheme;
     @param x0 is the initial unk;
     @param v0 is the initial velocity
     */
    void initialize_unk( const VectorType& x0, const VectorType& v0 );
   
     //! initialize the state vector
     /*!
     Initialize all the entries of the unknown vector to be derived with the vector x0, v0,w0 (duplicated).
     this class is virtual because used in Newamrk scheme;
     @param x0 is the initial solution;
     @param v0 is the initial velocity
     @param w0 is the initial accelerate
     */
    void initialize_unk(const VectorType& x0, const VectorType& v0, const VectorType& w0 );

  //! Initialize the state vector
    /*! Initialize all the entries of the unknown vector to be derived with a
     set of vectors x0
     note: this is taken as a copy (not a reference), since x0 is resized inside the method.
     */
    void initialize_unk( const vectorVector_Type& x0);
  
   //@}

  //!@name Get Methods
  //@{

   //!Return the \f$i\f$-th coefficient of the unk's extrapolation
  /*!
   @param i index of  extrapolation coefficient
   @returns beta
   */
   Real coeff_ext(const  UInt& i ) const;
  
   //! Return the \f$i\f$-th coefficient of the velocity's extrapolation
   /*!
   @param \f$i\f$ index of velocity's extrapolation  coefficient
   @returns \f$\beta^V_i\f$
   */
   Real coeff_extVelocity(const UInt& i ) const;

    //! Compute the polynomial extrapolation of solution
    /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
    */
    VectorType extrap() const;

    //! Compute the polynomial extrapolation of velocity
    /*!
    Compute the polynomial extrapolation approximation of order \f$n-1\f$ of
    \f$u^{n+1}\f$ defined by the n stored state vectors
    */
    VectorType extrapVelocity() const;
  
    //! Return the current velocity
    VectorType vnk()  const;

    //!Return the current accelerate
    VectorType wnk() const ;
  
  //@}

private:

    //! Coefficient of Newmark: \f$theta\f$
    Real M_theta;

    //! Coefficient of Newmark: \f$\gamma\f$
    Real M_gamma;

};

// ==================================================
// Constructors & Destructor
// ==================================================
template<typename VectorType>
Newmark <VectorType> ::Newmark():super()
  {
  }

// ===================================================
// Methods
// ===================================================
template<typename VectorType>
void Newmark <VectorType>::shift_right(const VectorType& solution)
{
    ASSERT (  this->M_dt != 0 ,  "M_dt must be different to 0");

    vectorVectorPtrIterator_Type it   =  this->M_unknowns.end();
    vectorVectorPtrIterator_Type itb1 =  this->M_unknowns.begin() +  static_cast<Real>(this->M_size)/2;
    vectorVectorPtrIterator_Type itb  =  this->M_unknowns.begin();

    for ( ; itb1 != it; ++itb1, ++itb)
        *itb1 = *itb;

    itb  =  this->M_unknowns.begin();

    // insert unk in unknowns[0];
    *itb = new VectorType( solution );

    itb++;

    //update velocity
    VectorType velocityTemp(solution);
    velocityTemp *=  this->M_alpha[0] / this-> M_dt;
    velocityTemp -= *this->M_rhs[0];

    // update unknows[1] with velocityTemp is current velocity
    *itb = new VectorType(velocityTemp);

    if ( this->M_orderDerivate == 2 )
    {
     itb++;

      //update accelerate
      VectorType accelerateTemp( solution );
      accelerateTemp *= this->M_xi[ 0 ] / ( this->M_dt * this->M_dt);
      accelerateTemp -= * this->M_rhs[ 1 ];

      *itb = new VectorType( accelerateTemp );
    }
    return;
}

template<typename VectorType>
VectorType
Newmark<VectorType>::time_der(const Real& dt )
{
    vectorVectorPtrIterator_Type it  =  this->M_rhs.begin();

    VectorType fv(* this->M_unknowns[0]);

    fv *=  this->M_alpha[ 1 ] / dt ;

    for (UInt i= 1; i  < this->M_firstOrderDerivateSize; ++i )
        fv += ( this->M_alpha[ i+1 ] * pow( dt, static_cast<Real>(i - 1 ) ) ) *  (* this->M_unknowns[ i ]);

    *it = new VectorType(fv);

    return fv;
}

template<typename VectorType>
VectorType
Newmark<VectorType>::time_derOrder2(const Real& dt )
{
    vectorVectorPtrIterator_Type it  =  this->M_rhs.end()-1;

    VectorType fw(* this->M_unknowns[0]);

    fw *=  this->M_xi[ 1 ] /(dt * dt) ;

    for ( UInt i = 1;  i < this->M_secondOrderDerivateSize; ++i )

      fw += ( this->M_xi[ i+1 ] * pow(dt, static_cast<Real>(i - 2) ) ) * ( *this->M_unknowns[ i ]);

    *it = new VectorType(fw);

    return fw;
}

template<typename VectorType>
void
Newmark<VectorType>::showMe() const
{
    std::cout << "*** NewmarkTime discretization maximum order of derivate "
                    << this->M_orderDerivate<< " ***"<< std::endl;
    std::cout <<" Coefficients : "      <<  std::endl;
    std::cout <<" theta :        "      << M_theta<<"\n"
              <<" gamma :  "            <<  M_gamma<<"\n"
              <<" size unknowns :"      << this->M_size<<"\n";

    for ( UInt i = 0; i <  this->M_alpha.size(); ++i )
        std::cout << "  alpha(" << i << ") = " <<  this->M_alpha[ i ]<< std::endl;

    if (this->M_orderDerivate == 2)
    {
        for ( UInt i = 0; i <  this->M_xi.size(); ++i )
            std::cout << "       xi(" << i << ") = " <<  this->M_xi[ i ] << std::endl;
    }

    std::cout <<"Delta Time : "<< this->M_dt<<"\n";
    std::cout <<"*************************************\n";
}

// ===================================================
// Set Methods
// ===================================================

template<typename VectorType>
void
Newmark<VectorType>::setup (const UInt& /*order*/, const  UInt& /*orderDerivate*/)
{
    ERROR_MSG("use setup for BDF but the time advance scheme is Newmark or  theta-method");
}

template<typename VectorType>
void
Newmark<VectorType>::setup(const std::vector<Real>& coefficients, const  UInt& orderDerivate)
{
    //initialize theta
    M_theta = coefficients[0];

    //initilialize gamma
    M_gamma = coefficients[1];

    //initialize Order Derivate
    this->M_orderDerivate= orderDerivate;

    // If theta equal 0, explicit meta method
    if (M_theta == 0)
    {
        ASSERT (this->M_orderDerivate == 2,  "theta is 0 must be different from 0 in Newmark");
        this->M_size = 4;
        this->M_alpha[ 0 ] =  1;
        this->M_alpha[ 1 ] =  1;
        this->M_alpha[ 2 ] =  1;
        this->M_order  = 1;
    }
    else
    {
        if (this->M_orderDerivate == 1 )  // Theta method
        {
            this->M_theta = 1;
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
        else     //NewmarkMethod
        {
            //unknown vector's dimension
            this->M_size = 6 ;
            this->M_alpha.resize(4);
            this->M_xi.resize(4);
            this->M_beta.resize(3);
            this->M_betaVelocity.resize(3);
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
            this->M_betaVelocity[ 0 ] = 0;
            this->M_betaVelocity[ 1 ] = 1;
            this->M_betaVelocity[ 2 ] = 1;

            this->M_coefficientsSize  = 4;
        }
        this->M_unknowns.reserve(this->M_size);
        this-> M_rhs.reserve(2);
        // check order  scheme
        if (this->M_alpha[0] == 0.5)
            this->M_order = 2;
        else
            this->M_order = 1;

        this->M_firstOrderDerivateSize  =  static_cast<Real>(this->M_size) / 2.0;
        this->M_secondOrderDerivateSize = static_cast<Real>(this->M_size) / 2.0;
    }
}

template<typename VectorType>
void Newmark<VectorType>::initialize_unk( const VectorType& x0)
{
    vectorVectorPtrIterator_Type iter     = this->M_unknowns.begin();
    vectorVectorPtrIterator_Type iter_end = this->M_unknowns.end();

    for ( ; iter != iter_end; ++iter )
    {
        delete *iter;
        *iter = new VectorType(x0.getMap(), Unique);
    }

    for ( UInt i(this->M_unknowns.size()) ; i < this->M_size; ++i )
        this->M_unknowns.push_back(new VectorType(x0));

    vectorVectorPtrIterator_Type iterRhs     = this->M_rhs.begin();
    vectorVectorPtrIterator_Type iterRhs_end = this->M_rhs.end();

    for ( ; iterRhs !=iterRhs_end ; ++iterRhs)
    {
        delete *iterRhs;
        *iterRhs = new VectorType(x0.getMap(), Unique);
    }
}

template<typename VectorType>
void Newmark<VectorType>::initialize_unk( const VectorType& x0, const VectorType& v0)
{
    vectorVectorPtrIterator_Type iter       = this->M_unknowns.begin();
    vectorVectorPtrIterator_Type iter_end   = this->M_unknowns.end();

    for ( ; iter != iter_end; ++iter )
    {
        delete *iter;
        *iter = new VectorType(x0.getMap(), Unique);
    }

    this->M_unknowns.push_back(new VectorType(x0));
    this->M_unknowns.push_back(new VectorType(v0));

    for (UInt i=0; i<4; ++i )
    {
        this->M_unknowns.push_back(new VectorType(x0.getMap(), Unique));
        *this->M_unknowns[i+2]  *=0;
    }
    for (UInt i=0; i<2; ++i )
    {
        this->M_rhs.push_back(new VectorType(x0.getMap(), Unique));
        *this->M_rhs[i]  *=0;
    }
}

template<typename VectorType>
void Newmark<VectorType>::initialize_unk( const VectorType& x0, const VectorType& v0, const VectorType& w0)
{
    vectorVectorPtrIterator_Type iter       = this->M_unknowns.begin();
    vectorVectorPtrIterator_Type iter_end   = this->M_unknowns.end();

    for ( ; iter != iter_end; ++iter )
    {
        delete *iter;
        *iter = new VectorType(x0.getMap(), Unique);
    }

    this->M_unknowns.push_back(new VectorType(x0));
    this->M_unknowns.push_back(new VectorType(v0));
    this->M_unknowns.push_back(new VectorType(w0));

    for (UInt i=0; i<3; ++i )
    {
        this->M_unknowns.push_back(new VectorType(x0.getMap(), Unique));
        *this->M_unknowns[i+3] *=0;
    }
    for (UInt i=0; i<2; ++i )
    {
        this->M_rhs.push_back(new VectorType(x0.getMap(), Unique));
        *this->M_rhs[i] *=0;
    }
}

template<typename VectorType>
void Newmark<VectorType>::initialize_unk( const vectorVector_Type& x0)
{
    const UInt n0 = x0.size();

    ASSERT( n0 != 0, "vector null " );

    vectorVectorPtrIterator_Type iter     = this->M_unknowns.begin();
    vectorVectorPtrIterator_Type iter_end = this->M_unknowns.end();

    UInt i(0);

    for ( ; iter != iter_end && i< n0 ; ++iter, ++i )
    {
        delete *iter;
        *iter = new VectorType(x0[i]);
    }

    for ( i = this->M_unknowns.size() ; i < this->M_size && i< n0; ++i )
        this->M_unknowns.push_back(new VectorType(x0[i]));

    for ( i = this->M_unknowns.size() ; i < this->M_size; i++ )
        this->M_unknowns.push_back(new VectorType( x0[ n0-1 ] ) );

    for (UInt i=0; i<2; i++ )
    {
        this->M_rhs.push_back(new VectorType(x0[0].getMap(), Unique));
        *this->M_rhs[i] *=0;
    }
}

// ===================================================
// Get Methods
// ===================================================

template<typename VectorType>
Real
Newmark<VectorType>::coeff_ext(const UInt& i) const
{
  ASSERT ( i <  3 ,  "coeff_der i must equal 0 or 1 because U^*= U^n + dt*V^n + dt^2 / 2 W^n");
  return  this->M_beta(i)*pow( this->M_dt, static_cast<Real>(i) );
}

template<typename VectorType>
Real
Newmark<VectorType>::coeff_extVelocity(const UInt& i ) const
{
ASSERT (i < 3,  "coeff_der i must equal 0 or 1 because V^* =  V^n + dt W^n");
 return  this->M_betaVelocity(i)*pow( this->M_dt, static_cast<Real>(i));
}

template<typename VectorType>
VectorType
Newmark<VectorType>::extrap()  const
{
    VectorType extrapolation(*this->M_unknowns[0]);
    extrapolation += this->M_dt * ( *this->M_unknowns[ 1 ]);
    if ( this->M_orderDerivate == 2 )
        extrapolation += ( this->M_dt * this->M_dt ) / 2.0 * ( *this->M_unknowns[2]);
    return extrapolation;
}

template<typename VectorType>
VectorType
Newmark<VectorType>::extrapVelocity() const
{
    VectorType extrapolation( *this->M_unknowns[1]);
    extrapolation += this->M_dt * ( *this->M_unknowns[ 2 ]);

    return extrapolation;
}

template<typename VectorType>
VectorType
Newmark<VectorType>::vnk() const
{
    VectorType velocity( *this->M_unknowns[1]);
    return velocity;
}

template<typename VectorType>
VectorType
Newmark<VectorType>::wnk() const
{
    VectorType accelerate( *this->M_unknowns[2]);
    return accelerate;
}

// ===================================================
// Macros
// ===================================================

//! define the Newmark factory
inline
TimeAdvance< EpetraVector >* createNewmark() { return new Newmark<EpetraVector>(); }

namespace
{
static bool registerNewmark= TimeAdvanceFactory::instance().registerProduct( "Newmark",  &createNewmark);
}

}
#endif  /*TIMEADVANCENEWMARK*/
