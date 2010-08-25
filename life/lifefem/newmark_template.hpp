/*
/ This file is part of the LifeV library

 Authors: Simone Deparis <simone.deparis@epfl.ch>

 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file bdf_template.hpp
  \author M. Pozzoli, F. Nobile, C .Vergara

  File containing a class for an easy handling of different order time
  discretizations/extrapolations Newmark based

*/
#ifndef _NEWMARK_TEMPLATE_H
#define _NEWMARK_TEMPLATE_H
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <life/lifefem/timeAdvance_template.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <math.h>

namespace LifeV
{
  
/*!
\class Newmark (templated)

This class can be used to appromate problems of the first order and the second order in time.
In the first case the temporal scheme is a theta-method, while in the second case is a Newmark scheme.

This class defines the state vector X^{n+1}, a suitible  extrapolation of vector X^* and  
opportune coefficients used to determinate X^{n+1} and X^*.

1) First order problem:

 \f$ M u' + A(u) = f \f$

the state vector is \f$X^{n+1} = (U^{n+1},V^{n+1}, U^{n}, V^{n})\f$. 
where U is an approximation of u and V  of u'.

We consider a parameter \theta, and we apply the following theta method: 

 \f$ U^{n+1} = U^{n} + \Delta t  \theta V^{n+1} + (1 − \theta) V^n  \f$,

so the approximation of velocity at timestep n+1 is:

\f$  V^{ n+1} = 1/ (\theta * \Delta t) (U^{n+1}-U^n)+( 1 -  1 / \theta) V^n\f$;

We can  linearize non-linear term A(u)  as A(U^*) U^{n+1},  where the U^* becomes: 

\f$ U^∗ = U^n + \Delta t V^n \f$.

The coefficients alpha, beta depend on theta, 
and we can use a explicit theta-method.

2) Second order Method:

 \f$ M u'' + D(u, u')+ A(u) = f \f$

the state vector is \f$X^{n+1} = (U^{n+1},V^{n+1}, W^{n+1}, U^{n}, V^{n}, W^{n})\f$. 
where U is an approximation of u and V  of u' and W of u'' .

We introduce the parameters (\gamma, \theta)  and we apply the following newmark method: 

 \f$ U^{n+1} = U^{n} + \Delta t  \theta V^{n+1} + (1 − \theta) V^n  \f$,

so the approximation of velocity at timestep n+1 is

\f$  V^{ n+1} = \gamma / (\theta * \Delta t) (U^{n+1}-U^n)+( 1 - gamma \theta) V^n+ (1- \gamma /(2\theta)) \Delta t W^n\f$;

and we determine W^{n+1} by

\f$ W^{n+1} =1/(\theta \Delta t^2) (U^{n+1}-U^{n}) - 1 / \Delta t V^n + (1-1/(2\theta)) W^n $\f.

We can  linearize non-linear term A(u)  as A(U^*) U^{n+1} and  D(u, u') as D(U^*, V^*) V^{n+1}, 
where U^* is given by

\f$ U^∗ = U^n + \Delta t V^n \f$.

and  V^* is

\f$  V^* =V^n+ \Delta t W^n  \f$

The coefficients \xi, \alpha, \beta, \beta2 depend on theta and \gamma.

*/

template<typename VectorType = EpetraVector >
class Newmark :
 public TimeAdvance< > 
{
public: 
  
  /** @name Typedefs
   */
  //@{
  typedef VectorType                      vector_raw_type;
  typedef std::vector< vector_raw_type* > vector_type;
  typedef typename vector_type::iterator  vector_type_iterator;
  
  // typedef SolverType                    solver_type;
  
  // typedef typename solver_type::matrix_type      matrix_type;

  //@}
  
  
  /*! Constructor
   *  @param n order of the Newmark
   */

  Newmark();

  // Newmark( const UInt orderDev,  std::vector<double> coefficients);
  
  ~Newmark() {}

  //! Initialize all the entries of the unknown vector to be derived with the
  //! vector u0 (duplicated)
  void initialize_unk( VectorType u0 );
  
  void initialize_unk( VectorType u0, VectorType v0 );
  
  void initialize_unk( VectorType u0, VectorType v0, VectorType const & w0 );
  
  //! Initialize all the entries of the unknown vector to be derived with a
  //! set of vectors uv0
  //! note: this is taken as a copy (not a reference), since uv0 is resized inside the method.
  void initialize_unk( std::vector<VectorType> uv0 );
  
  //! initialize parameters of time advance scheme
  void setup (const  std::vector<double>   coefficients, const  UInt orderDev);
  
  //! initialize parameters of time advance scheme
  void setup ( const UInt order, const  UInt orderDev);

  /*! Update the vectors of the previous time steps by shifting on the right
   *  the old values.
   *  @param u_curr current (new) value of the state vector
   */
  void shift_right(const VectorType & u_curr );
  
  //! Returns the right hand side \f$ f_V \f$ associated to the discretitation      of the first derivative
   VectorType time_der( Real dt = 1 ) /*const*/ ;
  
  //! Returns the right hand side \f$ f_W \f$ associated to the discretitation      of the second derivative
  VectorType time_derOrder2( Real dt = 1 ) const ;
  
  //! Compute the polynomial extrapolation approximation of order n-1 of
  //! u^{n+1} defined by the n stored state vectors
  VectorType extrap() const;
  
 //! Compute the polynomial extrapolation approximation of order n-1 of
  //! u^{n+1} defined by the n stored state vectors
  VectorType extrapVelocity() const;

  //! Return the i-th coefficient of the time derivative alpha_i
  double coeff_der( UInt i )const;
  
  //! Return the i-th coefficient of the time derivative xi_i
  double coeff_derOrder2( UInt i ) const;
  
  //! Return the i-th coefficient of the time extrapolation beta_i
  double coeff_ext( UInt i ) const ;

  //! Return the i-th coefficient of the time extrapolation betaV_i
  double coeff_extVelocity( UInt i ) const;

  // const vector_type& unk() const;
  // const VectorType unk(UInt i) const;

  VectorType vnk() const;

  VectorType wnk() const ;

  void showMe() const;

  private:
  
  //Coefficients of Newmark
  Real _M_theta;

  Real _M_gamma;

};
  

///
// template implementations
//

template<typename VectorType>
Newmark<VectorType>::Newmark() : 
  super()
{
  CONSTRUCTOR("Newmark");
}

template<typename VectorType>
void
Newmark<VectorType>::setup (const UInt /*order*/, const  UInt /*orderDev*/) 
{
 ERROR_MSG("use setup for BDF but the time advance scheme is Newmark or theta-method");
}

template<typename VectorType>
void
Newmark<VectorType>:: setup (const   std::vector<double> coefficients, const  UInt orderDev)
{
    //initialize theta
  _M_theta = coefficients[0];
  //initilialize gamma
  _M_gamma = coefficients[1];
  //initialize Order Derivate
  this->_M_orderDev = orderDev;
  
  // If theta equal 0, explicit meta method
  if (_M_theta == 0)
    {
      ASSERT (this->_M_orderDev==2,  "theta is 0 must be different from 0 in Newmark");
      this->_M_size = 4;
      this->_M_alpha[ 0 ] =  1;
      this->_M_alpha[ 1 ] =  1;
      this->_M_alpha[ 2 ] =  1;  
      this->_M_order  = 1;
    }
  else
    {  
      if (this->_M_orderDev == 1 )  // Theta method
	{ 
	  this->_M_gamma = 1; 
	  //  unknown vector's  dimension;
	  this->_M_size = 4;
	  this->_M_alpha.resize(3);
	  this->_M_xi.resize(3);
	  this->_M_beta.resize(3);
	  this->_M_alpha[ 0 ] =  _M_gamma/ _M_theta;
	  this->_M_alpha[ 1 ] =   _M_gamma/ _M_theta;
	  this->_M_alpha[ 2 ] =   _M_gamma/_M_theta - 1.0;  
	  this->_M_beta[0] =1;
	  this->_M_beta[1]=1;
	  this->_M_beta[2]=0.5;
	  this->_M_xi[0] =0;
	  this->_M_xi[1]=0;
	  this->_M_xi[2]=0;
	  this->_M_sizeCoefficients=3;
	}
      else     //Newmark Method
	{
	  //unknown vector's dimension
	  this->_M_size = 6 ;
	  this->_M_alpha.resize(4);
	  this->_M_xi.resize(4);
	  this->_M_beta.resize(3);
	  this->_M_beta2.resize(3);
	  //initialitation alpha coefficients
	 this->_M_alpha[ 0 ] =  _M_gamma/_M_theta;
	 this->_M_alpha[ 1 ] =  _M_gamma/_M_theta;
	 this->_M_alpha[ 2 ] =  _M_gamma/_M_theta - 1.0; 
	 this->_M_alpha[ 3 ] = _M_gamma/(2.0 * _M_theta) -1.0;
	 
	 //initialitation xi coefficients
	 this->_M_xi[ 0 ] =  1. /_M_theta;
	 this->_M_xi[ 1 ] =  1. /_M_theta;
	 this->_M_xi[ 2 ] =  1. /_M_theta; 
	 this->_M_xi[ 3 ] =  1. / (2.*_M_theta)-1.0; 
	 

	  //initialitation extrapolation coefficients
	  this->_M_beta[0] =1;
	  this->_M_beta[1]=1;
	  this->_M_beta[2]=0.5;
	  this->_M_beta2[0]=0;
	  this->_M_beta2[1]=1;
	  this->_M_beta2[2]=1;
	  
	 this->_M_sizeCoefficients=4;
	}  
      this->_M_unknowns.reserve(this->_M_size);
      this-> _M_rhs.reserve(2);
      // check order  scheme
      if(this->_M_alpha[0] == 0.5) 
	this->_M_order = 2;
      else
	this->_M_order = 1;

      this->_M_sizeTimeDer=this->_M_size / 2.0;
      this->_M_sizeTimeDer2=this->_M_size / 2.0;

    }
}

template<typename VectorType>
void Newmark<VectorType>::initialize_unk( VectorType u0 )
{
    vector_type_iterator iter     = this->_M_unknowns.begin();
    vector_type_iterator iter_end = this->_M_unknowns.end();

    for ( ; iter != iter_end; iter++ )
      {
	delete *iter;
	*iter = new VectorType(u0.getMap(), Unique);
      }

    for ( UInt i(this->_M_unknowns.size()) ; i < _M_size; i++ )
        this->_M_unknowns.push_back(new VectorType(u0));

    vector_type_iterator iterR     = this->_M_rhs.begin();
    vector_type_iterator iterR_end = this->_M_rhs.end();
    
    for ( ; iterR !=iterR_end ; iterR++ )
      {
	delete *iterR;
	*iterR = new VectorType(u0.getMap(), Unique);
      }
     return ;
}


template<typename VectorType>
void Newmark<VectorType>::initialize_unk( VectorType u0, VectorType v0 )
{
  vector_type_iterator iter       = _M_unknowns.begin();
  vector_type_iterator iter_end   = _M_unknowns.end();
  
  for ( ; iter != iter_end; iter++ )
    {
      delete *iter;
      *iter = new VectorType(u0.getMap(), Unique);
    }
  
  this->_M_unknowns.push_back(new VectorType(u0));
  this->_M_unknowns.push_back(new VectorType(v0));
  
  for (UInt ii=0; ii<4; ii++ )
    {    
      this->_M_unknowns.push_back(new VectorType(u0.getMap(), Unique));
      *this->_M_unknowns[ii+2]  *=0;
    }
  for (UInt i=0; i<2; i++ ) 
    {
      this->_M_rhs.push_back(new VectorType(u0.getMap(), Unique));
      *this->_M_rhs[i]  *=0; 
    }
  return ;
}


template<typename VectorType>
void Newmark<VectorType>::initialize_unk( VectorType u0, VectorType v0, VectorType const & w0 )
{
  vector_type_iterator iter       = this->_M_unknowns.begin();
  vector_type_iterator iter_end   = this->_M_unknowns.end();
  
  for ( ; iter != iter_end; iter++ )
    {
      delete *iter;
      *iter = new VectorType(u0.getMap(), Unique);
    }
  
  this->_M_unknowns.push_back(new VectorType(u0));
  this->_M_unknowns.push_back(new VectorType(v0));
  this->_M_unknowns.push_back(new VectorType(w0));

  for (UInt ii=0; ii<3; ii++ )
    {
      this->_M_unknowns.push_back(new VectorType(u0.getMap(), Unique));
     *this->_M_unknowns[ii+3] *=0;
    }
  for (UInt i=0; i<2; i++ ) 
    {
      this->_M_rhs.push_back(new VectorType(u0.getMap(), Unique));
      *this->_M_rhs[i] *=0; 
    }
  
  return ;
}

template<typename VectorType>
void Newmark<VectorType>::initialize_unk(const  std::vector<VectorType> uv0 )
{
  UInt n0 = uv0.size();
  
  ASSERT( n0 != 0, "vector null " );
  
  vector_type_iterator iter     = this->_M_unknowns.begin();
  vector_type_iterator iter_end =this->_M_unknowns.end();
  
  int i(0);
    
  for ( ; iter != iter_end && i< n0 ; ++iter, ++i )
    {
      delete *iter;
      *iter = new VectorType(uv0[i]);
    }
    
  for ( i = this->_M_unknowns.size() ; i < this->_M_size && i< n0; ++i )
    this->_M_unknowns.push_back(new VectorType(uv0[i]));
  
  for ( i = this->_M_unknowns.size() ; i < this->_M_size; i++ )
    this->_M_unknowns.push_back(new VectorType(uv0[n0-1]));

   for (UInt i=0; i<2; i++ ) 
    {
      this->_M_rhs.push_back(new VectorType(uv0[0].getMap(), Unique));
      *this->_M_rhs[i] *=0; 
    }
  

return ;
}

template<typename VectorType>
void
Newmark<VectorType>::showMe()  const
{
  std::cout << "*** Newmark Time discretization maximum order of derivate "<< this->_M_orderDev << " ***"<< std::endl;
  std::cout <<" Coefficients : "  << std::endl;
  std::cout <<" theta :        " <<_M_theta<<"\n"
	    <<" gamma :        " <<_M_gamma<<"\n"
	    <<" size unknowns :" << this->_M_size<<"\n";
  
    for ( UInt i = 0; i <  this->_M_alpha.size(); ++i )
       	 std::cout << "       alpha(" << i << ") = " <<  this->_M_alpha[ i ]
		   << std::endl;

	 if (_M_orderDev==2)
	   {
	     for ( UInt i = 0; i <  this->_M_xi.size(); ++i )
	       std::cout << "       xi(" << i << ") = " <<  this->_M_xi[ i ] << std::endl;
	   }

  std::cout <<"Delta T :"<< this->_M_dt<<"\n";
  std::cout <<"*************************************\n";

	 return ;
}

template<typename VectorType>
void
Newmark<VectorType>::
shift_right(const VectorType & unk_curr) 
{
  ASSERT ( _M_dt != 0 ,  "_M_dt must be different to 0");
  
   vector_type_iterator it   =  this->_M_unknowns.end();
  vector_type_iterator itb1 =  this->_M_unknowns.begin() +  this->_M_size/2;
  vector_type_iterator itb  =  this->_M_unknowns.begin();
 
  VectorType u_temp(unk_curr);
 
  for ( ; itb1 != it; itb1++, itb++ )
    *itb1 = *itb;  
  
  itb  =  this->_M_unknowns.begin();
 
  // insert unk_curr in unknows[0];
  *itb = new VectorType( unk_curr);
  
  itb++;
  
  //update vt;
   
   VectorType v_temp(unk_curr);
   v_temp *= _M_alpha[0] / _M_dt;
   v_temp -= *_M_rhs[0];

   // update unknows[1] with u_temp 
   // where u_temp is current velocity
  *itb = new VectorType(v_temp);
  
  if ( this->_M_orderDev == 2 ) 
    {
      itb++;
      u_temp *= 0;
      u_temp += unk_curr;
      
      //update wt;
      u_temp *= _M_xi[0] / ( this->_M_dt* this->_M_dt);
      u_temp -= * this->_M_rhs[1];   

      *itb = new VectorType(u_temp);
    }
  return;
}

template<typename VectorType>
VectorType
Newmark<VectorType>::time_der( Real dt ) /*const*/
{
  vector_type_iterator itb  =  this->_M_rhs.begin();

  VectorType fv(* this->_M_unknowns[0]);

  fv *=  this->_M_alpha[ 1 ] / dt ;
  
  for (int i= 1; i  < this->_M_size / 2.0 ; i++ )
      fv += ( this->_M_alpha[ i+1 ] * pow( dt, i - 1 )) *  (* this->_M_unknowns[ i ]);
   
  *it = new VectorType(fv);

   return fv;
}

template<typename VectorType>
VectorType
Newmark<VectorType>::time_derOrder2( Real dt )  /*const*/
{
  
  vector_type_iterator it  =  this->_M_rhs.end()-1;

  VectorType fw(* this->_M_unknowns[0]);
 
  fw *=  this->_M_xi[ 1 ] /(dt * dt) ;
  
  for ( int i = 1;  i < _M_size / 2.0; ++i )
    
      fw += ( this->_M_xi[ i+1 ] * pow(dt, i - 2 ) ) * ( *this->_M_unknowns[ i ]);
    
   *it = new VectorType(fw);
  
  return fw;
}

template<typename VectorType>
VectorType 
Newmark<VectorType>::vnk()  const
{
  VectorType v(*_M_unknowns[1]);
  return v;
}

template<typename VectorType>
VectorType 
Newmark<VectorType>::wnk()  const
{
  VectorType w(*_M_unknowns[2]);
  return w;
}


template<typename VectorType>
VectorType 
Newmark<VectorType>::extrap()  const
{
  VectorType ue(*_M_unknowns[0]);
    ue += this->_M_dt * (* this->_M_unknowns[ 1 ]);
    if ( this->_M_orderDev ==2 )
      ue += ( this->_M_dt* this->_M_dt) *(* this->_M_unknowns[2]);
  return ue;
}

template<typename VectorType>
VectorType 
Newmark<VectorType>::extrapVelocity()  const
{
  VectorType ve(* this->_M_unknowns[1]);
    ve += this->_M_dt * (* this->_M_unknowns[ 2 ]);
}

template<typename VectorType>
double
Newmark<VectorType>::coeff_ext(UInt i) const
{
  ASSERT ( i <  3 ,  "coeff_der i must equal 0 or 1 because U^*= U^n + dt*V^n + dt^2 / 2 W^n");
  return  this->_M_beta(i)*pow( this->_M_dt, (double) i);
}

template<typename VectorType>
double
Newmark<VectorType>::coeff_extVelocity(UInt i)  const
{
ASSERT (i < 3,  "coeff_der i must equal 0 or 1 because V^* =  V^n + dt W^n");
 return  this->_M_beta2(i)*pow( this->_M_dt, (double) i);
}

inline
TimeAdvance< >* createNewmark() { return new Newmark <>(); }
  namespace
  {
    static bool registerNewmark = TimeAdvanceFactory::instance().registerProduct( "Newmark",  &createNewmark );
  }



}  
#endif
