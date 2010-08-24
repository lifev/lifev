/*
 This file is part of the LifeV library

 Authors: Simone Deparis <simone.deparis@epfl.ch>

 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
  \author S. Deparis

  File containing a class for an easy handling of different order time
  discretizations/extrapolations BDF based

*/
#ifndef _BDF_TEMPLATE_H
#define _BDF_TEMPLATE_H
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <life/lifearray/EpetraVector.hpp>

#include <life/lifefem/timeAdvance_template.hpp>

namespace LifeV
{
const UInt BDFT_MAX_ORDER = 4;

//typedef Real ( *Funct ) ( const Real&, const Real&, const Real&, const Real&,
//                          const ID& );

/*!
  \class BdfT (templated)
  \brief Backward differencing formula time discretization for the first 
    and the second order problem in time.

 1) First order problem

  A differential equation of the form

  \f$ M u' + A u = f \f$

  is discretized in time as

  \f$ M V_{k+1}  + A U_{k+1} = f_{k+1} \f$

  where V denotes the polynomial of order n in t that interpolates
  (t_i,u_i) for i = k-n+1,...,k+1.

  The approximative time derivative \f$ V \f$ is a linear
  combination of state vectors u_i:

  \f$ V_{k+1} = \frac{1}{\Delta t} (\alpha_0 U_{k+1} - \sum_{i=0}^n \alpha_i U_{k+1-i} )\f$

  Thus we have

  \f$ \frac{\alpha_0}{\Delta t} M U_{k+1} = A U_{k+1} + f + M f_V$

  with

  \f$ f_V= \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i U_{k+1-i} \f$

  This class stores the n last state vectors in order to be able to
  calculate \f$ f_V \f$. It also provides alpha_i
  and can extrapolate the the new state from the n last states with a
  polynomial of order n-1:

  \f$ U_{k+1} \approx \sum_{i=0}^{n-1} \beta_i U_{k-i} \f$

2)  Second order problem

  A differential equation of the form

  \f$ M u'' = D( u, u') + A(u) + f \f$

  is discretized in time as

  \f$ M W_{k+1} = A(U^*) U_{k+1} +D(U^*, V^*) V_{k+1} + f_{k+1} \f$

  where W and V denotes the polynomial of order n+1 and order n, 
  while U^* and V^* are suitable extrapolations.  

  The velocity vector, as for first order  problem, is:

  \f$ V_{k+1} = \frac{1}{\Delta t} (\alpha_0 U_{k+1} - \sum_{i=0}^n \alpha_i U_{k+1-i} )\f$

  while the acceleration vector W^{n+1} is:

  \f$ W_{k+1} = \frac{1}{\Delta t^2} (\xi_0 U_{k+1} - \sum_{i=0}^{n+1} \xi_i U_{k+1-i} )\f$
 
  Thus we have

  \f$ \frac{\xi_0}{\Delta t^2} M U_{k+1} + \frac{\alpha_0}{\Delta t } D U_{k+1} + A U_{k+1} + f + M f_W + D f_V $

  with

  \f$ f_V= \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i U_{k+1-i} \f$

  and

  \f$  f_W =\frac{1}{\Delta t^2} \sum_{i=1}^{n+1} \xi_i U_{ k+1 - i}  \f$ 

  It also provides alpha_i  and can extrapolate the the new state from the n last states with a
  polynomial of order n-1:

  \f$ U_{k+1} \approx  U^*= \sum_{i=0}^{n-1} \beta_i U_{k-i} \f$

  and  V^* in following way:

  \f$  V_{k+1} \approx V^*=\sum_{i=0}^p \frac{\beta_i^V}{\Delta t}\Ubf^{n-i} = \Wbf^{n+1}+O(\Delta t^p),  \f$
*/

template<typename VectorType = EpetraVector >
class BdfT:  // T means template
    public  TimeAdvance < >
{
public:
  
    /** @name Typedefs
     */
    //@{
   typedef VectorType                      vector_raw_type;
   typedef std::vector< vector_raw_type* > vector_type;
   typedef typename vector_type::iterator  vector_type_iterator;
    //@}

    /*! Constructor
     *  @param n order of the BDF
     */
    BdfT();
   
    BdfT( const UInt order);

    BdfT( const UInt order, const  UInt orderDev );

    ~BdfT() {} 

    //! Initialize all the entries of the unknown vector to be derived with the
    //! vector u0 (duplicated)
    void initialize_unk( VectorType u0 );

    //! Initialize all the entries of the unknown vector to be derived with a
    //! set of vectors uv0
    //! note: this is taken as a copy (not a reference), since uv0 is resized inside the method.
    void initialize_unk( std::vector<VectorType> uv0 );

    void initialize_unk( VectorType u0, VectorType v0 );
  
    void initialize_unk( VectorType u0, VectorType v0, VectorType const & w0 );

    /* Initialize all the entries of the unknonwn vectors with a given function
        The array of initial conditions needed by the selected BDF is
        initialized as follows: _unk=[ u0(t0), u0(t0-dt), u0(t0-2*dt), ...]
        For the space dependence of the initial conditions we need informations
        on:
        -# the mesh (coordinates of points)
        -# the reference and the current FE (basis functions)
        -# is it a vector or a scalar problem ? bdf doesn't know it
        -# which is the initial time (t0) and the time step (for solutions
           before the initial instant)
        Based on the NavierStokesHandler::initialize by M. Fernandez
    */
/*    template <typename Mesh, typename RefFE, typename CurrFE, typename Dof>
    void initialize_unk( const Funct& u0, Mesh& mesh, RefFE& refFE, CurrFE& currFE,
                         Dof& dof, Real t0, Real dt, UInt nbComp );
*/

     //! initialize parameters of time advance scheme
     void setup ( const UInt order, const  UInt orderDev);

     void setup ( const  std::vector<double>  coefficients, const  UInt orderDev);

    /*! Update the vectors of the previous time steps by shifting on the right
     *  the old values.
     *  @param u_curr current (new) value of the state vector
     */
     void shift_right( VectorType const& u_curr );

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula
    VectorType time_der( Real dt = 1 ) const;

    //! Returns the right hand side \f_v$ \bar{p} \f$ of the time derivative
    //! formula
    VectorType time_derOrder2( Real dt = 1 ) const ;

    //! Returns the right hand side \f$ \bar{p} \Delta t \f$ of the time
    //! derivative formula. The timestep is taken into account elsewhere,
    //! e. g. in the mass matrix.
    //    VectorType time_der() const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    VectorType extrap() const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! v^{n+1} defined by the n stored state vectors
    VectorType extrapVelocity() const;
  
    //! Return the i-th coefficient of the time derivative alpha_i
  // double coeff_der( UInt i ) const;

    //! Return the i-th coefficient of the time derivative xi_i
    //double coeff_derOrder2( UInt i ) const;
    
  //! Return the i-th coefficient of the time extrapolation beta_i
    double coeff_ext( UInt i ) const;

    //! Return the i-th coefficient of the time derivative xi_i
    double coeff_extVelocity( UInt i ) const;

    //! Return a vector with the last n state vectors
    const vector_type& unk() const;

    //! Return a last velocity 
    VectorType vnk() const;
  
    //! Return a last acceleration 
    VectorType wnk() const;

    void showMe() const;

};

///
// template implementations
//

template<typename VectorType>
BdfT <VectorType> :: BdfT() :
//:  TimeAdvance()
  TimeAdvance<>::TimeAdvance( )
{
    CONSTRUCTOR( "BDF" );
}

 
template<typename VectorType>
BdfT<VectorType>::
BdfT( const UInt order )
  :
  TimeAdvance<> ::TimeAdvance()
{
    if ( order <= 0 || order > BDFT_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
        << " you want to use BDF order " << order << "\n"
        << " we support BDF order from 1 to " << BDFT_MAX_ORDER << "\n";
        throw std::invalid_argument( __ex.str() );
    }
    _M_alpha.resize( order +1 );
    _M_beta.resize( order );
    _M_size = order ;
    _M_sizeTimeDer = order;
   _M_sizeTimeDer2 = order ;
   _M_sizeCoefficients = order +1 ;
    switch ( order )
    {
    case 1:
        _M_alpha[ 0 ] = 1.; // Backward Euler
        _M_alpha[ 1 ] = 1.;
        _M_beta[ 0 ] = 1.; // u^{n+1} \approx u^n
        break;
    case 2:
        _M_alpha[ 0 ] = 3. / 2.;
        _M_alpha[ 1 ] = 2.;
        _M_alpha[ 2 ] = -1. / 2.;
        _M_beta[ 0 ] = 2.;
        _M_beta[ 1 ] = -1.;
        break;
    case 3:
        _M_alpha[ 0 ] = 11. / 6.;
        _M_alpha[ 1 ] = 3.;
        _M_alpha[ 2 ] = -3. / 2.;
        _M_alpha[ 3 ] = 1. / 3.;
        _M_beta[ 0 ] = 3.;
        _M_beta[ 1 ] = -3.;
        _M_beta[ 2 ] = 1.;
        break;
 case 4:
        _M_alpha[ 0 ] = 25. / 12.;
        _M_alpha[ 1 ] = 4.;
        _M_alpha[ 2 ] = -3. ;
        _M_alpha[ 3 ] = 4. / 3.;
	_M_alpha[ 4 ] = - 1 / 4.;
        _M_beta[ 0 ] = 4.;
        _M_beta[ 1 ] = -6.;
        _M_beta[ 2 ] = 4.;
	_M_beta[ 3 ] = 4.;
    }
    _M_unknowns.reserve(_M_size );
 }

template<typename VectorType>
BdfT <VectorType> ::
BdfT(const  UInt order, const UInt orderDev )
  :
  BdfT( order )
{
  _M_sizeCoefficients = order +orderDev ;
  if(_M_orderDev ==2)
    {
      _M_xi.resize( order + 2 );
      _M_beta2.resize( order + 1 );
    
      switch ( order )
	{
	case 1:
	  _M_xi[ 0 ] = 1.0;
	  _M_xi[ 1 ] = 2.0;
	  _M_xi[ 2 ] = -1.0;
	  _M_beta2[0]=2.0;
	  _M_beta2[1]=-1.0;
	  break;
	case 2:
	  _M_xi[ 0 ] =  2.0;
	  _M_xi[ 1 ] = 5.0;
	  _M_xi[ 2 ] = -4.0;
	  _M_xi[ 3 ] = 1.0;
	  _M_beta2[0] = 3.0;
	  _M_beta2[ 1 ] = -3.;
	  _M_beta2[ 2 ] = 1.;
	  break;
	case 3:
	  _M_beta[ 0 ] = 3.;
	  _M_beta[ 1 ] = -3.;
	  _M_beta[ 2 ] = 1.;
	  _M_xi[ 0 ] =  35./12.;
	  _M_xi[ 1 ] =  26./3.;
	  _M_xi[ 2 ] =  -19./2;
	  _M_xi[ 3 ] =  14./3.;
	  _M_xi[ 4 ] =  -11./12.;
	  _M_beta2[ 0 ] = 4.;
	  _M_beta2[ 1 ] = -6.;
	  _M_beta2[ 2 ] = 4.;
	  _M_beta2[ 3 ] = -1;
	  break;
	}
    }
  _M_unknowns.reserve(_M_size);
} 
  

template<typename VectorType>
void
BdfT<VectorType>::setup( const UInt order, const UInt orderDev)
{
  if ( order <= 0 || order > BDFT_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
        << " you want to use BDF order " << order << "\n"
        << " we support BDF order from 1 to " << BDFT_MAX_ORDER << "\n";
        throw std::invalid_argument( __ex.str() );
    }

  _M_order = order ;
  _M_orderDev = orderDev ;
  _M_size = order ;
  _M_alpha.resize( order + 1 );
  _M_xi.resize( order + 2 );
  _M_beta.resize( order );
  _M_beta2.resize( order + 1 );
  _M_sizeCoefficients = order + orderDev ; 

    switch ( order )
      {  
      case 1:
	_M_alpha[ 0 ] = 1.; // Backward Euler
	_M_alpha[ 1 ] = 1.;
	_M_beta[ 0 ] = 1.; // u^{n+1} \approx u^n
	_M_xi[ 0 ] = 1.;
	_M_xi[ 1 ] = 2;
	_M_xi[ 2 ] = -1;
	_M_beta2[0]=2.;
	_M_beta2[1]=-1;
	break;
      case 2:
	_M_alpha[ 0 ] = 3. / 2.;
	_M_alpha[ 1 ] = 2.;
	_M_alpha[ 2 ] = -1. / 2.;
	_M_beta[ 0 ] = 2.;
	_M_beta[ 1 ] = -1.;
	_M_xi[ 0 ] =  2.;
	_M_xi[ 1 ] = 5.;
	_M_xi[ 2 ] =  -4.;
	_M_xi[ 3 ] = 1.;
	_M_beta2[0]=3.;
	_M_beta2[ 1 ] = -3.;
	_M_beta2[ 2 ] = 1.;
	break;
      case 3:
	_M_alpha[ 0 ] = 11. / 6.;
	_M_alpha[ 1 ] = 3.;
	_M_alpha[ 2 ] = -3. / 2.;
	_M_alpha[ 3 ] = 1. / 3.;
	_M_beta[ 0 ] = 3.;
	_M_beta[ 1 ] = -3.;
	_M_beta[ 2 ] = 1.;
	_M_xi[ 0 ] =  35./12.;
	_M_xi[ 1 ] =  26./3.;
	_M_xi[ 2 ] =  -19./2;
	_M_xi[ 3 ] =  14./3.;
	_M_xi[ 4 ] =  -11./12.;
	_M_beta2[ 0 ] = 4.;
	_M_beta2[ 1 ] = -6.;
	_M_beta2[ 2 ] = 4.;
	_M_beta2[ 3 ] = -1;
	break;
      case 4:
	_M_alpha[ 0 ] = 25. / 12.;
	_M_alpha[ 1 ] = 4.;
	_M_alpha[ 2 ] = -3. ;
	_M_alpha[ 3 ] = 4. / 3.;
	_M_alpha[ 4 ] = - 1 / 4.;
	_M_beta[ 0 ] = 4.;
	_M_beta[ 1 ] = -6.;
	_M_beta[ 2 ] = 4.;
	_M_beta[ 3 ] = -1;
	break;
      }
  
    if( _M_orderDev == 2 )  
      _M_size++;  
    _M_unknowns.reserve( _M_size);
    _M_rhs.reserve(2);
}

template<typename VectorType>
void
BdfT<VectorType>:: setup ( const  std::vector<double>  /*coefficients*/,  const  UInt /*orderDev*/)
{
  ERROR_MSG("use setup for Newmark but the time advance scheme is BDF");
}

template<typename VectorType>
void BdfT<VectorType>::initialize_unk( VectorType u0 )
{
    vector_type_iterator iter     = _M_unknowns.begin();
    vector_type_iterator iter_end = _M_unknowns.end();

    for ( ; iter != iter_end; iter++ )
      {
	delete *iter;
	*iter = new VectorType(u0);
      }

    for ( UInt i(_M_unknowns.size()) ; i < _M_order; i++ )
        _M_unknowns.push_back(new VectorType(u0));

    for (UInt i=0; i<2; i++ ) 
      {
	_M_rhs.push_back(new VectorType(u0.getMap(), Unique));
	*_M_rhs[i] *=0;
      }	
    ASSERT ( _M_unknowns.size() == _M_order,
	     "_M_unknowns.size() and  _M_order must be equal" );

    return ;
}

template<typename VectorType>
void  BdfT<VectorType>::initialize_unk( VectorType /*u0*/, VectorType /*v0*/ )
{
  ERROR_MSG( "this method  is not yet implemented" );
  return;
}

template<typename VectorType>
void BdfT<VectorType>::initialize_unk( VectorType /*u0*/, VectorType /*v0*/, VectorType const & /*w0*/ )
{
  ERROR_MSG( "this method  is not yet implemented" );
  return;
}

template<typename VectorType>
void BdfT<VectorType>::initialize_unk(const  std::vector<VectorType> uv0 )
{
    UInt n0 = uv0.size();

    ASSERT( n0 != 0, "Initial data are not enough for the selected BDF" );

    vector_type_iterator iter     = _M_unknowns.begin();
    vector_type_iterator iter_end = _M_unknowns.end();

    UInt i(0);

    for ( ; iter != iter_end && i< n0 ; ++iter, ++i )
    {
    delete *iter;
    *iter = new VectorType(uv0[i]);
    }

    for ( i = _M_unknowns.size() ; i < _M_order && i< n0; ++i )
       _M_unknowns.push_back(new VectorType(uv0[i]));

    if(_M_orderDev==1)
	{
	for ( i = _M_unknowns.size() ; i < _M_order; i++ )
	  _M_unknowns.push_back(new VectorType(uv0[n0-1]));

	ASSERT ( _M_unknowns.size() == _M_order,
		 "_M_unknowns.size() and  _M_order must be equal" );
      }
    if(_M_orderDev==2)
      {
	for ( i = _M_unknowns.size() ; i < _M_order+1; i++ )
	  _M_unknowns.push_back(new VectorType(uv0[n0-1]));
	
	ASSERT ( _M_unknowns.size() == _M_order+1,
		 "_M_unknowns.size() and  _M_order must be equal" );
      }
    
    for (UInt i=0; i<2; i++ ) 
      {
	_M_rhs.push_back(new VectorType(uv0[0].getMap(), Unique));
	*_M_rhs[i] *=0;
      }
    return ;
}

//
//
//


template<typename VectorType>
const
typename BdfT<VectorType>::vector_type& BdfT<VectorType>::unk() const
{
    return _M_unknowns;
}

template<typename VectorType>
VectorType
BdfT<VectorType>::vnk()  const
{
  VectorType v(*_M_unknowns[0]);
   v  *= _M_alpha[ 0 ] / _M_dt;
   v  -= (*this->_M_rhs[0]);
   return v;
}
 
template<typename VectorType>
VectorType
BdfT<VectorType>::wnk()  const
{
  VectorType w(*_M_unknowns[0]);
  w  *= _M_xi[ 0 ]  /  (_M_dt*_M_dt);
   w  -=  ( *this->_M_rhs[1] );
   return w;
}

 /*
template<typename VectorType>
double
BdfT<VectorType>::coeff_der( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i < _M_order + 1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_alpha[ i ];
}
*/
 
template<typename VectorType>
double
BdfT<VectorType>::coeff_ext( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 & i < _M_order,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_beta[ i ];
}
  
template<typename VectorType>
double
BdfT<VectorType>::coeff_extVelocity( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 & i < _M_order + 1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_beta2[ i+1 ];
}


template<typename VectorType>
void
BdfT<VectorType>::showMe() const
{
  std::cout << "*** BDF Time discretization of order " << _M_order << " maximum order of derivate "<< _M_orderDev << " ***"
              << std::endl;
    std::cout << "    Coefficients: " << std::endl;
    for ( UInt i = 0;i < _M_order + 1;++i )
        std::cout << "       alpha(" << i << ") = " << _M_alpha[ i ]
                  << std::endl;
    for ( UInt i = 0;i < _M_order;++i )
        std::cout << "       beta (" << i << ") = " << _M_beta[ i ]
                  << std::endl;
    if (_M_orderDev ==2)
      {
	for ( UInt i = 0;i < _M_order + _M_orderDev; ++i )
	  std::cout << "     xi(" << i << ") = " << _M_xi[ i ]  << std::endl;
	for ( UInt i = 0;  i < _M_order+1; ++i  )
	  std::cout << "       beta2 (" << i << ") = " << _M_beta2[ i ]
		    << std::endl;
      }
    std::cout <<"Delta T :"<<_M_dt<<"\n";
    std::cout <<"*************************************\n";
    return ;
}

template<typename VectorType>
void
BdfT<VectorType>::shift_right( VectorType const& unk_curr )
{
    ASSERT ( _M_unknowns.size() == _M_size,
	     "_M_unknowns.size() and  _M_size must be equal" );

    vector_type_iterator it   = _M_unknowns.end() - 1;
    vector_type_iterator itm1 = _M_unknowns.end() - 1;
    vector_type_iterator itb  = _M_unknowns.begin();

    delete *itm1;

    for ( ; it != itb; --it )
    {
        itm1--;
        *it = *itm1;
    }

    *itb = new VectorType(unk_curr);

}

 
template<typename VectorType>
VectorType
BdfT<VectorType>::time_der( Real dt ) const
{
   * _M_rhs[0] *=0;
   VectorType ut(*_M_unknowns[ 0 ]);

   ut *= _M_alpha[ 1 ] / dt;
   
    for ( UInt i = 1; i < _M_order; ++i )
      {
	ut += (_M_alpha[ i + 1 ] / dt) * *_M_unknowns[ i ];
      }

    *_M_rhs[ 0 ]  = ut; 
    return ut;
}
  
// template<typename VectorType>
// VectorType
// BdfT<VectorType>::time_der( Real dt = 1.) const
// {
//     VectorType ut(_M_unknowns[ 0 ]);
//     ut *= _M_alpha[ 1 ];

//     for ( UInt i = 1;i < _M_order;++i )
//     {
//             ut += _M_alpha[ i + 1 ] * _M_unknowns[ i ];
//     }

//     return ut;

// }

template<typename VectorType>
VectorType
BdfT<VectorType>::time_derOrder2( Real dt )  const
{
  ASSERT ( _M_orderDev == 2 ,
	     " _M_orderDev must be equal two" );
  
  *_M_rhs[1] *=0;
  
  VectorType ut(*_M_unknowns[ 0 ]);
  
  ut *= _M_xi[ 1 ] / (dt*dt);
  
  for ( UInt i = 1; i < _M_order + 1; ++i )
    ut += ( _M_xi[ i + 1 ] / (dt*dt) ) * *_M_unknowns[ i ];
  
  *_M_rhs[1] = ut ;

  return ut;
}



template<typename VectorType>
VectorType
BdfT<VectorType>::extrap() const
{
    VectorType ue(*_M_unknowns[ 0 ]);
    ue *= _M_beta[ 0 ];

    for ( UInt i = 1;i < _M_order;++i )
    {
            ue += _M_beta[ i ] * *_M_unknowns[ i ];
    }

    return ue;
}

template<typename VectorType>
VectorType
BdfT<VectorType>::extrapVelocity() const
{
    VectorType ve(*_M_unknowns[ 0 ]);
    ve *= _M_beta2[ 0 ];

    for ( UInt i = 1; i < _M_order+1; ++i )
              ve += _M_beta2[ i ] * *_M_unknowns[ i ];

    return ve;
}

inline
TimeAdvance<>* createBDF() { return new BdfT<>(); }
  namespace
  {
    static bool registerBDF = TimeAdvanceFactory::instance().registerProduct( "BDF",  &createBDF );
  }
  
}

#endif
