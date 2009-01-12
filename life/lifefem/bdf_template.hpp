/*
 This file is part of the LifeV library

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

namespace LifeV
{
const UInt BDFT_MAX_ORDER = 3;

//typedef Real ( *Funct ) ( const Real&, const Real&, const Real&, const Real&,
//                          const ID& );

/*!
  \class BdfT (templated)
  \brief Backward differencing formula time discretization

  A differential equation of the form

  \f$ M u' = A u + f \f$

  is discretized in time as

  \f$ M p'(t_{k+1}) = A u_{k+1} + f_{k+1} \f$

  where p denotes the polynomial of order n in t that interpolates
  (t_i,u_i) for i = k-n+1,...,k+1.

  The approximative time derivative \f$ p'(t_{k+1}) \f$ is a linear
  combination of state vectors u_i:

  \f$ p'(t_{k+1}) = \frac{1}{\Delta t} (\alpha_0 u_{k+1} - \sum_{i=0}^n \alpha_i u_{k+1-i} )\f$

  Thus we have

  \f$ \frac{\alpha_0}{\Delta t} M u_{k+1} = A u_{k+1} + f + M \bar{p} \f$

  with

  \f$ \bar{p} = \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i u_{k+1-i} \f$

  This class stores the n last state vectors in order to be able to
  calculate \f$ \bar{p} \f$. It also provides alpha_i
  and can extrapolate the the new state from the n last states with a
  polynomial of order n-1:

  \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$
*/
template<typename VectorType = EpetraVector >
class BdfT // T means template
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
    BdfT( const UInt n );

    ~BdfT();

    //! Initialize all the entries of the unknown vector to be derived with the
    //! vector u0 (duplicated)
    void initialize_unk( VectorType u0 );

    //! Initialize all the entries of the unknown vector to be derived with a
    //! set of vectors uv0
    //! note: this is taken as a copy (not a reference), since uv0 is resized inside the method.
    void initialize_unk( std::vector<VectorType> uv0 );


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
    /*! Update the vectors of the previous time steps by shifting on the right
     *  the old values.
     *  @param u_curr current (new) value of the state vector
     */
    void shift_right( VectorType const& u_curr );

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula
    VectorType time_der( Real dt = 1 ) const;

    //! Returns the right hand side \f$ \bar{p} \Delta t \f$ of the time
    //! derivative formula. The timestep is taken into account elsewhere,
    //! e. g. in the mass matrix.
//    VectorType time_der() const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    VectorType extrap() const;

    //! Return the i-th coefficient of the time derivative alpha_i
    double coeff_der( UInt i ) const;

    //! Return the i-th coefficient of the time extrapolation beta_i
    double coeff_ext( UInt i ) const;

    //! Return a vector with the last n state vectors
    const vector_type& unk() const;

    void showMe() const;

    UInt order() const {return _M_order;}

private:
    //! Order of the BDF derivative/extrapolation: the time-derivative
    //! coefficients vector has size n+1, the extrapolation vector has size n
    UInt _M_order;

    // Size of the unknown vector
    // UInt _M_size;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    Vector _M_alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    Vector _M_beta;

    //! Last n state vectors
    vector_type _M_unknowns;
};


///
// template implementations
//

template<typename VectorType>
BdfT<VectorType>::
BdfT( const UInt n )
    :
    _M_order( n ),
    //    _M_size( 0 ),
    _M_alpha( n + 1 ),
    _M_beta( n ),
    _M_unknowns()
{
    if ( n <= 0 || n > BDFT_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
        << " you want to use BDF order " << n << "\n"
        << " we support BDF order from 1 to " << BDFT_MAX_ORDER << "\n";
        throw std::invalid_argument( __ex.str() );
    }
    switch ( n )
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
    }
    _M_unknowns.reserve( n );
}



template<typename VectorType>
BdfT<VectorType>::~BdfT()
{
    vector_type_iterator iter     = _M_unknowns.begin();
    vector_type_iterator iter_end = _M_unknowns.end();

    for ( ; iter != iter_end; iter++ )
        delete *iter;
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

    ASSERT ( _M_unknowns.size() == _M_order,
	     "_M_unknowns.size() and  _M_order must be equal" );

    return ;
}


template<typename VectorType>
void BdfT<VectorType>::initialize_unk( std::vector<VectorType> uv0 )
{
    UInt n0 = uv0.size();

    ASSERT( n0 != 0, "Initial data are not enough for the selected BDF" );

    vector_type_iterator iter     = _M_unknowns.begin();
    vector_type_iterator iter_end = _M_unknowns.end();

    int i(0);

    for ( ; iter != iter_end && i< n0 ; ++iter, ++i )
      {
	delete *iter;
	*iter = new VectorType(uv0[i]);
      }

    for ( i = _M_unknowns.size() ; i < _M_order && i< n0; ++i )
        _M_unknowns.push_back(new VectorType(uv0[i]));


    for ( i = _M_unknowns.size() ; i < _M_order; i++ )
        _M_unknowns.push_back(new VectorType(uv0[n0-1]));

    ASSERT ( _M_unknowns.size() == _M_order,
	     "_M_unknowns.size() and  _M_order must be equal" );

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
double
BdfT<VectorType>::coeff_der( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i < _M_order + 1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_alpha[ i ];
}

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
void
BdfT<VectorType>::showMe() const
{
    std::cout << "*** BDF Time discretization of order " << _M_order << " ***"
              << std::endl;
    std::cout << "    Coefficients: " << std::endl;
    for ( UInt i = 0;i < _M_order + 1;++i )
        std::cout << "       alpha(" << i << ") = " << _M_alpha[ i ]
                  << std::endl;
    for ( UInt i = 0;i < _M_order;++i )
        std::cout << "       beta (" << i << ") = " << _M_beta[ i ]
                  << std::endl;
    return ;
}

template<typename VectorType>
void
BdfT<VectorType>::shift_right( VectorType const& unk_curr )
{
    ASSERT ( _M_unknowns.size() == _M_order,
	     "_M_unknowns.size() and  _M_order must be equal" );

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
    VectorType ut(*_M_unknowns[ 0 ]);
    ut *= _M_alpha[ 1 ] / dt;

    for ( UInt i = 1;i < _M_order;++i )
    {
            ut += (_M_alpha[ i + 1 ] / dt) * *_M_unknowns[ i ];
    }

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
}

#endif
