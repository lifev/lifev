/*
 This file is part of the LifeV library
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
/* -------------------------------------------------------------------------*/
/*!
  \file vecUnknown.h

  Vector classes

  #purpose: provides vector classes handler useful in case of solving
   vector problem.
   Provides Vector and VectorBlock class necessary for IML++ library
   Alain Gauthier.
*/
#ifndef _VEC_UNKNOWN_HH
#define _VEC_UNKNOWN_HH

#include <life.hpp>
#include <vector>

#include "tab.hpp"

namespace LifeV
{

/*!\class VectorBlock
  Version for block matrices handling.
  19/11/01
*/
class VectorBlock
{
    //std::vector<Tab1d> _v; //!< container of block vector
    std::vector<Vector> _v; //!< container of block vector
public:
    typedef Vector sub_array_type;
    // Default constructor we need for IML++
    VectorBlock()
        {}
    explicit VectorBlock( int n, int blsize = 1 )  //!< default size block =1 !
        {
            sub_array_type bloc( blsize );
            bloc = ZeroVector( blsize );
            _v.resize( n, bloc );
        }
    VectorBlock( VectorBlock const &ex_v )
        {
            *this = ex_v;
        }
    std::vector<sub_array_type> v() const
        {
            return _v;
        }
    UInt size() const
        {
            return _v.size();
        }
    //operators
    VectorBlock & operator=( VectorBlock const &ex_v )
        {
            if ( &ex_v != this )
            {
                _v.clear();
                _v.reserve( ex_v.v().size() );
                for ( UInt ir = 0; ir < ex_v.v().size(); ir++ )
                    _v.push_back( ex_v.v() [ ir ] );
            }
            return *this;
        }
    VectorBlock operator=( double const val )
        {
            for ( std::vector<sub_array_type>::iterator ip = _v.begin(); ip < _v.end(); ip++ )
                *ip = ScalarVector( ip->size(), val );
            return *this;
        }
    VectorBlock operator-=( VectorBlock const &ex_v )
        {
            for ( UInt ir = 0; ir < ex_v.v().size(); ir++ )
                _v[ ir ] -= ex_v.numBlock( ir );
            return *this;
        }
    VectorBlock operator+=( VectorBlock const &ex_v )
        {
            for ( UInt ir = 0; ir < ex_v.v().size(); ir++ )
                _v[ ir ] += ex_v.numBlock( ir );
            return *this;
        }
    double operator() ( int i ) const
        {
            int blsize=_v[ 0 ].size();
            int bloc_i = i / blsize;
            int loc_i = i % blsize;
            return _v[ bloc_i ] ( loc_i );
        }
    double& operator[] ( int i )
        {
            int blsize=_v[ 0 ].size();
            int bloc_i = i / blsize;
            int loc_i = i % blsize;
            return _v[ bloc_i ] ( loc_i );
        }
    sub_array_type & numBlock( int i )
        {
            return _v[ i ];
        }
    const sub_array_type & numBlock( int i ) const
        {
            return _v[ i ];
        }

    friend VectorBlock operator+( VectorBlock const &ex_v1, VectorBlock
                                  const &ex_v2 );
    friend VectorBlock operator-( VectorBlock const &ex_v1, VectorBlock
                                  const &ex_v2 );
    friend VectorBlock operator*( VectorBlock const &ex_v, double const val );
    friend VectorBlock operator*( double const val, VectorBlock const &ex_v );

    //! vector inner product
    friend double dot( VectorBlock const &ex_v1, VectorBlock const &ex_v2 );
    //! norm derived from dot:
    friend inline double norm( VectorBlock const &ex_v )
        {
            return sqrt( dot( ex_v, ex_v ) );
        }

};

/*! \class PhysVectUnknown
  vector unknown which has the dimension of the physical domain
  The type VectorType could be a Vector, or a VectorBlock or a
  std::vectordouble> depending on the choice of the linear system solver and
  on the case of scalar or vectorial problem
*/
template <typename VectorType>
class PhysVectUnknown
    :
        public VectorType

{
    static const UInt _S_nbcomp = nDimensions;
public:

    typedef VectorType super;

    PhysVectUnknown()
        :
        super(),
        _M_name("vector_unknown")
        {}

    //  PhysVectUnknown(){}
    explicit PhysVectUnknown( std::string const& __s )
        :
        super(),
        _M_name( __s )
        {}

    //  PhysVectUnknown(){}
    explicit PhysVectUnknown( UInt const Ndof )
        :
        super( nDimensions*Ndof ),
        _M_name( "vector_unknown" )
        {}

    PhysVectUnknown( PhysVectUnknown<VectorType> const &RhPhysVectUnknown )
        :
        super( RhPhysVectUnknown ),
        _M_name( RhPhysVectUnknown._M_name )
        {}

    PhysVectUnknown& operator=( PhysVectUnknown const& __v )
        {
            if ( this == &__v )
                return * this;
            super::operator=( ( super const& )__v );
            _M_name = __v._M_name;
            return *this;
        }
    template<typename VectorExpr>
    PhysVectUnknown& operator=( VectorExpr const& __v )
        {
            super::operator=( __v );
            return *this;
        }
    //! gives the front of the vector
    Real * giveVec()
        {
            return & ( ( *this ) [ 0 ] );
        }
    //! gives the front of the vector
    Real const* giveVec() const
        {
            return & ( ( *this ) [ 0 ] );
        }
    UInt size() const
        {
            return super::size();
        }
    void resize( size_t __s, bool preserve = false )
        {
            super::resize( __s*_S_nbcomp, preserve );
        }
    static UInt nbcomp()
        {
            return _S_nbcomp;
        }

    void setName( std::string const& __name )
        {
            _M_name = __name;
        }
    std::string const& name() const
        {
            return _M_name;
        }
private:
    std::string _M_name;
};

//! the case of VectorBlock type
template <>
PhysVectUnknown<VectorBlock>::PhysVectUnknown( UInt const Ndof );


/*! \class ScalUnknown
  scalar unknown of dimension=1
  The type VectorType could be a Vector, or a VectorBlock or a
  std::vectordouble> depending on the choice of the linear system solver and
  on the case of scalar or vectorial problem
*/
template <typename VectorType>
class ScalUnknown
    :
    public VectorType
{
    static const UInt _S_nbcomp = 1;
public:

    typedef VectorType super;

    ScalUnknown()
        :
        super(),
        _M_name( "scalar_unknown" )

        {}
    ScalUnknown( std::string const& __s )
        :
        super(),
        _M_name( __s )

        {}
    explicit ScalUnknown( UInt const Ndof )
        :
        super( Ndof ),
        _M_name( "scalar_unknown" )
        {}

    ScalUnknown( const ScalUnknown<VectorType> &RhScalUnknown )
        :
        super( RhScalUnknown ),
        _M_name( RhScalUnknown._M_name )
        {}

    ScalUnknown& operator=( ScalUnknown const& __v )
        {
            if ( this == &__v )
                return * this;
            super::operator=( ( super const& )__v );
            _M_name = __v._M_name;
            return *this;
        }
    template<typename VectorExpr>
    ScalUnknown& operator=( VectorExpr const& __v )
        {
            super::operator=( __v );
            return *this;
        }
    //! gives the front of the vector
    Real * giveVec()
        {
            return & ( ( *this ) [ 0 ] );
        }
    //! gives the front of the vector
    Real const* giveVec() const
        {
            return & ( ( *this ) [ 0 ] );
        }
    UInt size() const
        {
            return super::size();
        }
    static UInt nbcomp()
        {
            return _S_nbcomp;
        }
    void setName( std::string const& __name )
        {
            _M_name = __name;
        }
    std::string const& name() const
        {
            return _M_name;
        }
private:
    std::string _M_name;
};

//---------------------------------------------------------------//

/*! \class GenericVecHdl
  vector problem handler
  The type VectorType could be a Vector, or a VectorBlock or a
  std::vectordouble> depending on the choice of the linear system solver and
  on the case of scalar or vectorial problem
*/
template <typename VectorType>
class GenericVecHdl
    :
    public VectorType
{
    UInt _size;
    UInt _nbcomp;
public:

    typedef VectorType super;

    GenericVecHdl()
        {}

    GenericVecHdl( UInt ex_size, UInt ex_nbcomp )
        :
        super( ex_size ),
        _size( ex_size ),
        _nbcomp( ex_nbcomp )
        {}
    ;

    GenericVecHdl( const GenericVecHdl<VectorType> &RhGenVec );
    //! construction from a scalar unknown:
    GenericVecHdl( const ScalUnknown<VectorType> &RhScalUnknown );
    //! construction from two scalar unknowns:
    GenericVecHdl( const ScalUnknown<VectorType> &RhScU1,
                   const ScalUnknown<VectorType> &RhScU2 );
    //! construction from a physical vectorial unknown:
    GenericVecHdl( const PhysVectUnknown<VectorType> &RhPhVU );
    //! construction from a physical vectorial unknown and a scalar unknown:
    GenericVecHdl( const PhysVectUnknown<VectorType> &RhPhVU,
                   const ScalUnknown<VectorType> &RhScU );
    GenericVecHdl( const ScalUnknown<VectorType> &RhScU,
                   const PhysVectUnknown<VectorType> &RhPhVU );
    //! construction from a physical vectorial unknown and a generic unknown:
    GenericVecHdl( const PhysVectUnknown<VectorType> &RhPhVU,
                   const GenericVecHdl<VectorType> &RhGenVec );
    GenericVecHdl( const GenericVecHdl<VectorType> &RhGenVec,
                   const PhysVectUnknown<VectorType> &RhPhVU );


    GenericVecHdl& operator=( double const __val )
        {
            super::operator=( __val );
            return *this;
        }

    GenericVecHdl& operator=( VectorType const& __v )
        {
            if ( this == &__v )
                return * this;
            super::operator=( ( super const& )__v );
            return *this;

        }
    //! gives the front of the vector
    inline Real * giveVec()
        {
            return & ( ( *this ) [ 0 ] );
        }
    inline UInt size() const
        {
            return _size;
        }
    //! gives the size of one block
    inline UInt nbcomp() const
        {
            return _nbcomp;
        }


};

///////////////////////////////////////////////////////////
inline VectorBlock
operator+( VectorBlock const &ex_v1, VectorBlock
           const &ex_v2 )
{
    ASSERT( ex_v1._v.size() == ex_v2._v.size(), " wrong dimension of vector" );
    VectorBlock ans( ex_v1._v.size(), ex_v1.numBlock( 0 ).size() );
    ans = 0.0;
    ans = ex_v1;
    ans += ex_v2;
    return ans;
}
inline VectorBlock
operator-( VectorBlock const &ex_v1, VectorBlock
           const &ex_v2 )
{
    ASSERT( ex_v1._v.size() == ex_v2._v.size(), " wrong dimension of vector" );
    VectorBlock ans( ex_v1._v.size(), ex_v1.numBlock( 0 ).size() );
    ans = 0.0;
    ans = ex_v1;
    ans -= ex_v2;
    return ans;
}
inline VectorBlock
operator*( VectorBlock const &ex_v, double const val )
{
    VectorBlock ans;
    ans = ex_v;
    for ( UInt i = 0; i < ex_v._v.size(); i++ )
        ans.numBlock( i ) *= val;
    return ans;
}
inline VectorBlock
operator*( double const val, VectorBlock const &ex_v )
{
    return ex_v * val;
}
inline double
dot( VectorBlock const &ex_v1, VectorBlock const &ex_v2 )
{
    double ans = 0;
    for ( UInt i = 0; i < ex_v1._v.size(); i++ )
        ans += dot( ex_v1.numBlock( i ), ex_v2.numBlock( i ) );
    return ans;
}

//! norm derived from dot:
inline double norm_2( VectorBlock const &ex_v )
    {
        return sqrt( dot( ex_v, ex_v ) );
    }


//---------------------------------------------------------------//
// Useful function

//! assign the values of a pointer to an existing vector
//! No check of the dimensions !!!
template <typename VectorType>
void
point2Vector( double const * point, VectorType & v )
{
    VectorType ans( v.size() );
    for ( UInt i = 0; i < ans.size(); i++ )
    {
        ans[ i ] = *point;
        point++;
    }
    v = ans;
}


/*-------------------------------------------------------/
  IMPLEMENTATION
  /-------------------------------------------------------*/


//////////////////////////////
// class GenericVecHdl
/////////////////////////////

template <typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl( const GenericVecHdl<VectorType> &RhGenVec )
    :
    super( RhGenVec ),
    _size( RhGenVec.size() ),
    _nbcomp( RhGenVec.nbcomp() )
{}

//! construction from a scalar unknown:
template <typename VectorType>
GenericVecHdl<VectorType>::
GenericVecHdl( const ScalUnknown<VectorType> &RhScalUnknown )
    :
    super( RhScalUnknown ),
    _size( RhScalUnknown.size() ),
    _nbcomp( RhScalUnknown.nbcomp() )
{}

//! construction from two scalar unknowns:
template <typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl( const ScalUnknown<VectorType> &RhScU1, const ScalUnknown<VectorType> &RhScU2 )
    :
    super( RhScU1.size() + RhScU2.size() ),
    _size( RhScU1.size() + RhScU2.size() ),
    _nbcomp( RhScU1.nbcomp() + RhScU2.nbcomp() )
{
    for ( UInt i = 0; i < RhScU1.size(); ++i )
        ( *this ) [ i ] = RhScU1[ i ];
    for ( UInt i = 0; i < RhScU2.size(); ++i )
        ( *this ) [ i + RhScU1.size() ] = RhScU2[ i ];
}
//! the case of VectorBlock type
template <>
GenericVecHdl<VectorBlock>::GenericVecHdl( const ScalUnknown<VectorBlock> &RhScU1,
                                           const ScalUnknown<VectorBlock> &RhScU2 );


//! construction from a physical vectorial unknown:
template <typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl( const PhysVectUnknown<VectorType> &RhPhVU )
    :
    super( RhPhVU ),
    _size( RhPhVU.size() ),
    _nbcomp( RhPhVU.nbcomp() )
{}

//! construction from a physical vectorial unknown and a scalar unknown:
template <typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl( const PhysVectUnknown<VectorType> &RhPhVU,
                                          const ScalUnknown<VectorType> &RhScU )
    :
    super( RhPhVU.size() + RhScU.size() ),
    _size( RhPhVU.size() + RhScU.size() ),
    _nbcomp( RhPhVU.nbcomp() + RhScU.nbcomp() )
{
    for ( UInt i = 0; i < RhPhVU.size(); ++i )
        ( *this ) [ i ] = RhPhVU[ i ];
    for ( UInt i = 0; i < RhScU.size(); ++i )
        ( *this ) [ RhPhVU.size() + i ] = RhScU[ i ];
}
//! the case of VectorBlock type
template <>
GenericVecHdl<VectorBlock>::GenericVecHdl( const PhysVectUnknown<VectorBlock> &RhPhVU,
                                           const ScalUnknown<VectorBlock> &RhScU );


template <typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl( const ScalUnknown<VectorType> &RhScU,
                                          const PhysVectUnknown<VectorType> &RhPhVU )
    :
    super( RhPhVU.size() + RhScU.size() ),
    _size( RhPhVU.size() + RhScU.size() ),
    _nbcomp( RhPhVU.nbcomp() + RhScU.nbcomp() )
{
    for ( UInt i = 0; i < RhScU.size(); ++i )
        ( *this ) [ i ] = RhScU[ i ];
    for ( UInt i = 0; i < RhPhVU.size(); ++i )
        ( *this ) [ RhScU.size() + i ] = RhPhVU[ i ];
}
//! the case of VectorBlock type
template <>
GenericVecHdl<VectorBlock>::GenericVecHdl( const ScalUnknown<VectorBlock> &RhScU,
                                           const PhysVectUnknown<VectorBlock> &RhPhVU );


//! construction from a physical vectorial unknown and a generic unknown:
template <typename VectorType>
GenericVecHdl<VectorType>::GenericVecHdl( const PhysVectUnknown<VectorType> &RhPhVU,
                                          const GenericVecHdl<VectorType> &RhGenVec )
    :
    super( RhPhVU.size() + RhGenVec.size() ),
    _size( RhPhVU.size() + RhGenVec.size() ),
    _nbcomp( RhPhVU.nbcomp() + RhGenVec.nbcomp() )
{
    for ( UInt i = 0;i < RhPhVU.size();++i )
        ( *this ) [ i ] = RhPhVU[ i ];
    for ( UInt i = 0;i < RhGenVec.size();++i )
        ( *this ) [ RhPhVU.size() + i ] = RhGenVec[ i ];
}
//! the case of VectorBlock type
template <>
GenericVecHdl<VectorBlock>::GenericVecHdl( const PhysVectUnknown<VectorBlock> &RhPhVU,
                                           const GenericVecHdl<VectorBlock> &RhGenVec );


template <typename VectorType>
GenericVecHdl<VectorType>::
GenericVecHdl( const GenericVecHdl<VectorType> &RhGenVec,
               const PhysVectUnknown<VectorType> &RhPhVU ) :
    super( RhPhVU.size() + RhGenVec.size() ),
    _size( RhPhVU.size() + RhGenVec.size() ),
    _nbcomp( RhPhVU.nbcomp() + RhGenVec.nbcomp() )
{
    for ( UInt i = 0;i < RhGenVec.size();++i )
        ( *this ) [ i ] = RhGenVec[ i ];
    for ( UInt i = 0;i < RhPhVU.size();++i )
        ( *this ) [ i + RhGenVec.size() ] = RhPhVU[ i ];
}
//! the case of VectorBlock type
template <>
GenericVecHdl<VectorBlock>::GenericVecHdl( const GenericVecHdl<VectorBlock> &RhGenVec,
                                           const PhysVectUnknown<VectorBlock> &RhPhVU );

}
#endif
