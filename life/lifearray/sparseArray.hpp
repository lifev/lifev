/*
 This file is part of the LifeV library
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
/* --------------------------------------------------------------------------*
/                                                                            /
/      ...                                                                   /
/                                                                            /
/                                                                            /
/ MATRIX VALUES                                                              /
/                                                                            /
/ 21/06/2000 Alessandro Veneziani                                            /
/                                                                            /
/ #Purpose: Provides the basic classes for sparse matrix handling            /
/ including common formats (CSR,MSR),                                        /
/ assembling in the LifeV environment and utilities                          /
/                                                                            /
/ #Note: It will be suitably modified and linked with                        /
/ yAMLib (by Bertolazzi, Formaggia, Manzini, Veneziani)                      /
/ and its sparse environment                                                 /
/                                                                            /
/ #attempt to define the VBR matrix class in order to handle vectorial       /
/  problems. 26/10/01, Alain Gauthier.                                       /
/                                                                            /
/ #modification of CSR matrix: _value can be a vector of matrices, each      /
/  matrix being a RNM class. Attempt of 15/11/01.                            /
/                                                                            /
/ #class used in IML++: construction of a DiagPreconditioner class.          /
/  21/11/01.                                                                 /
/                                                                            /
/ #definition of MixedMatr class: holds the values of block matrix           /
/  stored under a MixedPattern format. 30/04/02.                             /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _SPARSE_ARRAY_HH
#define _SPARSE_ARRAY_HH

#include <fstream>
#include<set>
#include<algorithm>
#include<string>
#include<utility>

#include <boost/utility.hpp>



#ifndef OFFSET
#define OFFSET 0 // for the Fortran vs C numbering
#endif
#include <life/lifecore/life.hpp>

#include <life/lifearray/pattern.hpp>
#include <life/lifearray/vecUnknown.hpp>

namespace LifeV
{
typedef std::vector<INDEX_T> Container;
typedef Container::iterator ContIter;

template <typename DataType>
class MSRMatr;

template <typename PatternType, typename DataType>
class CSRMatr;

// hide zero to the external world
namespace{
double zero( double /*val*/ ) {return 0.0;}
}

}

#include <life/lifearray/MSRMatrix.hpp>
#include <life/lifearray/CSRMatrix.hpp>
#include <life/lifearray/VBRMatrix.hpp>
#include <life/lifearray/MixedMatrix.hpp>

namespace LifeV
{


////////////////////////////////////////////////////////////////
//
// DiagPreconditioner class for IML++ dense preconditioner matrix
//
///////////////////////////////////////////////////////////////
/* !\class DiagPreconditioner
   class useful for IML++ dense preconditioner matrix
*/
template <typename VectorType>
class DiagPreconditioner
{
    VectorType _diag;
public:
    DiagPreconditioner()
    {}
    ;
    //!for CSR or MSR normal pattern
    DiagPreconditioner( const CSRMatr<CSRPatt, double> &M );
    DiagPreconditioner( const MSRMatr<double> &M );
    //!for VBR pattern
    DiagPreconditioner( const VBRMatr<double> &M );
    //!for CSR block pattern
    DiagPreconditioner( const CSRMatr<CSRPatt, Tab2d> &M );

    Vector solve( const Vector &x ) const;
    VectorBlock solve( const VectorBlock &x ) const;
    VectorType trans_solve( const VectorType &x ) const
    {
        return solve( x );
    }

    const double & diag( UInt i ) const
    {
        return _diag( i );
    }
    double & diag( UInt i )
    {
        return _diag( i );
    }
    const Tab1d & diagBlock( UInt i ) const
    {
        return _diag.numBlock( i );
    }
    Tab1d & diagBlock( UInt i )
    {
        return _diag.numBlock( i );
    }
};

////////////////////////////////////////////////////////////////
//
// IDPreconditioner class for IML++ dense preconditioner matrix
// Here the preconditioner is simply the identity matrix
//
///////////////////////////////////////////////////////////////
/*! \class IDPreconditioner
  class useful for IML++ dense preconditioner matrix.
  Here the preconditioner is simply the identity matrix (no preconditioning).
 */
template <typename VectorType>
class IDPreconditioner
{
    VectorType _diag;
public:
    IDPreconditioner()
    {}
    ;
    //! for CSR or MSR normal pattern
    IDPreconditioner( const CSRMatr<CSRPatt, double> &M );
    IDPreconditioner( const MSRMatr<double> &M );
    //! for VBR pattern
    IDPreconditioner( const VBRMatr<double> &M );
    //! for CSR block pattern
    IDPreconditioner( const CSRMatr<CSRPatt, Tab2d> &M );

    Vector solve( const Vector &x ) const;
    VectorBlock solve( const VectorBlock &x ) const;
    VectorType trans_solve( const VectorType &x ) const
    {
        return solve( x );
    }

    const double & diag( UInt i ) const
    {
        return _diag( i );
    }
    double & diag( UInt i )
    {
        return _diag( i );
    }
    const Tab1d & diagBlock( UInt i ) const
    {
        return _diag.numBlock( i );
    }
    Tab1d & diagBlock( UInt i )
    {
        return _diag.numBlock( i );
    }
};






//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// DiagPreconditioner
//-----------------------------------------------------------------------

//for CSR or MSR normal pattern
template <>
DiagPreconditioner<Vector>::DiagPreconditioner( const CSRMatr<CSRPatt, double> &M );

template <>
DiagPreconditioner<Vector>::DiagPreconditioner( const MSRMatr<double> &M );



//for VBR pattern
template <>
DiagPreconditioner<Vector>::DiagPreconditioner( const VBRMatr<double> &M );


//for CSR block pattern
template <>
DiagPreconditioner<VectorBlock>::DiagPreconditioner( const CSRMatr<CSRPatt, Tab2d> &M );


//solve the diagonal system
template <>
Vector
DiagPreconditioner<Vector>::solve( const Vector &x ) const;


template <>
VectorBlock
DiagPreconditioner<VectorBlock>::solve( const VectorBlock &x ) const;


//-----------------------------------------------------------------------
// DiagPreconditioner
//-----------------------------------------------------------------------

//for CSR or MSR normal pattern
template <>
IDPreconditioner<Vector>::IDPreconditioner( const CSRMatr<CSRPatt, double> &M );


template <>
IDPreconditioner<Vector>::IDPreconditioner( const MSRMatr<double> &M );


//for VBR pattern
template <>
IDPreconditioner<Vector>::IDPreconditioner( const VBRMatr<double> &M );



//for CSR block pattern
template <>
IDPreconditioner<VectorBlock>::IDPreconditioner( const CSRMatr<CSRPatt, Tab2d> &M );


//solve the diagonal system
template <>
Vector
IDPreconditioner<Vector>::solve( const Vector &x ) const;

template <>
VectorBlock
IDPreconditioner<VectorBlock>::solve( const VectorBlock &x ) const;

} // namespace LifeV

#endif // _SPARSE_ARRAY_HH
