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
/*!
  \file bcVector.cc
  \brief Implementations for bcVector.h
  \version 1.0
  \author M.A. Fernandez
  \date 11/2002
 
  modified by Vincent Martin
  02/2003
 
*/

#include "bcVector.hpp"
namespace LifeV
{
//
// Implementation for BCVector_Base
//


//! Default Constructor (the user must call setBCVector(..))
BCVector_Base::BCVector_Base() :
        _MixteCoef( 0.0 ), _finalized( false ), _type( 0 )
{}

//! Constructor
BCVector_Base::BCVector_Base( Vector& vec, UInt nbTotalDof ) :
        _vec( &vec ), _nbTotalDof( nbTotalDof ),
        _MixteCoef( 0.0 ), _finalized( false ), _type( 0 )
{}

//! Constructor
BCVector_Base::BCVector_Base( Vector& vec, const UInt nbTotalDof, UInt type ) :
        _vec( &vec ), _nbTotalDof( nbTotalDof ),
        _MixteCoef( 0.0 ), _finalized( false ), _type( type )
{}


//! set the Mixte coefficient
void BCVector_Base::setMixteCoef( const Real& coef )
{
    _MixteCoef = coef;
}


//! Return the value of the Mixte coefficient
Real BCVector_Base::MixteCoef() const
{
    return _MixteCoef;
}

//! Return the value of type
UInt BCVector_Base::type() const
{
    return _type;
}



//
// Implementation for BCVector
//

//! Default Constructor (the user must call setBCVector(..))
BCVector::BCVector()
{}

//! Constructor
BCVector::BCVector( Vector& vec, UInt nbTotalDof ) :
        BCVector_Base( vec, nbTotalDof )
{
    _finalized = true;
}

//! Constructor
BCVector::BCVector( Vector& vec, UInt const nbTotalDof, UInt type ) :
        BCVector_Base( vec, nbTotalDof, type )
{
    _finalized = true;
}




//!set the BC vector (after default construction)
void BCVector::setvector( Vector& vec, UInt nbTotalDof )
{
    ASSERT_PRE( !_finalized, "BC Vector cannot be set twice." );
    _vec = &vec ;
    _nbTotalDof = nbTotalDof;
    _finalized = true;
}


//! This method returns the value to be imposed in the component iComp of the dof iDof
Real BCVector::operator() ( const ID& iDof, const ID& iComp ) const
{
    ASSERT_PRE( _finalized, "BC Vector should be finalized before being accessed." );
    return ( *_vec ) ( ( iComp - 1 ) * _nbTotalDof + iDof - 1 );
}

//! Assignment operator for BCVector
BCVector & BCVector::operator=( const BCVector& BCv )
{
    _vec = BCv._vec;
    _nbTotalDof = BCv._nbTotalDof;
    _MixteCoef = BCv._MixteCoef;
    _finalized = BCv._finalized;
    _type = BCv._type;
    return *this;
}


//! Output
std::ostream& BCVector::showMe( bool verbose, std::ostream & out ) const
{
    ASSERT_PRE( _finalized, "BC Vector should be finalized before being accessed." );
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    out << "BC Vector Interface: " << std::endl;
    out << "number of interface vector Dof : " << _nbTotalDof << std::endl;
    out << "==>Interface Dof :\n";
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    return out;
}


//
// Implementation for BCVector_Interface
//

//! Default Constructor (the user must call setBCVector(..))
BCVector_Interface::BCVector_Interface()
{}

//! Constructor
BCVector_Interface::BCVector_Interface( Vector& vec, UInt nbTotalDof,
                                        DofInterfaceBase& dofIn ) :
        BCVector_Base( vec, nbTotalDof ),
        _dofIn( &dofIn )
{
    _finalized = true;
}


//!set the BC vector (after default construction)
void BCVector_Interface::setvector( Vector& vec, UInt nbTotalDof, DofInterfaceBase& dofIn )
{
    ASSERT_PRE( !_finalized, "BC Vector cannot be set twice." );
    _vec = &vec ;
    _dofIn = &dofIn;
    _nbTotalDof = nbTotalDof;
    _finalized = true;
}


//! This method returns the value to be imposed in the component iComp of the dof iDof
Real BCVector_Interface::operator() ( const ID& iDof, const ID& iComp ) const
{
    ASSERT_PRE( _finalized, "BC Vector should be finalized before being accessed." );
    return ( *_vec ) ( ( iComp - 1 ) * _nbTotalDof + _dofIn->getInterfaceDof( iDof ) - 1 );
}

//! Assignment operator for BCVector_Interface
BCVector_Interface & BCVector_Interface::operator=( const BCVector_Interface & BCv )
{
    _vec = BCv._vec;
    _dofIn = BCv._dofIn;
    _nbTotalDof = BCv._nbTotalDof;
    _MixteCoef = BCv._MixteCoef;
    _finalized = BCv._finalized;
    _type = BCv._type;
    return *this;
}


//! Output
std::ostream& BCVector_Interface::showMe( bool verbose, std::ostream & out ) const
{
    ASSERT_PRE( _finalized, "BC Vector should be finalized before being accessed." );
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    out << "BC Vector Interface: " << std::endl;
    out << "number of interface vector Dof : " << _nbTotalDof << std::endl;
    out << "==>Interface Dof :\n";
    _dofIn->showMe( verbose, out );  // no showMe(..) in Miguel's DofInterface
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    return out;
}
}
