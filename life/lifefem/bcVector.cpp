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
  \brief Implementations for bcVector.hpp
  \version 1.0
  \author M.A. Fernandez
  \date 11/2002

  modified by Vincent Martin
  02/2003

*/

#include <life/lifefem/bcVector.hpp>
namespace LifeV
{
//
// Implementation for BCVectorBase
//


//! Default Constructor (the user must call setBCVector(..))
BCVectorBase::BCVectorBase()
    :
   _M_mixteCoef( 0.0 ),
   _M_type( 0 ),
   _M_finalized( false )
{}

//! Constructor
BCVectorBase::BCVectorBase( const EpetraVector& vec, const UInt nbTotalDof, UInt type )
    :
    _M_vec       ( &vec ),
    _M_nbTotalDof( nbTotalDof ),
    _M_mixteCoef ( 0.0 ),
    _M_type      ( type ),
    _M_finalized ( false )
{}

BCVectorBase&
BCVectorBase::operator=( BCVectorBase const& __bcv )
{
    if ( this == &__bcv )
        return *this;

    _M_vec        = __bcv._M_vec;
    _M_nbTotalDof = __bcv._M_nbTotalDof;
    _M_mixteCoef  = __bcv._M_mixteCoef;
    _M_type       = __bcv._M_type;
    _M_finalized  = __bcv._M_finalized;

    return *this;
}

//! This method returns the value to be imposed in the component iComp of the dof iDof
Real
BCVectorBase::operator() ( const ID& iDof, const ID& iComp ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *_M_vec ) ( ( iComp - 1 ) * _M_nbTotalDof + iDof );
}

//! This method returns the value of the mixte coefficient to be imposed in the component iComp of the dof iDof
Real
BCVectorBase::MixteVec ( const ID& iDof, const ID& iComp ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *_M_vec_mixte ) ( ( iComp - 1 ) * _M_nbTotalDof + iDof );
}


void
BCVectorBase::setVector( EpetraVector& __vec, UInt __nbTotalDof, UInt type )
{
    _M_vec = &__vec ;
    _M_nbTotalDof = __nbTotalDof;
    _M_type = type;
    _M_finalized = true;
}

//
// Implementation for BCVector
//

//! Default Constructor (the user must call setBCVector(..))
BCVector::BCVector()
{}

//! Constructor
BCVector::BCVector( EpetraVector& vec, UInt const nbTotalDof, UInt type )
    :
    BCVectorBase( vec, nbTotalDof, type )
{
    this->setFinalized( true );
}


//! Assignment operator for BCVector
BCVector&
BCVector::operator=( const BCVector& BCv )
{
    if ( this == &BCv )
        return *this;

    super::operator=( ( super& )BCv );

    return *this;
}


//! Output
std::ostream&
BCVector::showMe( bool /* verbose */, std::ostream & out ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    out << "BC Vector Interface: " << std::endl;
    out << "number of interface vector Dof : " << this->nbTotalDOF() << std::endl;
    out << "==>Interface Dof :\n";
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    return out;
}
BCVectorBase*
createBCVector( BCVectorBase const* __bc )
{
    return new BCVector( ( BCVector const& )*__bc );
}
// register BCFunctionBase in factory for cloning
const bool __bcvec = FactoryCloneBCVector::instance().registerProduct( typeid(BCVector), &createBCVector );


//
// Implementation for BCVectorInterface
//

//! Default Constructor (the user must call setBCVector(..))
BCVectorInterface::BCVectorInterface()
{}

//! Constructor
BCVectorInterface::BCVectorInterface( const EpetraVector& vec, UInt nbTotalDof,
                                      dof_interface_type dofIn, UInt type )
    :
    BCVectorBase( vec, nbTotalDof, type ),
    _M_dofIn( dofIn )
{
    this->setFinalized( true );
}

//! setup after default constructor

void BCVectorInterface::setup( const EpetraVector& vec, UInt nbTotalDof, dof_interface_type dofIn, UInt type )
{
    _M_vec        = &vec;
    _M_nbTotalDof = nbTotalDof;
    _M_mixteCoef  = 0.0;
    _M_type       = type;
    _M_dofIn      = dofIn;
    this->setFinalized( true );
}

//!set the BC vector (after default construction)
void
BCVectorInterface::setVector( EpetraVector& vec, UInt nbTotalDof, dof_interface_type dofIn, UInt type )
{
    ASSERT_PRE( !this->isFinalized(), "BC Vector cannot be set twice." );

    super::setVector( vec, nbTotalDof, type );

    _M_dofIn = dofIn;

}


//! This method returns the value to be imposed in the component iComp of the dof iDof
Real
BCVectorInterface::operator() ( const ID& iDof, const ID& iComp ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
//     std::cout << iDof << " " << std::flush
//                  << ( iComp - 1 ) * _M_nbTotalDof + _M_dofIn->getInterfaceDof( iDof ) << " "
//                  <<  ( *_M_vec ) ( ( iComp - 1 ) * _M_nbTotalDof + _M_dofIn->getInterfaceDof( iDof )  ) << std::endl;
    return ( *_M_vec ) (( iComp - 1 ) * _M_nbTotalDof + _M_dofIn->getInterfaceDof( iDof ));
}


//! This method returns the value of the mixte coefficient to be imposed in the component iComp of the dof iDof
Real
BCVectorInterface::MixteVec( const ID& iDof, const ID& iComp ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    return ( *_M_vec_mixte ) (( iComp - 1 ) * _M_nbTotalDof + _M_dofIn->getInterfaceDof( iDof ));
}


//! Assignment operator for BCVectorInterface
BCVectorInterface&
BCVectorInterface::operator=( const BCVectorInterface & BCv )
{
    if ( this == &BCv )
        return *this;

    super::operator=( ( super& ) BCv );

    _M_dofIn = BCv._M_dofIn;
    return *this;
}


//! Output
std::ostream&
BCVectorInterface::showMe( bool verbose, std::ostream & out ) const
{
    ASSERT_PRE( this->isFinalized(), "BC Vector should be finalized before being accessed." );
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    out << "BC Vector Interface: " << std::endl;
    out << "number of interface vector Dof : " << _M_nbTotalDof << std::endl;
    out << "==>Interface Dof :\n";
    _M_dofIn->showMe( verbose, out );  // no showMe(..) in Miguel's DofInterface
    out << "+++++++++++++++++++++++++++++++" << std::endl;
    return out;
}

BCVectorBase*
createBCVectorInterface( BCVectorBase const* __bc )
{
    return new BCVectorInterface( ( BCVectorInterface const& )*__bc );
}
// register BCFunctionBase in factory for cloning
const bool __bcint = FactoryCloneBCVector::instance().registerProduct( typeid(BCVectorInterface), &createBCVectorInterface );

}
