/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003 LifeV Team
  
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

//! Default Constructor (the user must call setBCVector(..))
BCVector_Interface::BCVector_Interface():
  _MixteCoef(0.0),_finalized(false)
{}

//! Constructor
BCVector_Interface::BCVector_Interface(Vector& vec, UInt nbTotalDof, 
				       DofInterfaceBase& dofIn):
  _vec(&vec),_dofIn(&dofIn),_nbTotalDof(nbTotalDof),
  _MixteCoef(0.0),_finalized(true)
{}


//!set the BC vector (after default construction)
void BCVector_Interface::setBCVector_Interface(Vector& vec, UInt nbTotalDof, 
					       DofInterfaceBase& dofIn)
{
  ASSERT_PRE( !_finalized, "BC Vector cannot be set twice."); //!remove this...
  _vec        = &vec ;
  _dofIn      = &dofIn;
  _nbTotalDof = nbTotalDof; 
  _finalized  = true;
}


//! set the Mixte coefficient
void BCVector_Interface::setBCVec_MixteCoef( const Real& coef ){
  _MixteCoef = coef;
}

//! This method returns the value to be imposed in the component iComp of the dof iDof
Real BCVector_Interface::operator()(const ID& iDof, const ID& iComp) const {
  ASSERT_PRE(_finalized, "BC Vector should be finalized before being accessed.");
  return (*_vec)( (iComp-1)* _nbTotalDof + _dofIn->getInterfaceDof(iDof) - 1 ); 
} 
  
//! Return the value of the Mixte coefficient
Real BCVector_Interface::MixteCoef() const{
  return _MixteCoef;
}


//! Assignment operator for BCVector_Interface 
BCVector_Interface & BCVector_Interface::operator=(const BCVector_Interface & BCv) {
  _vec        = BCv._vec;
  _dofIn      = BCv._dofIn;
  _nbTotalDof = BCv._nbTotalDof;
  _MixteCoef  = BCv._MixteCoef;
  _finalized  = BCv._finalized;
  
  return *this;
}


//! Output
ostream& BCVector_Interface::showMe(bool verbose, ostream & out) const {
  ASSERT_PRE(_finalized, "BC Vector should be finalized before being accessed.");
  out << "+++++++++++++++++++++++++++++++"<<endl;
  out << "BC Vector Interface: "  << endl;
  out << "number of interface vector Dof : " << _nbTotalDof << endl;
  out << "==>Interface Dof :\n";
  _dofIn->showMe( verbose, out );  // no showMe(..) in Miguel's DofInterface
  out << "+++++++++++++++++++++++++++++++" << endl;
  return out;
} 

