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
  \file bcVector.h
  \brief classes to handle data vectors for boundary conditions.
  \version 1.0
  \author M.A. Fernandez
  \date 11/2002

  \version 1.1
  \author V. Martin
  \date 02/2003

  This file contains the classes which may be used to store boundary
  conditions.
*/

#ifndef __BCVECTOR__
#define __BCVECTOR__

#include "lifeV.hpp"
#include "dofInterfaceBase.hpp"
#include "vecUnknown.hpp"


namespace LifeV
{

// ============ BCVecto_Base ================

/*!

 \class BCVector_Base

 Base class that holding data vector for boundary conditions

*/

class BCVector_Base {
 public:

  //! Constructor
  /*!
    \param vec data vector holding data
    \param nbTotalDof number of total dof in the vector of data
    \param dofIn dofInterfaceBase object holding the connections between the interface dofs of the
           data vector and those of the associated to the boundary conditions
  */
  BCVector_Base(Vector& vec, const UInt nbTotalDof);

  //! Default Constructor (the user must call setBCVector(..))
  BCVector_Base();

  //! Do nothing destructor
  virtual ~BCVector_Base() {}

  //! set the Mixte coefficient
  void setMixteCoef( const Real& coef );

  //! This method returns the value to be imposed in the component iComp of the dof iDof
  /*!
    \param iDof the number of the Dof
    \param iComp the number of the component
  */
   virtual Real operator()(const ID& iDof, const ID& iComp) const =0;

   //! Return the value of the Mixte coefficient
   Real MixteCoef() const;


   //! Output
   virtual ostream &  showMe(bool verbose=false, ostream & out=cout) const =0;

 protected:

   //! The data vector
   Vector* _vec;

   //! Number of total dof in the vector of data
   UInt _nbTotalDof;

   //! Coefficient for mixte boundary conditions (Robin)
   /*! For the moment, it is the same for all the entries of the data vector.
    */
   Real _MixteCoef;

   //! true when the BCVector is updated (and can be used)
   bool _finalized;

};



// ============ BCVector ================

/*!

 \class BCVector

 Class that holds a user data vector for boundary conditions

*/

class BCVector:
  public BCVector_Base {
 public:

  //! Constructor
  /*!
    \param vec data vector holding data
    \param nbTotalDof number of total dof in the vector of data
           data vector and those of the associated to the boundary conditions
  */
  BCVector( Vector& vec, UInt nbTotalDof);

  //! Default Constructor (the user must call setvector(..))
  BCVector();

  //! set the BC vector (after default construction)
  void setvector( Vector& vec, UInt nbTotalDof);

  //! This method returns the value to be imposed in the component iComp of the dof iDof
  /*!
    \param iDof the number of the Dof
    \param iComp the number of the component
  */
  Real operator()(const ID& iDof, const ID& iComp) const;


  //! Assignment operator for BCVector_Interface
  BCVector & operator=(const BCVector & BCv);

  //! Output
  ostream &  showMe(bool verbose=false, ostream & out=cout) const;
};

// ============ BCVector_Interface ================

/*!

 \class BCVector_Interface

 Class that holds a user data vector for boundary conditions on interfaces

 The data funcitions given by the user must have the following declaration
 Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
*/

class BCVector_Interface:
  public BCVector_Base {
 public:

  //! Constructor
  /*!
    \param vec data vector holding data
    \param nbTotalDof number of total dof in the vector of data
    \param dofIn dofInterfaceBase object holding the connections between the interface dofs of the
           data vector and those of the associated to the boundary conditions
  */
  BCVector_Interface( Vector& vec, UInt nbTotalDof, DofInterfaceBase& dofIn);

  //! Default Constructor (the user must call setBCVector(..))
  BCVector_Interface ();

  //! set the BC vector (after default construction)
  void setvector( Vector& vec, UInt nbTotalDof, DofInterfaceBase& dofIn );

  //! This method returns the value to be imposed in the component iComp of the dof iDof
  /*!
    \param iDof the number of the Dof
    \param iComp the number of the component
  */
   Real operator()(const ID& iDof, const ID& iComp) const;


   //! Assignment operator for BCVector_Interface
   BCVector_Interface & operator=(const BCVector_Interface & BCv);


   //! Output
   ostream &  showMe(bool verbose=false, ostream & out=cout) const;

 protected:

  //! DofInterfaceBase object holding the connections between the interface dofs
  DofInterfaceBase* _dofIn;

};
}
#endif
