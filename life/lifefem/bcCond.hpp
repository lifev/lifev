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
  \file bcCond.h
  \brief classes to handle boundary conditions.
  \version 1.0
  \author M.A. Fernandez
  \date 06/2002
        11/2002 Introduction of data vectors and method bdUpdate (moved from Dof class).

  This file contains the classes which may be used to store boundary
  conditions. A boundary condition objet will have the following
  elements:

  -# a name identifying a specific BC,

  -# a flag identifying a specific part of the mesh boundary,

  -# a type (Essential, Natural, Mixte),

  -# a mode of implementation (Scalar, Full, Component, Normal,
     Tangential),

  -# a functor holding the data function,

  -# a bool vector  describing the components involved in this boundary condition

  -# a list of pointers to identifiers allowing the user to know to
     which DOF the boundary condition applies.
*/

#ifndef __BCCOND_HH__
#define __BCCOND_HH__

#include <set>
#include <map>

#include <boost/shared_ptr.hpp>

#include "lifeV.hpp"
#include "identifier.hpp"
#include "markers.hpp"
#include "dof.hpp"
#include "currentFE.hpp"
#include "currentBdFE.hpp"

#include "bcVector.hpp"
#include "bcFunction.hpp"


namespace LifeV
{
/** Boundary condition basic types
    Essential, Natural or Mixte
*/
enum BCType{Essential, Natural, Mixte};
//	,UDepEssential,UDepNatural,UDepMixte};

/** Type for boundary conditions application modes

-# Scalar: for scalar problems

-# Full: for vector problems involving all components

-# Component: for vector problems not involving all compontents

-# Normal: for vector problems dealing with the normal component

-# Tangential: for vector problems dealing with tangential components
*/
enum BCMode{Scalar, Full, Component, Normal, Tangential};

/*!
  \class BCBase

  Base class which holds the boundary condition information

  For each boundary condtion the user must give

  -# a name,

  -# a mesh flag,

  -# a type,

  -# a mode,

  -# a data BCFuncion,

  -# three (or two in 2D) bools describing the components involved in
  this boundary condition.

  Finally the list of pointers to identifiers will be updated in the
  Dof class (\c BCHandler::bdUpdate method).

  \warning The idea is to not use inheritance from this class
*/


//!
class BCBase
{
public:

    friend class BCHandler;

    //! iterator type in the identifiers list
    typedef std::set<IdentifierBase*, identifierComp>::iterator IDIterator0;
    typedef std::vector<IdentifierBase*>::iterator IDIterator;


    /**
       Constructor for BC

       \param name the name of the boundary condition

       \param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies

       \param type the boundary condition type: Natural, Essential, Mixte

       \param mode the boundary condition mode: Scalar, Full,
       Component, Normal, Tangential

       \param bcf the function holding the user defined function
       involved in this boundary condition

       \param std::vector<ID> storing the list of components involved
       in this boundary condition
    */
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionBase& bcf,
            const std::vector<ID>& comp );

    /**
       Constructor for BC without components for Scalar, Tangential or
       Normal mode problems

       \param name the name of the boundary condition

       \param flag the mesh flag identifying the part of the mesh where
       the boundary condiion applies

       \param type the boundary condition type: Natural, Essential, Mixte

       \param mode the boundary condition mode: Scalar, Full, Normal, Tangential

       \param bcf the function holding the user defined function
       involved in this boundary condition
    */
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionBase& bcf );

    /**
       Constructor for BC without list of components for Full mode problems

       \param name the name of the boundary condition

       \param flag the mesh flag identifying the part of the mesh where
       the boundary condiion applies

       \param type the boundary condition type: Natural, Essential,
       Mixte

       \param mode the boundary condition mode: Scalar, Full, Normal, Tangential

       \param bcf the function holding the user defined function
       involved in this boundary condition

       \param nComp the number of componets involved in this boundary
       condition
    */
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionBase& bcf,
            const UInt& nComp );

    //! Constructor for BC with data vector
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condition applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
      \param bcv data vector
      \param std::vector<ID> storing the list of components involved in this boundary condition
    */
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCVectorBase& bcv,
            const std::vector<ID>& comp );

    //! Constructor for BC with data vector, without components for Scalar, Tangential or Normal  mode problems
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
      \param bcv data vector
    */
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCVectorBase& bcv );

    //! Constructor for BC with data vector, without list of components for Full mode problems
    /*!
      \param name the name of the boundary condition
      \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
      \param type the boundary condition type: Natural, Essential, Mixte
      \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
      \param bcv data vector
      \param nComp the number of componets involved in this boundary condition
    */
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCVectorBase& bcv,
            const UInt& nComp );

/* constructors for BCFunctionUDepBase ... */
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionUDepBase& bcf,
            const std::vector<ID>& comp );
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionUDepBase& bcf);
    BCBase( const std::string& name,
            const EntityFlag& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionUDepBase& bcf,
            const UInt& nComp );



    //! Copy constructor for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
    /*!
      \param BCb a boundary condition
    */
    BCBase( const BCBase& BCb );

    //! Assignment operator for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
    /*!
      \param BCb a boundary condition
    */
    BCBase & operator=( const BCBase& );

    //! Destructor (we have a vector of pointers to ID's and a pointer to user defined functions)
    ~BCBase();

    //! Returns the BC name
    std::string name() const;

    //! Returns the BC associated flag
    EntityFlag flag() const;

    //! Returns the BC type
    BCType type() const;

    //! Returns the BC mode
    BCMode mode() const;

    //! Returns the number of components involved in this boundary condition
    UInt numberOfComponents() const;

    /**
       Returns the global i-th component involved in the boundary
       condition

       \param i the specified "local" component (from 1 to
       numberOfComponents)

       \return true if the specified component component is involved in
       the BC
    */
    ID component( const ID i ) const;

    /**
       Returns wether the list is finalised and the vector of ID's is
       then accessible.
    */
    bool finalised() const;

    //! Overloading function operator by calling the _M_bcf() user specified function
    /*!
      \param t time
      \param x coordinate
      \param y coordinate
      \param z coordinate
      \param i component of the vector function
      \return i-component of the user defined fonction evaluted in (t,x,y,z)
    */
    Real operator() ( const Real& t, const Real& x, const Real& y,
                      const Real& z, const ID& i ) const;
 
    /* new overloading for BCFunctionUDepending */
    Real operator() ( const Real& t, const Real& x, const Real& y,
                      const Real& z, const ID& i, const Real& u ) const;


    //! Returns a pointer  to the user defined STL functor
    const BCFunctionBase* pointerToFunctor() const;

    //! Returns a pointer to the BCVector
    const BCVectorBase* pointerToBCVector() const;

    const BCFunctionUDepBase* pointerToFunctorUDep() const;

    //! True if a data vector has been provided
    bool dataVector() const;

    //! use vector boundary conditions
    void setBCVector( BCVectorBase& __v );

    //! use function boundary conditions
    void setBCFunction( BCFunctionBase& __f );

    //! use function boundary conditions
    void setBCFunction( BCFunctionUDepBase& __f );

    bool isUDep() const;


    //! Returns a pointer  to the i-th elements in the (finalised) list
    //! (counting from 1 ' a la FORTRAN')
    const IdentifierBase* operator() ( const ID& i ) const;


    //! Overloading function operator by calling the
    /*!
      \param iDpof global dof number
      \param iComp component number
    */
    Real operator() ( const ID& iDof, const ID& iComp ) const;


    //! Returns the value of the mixte coefficient (in BC Vector)
    Real mixteCoef() const;

    //! Returns the value of the mixte coefficient vector (in BC Vector)
    Real MixteVec( const ID& iDof, const ID& iComp ) const;

    //! Returns a pointer to the i-th elements in the (finalised) list
    //! (counting from 0 ' a la C')
    const IdentifierBase* operator[] ( const Index_t& i ) const;

    //! Add a new indentifier in the list
    void addIdentifier( IdentifierBase* );

    //! Returns the liste size
    UInt list_size() const;

    //! Output
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;

    //! overloaded operator allowing decreasing ordering operations
    friend bool operator<( const BCBase& a, const BCBase& b )
        {
            return ( a.type() > b.type() );
        }

    //! overloaded operator allowing finding operations
    friend bool operator==( const BCBase& a, const EntityFlag flag )
        {
            return a.flag() == flag;
        }

private:
    bool _M_isUDep;

    //! name identifying a specific BC
    std::string _M_name;

    //! flag identifying a specific part of the mesh boundary
    EntityFlag _M_flag;

    //! the boundary condition type
    BCType _M_type;

    //! the boundary condition mode of application
    BCMode _M_mode;

    //! Pointer to a user defined functor
    boost::shared_ptr<BCFunctionBase> _M_bcf;

    boost::shared_ptr<BCFunctionUDepBase> _M_bcfUDep;


    //! Pointer to a user given data vector
    boost::shared_ptr<BCVectorBase> _M_bcv;

    //! True is a data vector has been provided
    bool _M_dataVector;

    //! the list of involved in this BC
    std::vector<ID> _M_comp;

    //! set of pointers to identifiers allowing the user to get hold the DOF
    //! to which the BC applies
    std::set<boost::shared_ptr<IdentifierBase>, identifierComp> list0;

    //! container for id's when the list is finalised
    std::vector<boost::shared_ptr<IdentifierBase> > _M_idList;

    //! true, when idlist updated
    bool _M_finalised;

    //! Transfert between list and vector containers
    void finalise();
};





}

#endif
