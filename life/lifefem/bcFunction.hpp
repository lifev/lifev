/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-12

  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
/**
   \file bcFunction.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-12
 */
#ifndef __bcFunction_H
#define __bcFunction_H 1


namespace LifeV
{
// ============ BCFunctionBase ================

/*!

 \class BCFunctionBase

 Base class (STL functor) that holds the function used for imposing BC.

  The data functions given by the user must have the following declaration
  Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
  We can use inheritance to hold specific boundary condition data. See, for instance,
  Mixed boundary conditions.

*/

class BCFunctionBase
{
public:
    //! Type for a generic user defined  function
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );

    //! Default constructor
    /*!
      The user must supply a function by calling setFunction(..)
    */
    BCFunctionBase()
    {}
    ;

    //! Constructing from a user defined function
    /*!
      \param g the user defined function
    */
    BCFunctionBase( Function g );

    //! Constructing from a user defined functor
    /*!
      \param bcf user defined functor
    */
    BCFunctionBase( const BCFunctionBase& bcf );

    //! Set the function
    /*!
      \param g the user defined function
    */
    void setFunction( Function g );

    //! Overloading function operator by calling _g
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

protected:
    //! user defined function
    Function _g;
};


/*!

 \class BCFunctionMixte

 Class (STL functor) that holds the user defined fonctions for a Mixte BC.

  The data funcitions given by the user must have the following declaration
 Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
*/
class BCFunctionMixte
    :
    public BCFunctionBase
{
public:

    typedef BCFunctionBase::Function Function;

    //! Default constructor
    /*!
      The user must supply a function by calling setFunction(..)
    */
    BCFunctionMixte()
        {}


    //! Constructing from user defined functions
    /*!
      \param g user defined function
      \param coef user defined function
    */
    BCFunctionMixte( Function g, Function coef );

    //! Constructing from a user defined functor
    /*!
      \param bcf user defined functor
    */
    BCFunctionMixte( const BCFunctionMixte& bcf );


    //! Set the functions in the mixte case (beware : plural!)
    /*!
      \param g : the user defined function
      \param coef : user defined function
    */
    void setFunctions_Mixte( Function g, Function coef );


    //! Method to call the auxiliary user defined function
    /*!
      \param t time
      \param x coordinate
      \param y coordinate
      \param z coordinate
      \param i component of the vector function
      \return i-component of the user defined fonction evaluted in (t,x,y,z)
    */
    Real coef( const Real& t, const Real& x, const Real& y,
               const Real& z, const ID& i ) const;
private:
    //! user defined function
    Function _coef;
};

}

#endif /* __bcFunction_H */
