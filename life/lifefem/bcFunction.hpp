/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): M.A. Fernandez
             Christophe Prud'homme <christophe.prudhomme@epfl.ch>
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
   \author M.A Fernandez
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-12
 */
#ifndef __bcFunction_H
#define __bcFunction_H 1

#include <boost/function.hpp>

#include <life/lifecore/singleton.hpp>
#include <life/lifecore/factory.hpp>

namespace LifeV
{

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
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& )> function_type;

    //! Default constructor
    /*!
      The user must supply a function by calling setFunction(..)
    */
    BCFunctionBase() {}

    //! Constructing from a user defined function
    /*!
      \param g the user defined function
    */
    BCFunctionBase( function_type g );

    //! Constructing from a user defined functor
    /*!
      \param bcf user defined functor
    */
    BCFunctionBase( const BCFunctionBase& bcf );

    virtual ~BCFunctionBase() {}

    //! Set the function
    /*!
      \param g the user defined function
    */
    void setFunction( function_type g );

    //! Get the function
    /*!
      \return the user defined function
    */
    function_type& Function();

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
    function_type _M_g;
};

/*!

 \class BCFunctionMixte

 Class (STL functor) that holds the user defined fonctions for a Mixte BC.

  The data funcitions given by the user must have the following declaration
 Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
*/
class BCFunctionMixte : public BCFunctionBase
{
public:

    typedef BCFunctionBase::function_type function_type;

    //! Default constructor
    /*!
      The user must supply a function by calling setFunction(..)
    */
    BCFunctionMixte() {}

    //! Constructing from user defined functions
    /*!
      \param g user defined function
      \param coef user defined function
    */
    BCFunctionMixte( function_type g, function_type coef );

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
    void setFunctions_Mixte( function_type g, function_type coef );

    //! Get the function mixte
    /*!
      \return user defined function
    */
    function_type& Functions_Mixte();

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
    function_type _M_coef;
};

typedef LifeV::singleton< LifeV::factoryClone< BCFunctionBase > > FactoryCloneBCFunction;

/* much similar to BCFunctionBase but different func prototipe
   for bc, I derive BCFunctionBaseUDepending from BCFunctionBase
   only becouse BCBase and BCHandler can work with us with little change,
   really I don't like that class BCFunctionBase has fixed function_type
   prototype in the base class, but I conform
 */
class BCFunctionUDepBase
{
public:
    //Real g(t,x,y,z,ID,U)
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID&, const Real& )> function_type;

    BCFunctionUDepBase(function_type g );
    BCFunctionUDepBase(const BCFunctionUDepBase& bcf );

    void setFunction(function_type g);
    Real operator()(const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i, const Real& U ) const;

protected:
    function_type _M_g;
};

class BCFunctionUDepMixte: public BCFunctionUDepBase
{
public:
    typedef BCFunctionUDepBase::function_type function_type;

    BCFunctionUDepMixte(function_type g,function_type coef);
    BCFunctionUDepMixte(const BCFunctionUDepMixte& bcf);

    void setFunctions_Mixte(function_type g, function_type coef );

    Real coef(const Real& t, const Real& x, const Real& y,
              const Real& z, const ID& i, const Real& U ) const;
private:
    function_type _M_coef;
};

typedef LifeV::singleton< LifeV::factoryClone< BCFunctionUDepBase > > FactoryCloneBCFunctionUDep;

/*!

 \class BCFunctionDirectional

 Class (STL functor) that holds the user defined fonctions for a directional Dirichlet bc

  The data funcitions given by the user must have the following declaration
  Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
*/
class BCFunctionDirectional
        :
        public BCFunctionBase
{
public:

    typedef BCFunctionBase::function_type function_type;

    //! Default constructor
    /*!
      The user must supply a function by calling setFunction(..)
    */
    BCFunctionDirectional() {}


    //! Constructing from user defined functions
    /*!
      \param g user defined function
      \param vectFct user defined function to defined the direction where the Dirichlet condition is imposed
    */
    BCFunctionDirectional( function_type g, function_type vectFct );

    //! Constructing from a user defined functor
    /*!
      \param bcf user defined functor
    */
    BCFunctionDirectional( const BCFunctionDirectional& bcf );


    //! Set the functions in the mixte case (beware : plural!)
    /*!
      \param g : the user defined function
      \param vectFct user defined function to defined the direction where the Dirichlet condition is imposed
    */
    void setFunctions_Directional( function_type g, function_type vectFct );

    //! Set the functions in the mixte case (beware : plural!)
    /*!
      \return user defined function to defined the direction where the Dirichlet condition is imposed
    */
    function_type& Functions_Directional();

    //! Method to call the auxiliary user defined function
    /*!
      \param t time
      \param x coordinate
      \param y coordinate
      \param z coordinate
      \param i component of the vector function
      \return i-component of the user defined fonction evaluted in (t,x,y,z)
    */
    Real vectFct( const Real& t, const Real& x, const Real& y,
                  const Real& z, const ID& i ) const;
private:
    //! user defined function
    function_type _M_vectFct;
};

}//End of namespace LifeV

#endif
