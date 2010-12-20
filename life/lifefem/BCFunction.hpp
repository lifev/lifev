//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief File contains BCManageNormal class for handling normal essential boundary conditions

	@author Miguel Fernandez <miguel.fernandez@inria.fr>
    @contributor Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    @date 10-12-2004
 *///@HEADER

#ifndef BCFUNCTION_H
#define BCFUNCTION_H 1

#include <boost/function.hpp>

#include <life/lifecore/FactorySingleton.hpp>
#include <life/lifecore/Factory.hpp>

namespace LifeV
{
//! BCFunctionBase - class that holds the function used for prescribing boundary conditions.
/*!
   @author Miguel Fernandez

  This class holds the function used for prescribing Essential or Natural boundary conditions.
  The data functions given by the user must have the following signature.

  @verbatim
  Real f(const Real& time, const Real& x, const Real& y, const Real& z, const ID& component).
  @endverbatim

  In the case of Essential boundary condition, we want to prescribe
  @verbatim
  u =  f
  @endverbatim
  on part of the boundary,
  where @c u is the solution and @c f the user defined function.
  In the case of Natural boundary condition, we want to add the to the right hand side of the equation the following term:
  @verbatim
  ( f, phi )_bd
  @endverbatim

  Here @c ( . , . )_bd  denote the L^2 inner product on the  boundary, phi is the test function

  Functions f and is set using the correct constructor or using @c setFunction(f). <br>
  To get the function f use @c getFunction(), to evaluate it use the @c operator().

  This is the base class for other BCFunctionXXX classes.
  Inheritance is used to hold specific boundary condition data.

*/
class BCFunctionBase
{
public:

    //! @name Public Types
    //@{

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& )> function_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
      The user must supply a function by calling setFunction(..)
    */
    BCFunctionBase() {}

    //! Constructor for a user defined function
    /*!
      @param userDefinedFunction the user defined function
    */
    BCFunctionBase( function_Type userDefinedFunction );

    //! Copy Constructor
    /*!
      @param bcFunctionBase The BCFunctionBase
    */
    BCFunctionBase( const BCFunctionBase& bcFunctionBase );


    //! Destructor
    virtual ~BCFunctionBase() {}

    //@}

    //! @name Operators
    //@{

    //! Assignment Operator
    /*!
        @param bcFunctionBase The BCFunctionBase
        @return Reference to a new BCFunctionBase object with the same content of bcFunctionBase
     */
    virtual BCFunctionBase& operator=( const BCFunctionBase& bcFunctionBase );


    //! Overloading function operator by calling M_userDefinedFunction
    /*!
      @param t Time
      @param x Coordinate x
      @param y Coordinate y
      @param z Coordinate z
      @param component The Component of the vector function
      @return The selected component of the user defined function evaluated in (t,x,y,z)
    */
    inline Real operator() ( const Real& t, const Real& x, const Real& y,
                             const Real& z, const ID& component ) const
    { return M_userDefinedFunction( t, x, y, z, component ); }

    //@}

    //! @name Set Methods
    //@{


    //! Set the user defined function
    /*!
      @param userDefinedFunction the user defined function
    */
    inline void setFunction( function_Type userDefinedFunction ) {M_userDefinedFunction = userDefinedFunction;}

    //@}

    //! @name Get Methods
    //@{

    //! Get the function
    /*!
      @return Reference to M_userDefinedFunction
    */
    inline const function_Type& Function() const {return M_userDefinedFunction;}

    //@}



protected:
    //! user defined function
    function_Type M_userDefinedFunction;
};


typedef LifeV::FactorySingleton< LifeV::FactoryClone< BCFunctionBase > > FactoryCloneBCFunction;






//! BCFunctionRobin - class that holds the function used for prescribing Robin boundary conditions.
/*!
   @author Miguel Fernandez

  This class holds the functions used for prescribing Robin boundary conditions.
  This class is derived by @c BCFunctionBase class

  The data functions given by the user must have the following signature

  @verbatim
  Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& component)
  @endverbatim

  In order to prescribe Robin boundary condition we want to add to the system matrix the boundary mass term
  @verbatim
  (coeff u, phi)_bd
  @endverbatim
  and to add to the right hand side the term
  @verbatim
  (f, phi)_bd
  @endverbatim
  where u is the solution, phi the test function and @c ( , )_bd the L2 inner product on the boundary

  Functions f and coeff are set using the correct constructor or using setFunctionRobin(f, coeff).<br>
  To get the function f use @c getFunction(), to evaluate it use the @c operator().  <br>
  To get the function coeff use @c getFunction_Robin(), to evaluate it use the method @c coef(...).

*/

class BCFunctionRobin : public BCFunctionBase
{
public:

    //! @name Public Types
    //@{

    typedef BCFunctionBase::function_Type function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Default constructor
    /*!
      The user must supply the functions by calling setFunction_Robin(..)
    */
    BCFunctionRobin() {}

    //! Constructing from user defined functions
    /*!
      @param rightHandSideFunction The user defined function for @c f
      @param massTermFunction The user defined function for @c coeff
    */
    BCFunctionRobin( const function_Type& rightHandSideFunction, const function_Type& massTermFunction  );


    //! Copy Constructor
    /*!
      @param bcFunctionUDepRobin The BCFunctionRobin object
    */
    BCFunctionRobin( const BCFunctionRobin& bcFunctionUDepRobin );

    //! Destructor
    virtual ~BCFunctionRobin() {}


    //@}

    //! @name Operators
    //@{

    //! Assignment operator
    /*!
      @param bcFunctionUDepRobin The BCFunctionRobin object
      @return Reference to a new BCFunctionRobin object which is a copy of bcFunctionMixt
    */
    BCFunctionRobin&
    operator=( const BCFunctionRobin& bcFunctionUDepRobin );


    //@}

    //! @name Methods
    //@{

    //! evaluate the user defined function M_robinBoundaryMassCoeffFunction
    /*!
      @param t Time
      @param x Coordinate
      @param y Coordinate
      @param z Coordinate
      @param component component of the vector function
      @return selected component of the user defined for @ coeff evaluated in (t,x,y,z)
    */
    Real coef( const Real& t, const Real& x, const Real& y,
               const Real& z, const ID& component ) const
    { return M_robinBoundaryMassCoeffFunction(t, x, y, z, component); }

    //@}


    //! @name Set Methods
    //@{


    //! Set the functions
    /*!
      @param rightHandSideFunction The user defined function for @c f
      @param massTermFunction The user defined function for @c coeff
    */
    void setFunctions_Robin( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction );

    //@}


    //! @name Get Methods
    //@{


    //! Get the user defined function M_robinBoundaryMassCoeffFunction
    /*!
      @return user defined function for @c coeff
    */
    inline const function_Type& Functions_Robin() const {return M_robinBoundaryMassCoeffFunction;}

    //@}

private:

    //! user defined function for the boundary mass coefficient in Robin conditions
    function_Type M_robinBoundaryMassCoeffFunction;
};


//! BCFunctionUDepBase - class that holds the function used for prescribing boundary conditions.
/*!
   @author Miguel Fernandez

  This class holds the function used for prescribing Essential or Natural boundary conditions in the case in which the function depend on a FE vector
  (usually the solution at the previous iteration).

  The data functions given by the user must have the following signature

  @verbatim
  Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& component, const Real& feVectorEvaluatedInThisPoint)
  @endverbatim

  In the case of Essential boundary condition, we want to prescribe
  @verbatim
  u =  f
  @endverbatim
  on part of the boundary,
  where @c u is the solution and @c f the user defined function.
  In the case of Natural boundary condition, we want to add the to the right hand side of the equation the following term:
  @verbatim
  ( f, phi )_bd
  @endverbatim

  Functions f is set using the correct constructor or using @c setFunction(f). <br>
  To get the function f use @c getFunction(), to evaluate it use the @c operator().
  */
class BCFunctionUDepBase
{
public:

    //! @name Public Types
    //@{

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID&, const Real& )> function_Type;

    //@}

    //! Empty Constructor
    /*!
      The user must supply a function by calling setFunction(..)
    */
    BCFunctionUDepBase() {}

    //! Constructor for a user defined function
    /*!
      @param userDefinedFunction the user defined function
    */
    BCFunctionUDepBase(const function_Type& userDefinedFunction );


    //! Copy Constructor
    /*!
      @param bcFunctionBase The BCFunctionBase
    */
    BCFunctionUDepBase(const BCFunctionUDepBase& bcFunctionUDepBase );


    //! Destructor
    virtual ~BCFunctionUDepBase() {}

    //! @name Operators
    //@{


    //! Assignment Operator
    /*!
    	@param bcFunctionUDepBase The BCFunctionUDepBase
    	@return Reference to a new BCFunctionUDepBase object with the same content of bcFunctionUDepBase
     */
    virtual BCFunctionUDepBase& operator=( const BCFunctionUDepBase& bcFunctionUDepBase);


    //! Overloading function operator by calling M_userDefinedFunction
    /*!
      @param t Time
      @param x Coordinate x
      @param y Coordinate y
      @param z Coordinate z
      @param component The Component of the vector function
      @param feVectorEvaluatedInThisPoint The FE vector evaluated in the point (x, y, z) at time t.
      @return The selected component of the user defined function evaluated in (t,x,y,z)
    */
    inline Real operator()(const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& component, const Real& feVectorEvaluatedInThisPoint ) const
    { return M_userDefinedFunction( t, x, y, z, component, feVectorEvaluatedInThisPoint ); }

    //@}


    //! @name Set Methods
    //@{


    //! Set the user defined function
    /*!
      @param userDefinedFunction the user defined function
    */
    inline void setFunction(const function_Type& userDefinedFunction ) { M_userDefinedFunction = userDefinedFunction;}

    //@}


    //! @name Get Methods
    //@{


    //! Get the function
    /*!
      @return Reference to M_userDefinedFunction
    */
    inline const function_Type& Function() const {return M_userDefinedFunction;}

    //@}

protected:
    //! user defined function
    function_Type M_userDefinedFunction;
};


//! BCFunctionUDepRobin - class that holds the function used for prescribing Robin boundary conditions.
/*!
   @author Miguel Fernandez

  This class holds the functions used for prescribing Robin boundary conditions in the case in which the function depend on a FE vector
  (usually the solution at the previous iteration).
  This class is derived by BCFunctionUDepBase class

  In order to prescribe Robin boundary condition we want to add to the system matrix the boundary mass term
  @verbatim
  (coeff u, phi)_bd
  @endverbatim
  and to add to the right hand side the term
  @verbatim
  (f, phi)_bd
  @endverbatim
  where u is the solution, phi the test function and @c ( , )_bd the L2 inner product on the boundary

  Functions f and coeff are set using the correct constructor or using setFunctionRobin(f, coeff). <br>
  To get the function f use @c getFunction(), to evaluate it use the @c operator(). <br>
  To get the function coeff use @c getFunction_Robin(), to evaluate it use the method @c coef(...).

  Function f and coeff have the following signature:
  @verbatim
  Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& component, const Real& feVectorEvaluatedInThisPoint)
  @endverbatim


*/
class BCFunctionUDepRobin: public BCFunctionUDepBase
{
public:

    //! @name Public Types
    //@{

    typedef BCFunctionUDepBase::function_Type function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    /*!
      The user must supply the functions by calling setFunction_Robin(..)
    */
    BCFunctionUDepRobin() {}

    //! Constructing from user defined functions
    /*!
      @param rightHandSideFunction The user defined function for f
      @param massTermFunction The user defined function for coeff
    */
    BCFunctionUDepRobin( const function_Type& rightHandSideFunction, const function_Type& massTermFunction  );


    //! Copy Constructor
    /*!
      @param bcFunctionUDepRobin The BCFunctionUDepRobin object
    */
    BCFunctionUDepRobin( const BCFunctionUDepRobin& bcFunctionUDepRobin );


    //! Destructor
    virtual ~BCFunctionUDepRobin() {}


    //@}


    //! @name Operators
    //@{


    //! Assignment operator
    /*!
      @param bcFunctionUDepRobin The BCFunctionUDepRobin object
      @return Reference to a new BCFunctionUDepRobin object which is a copy of bcFunctionMixt
    */
    BCFunctionUDepRobin&
    operator=( const BCFunctionUDepRobin& bcFunctionUDepRobin );


    //@}


    //! @name Methods
    //@{


    //! evaluate the user defined function M_robinBoundaryMassCoeffFunction
    /*!
      @param t Time
      @param x Coordinate
      @param y Coordinate
      @param z Coordinate
      @param component component of the vector function
      @param feVectorEvaluatedInThisPoint The FE vector evaluated in the point (x, y, z) at time t.
      @return selected component of the user defined for @c coeff evaluated in (t,x,y,z)
    */
    Real coef( const Real& t, const Real& x, const Real& y,
               const Real& z, const ID& component, const Real& feVectorEvaluatedInThisPoint ) const
    { return M_robinBoundaryMassCoeffFunction(t, x, y, z, component, feVectorEvaluatedInThisPoint); }


    ///@}

    //! @name Set Methods
    //@{



    //! Set the functions
    /*!
      @param rightHandSideFunction The user defined function for f
      @param massTermFunction The user defined function for @c coeff
    */
    void setFunctions_Robin( const function_Type& rightHandSideFunction, const function_Type& massCoeffFunction );

    //@}


    //! @name Get Methods
    //@{


    //! Get the user defined function M_robinBoundaryMassCoeffFunction
    /*!
      @return user defined function for @c coeff
    */
    inline const function_Type& Functions_Robin() const {return M_robinBoundaryMassCoeffFunction;}

    //@}

private:
    //! user defined function for the boundary mass coefficient in Robin conditions
    function_Type M_robinBoundaryMassCoeffFunction;
};

typedef LifeV::FactorySingleton< LifeV::FactoryClone< BCFunctionUDepBase > > FactoryCloneBCFunctionUDep;

/*!

 \class BCFunctionDirectional

 Class (STL functor) that holds the user defined fonctions for a directional Dirichlet bc

  The data funcitions given by the user must have the following declaration
  Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
*/

//! BCFunctionUDepBase - class that holds the function used for prescribing boundary conditions.
/*!
   @author Miguel Fernandez

  This class holds the function used for prescribing Essential boundary conditions along a direction.

  The data functions given by the user must have the following signature

  @verbatim
  Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& component)
  @endverbatim

  We want to prescribe < u, d > = f on part of the boundary. Here < u, d > is the projection of the solution u along the versor d, which is a given function

  Functions f and d are set using the correct constructor or using @c setFunction_Directional(f,d).

  To get the function d use @c Function_Directional(), to evaluate it use the method @c vectFct(...). <br>
  To get the function f use @c getFunction(), to evaluate it use the @c operator().
  */

class BCFunctionDirectional
        :
        public BCFunctionBase
{
public:

    //! @name Public Types
    //@{

    typedef BCFunctionBase::function_Type function_Type;

    //@}


    //! @name Constructors & Destructor

    //! Default constructor
    /*!
      The user must supply a function by calling setFunction_Directional(..)
    */
    BCFunctionDirectional() {}


    //! Constructing from user defined functions
    /*!
      @param userDefinedFunctional The user defined function
      @param userDefinedVersorsFunction user defined function for returning versors along which the essential boundary condition will be prescribed
    */
    BCFunctionDirectional( const function_Type& userDefinedFunctional, const function_Type& userDefinedVersorsFunction );

    //! Copy Constructor
    /*!
      @param bcFunctionDirectional The BCFunctionDirectional
    */
    BCFunctionDirectional( const BCFunctionDirectional& bcFunctionDirectional );

    ~BCFunctionDirectional() {}

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
      @param bdFunctionDirectional The BCFunctionDirectional object
      @return Reference to a new BCFunctionDirectional object which is a copy of bcFunctionDirectional
    */
    BCFunctionDirectional&
    operator=( const BCFunctionDirectional& bcFunctionDirectional );


    //@}


    //! @name Methods
    //@{

    //! Evaluate the versors' function
    /*!
      \param t Time
      \param x Coordinate
      \param y Coordinate
      \param z Coordinate
      \param component The component of the vectors function
      \return The selected component of the versors' function evaluated in (t,x,y,z)
    */
    inline Real vectFct( const Real& t, const Real& x, const Real& y,
                         const Real& z, const ID& component ) const
    { return M_userDefinedVersorsFunction(t, x, y, z, component); }

    //@}



    //! @name Set Methods
    //@{


    //! Set the functions
    /*!
       @param userDefinedFunctional User defined function
       @param userDefinedVersorsFunction User defined function for returning versors along which the essential boundary condition will be prescribed
    */
    void setFunctions_Directional( const function_Type& userDefinedFunctional, const function_Type& userDefinedVersorsFunction );

    //@}


    //! @name Get Methods
    //@{


    //! Get the versors' function
    /*!
      @return User defined function for returning versors along which the essential boundary condition will be prescribed
    */
    inline const function_Type& Functions_Directional() const { return M_userDefinedVersorsFunction; }

    //@}


private:
    //! user defined function returning versors along which the essential boundary condition will be prescribed
    function_Type M_userDefinedVersorsFunction;
};

}//End of namespace LifeV

#endif
