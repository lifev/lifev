#ifndef _FEFUNCTION_
#define _FEFUNCTION_ 1

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
    @brief File containing the FEFct, FEScalarFct, FEVectorFct and FEMatrixFct classes

    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author M. Kern <michel.kern@inria.fr>
    @date 29-03-2011

*/

#include <life/lifefem/FEField.hpp>

namespace LifeV
{

//! FEFct - This class gives an abstract implementation of a finite element function on finite elements fields.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  This class represents the concept of a function associated to a set of finite element fields.
  <br>
  This class is implemented as a simple container which stores pointers to finite element fields, both
  scalar and vector.
  <br>
  <br>
  The access operator (), which is abstract and virtual, is the manipulator of the fields.
  <br>
  For example the scalar function \f$ f \f$ is a function of the scalar field \f$ S \f$ and the vector 
  field \f$ u \f$, that is \f$ f = f(S, u)\f$. And, for example, the function \f$ f \f$ is
  \f[
	f(S, u) = S^2 + \sin( u_0 )
  \f]
  then create a new class MyFun which inherit from FEScalarFct where the access operator (), using the methods scalarField and vectorField, can be written as
  @code
  const Real S = scalarField(0)( iElem, point, time );
  const Vector u = vectorField(0)( iElem, point, time );
  return std::pow(S, 2) + std::sin( u[0] );
  @endcode
  Partial specialization of the type returned by the access operator is scalar, vector and matrix, implemented
  in the derived classes. The implementation of user defined classes can be derived from one of the specialized 
  classes.
  <br>
  The evaluation of the field in a given point is given in the access operator (). The implementation 
  of this method requires also the id of the geometrical element in the mesh, which is less general but
  more efficient compared to the case where the id of the geometrical element is not given.

  @todo Add a method without the element id, less efficient but more flexible.
*/
template < typename Mesh, typename Map, typename FunctionType >
class FEFct
{
public:

    //! @name Public Types
    //@{

    typedef FEScalarField<Mesh, Map> FEScalarField_Type;
    typedef FEVectorField<Mesh, Map> FEVectorField_Type;

    typedef typename boost::shared_ptr<FEScalarField_Type> FEScalarFieldPtr_Type;
    typedef typename boost::shared_ptr<FEVectorField_Type> FEVectorFieldPtr_Type;

    typedef typename std::vector<FEScalarFieldPtr_Type> FEScalarFieldPtrContainer_Type;
    typedef typename std::vector<FEVectorFieldPtr_Type> FEVectorFieldPtrContainer_Type;

    typedef typename FEScalarField_Type::point_Type point_Type;

    //@}

    //! @name Opertors
    //@{

    //! Abstract virtual access opertor.
    /*!
      In the case of derived and specialized classes it evaluates the function on a given point 
      in a given element
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @return The type template of value of the function.
    */
    virtual FunctionType operator() ( const UInt& iElem, 
                                      const point_Type& P, 
                                      const Real& time = 0. ) const = 0;

    //@}

    //! @name Methods
    //@{

    //! Add an external scalar field to the function.
    /*!
      @param fieldPtr Pointer of a scalar field which as to be added.
    */
    inline void addScalarField ( const FEScalarFieldPtr_Type& fieldPtr )
    {
        M_scalarFields.push_back( fieldPtr );
    }

    //! Add an external vector field to the function.
    /*!
      @param fieldPtr Pointer of a vector field which as to be added.
    */
    inline void addVectorField ( const FEVectorFieldPtr_Type& fieldPtr )
    {
        M_vectorFields.push_back( fieldPtr );
    }

    //@}

    //! @name Get Methods
    //@{

    //! Return the i-th scalar field.
    /*!
      @param i Index of the scalar field.
      @return Constant FEScalarField_Type reference stored of index i.
    */
    inline const FEScalarField_Type& scalarField ( const UInt& i ) const
    {
        ASSERT ( i < M_scalarFields.size() , "Index out of range.");
        return *( M_scalarFields[ i ] );
    }

    //! Return the i-th scalar field.
    /*!
      @param i Index of the scalar field.
      @return FEScalarField_Type reference stored of index i.
    */
    inline FEScalarField_Type& scalarField ( const UInt& i )
    {
        ASSERT ( i < M_scalarFields.size() , "Index out of range.");
        return *( M_scalarFields[ i ] );
    }

    //! Return the i-th vector field.
    /*!
      @param i Index of the vector field.
      @return Constant FEVectorField_Type reference stored of index i.
    */
    inline const FEVectorField_Type& vectorField ( const UInt& i ) const
    {
        ASSERT ( i < M_vectorFields.size() , "Index out of range.");
        return *( M_vectorFields[ i ] );
    }

    //! Return the i-th vector field.
    /*!
      @param i Index of the vector field.
      @return FEVectorField_Type reference stored of index i.
    */
    inline FEVectorField_Type& vectorField ( const UInt& i )
    {
        ASSERT ( i < M_vectorFields.size() , "Index out of range.");
        return *( M_vectorFields[ i ] );
    }
   
    //@}

private:

    //! Vector of pointers of scalar fields.
    FEScalarFieldPtrContainer_Type M_scalarFields;

    //! Vector of pointers of vector fields.
    FEVectorFieldPtrContainer_Type M_vectorFields;

};

//! FEScalarFct - This class gives an abstract implementation of a scalar finite element function on finite elements fields.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  Template partial specialization of the class FEFct for the scalar type functions.
*/
template < typename Mesh, typename Map >
class FEScalarFct :
public FEFct < Mesh, Map, Real >
{
public:

    //! @name Public Types
    //@{

    typedef typename FEFct < Mesh, Map, Real >::point_Type point_Type;

    //@}
};

//! FEVectorFct - This class gives an abstract implementation of a vector finite element function on finite elements fields.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  Template partial specialization of the class FEFct for the vector type functions.
*/
template < typename Mesh, typename Map >
class FEVectorFct :
public FEFct < Mesh, Map, Vector >
{
public:

    //! @name Public Types
    //@{

    typedef typename FEFct < Mesh, Map, Vector >::point_Type point_Type;

    //@}
};

//! FEMatrixFct - This class gives an abstract implementation of a matrix finite element function on finite elements fields.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  Template partial specialization of the class FEFct for the matrix type functions.
*/
template < typename Mesh, typename Map >
class FEMatrixFct :
public FEFct < Mesh, Map, Matrix >
{
public:

    //! @name Public Types
    //@{

    typedef typename FEFct < Mesh, Map, Matrix >::point_Type point_Type;

    //@}
};

} // namespace LifeV

#endif
