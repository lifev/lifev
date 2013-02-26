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
    @brief File containing the FEFunction class

    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author M. Kern <michel.kern@inria.fr>
    @date 29-03-2011

*/

#include <lifev/core/fem/FEField.hpp>

namespace LifeV
{

//! FEFunction - This class gives an abstract implementation of a finite element function on finite elements fields.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  This class represents the concept of a function associated to a set of finite element fields.
  <br>
  This class is implemented as a simple container which stores pointers to finite element fields, both
  scalar and vector.
  <br>
  <br>
  The method eval, which is abstract and virtual, is the manipulator of the fields.
  <br>
  For example the scalar function \f$ f \f$ is a function of the scalar field \f$ S \f$ and the vector
  field \f$ u \f$, that is \f$ f = f(S, u)\f$. And, for example, the function \f$ f \f$ is
  \f[
    f(S, u) = S^2 + \sin( u_0 )
  \f]
  then create a new class MyFun which inherit from FEFunction<MeshType, MapType, Real> where the eval method,
  using the methods scalarField and vectorField, can be written as
  @code
  const Real S = scalarField(0).eval( iElem, point, time );
  const Vector u = vectorField(0).eval( iElem, point, time );
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
template < typename MeshType, typename MapType, typename ReturnType >
class FEFunction
{
public:

    //! @name Public Types
    //@{

    typedef MeshType mesh_Type;
    typedef MapType map_Type;
    typedef ReturnType return_Type;

    typedef FEField < mesh_Type, map_Type, return_Type > FEField_Type;
    typedef FEScalarField < mesh_Type, map_Type > FEScalarField_Type;
    typedef FEVectorField < mesh_Type, map_Type > FEVectorField_Type;

    typedef typename boost::shared_ptr < FEField_Type > FEFieldPtr_Type;
    typedef typename boost::shared_ptr < FEScalarField_Type > FEScalarFieldPtr_Type;
    typedef typename boost::shared_ptr < FEVectorField_Type > FEVectorFieldPtr_Type;

    typedef typename std::vector < FEScalarFieldPtr_Type > FEScalarFieldPtrContainer_Type;
    typedef typename std::vector < FEVectorFieldPtr_Type > FEVectorFieldPtrContainer_Type;

    typedef typename FEField_Type::point_Type point_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Empty constructor for the class.
    FEFunction ()
    {}

    //! Virtual destructor.
    virtual ~FEFunction ()
    {}

    //! @name Methods
    //@{

    //! Abstract virtual eval function.
    /*!
      In the case of derived and specialized classes it evaluates the function on a given point
      in a given element
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @param time Time in the evaluation.
      @return The type template of value of the function.
    */
    virtual return_Type eval ( const UInt& iElem,
                               const point_Type& P,
                               const Real& time = 0. ) const = 0;

    //! Interpolate to a field
    /*!
      Fill the input field with the value of the function. Useful for exporting pourpose.
      @param interpolatedFct Field where the function is interpolated.
      @param time Time in the evaluation.
      @note The method does not work for matrix value functions.
    */
    void interpolate ( FEField_Type& interpolatedFct,
                       const Real& time = 0. ) const
    {
        interpolatedFct.getFESpace().interpolate ( this, interpolatedFct.getVector(), time );
    }

    //! Add an external scalar field to the function.
    /*!
      @param fieldPtr Pointer of a scalar field which as to be added.
    */
    void addScalarField ( const FEScalarFieldPtr_Type& fieldPtr )
    {
        M_scalarFields.push_back ( fieldPtr );
    }

    //! Add an external vector field to the function.
    /*!
      @param fieldPtr Pointer of a vector field which as to be added.
    */
    void addVectorField ( const FEVectorFieldPtr_Type& fieldPtr )
    {
        M_vectorFields.push_back ( fieldPtr );
    }

    //@}

    //! @name Get Methods
    //@{

    //! Return the i-th scalar field.
    /*!
      @param i Index of the scalar field.
      @return Constant FEScalarFieldPtr_Type reference stored of index i.
    */
    const FEScalarFieldPtr_Type& scalarFieldPtr ( const UInt& i ) const
    {
        ASSERT ( i < M_scalarFields.size() , "Index out of range.");
        return M_scalarFields[ i ];
    }

    //! Return the i-th scalar field.
    /*!
      @param i Index of the scalar field.
      @return Constant FEScalarField_Type reference stored of index i.
    */
    const FEScalarField_Type& scalarField ( const UInt& i ) const
    {
        ASSERT ( i < M_scalarFields.size() , "Index out of range.");
        return * ( M_scalarFields[ i ] );
    }

    //! Return the i-th vector field.
    /*!
      @param i Index of the vector field.
      @return Constant FEVectorFieldPtr_Type reference stored of index i.
    */
    const FEVectorFieldPtr_Type& vectorFieldPtr ( const UInt& i ) const
    {
        ASSERT ( i < M_vectorFields.size() , "Index out of range.");
        return M_vectorFields[ i ];
    }

    //! Return the i-th vector field.
    /*!
      @param i Index of the vector field.
      @return Constant FEVectorField_Type reference stored of index i.
    */
    const FEVectorField_Type& vectorField ( const UInt& i ) const
    {
        ASSERT ( i < M_vectorFields.size() , "Index out of range.");
        return * ( M_vectorFields[ i ] );
    }

    //@}

private:

    //! Vector of pointers of scalar fields.
    FEScalarFieldPtrContainer_Type M_scalarFields;

    //! Vector of pointers of vector fields.
    FEVectorFieldPtrContainer_Type M_vectorFields;

};

} // namespace LifeV

#endif
