#ifndef _FEFIELD_
#define _FEFIELD_ 1
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
    @brief File containing the FEAbstractField, FEScalarField and FEVectorField classes

    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author M. Kern <michel.kern@inria.fr>
    @date 29-03-2011

*/

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifefem/FESpace.hpp>
#include <life/lifearray/VectorEpetra.hpp>

namespace LifeV
{

typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::matrix<Real> Matrix;

//! FEAbstractField - This class gives an abstract implementation of a finite element field.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  This class represents the concept of a field associated to a finite element space.
  <br>
  This class is implemented as a simple container which stores a reference to a finite element space
  and a pointer to a vector. The vector represents the value of the field at each degree of freedom.
  The type of the value, e.g. scalar, is a template parameter which is specialized in the derived classes 
  for the case of scalar fields and vector fields.
  <br>
  The evaluation of the field in a given point is given in the access operator (). The implementation 
  of this method requires also the id of the geometrical element in the mesh, which is less general but
  more efficient compared to the case where the id of the geometrical element is not given.

  @todo Add a method without the element id, less efficient but more flexible.
*/
template < typename Mesh, typename Map, typename FunctionType >
class FEAbstractField
{

public:

    //! @name Public Types
    //@{

    typedef FESpace < Mesh, Map > FESpace_Type;

    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr < vector_Type > vectorPtr_Type;

    typedef std::vector < Real > point_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Full constructor for the class.
    /*!
      @param fESpace Finite element space where the field is defined.
      @param vector vector witch represent the solution.
    */
    FEAbstractField ( FESpace_Type& fESpace, const vectorPtr_Type& vector ):
    M_FESpace       ( fESpace ),
    M_vector        ( vector )
    {}

    //! Constructor for the class without the vector.
    /*!
      Create the field with the input finite element space and a new vector.
      The type of the map of the vector is Repeated (default value) or Unique.
      @param fESpace Finite element space where the field is defined.
      @param mapType Specify wether the map is Unique or Repeated. Default value: Repeated.
    */
    FEAbstractField ( FESpace_Type& fESpace, const MapEpetraType& mapType = Repeated ):
    M_FESpace       ( fESpace ),
    M_vector        ( new vector_Type ( M_FESpace.map(), mapType ) )
    {}

    //! Copy constructor.
    /*!
      @param field Finite element field to be copied.
    */
    FEAbstractField ( const FEAbstractField& field ):
    M_FESpace       ( field.M_FESpace ),
    M_vector        ( field.M_vector )
    {}

    //! Virtual destructor.
    virtual ~FEAbstractField () {};

    //@}

    //! @name Opertors
    //@{

    //! Abstract virtual access opertor.
    /*!
      Evaluate the field on a given point in a given element.
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @return The template type of the value of the field.
    */
    virtual FunctionType operator() ( const UInt& iElem, 
                                      const point_Type& point, 
                                      const Real& time = 0.) const = 0;

    //@}    

    //! @name Get Methods
    //@{
    
    //! Return the finite element space.
    /*!
      @return Constant FESpace_Type reference of the finite element space.
    */
    inline const FESpace_Type& getFESpace () const
    {
        return M_FESpace;
    }

    //! Return the finite element space.
    /*!
      @return FESpace_Type reference of the finite element space.
    */
    inline FESpace_Type& getFESpace ()
    {
        return M_FESpace;
    }

    //! Return the vector where the value are stored.
    /*!
      @return Constant vectorPtr_Type reference of the vector.
    */
    inline const vectorPtr_Type& getVector () const
    {
        return M_vector;
    }

    //! Return the vector where the value are stored.
    /*!
      @return vectorPtr_Type reference of the vector.
    */
    inline vectorPtr_Type& getVector ()
    {
        return M_vector;
    }

    //@}

protected:

    //! Reference of the finite element space.
    FESpace_Type& M_FESpace;

    //! Pointer of a vector where the value are stored.
    vectorPtr_Type M_vector;

};

//! FEScalarField - This class gives an abstract implementation of a finite element scalar field.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  This class, derived from FEAbstractField, implements the concept of a scalar field associated to a 
  finite element space.
*/
template < typename Mesh, typename Map >
class FEScalarField :
public FEAbstractField < Mesh, Map, Real >
{
public:

    //! @name Public Types
    //@{

    typedef FEAbstractField < Mesh, Map, Real > FEAbstractField_Type;

    typedef typename FEAbstractField_Type::FESpace_Type FESpace_Type;
    typedef typename FEAbstractField_Type::vectorPtr_Type vectorPtr_Type;

    typedef typename FEAbstractField_Type::point_Type point_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Full constructor for the class.
    /*!
      @param fESpace Finite element space where the field is defined.
      @param vector vector witch represent the solution.
    */
    FEScalarField ( FESpace_Type& fESpace, const vectorPtr_Type& vector ):
    FEAbstractField_Type ( fESpace, vector )
    {}

    //! Constructor for the class without the vector.
    /*!
      Create the field with the input finite element space and a new vector.
      The type of the map of the vector is Repeated (default value) or Unique.
      @param fESpace Finite element space where the field is defined.
      @param mapType Specify wether the map is Unique or Repeated. Default value: Repeated.
    */
    FEScalarField ( FESpace_Type& fESpace, const MapEpetraType& mapType = Repeated ):
    FEAbstractField_Type ( fESpace, mapType )
    {}

    //! Copy constructor.
    /*!
      @param field Finite element field to be copied.
    */
    FEScalarField ( const FEScalarField& field ):
    FEAbstractField_Type ( field )
    {}

    //! Virtual destructor.
    virtual ~FEScalarField () {};

    //@}

    //! @name Opertors
    //@{

    //! Access opertor.
    /*!
      Evaluate the field on a given point in a given element.
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @return The scalar value of the field.
    */
    inline virtual Real operator() ( const UInt& iElem, 
                                     const point_Type& P, 
                                     const Real& time = 0. ) const
    {
        // Evaluate the field using a method implemented in FESpace
        return this->M_FESpace.feInterpolateValue( iElem, *(this->M_vector), P );
    }
   
    //!}

};

//! FEVectorField - This class gives an abstract implementation of a finite element vector field.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr> 

  This class, derived from FEAbstractField, implements the concept of a vector field associated to a 
  finite element space.
*/
template < typename Mesh, typename Map >
class FEVectorField :
public FEAbstractField < Mesh, Map, Vector >
{

public:

    //! @name Public Types
    //@{

    typedef FEAbstractField < Mesh, Map, Vector > FEAbstractField_Type;

    typedef typename FEAbstractField_Type::FESpace_Type   FESpace_Type;
    typedef typename FEAbstractField_Type::vectorPtr_Type vectorPtr_Type;

    typedef typename FEAbstractField_Type::point_Type point_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Full constructor for the class.
    /*!
      @param fESpace Finite element space where the field is defined.
      @param vector vector witch represent the solution.
    */
    FEVectorField ( FESpace_Type& fESpace, const vectorPtr_Type& vector ):
    FEAbstractField_Type ( fESpace, vector )
    {}

    //! Constructor for the class without the vector.
    /*!
      Create the field with the input finite element space and a new vector.
      The type of the map of the vector is Repeated (default value) or Unique.
      @param fESpace Finite element space where the field is defined.
      @param mapType Specify wether the map is Unique or Repeated. Default value: Repeated.
    */
    FEVectorField ( FESpace_Type& fESpace, const MapEpetraType& mapType = Repeated ):
    FEAbstractField_Type ( fESpace, mapType )
    {}

    //! Copy constructor.
    /*!
      @param field Finite element field to be copied.
    */
    FEVectorField ( const FEVectorField& field ):
    FEAbstractField_Type ( field )
    {}

    //! Virtual destructor.
    virtual ~FEVectorField () {};

    //@}

    //! @name Opertors
    //@{

    //! Access opertor.
    /*!
      Evaluate the field on a given point in a given element.
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @return The vector value of the field.
    */
    inline virtual Vector operator() ( const UInt& iElem, 
                                       const point_Type& P, 
                                       const Real& time = 0. ) const
    {
        Vector value( P.size() );
        // Evaluate each component of the vector field with a method in the class FESpace.
        for ( UInt i = 0; i < P.size(); ++i )
        {
            value(i) = this->M_FESpace.feInterpolateValue( iElem, *(this->M_vector), P, i );
        }
        return value;
    }

};

typedef FEScalarField < RegionMesh3D < LinearTetra >, MapEpetra > FEScalarFieldTetra;

typedef FEVectorField < RegionMesh3D < LinearTetra >, MapEpetra > FEVectorFieldTetra;

} // namespace LifeV

#endif
