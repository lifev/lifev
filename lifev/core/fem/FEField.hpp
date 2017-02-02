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
    @brief File containing the FEField, FEScalarField and FEVectorField classes

    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author M. Kern <michel.kern@inria.fr>
    @date 29-03-2011

*/


#include <boost/numeric/ublas/matrix.hpp>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorSmall.hpp>

namespace LifeV
{

typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::matrix<Real> Matrix;

//! FEField - This class gives an abstract implementation of a finite element field.
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
  @note The attribute which represent a point is an array of three elements. Change in future with a standard container.
  @note If you need an Raviart-Thomas field you need to create a scalar field and not a vector field, since the field
        contains the value of a Raviart-Thomas times the outward unit normal of the face element.
*/
template < typename MeshType, typename MapType, typename ReturnType >
class FEField
{

public:

    //! @name Public Types
    //@{

    typedef MeshType mesh_Type;
    typedef MapType map_Type;
    typedef ReturnType return_Type;

    typedef FESpace < mesh_Type, map_Type > FESpace_Type;
    typedef std::shared_ptr < FESpace_Type > FESpacePtr_Type;

    typedef VectorEpetra vector_Type;
    typedef std::shared_ptr < vector_Type > vectorPtr_Type;

    typedef Vector3D point_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Empty constructor for the class.
    FEField ( )
    {}

    //! Full constructor for the class.
    /*!
      @param fESpace Finite element space where the field is defined.
      @param vector vector witch represent the solution.
    */
    FEField ( const FESpacePtr_Type& fESpace, const vectorPtr_Type& vector ) :
        M_FESpace ( fESpace ),
        M_vector ( vector )
    {}

    //! Constructor for the class without the vector.
    /*!
      Create the field with the input finite element space and a new vector.
      The type of the map of the vector is Repeated (default value) or Unique.
      @param fESpace Finite element space where the field is defined.
      @param mapType Specify wether the map is Unique or Repeated. Default value: Repeated.
    */
    FEField ( const FESpacePtr_Type& fESpace, const MapEpetraType& mapType = Repeated ) :
        M_FESpace ( fESpace ),
        M_vector ( new vector_Type ( M_FESpace->map(), mapType ) )
    {}

    //! Copy constructor.
    /*!
      @param field Finite element field to be copied.
    */
    FEField ( const FEField& field ) :
        M_FESpace ( field.M_FESpace ),
        M_vector ( field.M_vector )
    {}

    //! Virtual destructor.
    virtual ~FEField () {};

    //@}

    //! @name Operators
    //@{

    FEField& operator = ( const FEField& field )
    {
        if ( this != & field )
        {
            M_FESpace = field.M_FESpace;
            M_vector = field.M_vector;
        }
        return *this;
    }

    //@}

    //! @name Methods
    //@{

    //! Abstract virtual eval function.
    /*!
      Evaluate the field on a given point in a given element.
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @return The template type of the value of the field.
    */
    virtual return_Type eval ( const UInt& iElem,
                               const point_Type& point,
                               const Real& time = 0.) const = 0;

    //@}

    //! @name Set Methods
    //@{

    //! Set the FESpace
    /*!
      @param fESpace Pointer of a FESpace.
      @param createVector True if the vector is created, false otherwise. Default value: false.
      @param mapType Specify wether the map is Unique or Repeated. Default value: Repeated.
    */
    void setFESpacePtr ( const FESpacePtr_Type& fESpace,
                         const bool createVector = false,
                         const MapEpetraType& mapType = Repeated )
    {
        M_FESpace = fESpace;

        if ( createVector )
        {
            M_vector.reset ( new vector_Type ( M_FESpace->map(), mapType ) );
        }
    }

    //! Set the epetra vector
    /*!
      @param vector Pointer of an epetra vector.
    */
    void setVectorPtr ( const vectorPtr_Type& vector )
    {
        M_vector = vector;
    }

    //! Set the epetra vector
    /*!
      @param vector The epetra vector.
    */
    void setVector ( const vector_Type& vector )
    {
        M_vector.reset ( new vector_Type ( vector ) );
    }

    //! Set both FESpace and vector
    /*!
      @param fESpace Pointer of a FESpace.
      @param vector Pointer of a vector.
    */
    void setUp ( const FESpacePtr_Type& fESpace, const vectorPtr_Type& vector )
    {
        this->setFESpacePtr ( fESpace, false );
        this->setVectorPtr ( vector );
    }

    //! Clean the field.
    void cleanField ()
    {
        this->getVector().zero();
    }

    //@}

    //! @name Get Methods
    //@{

    //! Return the finite element space.
    /*!
      @return Constant FESpace_Type reference of the finite element space.
    */
    const FESpace_Type& getFESpace () const
    {
        return *M_FESpace;
    }

    //! Return the pointer to the finite element space.
    /*!
      @return Constant FESpacePtr_Type reference of the finite element space.
    */
    const FESpacePtr_Type& getFESpacePtr () const
    {
        return M_FESpace;
    }

    //! Return the finite element space.
    /*!
      @return FESpace_Type reference of the finite element space.
    */
    FESpace_Type& getFESpace ()
    {
        return *M_FESpace;
    }

    //! Return the vector where the value are stored.
    /*!
      @return Constant vectorPtr_Type reference of the vector.
    */
    const vectorPtr_Type& getVectorPtr () const
    {
        return M_vector;
    }

    //! Return the vector where the value are stored.
    /*!
      @return vector_Type reference of the vector.
    */
    vector_Type& getVector ()
    {
        return *M_vector;
    }

    //! Return the vector where the value are stored.
    /*!
      @return Constant vector_Type reference of the vector.
    */
    const vector_Type& getVector () const
    {
        return *M_vector;
    }

    //@}

protected:

    //! Reference of the finite element space.
    FESpacePtr_Type M_FESpace;

    //! Pointer of a vector where the value are stored.
    vectorPtr_Type M_vector;

};

//! FEScalarField - This class gives an abstract implementation of a finite element scalar field.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  This class, derived from FEField, implements the concept of a scalar field associated to a
  finite element space.
*/
template < typename MeshType, typename MapType >
class FEScalarField :
    public FEField < MeshType, MapType, Real >
{
public:

    //! @name Public Types
    //@{

    typedef MeshType mesh_Type;
    typedef MapType map_Type;
    typedef Real return_Type;

    typedef FEField < mesh_Type, map_Type, return_Type > FEField_Type;

    typedef typename FEField_Type::FESpacePtr_Type FESpacePtr_Type;
    typedef typename FEField_Type::vectorPtr_Type vectorPtr_Type;

    typedef typename FEField_Type::point_Type point_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Empty constructor for the class.
    FEScalarField ( )
    {}

    //! Full constructor for the class.
    /*!
      @param fESpace Finite element space where the field is defined.
      @param vector vector witch represent the solution.
    */
    FEScalarField ( const FESpacePtr_Type& fESpace, const vectorPtr_Type& vector ) :
        FEField_Type ( fESpace, vector )
    {}

    //! Constructor for the class without the vector.
    /*!
      Create the field with the input finite element space and a new vector.
      The type of the map of the vector is Repeated (default value) or Unique.
      @param fESpace Finite element space where the field is defined.
      @param mapType Specify wether the map is Unique or Repeated. Default value: Repeated.
    */
    FEScalarField ( const FESpacePtr_Type& fESpace, const MapEpetraType& mapType = Repeated ) :
        FEField_Type ( fESpace, mapType )
    {}

    //! Copy constructor.
    /*!
      @param field Finite element field to be copied.
    */
    FEScalarField ( const FEScalarField& field ) :
        FEField_Type ( field )
    {}

    //! Virtual destructor.
    virtual ~FEScalarField () {};

    //@}

    //! @name Methods
    //@{

    //! Eval function.
    /*!
      Evaluate the field on a given point in a given element.
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @return The scalar value of the field.
    */
    virtual return_Type eval ( const UInt& iElem,
                               const point_Type& P,
                               const Real& time = 0. ) const;

    //!}

};

template < typename MeshType, typename MapType >
inline Real
FEScalarField < MeshType, MapType >::
eval ( const UInt& iElem, const point_Type& P, const Real& /*time*/ ) const
{
    // Evaluate the field using a method implemented in FESpace
    return this->M_FESpace->feInterpolateValue ( iElem, * (this->M_vector), P );
} // eval


//! FEVectorField - This class gives an abstract implementation of a finite element vector field.
/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author M. Kern <michel.kern@inria.fr>

  This class, derived from FEField, implements the concept of a vector field associated to a
  finite element space.
*/
template < typename MeshType, typename MapType >
class FEVectorField :
    public FEField < MeshType, MapType, Vector >
{

public:

    //! @name Public Types
    //@{

    typedef MeshType mesh_Type;
    typedef MapType map_Type;
    typedef Vector return_Type;

    typedef FEField < mesh_Type, map_Type, return_Type > FEField_Type;

    typedef typename FEField_Type::FESpacePtr_Type FESpacePtr_Type;
    typedef typename FEField_Type::vectorPtr_Type vectorPtr_Type;

    typedef typename FEField_Type::point_Type point_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Empty constructor for the class.
    FEVectorField ( )
    {}

    //! Full constructor for the class.
    /*!
      @param fESpace Finite element space where the field is defined.
      @param vector vector witch represent the solution.
    */
    FEVectorField ( const FESpacePtr_Type& fESpace, const vectorPtr_Type& vector ) :
        FEField_Type ( fESpace, vector )
    {}

    //! Constructor for the class without the vector.
    /*!
      Create the field with the input finite element space and a new vector.
      The type of the map of the vector is Repeated (default value) or Unique.
      @param fESpace Finite element space where the field is defined.
      @param mapType Specify wether the map is Unique or Repeated. Default value: Repeated.
    */
    FEVectorField ( const FESpacePtr_Type& fESpace, const MapEpetraType& mapType = Repeated ) :
        FEField_Type ( fESpace, mapType )
    {}

    //! Copy constructor.
    /*!
      @param field Finite element field to be copied.
    */
    FEVectorField ( const FEVectorField& field ) :
        FEField_Type ( field )
    {}

    //! Virtual destructor.
    virtual ~FEVectorField () {};

    //@}

    //! @name Opertors
    //@{

    //! Eval function.
    /*!
      Evaluate the field on a given point in a given element.
      @param iElem Element id in the mesh.
      @param point Point where the field is evaluated, vector format.
      @return The vector value of the field.
    */
    virtual return_Type eval ( const UInt& iElem,
                               const point_Type& P,
                               const Real& time = 0. ) const;

};

template < typename MeshType, typename MapType >
inline Vector
FEVectorField < MeshType, MapType >::
eval ( const UInt& iElem, const point_Type& P, const Real& /*time*/ ) const
{
    Vector value ( P.size() );

    // Evaluate each component of the vector field with a method in the class FESpace.
    for ( UInt i = 0; i < P.size(); ++i )
    {
        value (i) = this->M_FESpace->feInterpolateValue ( iElem, * (this->M_vector), P, i );
    }
    return value;

} // eval

typedef FEScalarField < RegionMesh < LinearTetra >, MapEpetra > FEScalarFieldTetra;

typedef FEVectorField < RegionMesh < LinearTetra >, MapEpetra > FEVectorFieldTetra;

} // namespace LifeV

#endif
