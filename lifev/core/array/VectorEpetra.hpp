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
    @brief VectorEpetra

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 04-10-2006
 */

#ifndef _EPETRAVECTOR_HPP_
#define _EPETRAVECTOR_HPP_


#include <Epetra_FEVector.h>
#include <Epetra_Export.h>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MapEpetra.hpp>

namespace LifeV
{

//! VectorEpetra - The Epetra Vector format Wrapper
/*!
    @author Gilles Fourestey, Simone Deparis, Cristiano Malossi

    The VectorEpetra class provides a general interface for the Epetra_Vector of Trilinos.

    Visit http://trilinos.sandia.gov for more informations about Epetra_Vector.
 */
class VectorEpetra
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_FEVector                  vector_type;
    typedef boost::shared_ptr< vector_type > Vector_PtrType;
    typedef Real                             data_type;
    typedef Epetra_CombineMode               combineMode_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
      @param mapType Specify whether the map is Unique or Repeated
     */
    VectorEpetra ( const MapEpetraType& mapType = Unique );

    //! Constructor - Using Maps
    /*!
      This constructor uses maps to build the vector
      @param map Map to be used to split the vector between the processors
      @param mapType Specify wether the map is Unique or Repeated
     */
    VectorEpetra ( const MapEpetra& map, const MapEpetraType& mapType = Unique );

    //! Constructor - Using Maps
    /*!
      This constructor uses maps to build the vector
      @param mapPtr Pointer to the map which has to be used to split the vector between the processors
      @param mapType Specify wether the map is Unique or Repeated
     */
    VectorEpetra ( const boost::shared_ptr< MapEpetra >& mapPtr,
                   const MapEpetraType& mapType = Unique );

    //! Copy constructor
    /*!
      @param vector Vector to be used to create the new vector
     */
    VectorEpetra ( const VectorEpetra& vector );

    //! Copy constructor
    /*!
      @param vector Vector to be used to create the new vector
      @param mapType Specify whether the map is Unique or Repeated
     */
    VectorEpetra ( const VectorEpetra& vector, const MapEpetraType& mapType );

    //! Copy constructor
    /*!
      @param vector Vector to be used to create the new vector
      @param mapType Specify wether the map is Unique or Repeated
      @param combineMode Parameter used during the copy, and not in subsequent calls.
    */
    VectorEpetra ( const VectorEpetra& vector,
                   const MapEpetraType& mapType,
                   const combineMode_Type& combineMode );

    //! Copy constructor
    /*!
      Copies vector to FEvector that comes as Multivector
      @param vector Vector to be used to create the new vector
      @param map Map to be used to split the vector between the processors
      @param mapType Specify wether the map is Unique or Repeated
     */
    VectorEpetra ( const Epetra_MultiVector& vector,
                   const boost::shared_ptr< MapEpetra > map,
                   const MapEpetraType& mapType );

    //! Copy constructor
    /*!
      Copies vector to a vector which resides only on the processor "reduceToProc".
      @param vector Vector to be used to create the new vector
      @param map Map to be used to split the vector between the processors
      @param mapType Specify wether the map is Unique or Repeated
     */
    VectorEpetra ( const VectorEpetra& vector, const Int& reduceToProc );

    //! Destructor
    virtual ~VectorEpetra() {}

    //@}


    //! @name Operators
    //@{
    //! Access operators
    /*!
      @param row Index of the entry to be accessed
     */
    data_type& operator[] ( const UInt row );

    //! Access operators
    /*!
      @param row Index of the entry to be accessed
     */
    const data_type& operator[] ( const UInt row ) const;

    //! Access operators
    /*!
      @param row Index of the entry to be accessed
     */
    data_type& operator() ( const UInt row );

    //! Access operators
    /*!
      @param row Index of the entry to be accessed
     */
    const data_type& operator() ( const UInt row ) const;

    //! Affectation operator
    /*!
      Copies the value of the given vector inside the vector.
      If the map is not the same, try to import the values.
      Calls Import with Add.
      @param vector Vector to be affected to the current vector
     */
    VectorEpetra& operator= ( const VectorEpetra& vector );

    //! Affectation operator
    /*!
      Copies the value of the given vector inside the vector.
      If the map is not the same, try to import the values.
      Calls Import with Add.
      @param vector Vector to be affected to the current vector
     */
    VectorEpetra& operator= ( const Epetra_MultiVector& vector );

    //! Affectation operator
    /*!
      Put the given scalar in each component of the vector
      @param scalar Scalar to be used to fill the current vector
     */
    VectorEpetra& operator= ( data_type scalar );

    //! Addition operator
    /*!
      Element by Element addition (if the map is not the same, try to import values)
      @param vector Vector to be added to the current vector
     */
    VectorEpetra& operator+= ( const VectorEpetra& vector );

    //! Subtraction operator
    /*!
      Element by Element subtraction (if the map is not the same, try to import values)
      @param vector Vector to be subtracted to the current vector
     */
    VectorEpetra& operator-= ( const VectorEpetra& vector );

    //! Multiplication operator
    /*!
      Element by Element multiplication (if the map is not the same, try to import values)
      @param vector Vector to be perform the multiplication
     */
    VectorEpetra& operator*= ( const VectorEpetra& vector );

    //! Division operator
    /*!
      Element by Element division (if the map is not the same, try to import values)
      @param vector Vector to be perform the division
     */
    VectorEpetra& operator/= ( const VectorEpetra& vector );

    //! Addition operator
    /*!
      Element by Element addition (do not modify the vector of the class)
      @param vector Vector to be added to the current vector
     */
    const VectorEpetra operator+ ( const VectorEpetra& vector ) const;

    //! Subtraction operator
    /*!
      Element by Element subtraction (do not modify the vector of the class)
      @param vector Vector to be subtracted to the current vector
     */
    const VectorEpetra operator- ( const VectorEpetra& vector ) const;

    //! Multiplication operator
    /*!
      Element by Element multiplication (do not modify the vector of the class)
      @param vector Vector to be perform the multiplication
     */
    const VectorEpetra operator* ( const VectorEpetra& vector ) const;

    //! Division operator
    /*!
      Element by Element division (do not modify the vector of the class)
      @param vector Vector to be perform the division
     */
    const VectorEpetra operator/ ( const VectorEpetra& vector ) const;

    //! Addition operator
    /*!
      Add a scalar value to the components of the current vector.
      (modify the vector of the class)
      @param scalar Value to be added
     */
    VectorEpetra& operator+= ( const data_type& scalar );

    //! Subtraction operator
    /*!
      Subtract a scalar value to the components of the current vector.
      (modify the vector of the class)
      @param scalar Value to be subtracted
     */
    VectorEpetra& operator-= ( const data_type& scalar );

    //! Multiplication operator
    /*!
      Multiply by a scalar value the components of the current vector.
      (modify the vector of the class)
      @param scalar Value for the multiplication
     */
    VectorEpetra& operator*= ( const data_type& scalar );

    //! Division operator
    /*!
      Division by a scalar value the components of the current vector.
      (modify the vector of the class)
      @param scalar Value for the division
     */
    VectorEpetra& operator/= ( const data_type& scalar );

    //! Operations with scalar values (do not modify the vector of the class)

    //! Addition operator
    /*!
      Add a scalar value to the components of the current vector.
      (do not modify the vector of the class)
      @param scalar Value to be added
     */
    const VectorEpetra operator+ ( const data_type& scalar ) const;

    //! Subtraction operator
    /*!
      Subtract a scalar value to the components of the current vector.
      (do not modify the vector of the class)
      @param scalar Value to be subtracted
     */
    const VectorEpetra operator- ( const data_type& scalar ) const;

    //! Multiplication operator
    /*!
      Multiply by a scalar value the components of the current vector.
      (do not modify the vector of the class)
      @param scalar Value for the multiplication
     */
    const VectorEpetra operator* ( const data_type& scalar ) const;

    //! Division operator
    /*!
      Division by a scalar value the components of the current vector.
      (do not modify the vector of the class)
      @param scalar Value for the division
     */
    const VectorEpetra operator/ ( const data_type& scalar ) const;

    //! Equality operator
    /*!
      Return a vector containing 1 where vector elements are == scalar
      @param scalar Value for the comparison.
     */
    VectorEpetra operator== ( const Real& scalar );

    //! Inequality operator
    /*!
      Return a vector containing 1 where vector elements are != scalar
      @param scalar Value for the comparison.
     */
    VectorEpetra operator!= ( const Real& scalar );

    //! Less than operator
    /*!
      Return a vector containing 1 where vector elements are < scalar
      @param scalar Value for the comparison.
     */
    VectorEpetra operator< ( const Real& scalar );

    //! Greater than operator
    /*!
      Return a vector containing 1 where vector elements are > scalar
      @param scalar Value for the comparison.
     */
    VectorEpetra operator> ( const Real& scalar );

    //! Less than or equal to operator
    /*!
      Return a vector containing 1 where vector elements are <= scalar
      @param scalar Value for the comparison.
     */
    VectorEpetra operator<= ( const Real& scalar );

    //! Greater than or equal to operator
    /*!
      Return a vector containing 1 where vector elements are >= scalar
      @param scalar Value for the comparison.
     */
    VectorEpetra operator>= ( const Real& scalar );

    //! Logical AND operator
    /*!
      Return a vector containing one where both elements are != zero
      @param vector Vector for the logical comparison.
     */
    VectorEpetra operator&& ( const VectorEpetra& vector );

    //! Logical OR operator
    /*!
      Return a vector containing one where one of the elements is != zero
      @param vector Vector for the logical comparison.
     */
    VectorEpetra operator|| ( const VectorEpetra& vector );

    //! Logical NOT operator
    /*!
      Return a vector containing one where the vector is equal to zero
     */
    VectorEpetra operator! ( void );

    //@}


    //! @name Methods
    //@{

    //! Access operators
    /**
     * It returns true if the element is present in the vector
     * @param row The element to test
     */
    bool isGlobalIDPresent (const UInt row) const;

    //! Assemble the vector
    /*!
       Gather any shared data into the non-overlapping partitioning defined by
       the Map that was passed to this vector at construction time.
       Data imported from other processors is stored on the owning processor
       with a the given operation
       @param mode Combining mode used to gather the data
    */
    Int globalAssemble ( combineMode_Type mode = Add )
    {
        return M_epetraVector->GlobalAssemble ( mode );
    }

    //! Return the local Id of a global row
    /*!
      @param row Global row Id
      <ol>
      <li> if row is mine returns the LID
      <li> if row is not mine and if the numCpus > 1, returns -1
      <li> if row is not mine and if the numCpus == 1, asserts
      </ol>
     */
    Int globalToLocalRowId ( const UInt row ) const;

    //! set zero in all the vector entries
    void zero()
    {
        M_epetraVector->PutScalar (0.);
    }

    //! Look for the given global row and set its value
    /*!
      @param row Global row Id
      <ol>
      <li> if row is mine sets this[row] = value and return true
      <li> if row is not mine and if the numCpus > 1, returns false
      <li> if row is not mine and if the numCpus == 1, asserts
      </ol>
      @param value Value to be set at the given row
      @param offset Offset used in the map numbering
     */
    bool setCoefficient ( const UInt row, const data_type& value, UInt offset = 0 );

    //! Set the row of the vector to the given value
    /*!
      @param rowsVector Vector containing the row Ids
      @param valuesVector Vector containing the values
     */
    Int setCoefficients ( std::vector< Int >& rowsVector, std::vector< Real >& valuesVector );

    //! insert a global value
    /*!
      After insertion, you will have to call global assemble
      @param GID Global Id of the row where the value should be inserted
      @param value Value to be inserted
     */
    Int sumIntoGlobalValues ( const Int GID, const Real value );

    //! Add a vector to the current vector with an offset
    /*!
      typically to do: (u,p) += p or (u,p) += u.
      Note: the nodes to add are taken by the map of vector, hence:
      <ol>
      <li> if this has a unique map: then vector should also (otherwise run time error)
      <li> if this has a repeated map: then vector should also. (otherwise wrong)
      </ol>
      @param vector Vector to be added
      @param offset Offset to shift the value
     */
    VectorEpetra& add ( const VectorEpetra& vector, const Int offset = 0 );

    //! Replace part of the vector with a given vector
    /*!
     * Typical examples are: (u,p) = p or (u,p) = u.
     * Note: the nodes to add are taken by the map of vector, hence:
     * <ol>
     *   <li> if this has a unique map: then vector should also (otherwise run time error)
     *   <li> if this has a repeated map: then vector should also. (otherwise wrong)
     * </ol>
     * @param vector given vector
     * @param offset identify the first element to be replaced
     */
    VectorEpetra& replace ( const VectorEpetra& vector, const Int& offset );

    //! Set the current vector to a subset of the given vector with an offset
    /*!
      typically to do: p = (u,p) or u = (u,p).
      Note: the nodes to add are taken by the map of this, hence:
      <ol>
      <li> if vector has a unique map: then this should also (otherwise run time error)
      <li> if vector has a repeated map: then this should also. (otherwise wrong)
      @param vector Vector of value to set the current vector
      @param offset Offset to shift the value
     */
    VectorEpetra& subset ( const VectorEpetra& vector, const UInt offset = 0 );

    //! Set the current vector to a subset of  vector with an offset.
    /*!
      similar to subset( const VectorEpetra& , const Int ), but with
      additional parameters:
      @param vector  vector from which to copy data
      @param map     map from which to select indeces to copy
      @param offset1 offset to apply to input vector
      @param offset2 offset to apply to this vector
     */
    VectorEpetra& subset ( const VectorEpetra& vector,
                           const MapEpetra& map,
                           const UInt offset1,
                           const UInt offset2 );

    //! Set the current vector to a subset of  vector with an offset.
    /*!
      similar to subset( const VectorEpetra& , const Int ), but with
      additional parameters:
      @param vector  Epetra_MultiVector, instead of VectorEpetra, from which to copy data
      @param map     map from which to select indeces to copy
      @param offset1 offset to apply to input vector
      @param offset2 offset to apply to this vector
      @param column  column of the multivector from which to extract the data
    */
    VectorEpetra& subset (const Epetra_MultiVector& vector,
                          const MapEpetra& map,
                          const UInt offset1,
                          const UInt offset2,
                          const UInt column = 0);

    //! Compute the mean value of the vector components and store it in the given variable
    /*!
      @param result Variable where the result should be stored
     */
    void meanValue ( Real* result ) const;

    //! Compute and return the norm 1
    Real norm1   () const;

    //! Compute and store the norm 1 in the given pointed variable
    /*!
      @param result Pointer on the variable  where the result should be stored
     */
    void norm1   ( Real* result ) const;

    //! Compute and store the norm 1 in the given variable
    /*!
      @param result Variable where the result should be stored
     */
    void norm1   ( Real& result ) const;

    //! Compute and return the norm 2
    Real norm2   () const;

    //! Compute and store the norm 2 in the given pointed variable
    /*!
      @param result Pointer on the variable  where the result should be stored
     */
    void norm2   ( Real* result ) const;

    //! Compute and store the norm 2 in the given variable
    /*!
      @param result Variable where the result should be stored
     */
    void norm2   ( Real& result ) const;

    //! Compute and return the norm inf
    Real normInf () const;

    //! Compute and store the norm inf in the given pointed variable
    /*!
      @param result Pointer on the variable  where the result should be stored
     */
    void normInf ( Real* result ) const;

    //! Compute and store the norm inf in the given variable
    /*!
      @param result Variable where the result should be stored
     */
    void normInf ( Real& result ) const;

    //! Compute and return the minimum value in the vector
    Real minValue() const;

    //! Compute and store the minimum value of the vector in the given pointed variable
    /*!
      @param result Pointer on the variable  where the result should be stored
     */
    void minValue ( Real* result ) const;

    //! Compute and store the minimum value of the vector in the given variable
    /*!
      @param result Variable where the result should be stored
     */
    void minValue ( Real& result ) const;

    //! Compute and return the maximum value in the vector
    Real maxValue() const;

    //! Compute and store the maximum value of the vector in the given pointed variable
    /*!
      @param result Pointer on the variable  where the result should be stored
     */
    void maxValue ( Real* result ) const;

    //! Compute and store the maximum value of the vector in the given variable
    /*!
      @param result Variable where the result should be stored
     */
    void maxValue ( Real& result ) const;

    //! Replace the vector with his absolute value
    void abs ( void );

    //! Compute the absolute value of a vector and store it in an other vector
    /*!
      @param vector Output vector to store the absolute value of the vector
     */
    void abs ( VectorEpetra& vector );

    //! Compute the scalar product of two vectors
    /*!
      @param vector Second vector for the scalar product
     */
    data_type dot ( const VectorEpetra& vector ) const;

    //! Compute the scalar product of two vectors and store the result in a given variable
    /*!
      @param vector Second vector for the scalar product
      @param scalarProduct Variable to store the result
     */
    void dot ( const VectorEpetra& vector, data_type& scalarProduct );

    //! Save the values of the matrix into a file
    /*!
      To read the file in Matlab type load filename;
      @param filename File where to save the vector
     */
    void matrixMarket ( std::string const& fileName, const bool headers = true ) const;

    //! Save the values of the matrix into a file
    /*!
      To read the file in Matlab type load filename;
      @param filename File where to save the vector
     */
    void spy ( std::string const& fileName ) const;

    //! Print the contents of the vector
    /*!
      @param output Stream where the informations must be printed
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Sets the combine mode for the import/export operations.
    /*!
      Most of the LifeV library is structured to use combine mode
      equal to Add. In some cases (cf test_filters) it is necessary
      to discard the data coming from other processors.
      @param combineMode combien mode to use for this vector from now on
     */
    void setCombineMode ( combineMode_Type combineMode );

    //! Sets the combine mode for the import/export operations to default.
    /*!
      Most of the LifeV library is structured to use combine mode
      equal to Add.
     */
    void setDefaultCombineMode();

    //! Sets the map to use for the epetra vector
    /*!
      This method can be used when building the VectorEpetra using
      empty constructor.
      @param map the map of the vector
     */
    void setMap ( const MapEpetra& map );

    //@}

    //! @name Get Methods
    //@{

    //! Return the communicator of the vector
    const Epetra_Comm& comm() const
    {
        return blockMap().Comm();
    }

    //! Return the VectorEpetra in the wrapper
    vector_type& epetraVector()
    {
        return *M_epetraVector;
    }

    //! Return the VectorEpetra in the wrapper
    const vector_type& epetraVector() const
    {
        return *M_epetraVector;
    }

    //! Return the shared pointer on the raw VectorEpetra
    const Vector_PtrType& epetraVectorPtr() const
    {
        return M_epetraVector;
    }

    //! Return the Epetra_BlockMap of the vector
    const Epetra_BlockMap& blockMap() const
    {
        return M_epetraVector->Map();
    }

    //! Return the map type of the vector (Unique or Repeated)
    MapEpetraType mapType() const
    {
        return M_mapType;
    }

    //! Return the MapEpetra of the vector
    const MapEpetra& map() const
    {
        return *M_epetraMap;
    }

    //! Return a shared pointer on the MapEpetra
    const boost::shared_ptr< MapEpetra > mapPtr() const
    {
        return M_epetraMap;
    }

    //! Return the MapEpetra of the vector
    const Epetra_Map& epetraMap() const
    {
        return * ( M_epetraMap->map ( M_mapType ) );
    }

    //! Return the size of the vector
    Int size() const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! Import the value of a vector
    /*!
      Copies the value of a vector u. If the map is not the same,
      try to import the values. Let you decide wether to add or replace shared nodes:
      note:
      <ol>
      <li> if the original source vector vector is not repeated : use Import
      <li> if the original source vector vector is repeated : use Export
      </ol>
      CombineMode Values:
      <ol>
      <li> Add - Components on the receiving processor will be added together.
      <li> Insert - Off-processor components will be inserted into locations on receiving processor replacing existing values.
      <li> InsertAdd - Off-processor components will be inserted into locations on receiving processor and added to existing values.
      <li> Average - Off-processor components will be averaged with existing components on the receiving processor.
      <li> (+ Zero and AbsMax, probably never useful)
      </ol>
      @param vector Vector to be imported
      @param combineMode Mode to be used to combine the vector
     */
    VectorEpetra& Import ( const Epetra_FEVector& vector, combineMode_Type combineMode );

    //! Export the value of a vector
    /*!
      Copies the value of this to a vector vector. If the map is not the same,
      try to import the values. Let you decide wether to add or replace shared nodes:
      note: tested only if the destination source vector vector is not repeated
      CombineMode Values:
      <ol>
      <li> Add - Components on the receiving processor will be added together.
      <li> Insert - Off-processor components will be inserted into locations on receiving processor replacing existing values.
      <li> InsertAdd - Off-processor components will be inserted into locations on receiving processor and added to existing values.
      <li> Average - Off-processor components will be averaged with existing components on the receiving processor.
      <li> (+ Zero and AbsMax, probably never useful)
      </ol>
      @param vector Vector where to store the exportation
      @param combineMode Mode to be used to combine the vector
     */
    VectorEpetra& Export ( const Epetra_FEVector& vector, combineMode_Type combineMode );

    //@}

    boost::shared_ptr< MapEpetra > M_epetraMap;
    MapEpetraType                  M_mapType;
    Vector_PtrType                 M_epetraVector;
    combineMode_Type               M_combineMode;
};

VectorEpetra operator- ( const VectorEpetra& vector );
VectorEpetra operator+ ( const VectorEpetra::data_type& scalar, const VectorEpetra& vector );
VectorEpetra operator- ( const VectorEpetra::data_type& scalar, const VectorEpetra& vector );
VectorEpetra operator* ( const VectorEpetra::data_type& scalar, const VectorEpetra& vector );

} // end namespace LifeV

#endif
