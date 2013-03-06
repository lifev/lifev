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
    @brief MatrixEpetra

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch> Cristiano Malossi <cristiano.malossi@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-11-2009
 */

#ifndef _EPETRAMATRIX_HPP_
#define _EPETRAMATRIX_HPP_

#include <lifev/core/LifeV.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_MpiComm.h>
#include <Epetra_FECrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <ml_epetra_utils.h>

#ifdef HAVE_HDF5
#include <EpetraExt_HDF5.h>
#endif

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/array/VectorEpetra.hpp>

//@@
//#define OFFSET 0

namespace LifeV
{

//! MatrixEpetra - The Epetra Matrix format Wrapper
/*!
  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

  The MatrixEpetra class provides a general interface for the Epetra_FECrsMatrix of Trilinos.

  Visit http://trilinos.sandia.gov for more informations about Epetra_FECrsMatrix.
 */
template <typename DataType>
class MatrixEpetra
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_FECrsMatrix             matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef VectorEpetra                   vector_type;
    typedef boost::shared_ptr<vector_type> vectorPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor from a graph
    /*!
      @param map Row map. The column map will be defined in MatrixEpetra<DataType>::GlobalAssemble(...,...)
      @param graph A sparse compressed row graph.
     */
    MatrixEpetra ( const MapEpetra& map, const Epetra_CrsGraph& graph, bool ignoreNonLocalValues = false );

    //! Constructor for square and rectangular matrices
    /*!
      @param map Row map. The column map will be defined in MatrixEpetra<DataType>::GlobalAssemble(...,...)
      @param numEntries The average number of entries for each row.
     */
    MatrixEpetra ( const MapEpetra& map, Int numEntries = 50 );

    //! Copy Constructor
    /*!
      @param matrix Matrix used to create the new occurence
     */
    MatrixEpetra ( const MatrixEpetra& matrix);

    //! Copies the matrix to a matrix which resides only on the processor
    /*!
      @param matrix Matrix where the content of the matrix should be stored
      @param reduceToProc Processor where the matrix should be stored
     */
    MatrixEpetra ( const MatrixEpetra& matrix, const UInt reduceToProc );

    //! Constructs an MatrixEpetra view of an Epetra_FECrsMatrix.
    /*!
      This constructor can be used when we need to modify an Epetra_FECrsMatrix
      using a method of the class MatrixEpetra
      @param crsMatrixPtr Pointer on a Epetra_FECrsMatrix of Trilinos
     */
    LIFEV_DEPRECATED ( MatrixEpetra ( matrix_ptrtype crsMatrixPtr) );

    //! Constructs an MatrixEpetra view of an Epetra_FECrsMatrix.
    /*!
      This constructor can be used when we need to modify an Epetra_FECrsMatrix
      using a method of the class MatrixEpetra
      @param map Row map. The column map will be defined in MatrixEpetra<DataType>::GlobalAssemble(...,...)
      @param crsMatrixPtr Pointer on a Epetra_FECrsMatrix of Trilinos
     */
    MatrixEpetra ( const MapEpetra& map, matrix_ptrtype crsMatrixPtr);

    //! Destructor
    ~MatrixEpetra() {};

    //@}


    //! @name Operators
    //@{

    //! Addition operator
    /*!
      @param matrix matrix to be added to the current matrix
     */
    MatrixEpetra& operator += ( const MatrixEpetra& matrix );

    //! Assignment operator
    /*!
      @param matrix matrix to be assigned to the current matrix
     */
    MatrixEpetra& operator= ( const MatrixEpetra& matrix );

    //! Matrix-Vector multiplication
    /*!
      @param vector Vector to be multiplied by the current matrix
     */
    vector_type   operator * ( const vector_type& vector ) const;

    //! Multiplication operator
    /*!
      Multiply by a scalar value the components of the current matrix.
      (modify the matrix of the class)
      @param scalar Value for the multiplication
     */
    MatrixEpetra& operator *= ( const DataType scalar );

    //! Multiplication operator
    /*!
      Multiply by a scalar value the components of the current matrix.
      (do not modify the matrix of the class)
      @param scalar Value for the multiplication
     */
    MatrixEpetra  operator *  ( const DataType scalar ) const;

    //@}


    //! @name Methods
    //@{

    //! If the matrix has been filled, this function will reopen the Matrix
    void openCrsMatrix();

    //! This function removes all the zeros in the matrix and add zero on the diagonal
    /*!
      The zeros added on the diagonal are used to impose boundary conditions
     */
    void removeZeros();

    //! Swap the given shared pointer with the one of the matrix
    /*!
      @param matrixPtr pointer on the matrix
     */
    void swapCrsMatrix ( matrix_ptrtype& matrixPtr );

    //! Swap the matrix with the one given as argument
    /*!
      @param matrix matrix which is swapped
     */
    void swapCrsMatrix ( MatrixEpetra<DataType>& matrix );

    //! Multiply the MatrixEpetra by the first given matrix and put the result in the second given matrix
    /*!
      @param transposeCurrent If true, it transposes the MatrixEpetra
      @param matrix Matrix that multiply the MatrixEpetra
      @param transposeMatrix if true, it transposes the matrix "matrix"
      @param result Matrix to store the result
      @param callFillCompleteOnResult If true, the matrix "result" will be filled (i.e. closed) after the multiplication
     */
    Int multiply ( bool transposeCurrent,
                   const MatrixEpetra<DataType>& matrix,
                   bool transposeMatrix,
                   MatrixEpetra<DataType>& result,
                   bool callFillCompleteOnResult = true ) const;

    //! Multiply the first VectorEpetra given as a parameter by the MatrixEpetra and put the result into the second given VectorEpetra
    /*!
      @param transposeCurrent If true, it transposes the MatrixEpetra
      @param vector1 Vector that will be multiply by the MatrixEpetra
      @param vector2 Vector that will store the result
     */
    Int multiply ( bool transposeCurrent, const vector_type& vector1, vector_type& vector2 ) const;

    //! Add to the matrix a dyadic product: \f[ C += A \otimes B + \f], where \f[ C \f] is the matrix.
    /*!
     *  NOTE: This method has been tested only for square matrices and unique vectors
     *  @param vector1 unique vector A
     *  @param vector2 unique vector B
     */
    void addDyadicProduct ( const vector_type& uniqueVector1, const vector_type& uniqueVector2 );

    //! Add a multiple of a given matrix:  *this += scalar*matrix
    /*!
      @param scalar Scalar which multiplies the matrix
      @param matrix Matrix to be added
     */
    void add ( const DataType scalar, const MatrixEpetra& matrix );

    //! Returns a pointer to a new matrix which contains the transpose of the current matrix
    boost::shared_ptr<MatrixEpetra<DataType> > transpose( );

    //! Set entries (rVec(i),rVec(i)) to coefficient and the rest of the row entries to zero
    /*!
      @param rVec Vector of the Id that should be set to "coefficient"
      @param coefficient Value to be set on the diagonal
      @param offset Offset used for the indices
     */
    void diagonalize ( std::vector<UInt> rVec, DataType const coefficient, UInt offset = 0 );

    //! Set entry (entryIndex,entryIndex) to coefficient and the rest of the row to zero
    /*!
      @param entryIndex Index of the row where we set the "coefficient"
      @param coefficient Value to be set on the diagonal
      @param offset Offset used for the indices
     */
    void diagonalize ( UInt const entryIndex, DataType const coefficient, UInt offset = 0 );

    //! apply constraint on all rows rVec
    /*!
      @param rVec vector of rows
      @param coefficient Value to set entry (r,r) at
      @param rhs Right hand side Vector of the system to be adapted accordingly
      @param datumVector vector of values to constrain entry r of the solution at
      @param offset Offset used for the indices
     */
    void diagonalize ( std::vector<UInt> rVec,
                       DataType const coefficient,
                       vector_type& rhs,
                       std::vector<DataType> datumVector,
                       UInt offset = 0 );

    //! apply constraint on row "row"
    /*!
      @param row row number
      @param coefficient value to set entry (row,row) at
      @param rhs Right hand side vector of the system to be adapted accordingly
      @param datumVector value to constrain entry r of the solution at
      @param offset Offset used for the indices
     */
    void diagonalize ( UInt const row, DataType const coefficient, vector_type& rhs,
                       DataType datum,
                       UInt offset = 0 );

    //! Save the matrix into a MatrixMarket (.mtx) file
    /*!
      @param filename file where the matrix will be saved
      @param headers boolean to write the MM headers or not
     */
    void matrixMarket ( std::string const& fileName, const bool headers = true );

    //! Save the matrix into a Matlab (.m) file
    /*!
      @param filename file where the matrix will be saved
     */
    void spy ( std::string const& fileName );

#ifdef HAVE_HDF5
    //! Save the matrix into a HDF5 (.h5) file
    /*!
      @param fileName Name of the file where the matrix will be saved, without extension (.h5)
      @param matrixName Name of the matrix in the HDF5 file
      @param truncate True if the file has to be truncated; False if the file already exist and should not be truncated
     */
    void exportToHDF5 ( std::string const& fileName, std::string const& matrixName = "matrix", bool const& truncate = true );

    //! Read a matrix from a HDF5 (.h5) file
    /*!
      @param fileName Name of the file where the matrix will be saved, without extension (.h5)
      @param matrixName Name of the matrix in the HDF5 file
     */
    void importFromHDF5 ( std::string const& fileName, std::string const& matrixName = "matrix" );
#endif

    //! Print the contents of the matrix
    /*!
      @param output Stream where the informations must be printed
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //! Global assemble of a square matrix with default domain and range map
    /*
      <ol>
      <li> Calls insertZeroDiagonal and then Epetra_FECsrMatrix::GlobalAssemble();
      <li> Set M_domainMap and M_rangeMap
      </ol>

      EpetraFECsrMatrix will assume that both the domain and range map are the same
      of the row map defined in the constructor.

      NOTE: domain and range map must be one-to-one and onto. (Unique map)
     */
    Int globalAssemble();

    //! Global assemble for rectangular matrices.
    /*
     <ol>
     <li> Calls Epetra_FECsrMatrix::GlobalAssemble(Epetra_Map & domainMap, Epetra_Map & rangeMap);
     <li> Set M_domainMap and M_rangeMap
     </ol>
     @param domainMap The domain map
     @param rangeMap The range map
     */
    Int globalAssemble ( const boost::shared_ptr<const MapEpetra>& domainMap,
                         const boost::shared_ptr<const MapEpetra>& rangeMap );

    //! insert the given value into the diagonal
    /*!
      Pay intention that this will add values to the diagonal,
      so for later added values with set_mat_inc, the one
      will be added.
      Inserts Value on the local diagonal for diagonal elements specified by the input MapEpetra;
      This methods works only if matrix is not closed.

      @param entry The entry that is inserted in the diagonal
      @param Map The MapEpetra
      @param offset An offset for the insertion of the diagonal entries
    */
    void insertValueDiagonal ( const DataType entry, const MapEpetra& Map, const UInt offset = 0 );

    //! insert the given value into the diagonal
    /*!
      Pay intention that this will add values to the diagonal,
      so for later added values with set_mat_inc, the one
      will be added

      Inserts Value on the local diagonal for diagonal elements >= from and < to;
      <ol>
      <li> If from > to, process all diagonal entries entries
      <li> If from = to, do nothing
      </ol>
      This methods works only if matrix is not closed.
      @param value Value to be inserted on the diagonal
      @param from Starting row
      @param to Ending row
    */
    void insertValueDiagonal ( const DataType& value, Int from = -1, Int to = -2 );

    //! insert ones into the diagonal to ensure the matrix' graph has a entry there
    /*
      Pay intention that this will add ones to the diagonal,
      so for later added values with set_mat_inc, the one
      will be added

      Inserts Value on the local diagonal for diagonal elements >= from and < to;
      <ol>
      <li> If from > to, process all diagonal entries entries
      <li> If from = to, do nothing
      </ol>
      This methods works only if matrix is not closed.
      @param from Starting row
      @param to Ending row
    */
    void insertOneDiagonal ( Int from = -1, Int to = -2 );

    //! insert zeros into the diagonal to ensure the matrix' graph has a entry there
    /*!
      This method does not remove non zero entries in the diagonal.

      Inserts Value on the local diagonal for diagonal elements >= from and < to;
      <ol>
      <li> If from > to, process all diagonal entries entries
      <li> If from = to, do nothing
      </ol>
      This methods works only if matrix is not closed.
      @param from Starting row
      @param to Ending row
    */
    void insertZeroDiagonal ( Int from = -1, Int to = -2 );

    //! Compute the norm 1 of the global matrix
    /*!
      @return norm 1
     */
    Real norm1() const;

    //! Compute the norm inf of the global matrix
    /*!
      @return norm inf
     */
    Real normInf() const;

    //@}


    //! @name Set Methods
    //@{

    //! set zero in all the matrix entries
    void zero()
    {
        M_epetraCrs->PutScalar (0.);
    }

    //! Set a set of values to the corresponding set of coefficient in the matrix
    /*!
      @param numRows Number of rows into the list given in "localValues"
      @param numColumns Number of columns into the list given in "localValues"
      @param rowIndices List of row indices
      @param columnIndices List of column indices
      @param localValues 2D array containing the coefficient related to "rowIndices" and "columnIndices"
      @param format Format of the matrix (Epetra_FECrsMatrix::COLUMN_MAJOR or Epetra_FECrsMatrix::ROW_MAJOR)
     */
    void setCoefficients ( Int const numRows, Int const numColumns,
                           std::vector<Int> const& rowIndices, std::vector<Int> const& columnIndices,
                           DataType* const* const localValues,
                           Int format = Epetra_FECrsMatrix::COLUMN_MAJOR );

    //! Set a coefficient of the matrix
    /*!
      @param row Row index of the coefficient
      @param column column index of the coefficient
     */
    void setCoefficient ( UInt row, UInt column, DataType localValue );

    //! Add a set of values to the corresponding set of coefficient in the matrix
    /*!
      @param numRows Number of rows into the list given in "localValues"
      @param numColumns Number of columns into the list given in "localValues"
      @param rowIndices List of row indices
      @param columnIndices List of column indices
      @param localValues 2D array containing the coefficient related to "rowIndices" and "columnIndices"
      @param format Format of the matrix (Epetra_FECrsMatrix::COLUMN_MAJOR or Epetra_FECrsMatrix::ROW_MAJOR)
     */
    void addToCoefficients ( Int const numRows, Int const numColumns,
                             std::vector<Int> const& rowIndices, std::vector<Int> const& columnIndices,
                             DataType* const* const localValues,
                             Int format = Epetra_FECrsMatrix::COLUMN_MAJOR );

    //! Add a value at a coefficient of the matrix
    /*!
      @param row Row index of the value to be added
      @param column Column index of the value to be added
      @param localValue Value to be added to the coefficient
     */
    void addToCoefficient ( UInt const row, UInt const column, DataType const localValue );

    //! If set true, transpose of this operator will be applied.
    /*!
      @param transpose flag to identify whether to transpose or not
     */
    void setUseTranspose ( const bool& transpose = false )
    {
        M_epetraCrs->SetUseTranspose ( transpose );
    }

    //@}


    //! @name Get Methods
    //@{

    //! Return the shared_pointer of the Epetra_FECrsMatrix
    matrix_ptrtype& matrixPtr()
    {
        return M_epetraCrs;
    }

    //! Return the shared_pointer of the Epetra_Map
    boost::shared_ptr<MapEpetra> mapPtr()
    {
        return M_map;
    }

    //! Return the const shared_pointer of the Epetra_FECrsMatrix
    const matrix_ptrtype& matrixPtr() const
    {
        return M_epetraCrs;
    }

    //! Return the mean number of entries in the matrix rows
    Int meanNumEntries() const ;

    //! Return the Id of the processor
    Int processorId();

    //! Return the row MapEpetra of the MatrixEpetra used in the assembling
    /*!
      Note: This method should be call when MapEpetra is still open.
     */
    const MapEpetra& map() const;

    //! Return the domain MapEpetra of the MatrixEpetra
    /*!
      This function should be called only after MatrixEpetra<DataType>::GlobalAssemble(...) has been called.
      If this is an open matrix that M_domainMap is an invalid pointer
     */
    const MapEpetra& domainMap() const;

    //! Return the range MapEpetra of the MatrixEpetra
    /*!
      This function should be called only after MatrixEpetra<DataType>::GlobalAssemble(...) has been called.
      If this is an open matrix that M_domainMap is an invalid pointer
     */
    const MapEpetra& rangeMap() const;

    //@}

    //! Friend Functions
    //@{
    //! RAP matrix matrix multiplication result = R * A * P
    //! User is responsible to wrap the row pointer returned by this method with his favorite pointer
    template <typename DType>
    friend MatrixEpetra<DType>* RAP (const MatrixEpetra<DType>& R, const MatrixEpetra<DType>& A, const MatrixEpetra<DType>& P);

    //! PtAP matrix matrix multiplication result = Pt * A * P
    //! User is responsible to wrap the row pointer returned by this method with his favorite pointer
    template <typename DType>
    friend MatrixEpetra<DType>* PtAP (const MatrixEpetra<DType>& A, const MatrixEpetra<DType>& P);
    //@}
private:


    // Shared pointer on the row MapEpetra used in the assembling
    boost::shared_ptr< MapEpetra > M_map;

    // Shared pointer on the domain MapEpetra.
    /*
      if y = this*x,
      then x.getMap() is the domain map.
      NOTE: Epetra assume the domain map to be 1-1 and onto (Unique)
      M_domainMap is a NULL pointer until MatrixEpetra<DataType> is called.
     */
    boost::shared_ptr< const MapEpetra > M_domainMap;

    //! Shared pointer on the range MapEpetra.
    /*
      if y = this*x,
      then y.getMap() is the range map.
      NOTE: Epetra assume the domain map to be 1-1 and onto (Unique)
      M_rangeMap is a NULL pointer until MatrixEpetra<DataType> is called.
     */
    boost::shared_ptr< const MapEpetra > M_rangeMap;

    // Pointer on a Epetra_FECrsMatrix
    matrix_ptrtype  M_epetraCrs;
};


// ===================================================
// Constructors & Destructor
// ===================================================
template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra ( const MapEpetra& map, const Epetra_CrsGraph& graph, bool ignoreNonLocalValues ) :
    M_map       ( new MapEpetra ( map ) ),
    M_epetraCrs ( new matrix_type ( Copy, graph, ignoreNonLocalValues ) )
{

}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra ( const MapEpetra& map, Int numEntries ) :
    M_map       ( new MapEpetra ( map ) ),
    M_epetraCrs ( new matrix_type ( Copy, *M_map->map ( Unique ), numEntries, false) )
{

}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra ( const MatrixEpetra& matrix ) :
    M_map      ( matrix.M_map ),
    M_domainMap ( matrix.M_domainMap ),
    M_rangeMap ( matrix.M_rangeMap ),
    M_epetraCrs ( new matrix_type ( *matrix.M_epetraCrs ) )
{

}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra ( const MatrixEpetra& matrix, const UInt reduceToProc ) :
    M_map      ( matrix.getMap().createRootMap ( reduceToProc ) ),
    M_epetraCrs ( new matrix_type ( Copy, *M_map->map ( Unique ),
                                    numEntries ( matrix.M_epetraCrs->Map().Comm().MyPID() == reduceToProc ) * 20,
                                    false ) )
{
    Epetra_Export reducedExport ( M_epetraCrs->Map(), matrix.M_epetraCrs->Map() );
    M_epetraCrs->Import ( *matrix.M_epetraCrs, reducedExport, Add );

    if (M_epetraCrs->Filled() )
    {
        M_domainMap = matrix.getDomainMap().createRootMap (reduceToProc);
        M_rangeMap  = matrix.getRangeMap().createRootMap (reduceToProc);
    }
}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra ( matrix_ptrtype CRSMatrixPtr ) :
    M_map( ),
    M_domainMap(),
    M_rangeMap()
{
    M_epetraCrs = CRSMatrixPtr;
}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra (const MapEpetra& map, matrix_ptrtype CRSMatrixPtr ) :
    M_map ( new MapEpetra (map) ),
    M_domainMap(),
    M_rangeMap()
{
    M_epetraCrs = CRSMatrixPtr;
}



// ===================================================
// Operators
// ===================================================
template <typename DataType>
MatrixEpetra<DataType>&
MatrixEpetra<DataType>::operator += ( const MatrixEpetra& matrix )
{
    EpetraExt::MatrixMatrix::Add ( *matrix.matrixPtr(), false, 1., *this->matrixPtr(), 1. );

    return *this;
}

template<typename DataType>
MatrixEpetra<DataType>&  MatrixEpetra<DataType>::operator= ( const MatrixEpetra& matrix )
{
    M_map        = matrix.M_map;
    M_domainMap  = matrix.M_domainMap;
    M_rangeMap   = matrix.M_rangeMap;
    *M_epetraCrs = * ( matrix.M_epetraCrs );

    return *this;
}

template<typename DataType>
typename MatrixEpetra<DataType>::vector_type
MatrixEpetra<DataType>::operator * ( const vector_type& vector ) const
{
    ASSERT_PRE ( M_epetraCrs->Filled(),
                 "MatrixEpetra::Operator*: globalAssemble(...) should be called first" );
    ASSERT_PRE ( vector.map().mapsAreSimilar (*M_domainMap),
                 "MatrixEpetra::Operator*: the map of vec is not the same of domainMap" );
    ASSERT_PRE ( M_rangeMap.get(),
                 "MatrixEpetra::Operator*: the rangeMap is not set" );

    vector_type result (*M_rangeMap);
    M_epetraCrs->Apply ( vector.epetraVector(), result.epetraVector() );

    return result;
}


template<typename DataType>
MatrixEpetra<DataType>&
MatrixEpetra<DataType>::operator *= ( const DataType value )
{
    M_epetraCrs->Scale ( value );
    return *this;
}

template<typename DataType>
MatrixEpetra<DataType>
MatrixEpetra<DataType>::operator * ( const DataType scalar ) const
{
    MatrixEpetra<DataType> matrix ( *this );
    return matrix *= scalar;
}

// ===================================================
// Methods
// ===================================================
template <typename DataType>
void MatrixEpetra<DataType>::openCrsMatrix()
{
    if ( M_epetraCrs->Filled() )
    {
        Int meanNumEntries = this->meanNumEntries();
        matrix_ptrtype tmp ( M_epetraCrs );
        M_epetraCrs.reset ( new matrix_type ( Copy, M_epetraCrs->RowMap(), meanNumEntries ) );

        EpetraExt::MatrixMatrix::Add ( *tmp, false, 1., *M_epetraCrs, 1. );

        M_domainMap.reset();
        M_rangeMap.reset();

    }
}

template <typename DataType>
void MatrixEpetra<DataType>::removeZeros()
{
    if ( M_epetraCrs->Filled() )
    {
        Int meanNumEntries = 1/*this->getMeanNumEntries()*/;
        matrix_ptrtype tmp( M_epetraCrs );
        M_epetraCrs.reset(new matrix_type( Copy, M_epetraCrs->RowMap(), meanNumEntries ) );

        //Variables to store the informations
        Int NumEntries;
        Real* Values;
        Int* Indices;
        Int row (0);

        for ( Int i (0); i < tmp->NumGlobalRows(); ++i )
        {
            row = tmp->LRID( i );
            // Check if the row belong to this process
            if(row==-1)
                continue;
            tmp->ExtractMyRowView( row, NumEntries, Values, Indices );

            std::vector<Int> Indices2 ( NumEntries );
            std::vector<Real> Values2 ( NumEntries );
            Int NumEntries2 (0);

            for (Int j (0); j < NumEntries; ++j)
            {
                if (Values[j] != 0.0)
                {
                    Indices2[NumEntries2] = tmp->GCID (Indices[j]);
                    Values2[NumEntries2]  = Values[j];
                    NumEntries2++;
                }
            }
            M_epetraCrs->InsertGlobalValues( i, NumEntries2, &Values2[0], &Indices2[0] );
        }
        insertZeroDiagonal();
        M_epetraCrs->GlobalAssemble();
    }
}

template <typename DataType>
void MatrixEpetra<DataType>::swapCrsMatrix ( matrix_ptrtype& matrixPtr )
{
    M_epetraCrs.swap ( matrixPtr );
}


template <typename DataType>
void MatrixEpetra<DataType>::swapCrsMatrix ( MatrixEpetra<DataType>& matrix )
{
    M_epetraCrs.swap ( matrix.M_epetraCrs );
}

template <typename DataType>
Int MatrixEpetra<DataType>::multiply ( bool transposeCurrent,
                                       const MatrixEpetra<DataType>& matrix, bool transposeMatrix,
                                       MatrixEpetra<DataType>& result, bool callFillCompleteOnResult ) const
{
    Int errCode = EpetraExt::MatrixMatrix::Multiply ( *M_epetraCrs, transposeCurrent,
                                                      *matrix.matrixPtr(), transposeMatrix,
                                                      *result.matrixPtr(), false );
    if (callFillCompleteOnResult)
    {
        boost::shared_ptr<const MapEpetra> domainMap, rangeMap;
        if (transposeCurrent)
        {
            rangeMap = M_domainMap;
        }
        else
        {
            rangeMap = M_rangeMap;
        }

        if (transposeMatrix)
        {
            domainMap = matrix.M_rangeMap;
        }
        else
        {
            domainMap  = matrix.M_domainMap;
        }

        result.globalAssemble (domainMap, rangeMap);
    }

    return errCode;
}

template <typename DataType>
Int MatrixEpetra<DataType>::multiply ( bool transposeCurrent, const vector_type& vector1, vector_type& vector2 ) const
{
    ASSERT_PRE ( M_epetraCrs->Filled(),
                 "MatrixEpetra<DataType>::Multiply: GlobalAssemble(...) must be called first" );
    ASSERT_PRE ( transposeCurrent ? vector1.map().mapsAreSimilar (*M_rangeMap) : vector1.map().mapsAreSimilar (*M_domainMap),
                 "MatrixEpetra<DataType>::Multiply: x map is different from M_domainMap (or M_rangeMap if transposeCurrent == true)" );
    ASSERT_PRE ( transposeCurrent ? vector2.map().mapsAreSimilar (*M_domainMap) : vector2.map().mapsAreSimilar (*M_rangeMap),
                 "MatrixEpetra<DataType>::Multiply: y map is different from M_rangeMap (or M_domainMap if transposeCurrent == true)" );


    return M_epetraCrs->Multiply ( transposeCurrent, vector1.epetraVector(), vector2.epetraVector() );
}

template <typename DataType>
void MatrixEpetra<DataType>::addDyadicProduct ( const vector_type& uniqueVector1, const vector_type& uniqueVector2 )
{
    // Check if the matrix is open
    if ( M_epetraCrs->Filled() )
    {
        if ( M_epetraCrs->Comm().MyPID() == 0 )
        {
            std::cout << "ERROR: Matrix already filled: cannot add dyadic product!" << std::endl;
        }
        return;
    }

    // Build a repeated list of globalElements
    std::vector<Int> myGlobalElements ( uniqueVector2.size() );
    for ( UInt i = 0 ; i < myGlobalElements.size() ; ++i )
    {
        myGlobalElements[i] = i;
    }

    // Build a repeated map
    MapEpetra repeatedMap ( -1, static_cast< Int > ( myGlobalElements.size() ), &myGlobalElements[0], uniqueVector2.map().commPtr() );

    // Create a fully repeated copy of uniqueVector2
    VectorEpetra repeatedVector2 ( repeatedMap, Repeated );
    repeatedVector2 = uniqueVector2;

    // Fill the matrix with the result of the dyadic Product
    for ( Int row (0); row < M_epetraCrs->NumMyRows(); ++row )
        for ( Int column (0); column < M_epetraCrs->NumGlobalCols(); ++column )
        {
            addToCoefficient ( M_epetraCrs->GRID ( row ), column, uniqueVector1[M_epetraCrs->GRID ( row )] * repeatedVector2[column] );
        }
}

template <typename DataType>
void MatrixEpetra<DataType>::add ( const DataType scalar, const MatrixEpetra& matrix )
{
    EpetraExt::MatrixMatrix::Add ( *matrix.matrixPtr(), false, scalar, *this->matrixPtr(), 1. );
}

template <typename DataType>
boost::shared_ptr<MatrixEpetra<DataType> > MatrixEpetra<DataType>::transpose( )
{
    ASSERT_PRE (M_epetraCrs->Filled(), "The transpose can be formed only if the matrix is already filled!");
    boost::shared_ptr<Epetra_FECrsMatrix> transposedFE;
    transposedFE.reset (new Epetra_FECrsMatrix (Copy, M_epetraCrs->OperatorDomainMap(), M_epetraCrs->OperatorRangeMap(), 0, false) );
    EpetraExt::RowMatrix_Transpose transposer;
    *dynamic_cast<Epetra_CrsMatrix*> (& (*transposedFE) ) = dynamic_cast<Epetra_CrsMatrix&> (transposer (*M_epetraCrs) );
    transposedFE->FillComplete();
    boost::shared_ptr<MatrixEpetra<DataType> > transposedMatrix (new MatrixEpetra<DataType> (*M_domainMap) );
    transposedMatrix->globalAssemble (M_rangeMap, M_domainMap);
    transposedMatrix->matrixPtr() = transposedFE;

    return transposedMatrix;
}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize ( std::vector<UInt> rVec, DataType const coefficient, UInt offset )
{

    const Epetra_Comm&  Comm ( M_epetraCrs->Comm() );
    Int numProcs ( Comm.NumProc() );
    Int MyPID   ( Comm.MyPID() );
    Int i;


    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    Int*  r;
    UInt* Ur;


    // loop on all proc
    for ( Int p (0); p < numProcs; p++ )
    {
        Int sizeVec ( rVec.size() );

        Comm.Broadcast ( &sizeVec, 1, p );

        if ( p == MyPID )
        {
            Ur = &rVec.front();
        }
        else
        {
            Ur = new UInt[sizeVec];
        }

        r = (Int*) Ur;

        Comm.Broadcast ( r, sizeVec, p );

        for ( i = 0; i < sizeVec; i++ )
        {
            diagonalize ( Ur[i], coefficient, offset );
        }

        if ( p != MyPID )
        {
            delete[] Ur;
        }

    }

}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize ( UInt const row,
                                           DataType const coefficient,
                                           UInt offset )
{

    if ( !M_epetraCrs->Filled() )
    {
        // if not filled, I do not know how to diagonalize.
        ERROR_MSG ( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap ( M_epetraCrs->RowMap() );
    const Epetra_Map& colMap ( M_epetraCrs->ColMap() );


    Int myCol = colMap.LID ( row + offset );

    // row: if r is mine, zero out values
    Int myRow = rowMap.LID ( row + offset );

    if ( myRow >= 0 )  // I have this row
    {
        Int    NumEntries;
        Real* Values;
        Int* Indices;

        M_epetraCrs->ExtractMyRowView ( myRow, NumEntries, Values, Indices );

        for (Int i (0); i <  NumEntries; i++)
        {
            Values[i] = 0;
        }

        DataType coeff ( coefficient );
        M_epetraCrs->ReplaceMyValues ( myRow, 1, &coeff, &myCol ); // A(r,r) = coefficient
    }

}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize ( std::vector<UInt> rVec,
                                           DataType const coefficient,
                                           vector_type& rhs,
                                           std::vector<DataType> datumVec,
                                           UInt offset )
{


    const Epetra_Comm&  Comm ( M_epetraCrs->Comm() );
    Int numProcs ( Comm.NumProc() );
    Int MyPID   ( Comm.MyPID() );
    Int i;

    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    Int*     r;
    UInt*    Ur;
    DataType* datum;

#if 1
    // 2 arrays to store le local and remote IDs

    Int me = Comm.MyPID();

    std::vector<Int> localIDs;
    std::vector<Int> remoteIDs;

    std::vector<Real> localData;
    std::vector<Real> remoteData;

    std::map<Int, Real>           localBC;
    std::map<Int, Real>::iterator im;

    const Epetra_Map& rowMap ( M_epetraCrs->RowMap() );


    //Comm.Barrier();
    // we want to know which IDs are our or not

    for ( Int ii = 0; ii < (Int) rVec.size(); ++ii )
    {
        Int lID = rowMap.LID (rVec[ii]);
        if ( ! ( lID < 0 ) )
        {

            localIDs.push_back ( rVec[ii] );
            localData.push_back ( datumVec[ii] );
            localBC.insert ( std::pair<Int, Real> ( rVec[ii], datumVec[ii] ) );
        }
        else
        {
            remoteIDs.push_back ( rVec[ii] );
            remoteData.push_back ( datumVec[ii] );
        }
    }

    // now, we have to fill our localIDs with IDs from other processors
    // first, we have to build the map of all the remoteIDs and their processor owner


    Int numIDs = remoteIDs.size();

    Int* PIDList = new Int[numIDs];
    Int* LIDList = new Int[numIDs];

    rowMap.RemoteIDList ( numIDs,
                          &remoteIDs[0],
                          PIDList,
                          LIDList );

    std::vector< std::vector<Int> > procToID  ( Comm.NumProc() );
    std::vector< std::vector<Real> > procToData ( Comm.NumProc() );

    for ( Int ii = 0; ii < numIDs; ++ii )
    {
        Int pi = PIDList[ii];
        procToID[pi].push_back ( remoteIDs[ii] );
        procToData[pi].push_back ( remoteData[ii] );
    }

    // then, we send all the nodes where they belong

    const Epetra_MpiComm* comm = dynamic_cast<Epetra_MpiComm const*> ( &Comm );

    assert ( comm != 0 );

    for ( Int ii = 0; ii < (Int) procToID.size(); ++ii )
    {
        if ( ii != me )
        {
            Int length;
            length = procToID[ii].size();
            MPI_Send ( &length, 1, MPI_INT, ii, 666, comm->Comm() );
            if ( length > 0 )
            {
                MPI_Send ( &procToID[ii][0], length, MPI_INT, ii, 667, comm->Comm() );
                MPI_Send ( &procToData[ii][0], length, MPI_DOUBLE, ii, 668, comm->Comm() );
            }
        }

    }

    for ( Int ii = 0; ii < (Int) procToID.size(); ++ii )
    {
        if ( ii != me )
        {
            Int length;
            MPI_Status status;
            MPI_Recv ( &length, 1, MPI_INT, ii, 666, comm->Comm(), &status );


            if ( length > 0 )
            {
                Int* bufferID = new Int[length];
                Int* ptrID (0); //    = new Int[length];

                MPI_Recv ( bufferID, length, MPI_INT, ii, 667, comm->Comm(), &status );

                ptrID = bufferID;

                Real* bufferData = new Real[length];
                Real* ptrData (0);

                MPI_Recv ( bufferData, length, MPI_DOUBLE, ii, 668, comm->Comm(), &status );
                ptrData = bufferData;

                for ( Int ii = 0; ii < length; ++ii, ++ptrID, ++ptrData )
                {
                    localBC.insert ( std::pair<Int, Real>
                                     ( *ptrID, *ptrData ) );
                }

                delete[] bufferID;
                delete[] bufferData;

            }

        }
    }

    delete[] PIDList;
    delete[] LIDList;

    for ( im = localBC.begin(); im != localBC.end(); ++im )
    {
        diagonalize ( im->first, coefficient, rhs, im->second, offset );
    }

    return;
#endif


    // loop on all proc
    for ( Int p (0); p < numProcs; p++ )
    {
        Int sizeVec (rVec.size() );
        if ( sizeVec != Int (datumVec.size() ) )
        {
            // vectors must be of the same size
            ERROR_MSG ( "diagonalize: vectors must be of the same size\n" );
        }

        Comm.Broadcast ( &sizeVec, 1, p );

        if ( p == MyPID )
        {
            Ur    =  &rVec.front();
            datum = &datumVec.front();
        }
        else
        {
            Ur    = new UInt    [sizeVec];
            datum = new DataType[sizeVec];
        }

        r = (Int*) Ur;

        Comm.Broadcast ( r,    sizeVec, p );
        Comm.Broadcast ( datum, sizeVec, p );

        for ( i = 0; i < sizeVec; i++ )
        {
            diagonalize ( Ur[i], coefficient, rhs, datum[i], offset );
        }

        if ( p != MyPID )
        {
            delete[] Ur;
            delete[] datum;
        }

    }

}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize ( UInt const row,
                                           DataType const coefficient,
                                           vector_type& rhs,
                                           DataType datum,
                                           UInt offset )
{

    if ( !M_epetraCrs->Filled() )
    {
        // if not filled, I do not know how to diagonalize.
        ERROR_MSG ( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap ( M_epetraCrs->RowMap() );
    const Epetra_Map& colMap ( M_epetraCrs->ColMap() );


    Int myCol = colMap.LID ( row + offset );

#ifdef EPETRAMATRIX_SYMMETRIC_DIAGONALIZE
    if ( myCol >= 0 )  // I have this column
    {
        Real zero (0);
        for ( Int i (0); i < rowMap.NumMyElements(); i++ )
            // Note that if a value is not already present for the specified location in the matrix,
            // the input value will be ignored and a positive warning code will be returned.
        {
            M_epetraCrs->ReplaceMyValues (i, 1, &zero, &myCol);
        }
    }
#endif

    // row: if r is mine, zero out values
    Int myRow = rowMap.LID ( row + offset );

    if ( myRow >= 0 )  // I have this row
    {
        Int    NumEntries;
        Real* Values;
        Int* Indices;

        M_epetraCrs->ExtractMyRowView ( myRow, NumEntries, Values, Indices );

        for ( Int i (0); i <  NumEntries; i++ )
        {
            Values[i] = 0;
        }

        DataType coeff ( coefficient );

        M_epetraCrs->ReplaceMyValues ( myRow, 1, &coeff, &myCol ); // A(r,r) = coeff
        rhs[row + offset] = coefficient * datum; // correct right hand side for row r

    }

}

template <typename DataType>
void MatrixEpetra<DataType>::matrixMarket ( std::string const& fileName, const bool headers )
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName;
    std::string desc = "Created by LifeV";

    name = fileName + ".mtx";

    EpetraExt::RowMatrixToMatrixMarketFile ( name.c_str(),
                                             *M_epetraCrs,
                                             name.c_str(),
                                             desc.c_str(),
                                             headers
                                           );
}

template <typename DataType>
void MatrixEpetra<DataType>::spy ( std::string const& fileName )
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName, uti = " , ";

    Int  me = M_epetraCrs->Comm().MyPID();
    std::ostringstream myStream;
    myStream << me;
    name = fileName + ".m";

    EpetraExt::RowMatrixToMatlabFile ( name.c_str(), *M_epetraCrs );

}

#ifdef HAVE_HDF5
template <typename DataType>
void MatrixEpetra<DataType>::exportToHDF5 ( std::string const& fileName, std::string const& matrixName, bool const& truncate )
{
    EpetraExt::HDF5 HDF5 ( M_epetraCrs->Comm() );

    if ( truncate )
    {
        // Create and open the file / Truncate and open the file
        HDF5.Create ( ( fileName + ".h5" ).data() );
    }
    else
    {
        // Open an existing file without truncating it
        HDF5.Open ( ( fileName + ".h5" ).data() );
    }

    // Check if the file is created
    if ( !HDF5.IsOpen () )
    {
        std::cerr << "Unable to create " + fileName + ".h5";
        abort();
    }

    // Save the matrix into the file
    HDF5.Write ( matrixName.data(), *M_epetraCrs );

    // Close the file
    HDF5.Close();

} // exportToHDF5

template <typename DataType>
void MatrixEpetra<DataType>::importFromHDF5 ( std::string const& fileName, std::string const& matrixName )
{
    EpetraExt::HDF5 HDF5 ( M_epetraCrs->Comm() );

    // Open an existing file
    HDF5.Open ( ( fileName + ".h5" ).data() );

    // Check if the file is created
    if ( !HDF5.IsOpen () )
    {
        std::cerr << "Unable to open " + fileName + ".h5";
        abort();
    }

    // Read the matrix from the file
    Epetra_CrsMatrix* importedMatrix (0);
    HDF5.Read ( matrixName.data(), M_epetraCrs->DomainMap(), M_epetraCrs->RangeMap(), importedMatrix );

    // Copy the loaded matrix to the member object
    M_epetraCrs.reset ( new matrix_type ( *dynamic_cast< Epetra_FECrsMatrix* > ( importedMatrix ) ) );

    // Close the file
    HDF5.Close();

} // importFromHDF5
#endif

template <typename DataType>
void MatrixEpetra<DataType>::showMe ( std::ostream& output ) const
{
    output << "showMe must be implemented for the MatrixEpetra class" << std::endl;
}

template <typename DataType>
Int MatrixEpetra<DataType>::globalAssemble()
{
    if ( !M_epetraCrs->Filled() )
    {
        insertZeroDiagonal();
    }
    M_domainMap = M_map;
    M_rangeMap  = M_map;
    return  M_epetraCrs->GlobalAssemble();
}

template <typename DataType>
Int MatrixEpetra<DataType>::globalAssemble ( const boost::shared_ptr<const MapEpetra>& domainMap,
                                             const boost::shared_ptr<const MapEpetra>& rangeMap )
{

    if ( !M_epetraCrs->Filled() && domainMap->mapsAreSimilar ( *rangeMap) )
    {
        insertZeroDiagonal();
    }


    M_domainMap = domainMap;
    M_rangeMap  = rangeMap;
    return  M_epetraCrs->GlobalAssemble ( *domainMap->map (Unique), *rangeMap->map (Unique) );
}

template <typename DataType>
void
MatrixEpetra<DataType>::insertValueDiagonal ( const DataType entry, const MapEpetra& Map, const UInt offset )
{
    for ( Int i = 0 ; i < Map.map (Unique)->NumMyElements(); ++i )
    {
        addToCoefficient ( offset + Map.map (Unique)->GID (i) , offset + Map.map (Unique)->GID (i), entry );
    }
}

template <typename DataType>
void MatrixEpetra<DataType>::insertValueDiagonal ( const DataType& value, Int from, Int to )
{
    if ( M_epetraCrs->Filled() )
    {
        if ( M_epetraCrs->Comm().MyPID() == 0 )
        {
            std::cout << "Matrix is already filled, it is impossible to insert the diagonal now" << std::endl;
        }
        return;
    }

    if ( to == from )
    {
        return;
    }

    if ( to < from ) // do all entries
    {
        from = M_epetraCrs->RowMap().MinMyGID ();
        to = M_epetraCrs->RowMap().MaxMyGID () + 1;
    }

    Int* p =  M_epetraCrs->RowMap().MyGlobalElements();
    Int ierr;

    for ( Int i (0); i <  M_epetraCrs->RowMap().NumMyElements(); ++i, ++p )
    {
        if ( *p < from || *p >= to )
        {
            continue;
        }

        ierr = M_epetraCrs->InsertGlobalValues ( 1, p, 1, p, &value );

        if ( ierr < 0 )
        {
            std::cerr << " error in matrix insertion " << ierr << std::endl;
        }
    }
}

template <typename DataType>
void MatrixEpetra<DataType>::insertOneDiagonal ( Int from, Int to )
{
    insertValueDiagonal ( 1.0, from, to );
}

template <typename DataType>
void MatrixEpetra<DataType>::insertZeroDiagonal ( Int from, Int to )
{
    insertValueDiagonal ( 0.0, from, to );
}

template <typename DataType>
Real MatrixEpetra<DataType>::norm1() const
{
    return M_epetraCrs->NormOne();
}

template <typename DataType>
Real MatrixEpetra<DataType>::normInf() const
{
    return M_epetraCrs->NormInf();
}


// ===================================================
// Set Methods
// ===================================================
template <typename DataType>
void MatrixEpetra<DataType>::
setCoefficients ( Int const numRows, Int const numColumns,
                  std::vector<Int> const& rowIndices, std::vector<Int> const& columnIndices,
                  DataType* const* const localValues,
                  Int format )
{
    Int ierr;
    ierr = M_epetraCrs->ReplaceGlobalValues ( numRows, &rowIndices[0], numColumns, &columnIndices[0], localValues, format );

    if ( ierr < 0 ) std::cout << " error in matrix insertion [setCoefficients] " << ierr
                                  << " when inserting in (" << rowIndices[0] << ", " << columnIndices[0] << ")" << std::endl;
}



template <typename DataType>
void MatrixEpetra<DataType>::
setCoefficient ( UInt row, UInt column, DataType localValue )
{
    Int irow (    row );
    Int icol ( column );

    Int ierr = M_epetraCrs->ReplaceGlobalValues ( 1, &irow, 1, &icol, &localValue );
    if (ierr != 0)
    {
        std::cerr << " error in matrix replacement " << ierr << std::endl;
    }

}

template <typename DataType>
void MatrixEpetra<DataType>::
addToCoefficient ( UInt const row, UInt const column, DataType const localValue )
{

    Int irow = static_cast<Int> (row);
    Int icol = static_cast<Int> (column);

    Int ierr = ( M_epetraCrs->Filled() )   ?
               M_epetraCrs->SumIntoGlobalValues ( 1, &irow, 1, &icol, &localValue )
               :
               M_epetraCrs->InsertGlobalValues ( 1, &irow, 1, &icol, &localValue );

    std::stringstream errorMessage;
    errorMessage << " error in matrix insertion [addToCoefficient] " << ierr
                 << " when inserting " << localValue << " in (" << irow << ", " << icol << ")" << std::endl;
    ASSERT ( ierr >= 0, errorMessage.str() );

}

template <typename DataType>
void MatrixEpetra<DataType>::
addToCoefficients ( Int const numRows, Int const numColumns,
                    std::vector<Int> const& rowIndices, std::vector<Int> const& columnIndices,
                    DataType* const* const localValues,
                    Int format )
{
    Int ierr = ( M_epetraCrs->Filled() ) ?
               M_epetraCrs->SumIntoGlobalValues ( numRows, &rowIndices[0], numColumns,
                                                  &columnIndices[0], localValues, format )
               :
               M_epetraCrs->InsertGlobalValues ( numRows, &rowIndices[0], numColumns,
                                                 &columnIndices[0], localValues, format );

    std::stringstream errorMessage;
    errorMessage << " error in matrix insertion [addToCoefficients] " << ierr
                 << " when inserting in (" << rowIndices[0] << ", " << columnIndices[0] << ")" << std::endl;
    ASSERT ( ierr >= 0, errorMessage.str() );

}

// ===================================================
// Get Methods
// ===================================================
template <typename DataType>
Int MatrixEpetra<DataType>::meanNumEntries() const
{
    const Int minEntries = M_epetraCrs->MaxNumEntries() / 2;
    if ( !M_epetraCrs->NumMyRows() )
    {
        return minEntries;
    }

    Int meanNumEntries = M_epetraCrs->NumMyNonzeros() / M_epetraCrs->NumMyRows();
    if ( meanNumEntries < minEntries || meanNumEntries > 2 * minEntries )
    {
        return minEntries;
    }
    return meanNumEntries;
}

template <typename DataType>
Int MatrixEpetra<DataType>::processorId()
{
    return  M_epetraCrs->Comm().MyPID();
}

template <typename DataType>
const MapEpetra& MatrixEpetra<DataType>::map() const
{
    ASSERT ( M_map.get() != 0, "MatrixEpetra::getMap: Error: M_map pointer is null" );
    return *M_map;
}

template <typename DataType>
const MapEpetra& MatrixEpetra<DataType>::domainMap() const
{
    ASSERT ( M_domainMap.get() != 0, "MatrixEpetra::getdomainMap: Error: M_domainMap pointer is null" );
    return *M_domainMap;
}

template <typename DataType>
const MapEpetra& MatrixEpetra<DataType>::rangeMap() const
{
    ASSERT ( M_rangeMap.get() != 0, "MatrixEpetra::getRangeMap: Error: M_rangeMap pointer is null" );
    return *M_rangeMap;
}

template <typename DType>
MatrixEpetra<DType>* RAP (const MatrixEpetra<DType>& R, const MatrixEpetra<DType>& A, const MatrixEpetra<DType>& P)
{
    //Optimized implementation requires Trilinos 10.8
    /*
        typename MatrixEpetra<DType>::matrix_type * result = NULL;
        ML_Epetra::ML_Epetra_RAP (*A.matrixPtr(), *P.matrixPtr(), *R.matrixPtr(), result, true);

        MatrixEpetra<DType> * matrix(new MatrixEpetra<DType>(P.map()));
        matrix->M_epetraCrs.reset(result);
        matrix->M_domainMap = P.M_domainMap;
        matrix->M_rangeMap = R.M_rangeMap;

        return result;

    */
    // Slower implementation (no prerequisites on Trilinos)
    MatrixEpetra<DType>* result (new MatrixEpetra<DType> ( *R.M_rangeMap ) );
    MatrixEpetra<DType> tmp (A.map() );
    EpetraExt::MatrixMatrix::Multiply (*A.matrixPtr(), false, *P.matrixPtr(), false, *tmp.matrixPtr(), true);
    EpetraExt::MatrixMatrix::Multiply (*R.matrixPtr(), false, *tmp.matrixPtr(), false, * (result->matrixPtr() ), true);

    result->M_domainMap = P.M_domainMap;
    result->M_rangeMap = R.M_rangeMap;

    return result;
}

template <typename DType>
MatrixEpetra<DType>* PtAP (const MatrixEpetra<DType>& A, const MatrixEpetra<DType>& P)
{
    typename MatrixEpetra<DType>::matrix_type* result = NULL;
    ML_Epetra::ML_Epetra_PtAP (*A.matrixPtr(), *P.matrixPtr(), result, true);
    MatrixEpetra<DType>* matrix (new MatrixEpetra<DType> (P.M_domainMap) );
    matrix->M_epetraCrs.reset (result);
    matrix->M_domainMap = P.M_domainMap;
    matrix->M_rangeMap = P.M_domainMap;

    return matrix;
}


} // end namespace LifeV
//@@
//#undef OFFSET

#endif
