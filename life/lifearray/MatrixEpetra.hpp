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

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_MpiComm.h>
#include <Epetra_FECrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>
#include <life/lifearray/VectorEpetra.hpp>

#include <vector>

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

    //! Constructor for square and rectangular matrices
    /*!
      @param map Row map. The column map will be defined in MatrixEpetra<DataType>::GlobalAssemble(...,...)
      @param numEntries The average number of entries for each row.
     */
    MatrixEpetra( const MapEpetra& map, Int numEntries = 50 );

    //! Copy Constructor
    /*!
      @param matrix Matrix used to create the new occurence
     */
    MatrixEpetra( const MatrixEpetra& matrix);

    //! Copies the matrix to a matrix which resides only on the processor
    /*!
      @param matrix Matrix where the content of the matrix should be stored
      @param reduceToProc Processor where the matrix should be stored
     */
    MatrixEpetra( const MatrixEpetra& matrix, const UInt reduceToProc );

    //! Constructs an MatrixEpetra view of an Epetra_FECrsMatrix.
    /*!
      This constructor can be used when we need to modify an Epetra_FECrsMatrix
      using a method of the class MatrixEpetra
      @param crsMatrixPtr Pointer on a Epetra_FECrsMatrix of Trilinos
     */
    MatrixEpetra( matrix_ptrtype crsMatrixPtr );

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
    void swapCrsMatrix( matrix_ptrtype& matrixPtr );

    //! Swap the matrix with the one given as argument
    /*!
      @param matrix matrix which is swapped
     */
    void swapCrsMatrix( MatrixEpetra<DataType>& matrix );

    //! Multiply the MatrixEpetra by the first given matrix and put the result in the second given matrix
    /*!
      @param transposeCurrent If true, it transposes the MatrixEpetra
      @param matrix Matrix that multiply the MatrixEpetra
      @param transposeMatrix if true, it transposes the matrix "matrix"
      @param result Matrix to store the result
      @param callFillCompleteOnResult If true, the matrix "result" will be filled (i.e. closed) after the multiplication
     */
    Int multiply( bool transposeCurrent,
                  const MatrixEpetra<DataType> &matrix,
                  bool transposeMatrix,
                  MatrixEpetra<DataType> &result,
                  bool callFillCompleteOnResult=true ) const;

    //! Multiply the first VectorEpetra given as a parameter by the MatrixEpetra and put the result into the second given VectorEpetra
    /*!
      @param transposeCurrent If true, it transposes the MatrixEpetra
      @param vector1 Vector that will be multiply by the MatrixEpetra
      @param vector2 Vector that will store the result
     */
    Int multiply( bool transposeCurrent, const vector_type& vector1, vector_type& vector2 ) const;

    //! Add to the matrix a dyadic product: \f[ C += A \otimes B + \f], where \f[ C \f] is the matrix.
    /*!
     *  NOTE: This method has been tested only for square matrices.
      @param vector1 Vector A
      @param vector2 Vector B
     */
    void addDyadicProduct( const vector_type& vector1, const vector_type& vector2 );

    //! Add a multiple of a given matrix:  *this += scalar*matrix
    /*!
      @param scalar Scalar which multiplies the matrix
      @param matrix Matrix to be added
     */
    void add( const DataType scalar, const MatrixEpetra& matrix );

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
    void diagonalize( std::vector<UInt> rVec,
                      DataType const coefficient,
                      vector_type &rhs,
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
    void diagonalize( UInt const row, DataType const coefficient, vector_type &rhs,
                      DataType datum,
                      UInt offset = 0 );

    //! Save the matrix into a MatrixMarket (.mtx) file
    /*!
      @param filename file where the matrix will be saved
      @param headers boolean to write the MM headers or not
     */
    void matrixMarket( std::string const &fileName, const bool headers = true );

    //! Save the matrix into a Matlab (.m) file
    /*!
      @param filename file where the matrix will be saved
     */
    void spy( std::string const &fileName );

    //! Print the contents of the matrix
    /*!
      @param output Stream where the informations must be printed
     */
    void showMe( std::ostream& output = std::cout ) const;

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
    Int globalAssemble( const boost::shared_ptr<const MapEpetra> & domainMap,
                        const boost::shared_ptr<const MapEpetra> & rangeMap );

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
    void insertValueDiagonal( const DataType entry, const MapEpetra& Map, const UInt offset = 0 );

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
    void insertValueDiagonal( const DataType& value, Int from = -1, Int to = -2 );

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
    void insertOneDiagonal( Int from = -1, Int to = -2 );

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
    void insertZeroDiagonal( Int from = -1, Int to = -2 );

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

    //! Set a coefficient of the matrix
    /*!
      @param row Row index of the coefficient
      @param column column index of the coefficient
     */
    void setCoefficient( UInt row, UInt column, DataType localValue );

    //! Add a set of values to the corresponding set of coefficient in the matrix
    /*!
      @param numRows Number of rows into the list given in "localValues"
      @param numColumns Number of columns into the list given in "localValues"
      @param rowIndices List of row indices
      @param columnIndices List of column indices
      @param localValues 2D array containing the coefficient related to "rowIndices" and "columnIndices"
      @param format Format of the matrix (Epetra_FECrsMatrix::COLUMN_MAJOR or Epetra_FECrsMatrix::ROW_MAJOR)
     */
    void addToCoefficients( Int const numRows, Int const numColumns,
                            std::vector<Int> const& rowIndices, std::vector<Int> const& columnIndices,
                            DataType* const* const localValues,
                            Int format = Epetra_FECrsMatrix::COLUMN_MAJOR );

    //! Add a value at a coefficient of the matrix
    /*!
      @param row Row index of the value to be added
      @param column Column index of the value to be added
      @param localValue Value to be added to the coefficient
     */
    void addToCoefficient( UInt row, UInt column, DataType localValue );

    //! If set true, transpose of this operator will be applied.
    /*!
      @param transpose flag to identify whether to transpose or not
     */
    void setUseTranspose( const bool& transpose = false ) { M_epetraCrs->SetUseTranspose( transpose ); }

    //@}


    //! @name Get Methods
    //@{

    //! Return the shared_pointer of the Epetra_FECrsMatrix
    matrix_ptrtype& matrixPtr(){ return M_epetraCrs; }

    //! Return the const shared_pointer of the Epetra_FECrsMatrix
    const matrix_ptrtype& matrixPtr() const{ return M_epetraCrs; }

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
MatrixEpetra<DataType>::MatrixEpetra( const MapEpetra& map, Int numEntries ) :
    M_map       ( new MapEpetra( map ) ),
    M_epetraCrs ( new matrix_type( Copy, *M_map->map( Unique ), numEntries, false) )
{

}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra( const MatrixEpetra& matrix ) :
        M_map      ( matrix.M_map ),
        M_domainMap( matrix.M_domainMap ),
        M_rangeMap ( matrix.M_rangeMap ),
        M_epetraCrs( new matrix_type( *matrix.M_epetraCrs ) )
{

}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra( const MatrixEpetra& matrix, const UInt reduceToProc ) :
    M_map      ( matrix.getMap().createRootMap( reduceToProc ) ),
    M_epetraCrs( new matrix_type( Copy, *M_map->map( Unique ),
                                  numEntries( matrix.M_epetraCrs->Map().Comm().MyPID() == reduceToProc ) * 20,
                                  false ) )
{
    Epetra_Export reducedExport( M_epetraCrs->Map(), matrix.M_epetraCrs->Map() );
    M_epetraCrs->Import( *matrix.M_epetraCrs, reducedExport, Add );

    if (M_epetraCrs->Filled())
    {
        M_domainMap = matrix.getDomainMap().createRootMap(reduceToProc);
        M_rangeMap  = matrix.getRangeMap().createRootMap(reduceToProc);
    }
}

template <typename DataType>
MatrixEpetra<DataType>::MatrixEpetra( matrix_ptrtype CRSMatrixPtr ):
    M_map(),
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
#ifdef HAVE_TRILINOS_EPETRAEXT_31 // trilinos6
    EpetraExt::MatrixMatrix::Add( *matrix.matrixPtr(), false, 1., *this->matrixPtr(), 1., false );
#elif defined HAVE_TRILINOS_EPETRAEXT // trilinos8
    EpetraExt::MatrixMatrix::Add( *matrix.matrixPtr(), false, 1., *this->matrixPtr(), 1. );
#else
#error error: do not have nor EpetraExt 6 nor 7 or 8
#endif

    return *this;
}

template<typename DataType>
MatrixEpetra<DataType>&  MatrixEpetra<DataType>::operator=( const MatrixEpetra& matrix )
{
    M_map        = matrix.M_map;
    M_domainMap  = matrix.M_domainMap;
    M_rangeMap   = matrix.M_rangeMap;
    *M_epetraCrs = *( matrix.M_epetraCrs );

    return *this;
}

template<typename DataType>
typename MatrixEpetra<DataType>::vector_type
MatrixEpetra<DataType>::operator * ( const vector_type& vector ) const
{
    ASSERT_PRE( M_epetraCrs->Filled(),
                "MatrixEpetra::Operator*: globalAssemble(...) should be called first" );
    ASSERT_PRE( vector.map().mapsAreSimilar(*M_domainMap),
                "MatrixEpetra::Operator*: the map of vec is not the same of domainMap" );
    vector_type result(vector);

    M_epetraCrs->Apply( vector.epetraVector(), result.epetraVector() );

    return result;
}


template<typename DataType>
MatrixEpetra<DataType>&
MatrixEpetra<DataType>::operator *= ( const DataType value )
{
    M_epetraCrs->Scale( value );
    return *this;
}

template<typename DataType>
MatrixEpetra<DataType>
MatrixEpetra<DataType>::operator * ( const DataType scalar ) const
{
    MatrixEpetra<DataType> matrix( *this );
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
        matrix_ptrtype tmp( M_epetraCrs );
        M_epetraCrs.reset( new matrix_type( Copy, M_epetraCrs->RowMap(), meanNumEntries ) );

#ifdef HAVE_TRILINOS_EPETRAEXT_31 // trilinos6
        EpetraExt::MatrixMatrix::Add( *tmp, false, 1., *M_epetraCrs, 1., false );
#elif defined HAVE_TRILINOS_EPETRAEXT // trilinos8
        EpetraExt::MatrixMatrix::Add( *tmp, false, 1., *M_epetraCrs, 1. );
#else
#error error: do not have nor EpetraExt 6 nor 7 or 8
#endif
        M_domainMap.reset();
        M_rangeMap.reset();

    }
}

template <typename DataType>
void MatrixEpetra<DataType>::removeZeros()
{
    if ( M_epetraCrs->Filled() )
    {
        Int meanNumEntries = this->getMeanNumEntries();
        matrix_ptrtype tmp( M_epetraCrs );
        M_epetraCrs.reset(new matrix_type( Copy, M_epetraCrs->RowMap(), meanNumEntries ) );

        //Variables to store the informations
        Int NumEntries;
        Real* Values;
        Int* Indices;
        Int row(0);

        for ( Int i(0); i<tmp->NumGlobalRows(); ++i )
        {
            row = tmp->LRID( i );
            tmp->ExtractMyRowView( row, NumEntries, Values, Indices );

            Int Indices2[NumEntries];
            Real Values2[NumEntries];
            Int NumEntries2(0);

            for (Int j(0); j<NumEntries; ++j)
            {
                if (Values[j] != 0.0)
                {
                    Indices2[NumEntries2] = tmp->GCID(Indices[j]);
                    Values2[NumEntries2]  = Values[j];
                    NumEntries2++;
                }
            }
            M_epetraCrs->InsertGlobalValues( row, NumEntries2, Values2, Indices2 );
        }
        insertZeroDiagonal();
        M_epetraCrs->GlobalAssemble();
    }
}

template <typename DataType>
void MatrixEpetra<DataType>::swapCrsMatrix( matrix_ptrtype& matrixPtr )
{
    M_epetraCrs.swap( matrixPtr );
}


template <typename DataType>
void MatrixEpetra<DataType>::swapCrsMatrix( MatrixEpetra<DataType>& matrix )
{
    M_epetraCrs.swap( matrix.M_epetraCrs );
}

template <typename DataType>
Int MatrixEpetra<DataType>::multiply( bool transposeCurrent,
                                      const MatrixEpetra<DataType> &matrix, bool transposeMatrix,
                                      MatrixEpetra<DataType> &result, bool callFillCompleteOnResult ) const
{
    Int errCode = EpetraExt::MatrixMatrix::Multiply( *M_epetraCrs, transposeCurrent,
                                                     *matrix.matrixPtr(), transposeMatrix,
                                                     *result.matrixPtr(), false );
    if (callFillCompleteOnResult)
        result.globalAssemble();

    return errCode;
}

template <typename DataType>
Int MatrixEpetra<DataType>::multiply( bool transposeCurrent, const vector_type& vector1, vector_type &vector2 ) const
{
    ASSERT_PRE( M_epetraCrs->Filled(),
                "MatrixEpetra<DataType>::Multiply: GlobalAssemble(...) must be called first" );
    ASSERT_PRE( vector1.map().mapsAreSimilar(*M_domainMap),
                "MatrixEpetra<DataType>::Multiply: x map is different from M_domainMap" );
    ASSERT_PRE( vector2.map().mapsAreSimilar(*M_rangeMap),
                "MatrixEpetra<DataType>::Multiply: y map is different from M_rangeMap" );


    return M_epetraCrs->Multiply( transposeCurrent, vector1.epetraVector(), vector2.epetraVector() );
}

template <typename DataType>
void MatrixEpetra<DataType>::addDyadicProduct( const vector_type& vector1, const vector_type& vector2 )
{
    // Check if the matrix is open
    if ( M_epetraCrs->Filled() )
    {
        if ( M_epetraCrs->Comm().MyPID() == 0 )
            std::cout << "ERROR: Matrix already filled: cannot add dyadic product!" << std::endl;
        return;
    }

    // Creating a fully repeated copy of vector2
    VectorEpetra repeatedVector2( vector2, Repeated );
    //Epetra_Map repeatedMap( Epetra_Util::Create_OneToOne_Map ( vector2.epetraMap(), -1 ) );
    //boost::VectorEpetra repeatedVector2( repeatedMap, Repeated );
    //repeatedVector2 = vector2;

    // Add the new coefficients to the matrix
    for( Int row(0); row < M_epetraCrs->NumMyRows(); ++row )
        for( Int column(0); column < M_epetraCrs->NumGlobalCols(); ++column )
            addToCoefficient( M_epetraCrs->GRID(row), column, vector1[row]*repeatedVector2[column] );
}

template <typename DataType>
void MatrixEpetra<DataType>::add ( const DataType scalar, const MatrixEpetra& matrix )
{
#if defined HAVE_TRILINOS_EPETRAEXT // trilinos8
    EpetraExt::MatrixMatrix::Add( *matrix.matrixPtr(), false, scalar, *this->matrixPtr(), 1. );
#else
#error error: do not have nor EpetraExt  8+
#endif
}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize( std::vector<UInt> rVec, DataType const coefficient, UInt offset )
{

    const Epetra_Comm&  Comm( M_epetraCrs->Comm() );
    Int numProcs( Comm.NumProc() );
    Int MyPID   ( Comm.MyPID() );
    Int i;


    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    Int*  r;
    UInt* Ur;


    // loop on all proc
    for ( Int p(0); p < numProcs; p++ )
    {
        Int sizeVec( rVec.size() );

        Comm.Broadcast( &sizeVec, 1, p );

        if ( p == MyPID )
        {
            Ur = &rVec.front();
        }
        else
        {
            Ur = new UInt[sizeVec];
        }

        r = (Int*) Ur;

        Comm.Broadcast( r, sizeVec, p );

        for ( i=0; i < sizeVec; i++ )
            diagonalize( Ur[i], coefficient, offset );

        if ( p != MyPID )
        {
            delete[] Ur;
        }

    }

}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize( UInt const row,
                                          DataType const coefficient,
                                          UInt offset )
{

    if ( !M_epetraCrs->Filled() )
    { // if not filled, I do not know how to diagonalize.
        ERROR_MSG( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap( M_epetraCrs->RowMap() );
    const Epetra_Map& colMap( M_epetraCrs->ColMap() );


    Int myCol = colMap.LID( row + offset );

    // row: if r is mine, zero out values
    Int myRow = rowMap.LID( row + offset );

    if ( myRow >= 0 )  // I have this row
    {
        Int    NumEntries;
        Real* Values;
        Int* Indices;

        M_epetraCrs->ExtractMyRowView( myRow, NumEntries, Values, Indices );

        for (Int i(0); i <  NumEntries; i++)
        {
            Values[i] = 0;
        }

        DataType coeff( coefficient );
        M_epetraCrs->ReplaceMyValues( myRow, 1, &coeff, &myCol ); // A(r,r) = coefficient
    }

}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize( std::vector<UInt> rVec,
                                          DataType const coefficient,
                                          vector_type &rhs,
                                          std::vector<DataType> datumVec,
                                          UInt offset )
{


    const Epetra_Comm&  Comm( M_epetraCrs->Comm() );
    Int numProcs( Comm.NumProc() );
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

    const Epetra_Map& rowMap( M_epetraCrs->RowMap() );


    //Comm.Barrier();
    // we want to know which IDs are our or not

    for ( Int ii = 0; ii < (Int)rVec.size(); ++ii )
    {
        Int lID = rowMap.LID(rVec[ii]);
        if ( !( lID < 0 ) )
        {

            localIDs.push_back( rVec[ii] );
            localData.push_back( datumVec[ii] );
            localBC.insert( std::pair<Int, Real>( rVec[ii], datumVec[ii] ) );
        }
        else
        {
            remoteIDs.push_back( rVec[ii] );
            remoteData.push_back( datumVec[ii] );
        }
    }

    // now, we have to fill our localIDs with IDs from other processors
    // first, we have to build the map of all the remoteIDs and their processor owner


    Int numIDs = remoteIDs.size();

    Int* PIDList = new Int[numIDs];
    Int* LIDList = new Int[numIDs];

    rowMap.RemoteIDList( numIDs,
                         &remoteIDs[0],
                         PIDList,
                         LIDList );

    std::vector< std::vector<Int> > procToID  ( Comm.NumProc() );
    std::vector< std::vector<Real> > procToData( Comm.NumProc() );

    for ( Int ii = 0; ii < numIDs; ++ii )
    {
        Int pi = PIDList[ii];
        procToID[pi].push_back( remoteIDs[ii] );
        procToData[pi].push_back( remoteData[ii] );
    }

    // then, we send all the nodes where they belong

    const Epetra_MpiComm* comm = dynamic_cast<Epetra_MpiComm const*>( &Comm );

    assert( comm != 0 );

    for ( Int ii = 0; ii < (Int)procToID.size(); ++ii )
    {
        if ( ii != me )
        {
            Int length;
            length = procToID[ii].size();
            MPI_Send( &length, 1, MPI_INT, ii, 666, comm->Comm() );
            if ( length > 0 )
            {
                MPI_Send( &procToID[ii][0], length, MPI_INT, ii, 667, comm->Comm() );
                MPI_Send( &procToData[ii][0], length, MPI_DOUBLE, ii, 668, comm->Comm() );
            }
        }

    }

    for ( Int ii = 0; ii < (Int)procToID.size(); ++ii )
    {
        if ( ii != me )
        {
            Int length;
            MPI_Status status;
            MPI_Recv( &length, 1, MPI_INT, ii, 666, comm->Comm(), &status );


            if ( length > 0 )
            {
                Int* bufferID = new Int[length];
                Int* ptrID(0);//    = new Int[length];

                MPI_Recv( bufferID, length, MPI_INT, ii, 667, comm->Comm(), &status );

                ptrID = bufferID;

                Real* bufferData = new Real[length];
                Real* ptrData(0);

                MPI_Recv( bufferData, length, MPI_DOUBLE, ii, 668, comm->Comm(), &status );
                ptrData = bufferData;

                for ( Int ii = 0; ii < length; ++ii, ++ptrID, ++ptrData )
                {
                    localBC.insert( std::pair<Int, Real>
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
        diagonalize( im->first, coefficient, rhs, im->second, offset );
    }

    return;
#endif


    // loop on all proc
    for ( Int p(0); p < numProcs; p++ )
    {
        Int sizeVec(rVec.size());
        if ( sizeVec != Int(datumVec.size()))
        {
            // vectors must be of the same size
            ERROR_MSG( "diagonalize: vectors must be of the same size\n" );
        }

        Comm.Broadcast( &sizeVec, 1, p );

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

        Comm.Broadcast( r,    sizeVec, p );
        Comm.Broadcast( datum,sizeVec, p );

        for( i=0; i < sizeVec; i++ )
            diagonalize( Ur[i], coefficient, rhs, datum[i], offset );

        if ( p != MyPID )
        {
            delete[] Ur;
            delete[] datum;
        }

    }

}

template <typename DataType>
void MatrixEpetra<DataType>::diagonalize( UInt const row,
                                          DataType const coefficient,
                                          vector_type &rhs,
                                          DataType datum,
                                          UInt offset )
{

    if ( !M_epetraCrs->Filled() )
    {
        // if not filled, I do not know how to diagonalize.
        ERROR_MSG( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap( M_epetraCrs->RowMap() );
    const Epetra_Map& colMap( M_epetraCrs->ColMap() );


    Int myCol = colMap.LID( row + offset );

#ifdef EPETRAMATRIX_SYMMETRIC_DIAGONALIZE
    if ( myCol >= 0 )  // I have this column
    {
        Real zero(0);
        for ( Int i(0); i < rowMap.NumMyElements(); i++ )
            // Note that if a value is not already present for the specified location in the matrix,
            // the input value will be ignored and a positive warning code will be returned.
            M_epetraCrs->ReplaceMyValues(i,1, &zero, &myCol);
    }
#endif

    // row: if r is mine, zero out values
    Int myRow = rowMap.LID( row + offset );

    if ( myRow >= 0 )  // I have this row
    {
        Int    NumEntries;
        Real* Values;
        Int* Indices;

        M_epetraCrs->ExtractMyRowView( myRow, NumEntries, Values, Indices );

        for ( Int i(0); i <  NumEntries; i++ )
        {
            Values[i] = 0;
        }

        DataType coeff( coefficient );

        M_epetraCrs->ReplaceMyValues( myRow, 1, &coeff, &myCol ); // A(r,r) = coeff
        rhs[row + offset] = coefficient * datum; // correct right hand side for row r

    }

}

template <typename DataType>
void MatrixEpetra<DataType>::matrixMarket( std::string const &fileName, const bool headers )
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName;
    std::string desc = "Created by LifeV";

    name = fileName + ".mtx";

    EpetraExt::RowMatrixToMatrixMarketFile( name.c_str(),
                                            *M_epetraCrs,
                                            name.c_str(),
                                            desc.c_str(),
                                            headers
                                          );
}

template <typename DataType>
void MatrixEpetra<DataType>::spy( std::string const &fileName )
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName, uti = " , ";

    Int  me = M_epetraCrs->Comm().MyPID();
    std::ostringstream myStream;
    myStream << me;
    name = fileName + ".m";

    EpetraExt::RowMatrixToMatlabFile( name.c_str(), *M_epetraCrs );

}

template <typename DataType>
void MatrixEpetra<DataType>::showMe( std::ostream& output ) const
{
    output << "showMe must be implemented for the MatrixEpetra class" << std::endl;
}

template <typename DataType>
Int MatrixEpetra<DataType>::globalAssemble()
{
    if ( M_epetraCrs->Filled() )
    {
        return -1;
    }

    insertZeroDiagonal();
    M_domainMap = M_map;
    M_rangeMap  = M_map;
    return  M_epetraCrs->GlobalAssemble();
}

template <typename DataType>
Int MatrixEpetra<DataType>::globalAssemble( const boost::shared_ptr<const MapEpetra> & domainMap,
                                            const boost::shared_ptr<const MapEpetra> & rangeMap )
{
    if ( M_epetraCrs->Filled() )
    {
        return -1;
    }

    M_domainMap = domainMap;
    M_rangeMap  = rangeMap;
    return  M_epetraCrs->GlobalAssemble( *domainMap->map(Unique), *rangeMap->map(Unique) );
}

template <typename DataType>
void
MatrixEpetra<DataType>::insertValueDiagonal( const DataType entry, const MapEpetra& Map, const UInt offset )
{
    for ( Int i=0 ; i<Map.map(Unique)->NumMyElements(); ++i )
    {
        addToCoefficient( offset + Map.map(Unique)->GID(i) , offset + Map.map(Unique)->GID(i), entry );
    }
}

template <typename DataType>
void MatrixEpetra<DataType>::insertValueDiagonal( const DataType& value, Int from, Int to )
{
    if ( M_epetraCrs->Filled() )
    {
        if ( M_epetraCrs->Comm().MyPID() == 0 )
            std::cout << "Matrix is already filled, it is impossible to insert the diagonal now" << std::endl;
        return;
    }

    if ( to == from ) return;

    if ( to < from ) // do all entries
    {
        from = M_epetraCrs->RowMap().MinMyGID ();
        to = M_epetraCrs->RowMap().MaxMyGID () + 1;
    }

    Int* p =  M_epetraCrs->RowMap().MyGlobalElements();
    Int ierr;

    for ( Int i(0); i <  M_epetraCrs->RowMap().NumMyElements(); ++i, ++p )
    {
        if ( *p < from || *p >= to ) continue;

        ierr = M_epetraCrs->InsertGlobalValues( 1, p, 1, p, &value );

        if( ierr < 0 ) std::cout << " error in matrix insertion " << ierr << std::endl;
    }
}

template <typename DataType>
void MatrixEpetra<DataType>::insertOneDiagonal( Int from, Int to )
{
    insertValueDiagonal( 1.0, from, to );
}

template <typename DataType>
void MatrixEpetra<DataType>::insertZeroDiagonal( Int from, Int to )
{
    insertValueDiagonal( 0.0, from, to );
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
setCoefficient( UInt row, UInt column, DataType localValue )
{
    Int irow(    row );
    Int icol( column );

    Int ierr=M_epetraCrs->ReplaceGlobalValues( 1, &irow, 1, &icol, &localValue );
    if (ierr!=0)
        { std::cout << " error in matrix replacement " << ierr << std::endl; }

}

template <typename DataType>
void MatrixEpetra<DataType>::
addToCoefficient( UInt row, UInt column, DataType localValue )
{

    Int irow(    row );
    Int icol( column );

    Int ierr = M_epetraCrs->InsertGlobalValues( 1, &irow, 1, &icol, &localValue );

    if ( ierr < 0 ) std::cout << " error in matrix insertion " << ierr << std::endl;
}

template <typename DataType>
void MatrixEpetra<DataType>::
addToCoefficients( Int const numRows, Int const numColumns,
                   std::vector<Int> const& rowIndices, std::vector<Int> const& columnIndices,
                   DataType* const* const localValues,
                   Int format )
{
    Int ierr = M_epetraCrs->InsertGlobalValues( numRows, &rowIndices[0], numColumns,
                                                &columnIndices[0], localValues, format );

    if ( ierr < 0 ) std::cout << " error in matrix insertion " << ierr << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
template <typename DataType>
Int MatrixEpetra<DataType>::meanNumEntries() const
{
    const Int minEntries = M_epetraCrs->MaxNumEntries ()/2;
    if ( M_epetraCrs->NumMyRows() )
        return minEntries;

    Int meanNumEntries = M_epetraCrs->NumMyNonzeros()/M_epetraCrs->NumMyRows();
    if ( meanNumEntries < minEntries || meanNumEntries > 2*minEntries )
        return minEntries;
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
    ASSERT( M_map.get() != 0, "MatrixEpetra::getMap: Error: M_map pointer is null" );
    return *M_map;
}

template <typename DataType>
const MapEpetra& MatrixEpetra<DataType>::domainMap() const
{
    ASSERT( M_domainMap.get() != 0, "MatrixEpetra::getdomainMap: Error: M_domainMap pointer is null" );
    return *M_domainMap;
}

template <typename DataType>
const MapEpetra& MatrixEpetra<DataType>::rangeMap() const
{
    ASSERT( M_rangeMap.get() != 0, "MatrixEpetra::getRangeMap: Error: M_rangeMap pointer is null" );
    return *M_rangeMap;
}

} // end namespace LifeV
//@@
//#undef OFFSET

#endif
