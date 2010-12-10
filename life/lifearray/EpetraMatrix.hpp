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
    @brief EpetraMatrix

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
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

#include <lifeconfig.h>
#include <life/lifecore/life.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <vector>

//@@
//#define OFFSET 0

namespace LifeV
{

//! EpetraMatrix - The Epetra Matrix format Wrapper
/*!
 *  @author Gilles Fourestey, Simone Deparis, Gwenol Grandperrin
 *
 *  The EpetraMatrix class provides a general interface for the Epetra_FECrsMatrix of Trilinos.
 *
 *  TODO Reorder the class and to add more doxygen comments.
 */
template <typename DataType>
class EpetraMatrix
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_FECrsMatrix             matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef EpetraVector                   vector_type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor for square and rectangular matrices
    /*!
     * Constructor for square and rectangular matrices
     * @param _map: the row map. The column map will be defined in EpetraMatrix<DataType>::GlobalAssemble(...,...)
     * @param numEntries: the average number of entries for each row.
     * @param indexBase: the base index to address entries in the matrix (Usually 0 o 1)
     */
    EpetraMatrix( const EpetraMap& map, Int numEntries = 50, Int indexBase = 1 );

    //! Copy Constructor
    EpetraMatrix( const EpetraMatrix& matrix);

    //! Copies _matrix to a matrix which resides only on the processor "reduceToProc"
    EpetraMatrix( const EpetraMatrix& matrix, const UInt reduceToProc );

    /*! Constructs an EpetraMatrix view of an Epetra_FECrsMatrix. This constructor can be used
     *  when we need to modify an Epetra_FECrsMatrix using a method of the class EpetraMatrix.
     */
    EpetraMatrix( matrix_ptrtype CRSMatrixPtr );

    ~EpetraMatrix() {};

    //@}


    //! @name Operators
    //@{

    //! Unary operator+
    EpetraMatrix& operator += ( const EpetraMatrix& matrix );

    //! Assignment operator
    EpetraMatrix& operator= ( const EpetraMatrix& matrix );

    //! Matrix-Vector multiplication
    vector_type   operator * ( const vector_type& vector ) const;

    //! Unary operator* (it scales the matrix by the factor val)
    EpetraMatrix& operator *= ( const DataType scalar );

    //! const operator* (it returns a rescaled matrix this)
    EpetraMatrix  operator *  ( const DataType scalar ) const;

    //@}


    //! @name Methods
    //@{

    //! If the matrix has been filled, this function will reopen the Matrix
    void openCrsMatrix();

    //! This function removes all the zeros in the matrix and add zero on the diagonal
    void removeZeros();

    /*! Swap the given shared pointer with the one of the matrix
     *  @param p pointer on the matrix
     */
    void swapCrsMatrix( matrix_ptrtype& matrixPtr );

    /*! Swap the matrix with the one given as argument
     *  @param matrix matrix which is swapped
     */
    void swapCrsMatrix( EpetraMatrix<DataType>& matrix );

    /*! Multiply the EpetraMatrix by the first given matrix and put the result in the second given matrix
     *  @param transposeA if true, it transposes the EpetraMatrix
     *  @param B matrix that multiply the EpetraMatrix
     *  @param transposeB if true, it transposes the matrix B
     *  @param C matrix to store the result
     *  @param call_FillComplete_on_result if true, the matrix C will be filled (i.e. closed) after the multiplication
     */
    Int Multiply( bool transposeCurrent,
                  const EpetraMatrix<DataType> &matrix,
                  bool transposeMatrix,
                  EpetraMatrix<DataType> &result,
                  bool callFillCompleteOnResult=true ) const;

    /*! Multiply the first EpetraVector given as a parameter by the EpetraMatrix and put the result into the second given EpetraVector
     *  @param transposeA if true, it transposes the EpetraMatrix
     *  @param x vector that will be multiply by the EpetraMatrix
     *  @param y vector that will store the result
     */
    Int Multiply( bool transposeCurrent, const vector_type& vector1, vector_type &vector2 ) const;

    //! Add a multiple of a given matrix:  *this += val*_matrix
    void add( const DataType scalar, const EpetraMatrix& matrix );

    //! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
    void diagonalize ( std::vector<UInt> rVec, DataType const coefficient, UInt offset = 0 );

    //! set entry (r,r) to coeff and rest of row r to zero
    void diagonalize ( UInt const entryIndex, DataType const coefficient, UInt offset = 0 );

    /*! apply constraint on all rows rVec
     *  @param rVec vector of rows
     *  @param coeff value to set entry (r,r) at
     *  @param b right hand side Vector of the system to be adapted accordingly
     *  @param datumVec vector of values to constrain entry r of the solution at
     */
    void diagonalize( std::vector<UInt> rVec,
                      DataType const coefficient,
                      vector_type &rhs,
                      std::vector<DataType> datumVector,
                      UInt offset = 0 );

    /*! apply constraint on row r
     *  @param row row number
     *  @param coeff value to set entry (r,r) at
     *  @param b right hand side vector of the system to be adapted accordingly
     *  @param datum value to constrain entry r of the solution at
     */
    void diagonalize( UInt const row, DataType const coefficient, vector_type &rhs,
                      DataType datum,
                      UInt offset = 0 );

    /*! Save the matrix into a MatrixMarket (.mtx) file
     *  @param filename file where the matrix will be saved
     *  @param headers boolean to write the MM headers or not
     */

    void matrixMarket( std::string const &fileName, const bool headers = true );

    /*! Save the matrix into a Matlab (.m) file
     *  @param filename file where the matrix will be saved
     */
    void spy( std::string const &fileName );

    void showMe( std::ostream& output = std::cout ) const;

    //! Global assemble of a square matrix with default domain and range map
    /*
     * !
     * 1) Calls insertZeroDiagonal and then Epetra_FECsrMatrix::GlobalAssemble();
     * 2) Set M_domainMap and M_rangeMap
     *
     * EpetraFECsrMatrix will assume that both the domain and range map are the same
     * of the row map defined in the constructor.
     *
     * NOTE: domain and range map must be one-to-one and onto. (Unique map)
     *
     */
    Int GlobalAssemble();

    //! Global assemble for rectangular matrices.
    /*
     * !
     * 1) Calls Epetra_FECsrMatrix::GlobalAssemble(Epetra_Map & domainMap, Epetra_Map & rangeMap);
     * 2) Set M_domainMap and M_rangeMap
     * Input:
     * @param domainMap the domain map
     * @param rangeMap the range map
     */
    Int GlobalAssemble( const boost::shared_ptr<const EpetraMap> & domainMap,
                        const boost::shared_ptr<const EpetraMap> & rangeMap );

    //! insert the given value into the diagonal
    /*! Pay intention that this will add values to the diagonal,
        so for later added values with set_mat_inc, the one
        will be added
        Inserts Value on the local diagonal for diagonal elements specified by the input EpetraMap;
        This methods works only if matrix is not closed.

        \param entry: the entry that is inserted in the diagonal
        \param Map: the EpetraMap
        \param offset: an offset for the insertion of the diagonal entries
        \param replace: if the matrix is filled and the diagonal entries are present this bool should be set to true.
        Otherwise is the matrix is not filled it should be set to false.
    */
    void insertValueDiagonal( const DataType entry, const EpetraMap& Map, const UInt offset = 0 );

    //! insert the given value into the diagonal
    /*! Pay intention that this will add values to the diagonal,
        so for later added values with set_mat_inc, the one
    will be added
    Inserts Value on the local diagonal for diagonal elements >= from and < to;
    If from > to, process all diagonal entries entries
    If from = to, do nothing
    This methods works only if matrix is not closed.
    */
    void insertValueDiagonal( const DataType& value, Int from = -1, Int to = -2 );

    //! insert ones into the diagonal to ensure the matrix' graph has a entry there
    //! Pay intention that this will add ones to the diagonal,
    //! so for later added values with set_mat_inc, the one
    //! will be added
    /** Inserts Zero on the local diagonal for diagonal elements >= from and < to;
    If from > to, process all diagonal entries entries
    If from = to, do nothing
    This methods works only if matrix is not closed.
    */
    void insertOneDiagonal( Int from = -1, Int to = -2 );

    //! insert zeros into the diagonal to ensure the matrix' graph has a entry there
    //! This method does not remove non zero entries in the diagonal.
    /** Inserts Zero on the local diagonal for diagonal elements >= from and < to;
    If from > to, process all diagonal entries entries
    If from = to, do nothing
    This methods works only if matrix is not closed.
    */
    void insertZeroDiagonal( Int from = -1, Int to = -2 );

    //! Compute the norm 1 of the global matrix
    /*!
     * @return norm 1
     */
    Real NormOne() const;

    //! Compute the norm inf of the global matrix
    /*!
     * @return norm inf
     */
    Real NormInf() const;

    //@}


    //! @name Set Methods
    //@{

    void set_mat( UInt row, UInt column, DataType localValue );

    void set_mat_inc( Int const numRows, Int const numColumns,
                      std::vector<Int> const row, std::vector<Int> const column,
                      DataType* const* const localValue,
                      Int format = Epetra_FECrsMatrix::COLUMN_MAJOR );

    void set_mat_inc( UInt row, UInt column, DataType localValue );

    //@}


    //! @name Get Methods
    //@{

    //! Return the shared_pointer of the Epetra_FECrsMatrix
    matrix_ptrtype& getMatrixPtr(){ return M_epetraCrs; }

    //! Return the const shared_pointer of the Epetra_FECrsMatrix
    const matrix_ptrtype& getMatrixPtr() const{ return M_epetraCrs; }

    Int getMeanNumEntries() const ;

    Int MyPID();

    //! Return the row EpetraMap of the EpetraMatrix used in the assembling
    /*!
     * This method should be call when EpetraMap is still open.
     */
    const EpetraMap& getMap() const;

    //! Return the domain EpetraMap of the EpetraMatrix
    /*!
     * This function should be called only after EpetraMatrix<DataType>::GlobalAssemble(...) has been called.
     * If this is an open matrix that M_domainMap is an invalid pointer
     */
    const EpetraMap& getDomainMap() const;

    //! Return the range EpetraMap of the EpetraMatrix
    /*!
     * This function should be called only after EpetraMatrix<DataType>::GlobalAssemble(...) has been called.
     * If this is an open matrix that M_domainMap is an invalid pointer
     */
    const EpetraMap& getRangeMap() const;

    //@}

private:


    //! Shared pointer on the row EpetraMap used in the assembling
    boost::shared_ptr< EpetraMap > M_map;

    //! Shared pointer on the domain EpetraMap.
    /*
     * !
     * if y = this*x,
     * then x.getMap() is the domain map.
     * NOTE: Epetra assume the domain map to be 1-1 and onto (Unique)
     * M_domainMap is a NULL pointer until EpetraMatrix<DataType> is called.
     */
    boost::shared_ptr< const EpetraMap > M_domainMap;

    //! Shared pointer on the range EpetraMap.
    //! Shared pointer on the domain EpetraMap.
    /*
     * !
     * if y = this*x,
     * then y.getMap() is the range map.
     * NOTE: Epetra assume the domain map to be 1-1 and onto (Unique)
     * M_rangeMap is a NULL pointer until EpetraMatrix<DataType> is called.
     */
    boost::shared_ptr< const EpetraMap > M_rangeMap;


    //!Pointer on a Epetra_FECrsMatrix
    matrix_ptrtype  M_epetraCrs;

    //! The index base for the matrix
    Int             M_indexBase;

};


// ===================================================
// Constructors & Destructor
// ===================================================
template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMap& map, Int numEntries, Int indexBase ) :
    M_map       ( new EpetraMap( map ) ),
    M_epetraCrs ( new matrix_type( Copy, *M_map->getMap( Unique ), numEntries, false) ),
    M_indexBase ( indexBase )
{

}

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& matrix ) :
        M_map      ( matrix.M_map ),
        M_domainMap( matrix.M_domainMap ),
        M_rangeMap ( matrix.M_rangeMap ),
        M_epetraCrs( new matrix_type( *matrix.M_epetraCrs ) ),
        M_indexBase( matrix.M_indexBase )
{

}

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& matrix, const UInt reduceToProc ) :
    M_map      ( matrix.getMap().createRootMap( reduceToProc ) ),
    M_epetraCrs( new matrix_type( Copy, *M_map->getMap( Unique ),
                                  numEntries( matrix.M_epetraCrs->Map().Comm().MyPID() == reduceToProc ) * 20,
                                  false ) ),
    M_indexBase( matrix.M_indexBase )
{
    Int  me = M_epetraCrs->Comm().MyPID();
    if (!me)
        std::cout << "matrix.M_epetraCrs->Map().IndexBase() = "
                  << matrix.M_epetraCrs->Map().IndexBase()
                  << std::endl
                  << "M_epetraCrs->Map().IndexBase() = "
                  << M_epetraCrs->Map().IndexBase()
                  << std::endl;

    Epetra_Export reducedExport( M_epetraCrs->Map(), matrix.M_epetraCrs->Map() );
    M_epetraCrs->Import( *matrix.M_epetraCrs, reducedExport, Add );

    if (M_epetraCrs->Filled())
    {
        M_domainMap = matrix.getDomainMap().createRootMap(reduceToProc);
        M_rangeMap  = matrix.getRangeMap().createRootMap(reduceToProc);
    }
}

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( matrix_ptrtype CRSMatrixPtr ):
    M_map(),
    M_domainMap(),
    M_rangeMap(),
    M_indexBase( CRSMatrixPtr->Map().IndexBase() )
{
    M_epetraCrs = CRSMatrixPtr;
}


// ===================================================
// Operators
// ===================================================
template <typename DataType>
EpetraMatrix<DataType>&
EpetraMatrix<DataType>::operator += ( const EpetraMatrix& matrix )
{
#ifdef HAVE_TRILINOS_EPETRAEXT_31 // trilinos6
    EpetraExt::MatrixMatrix::Add( *matrix.getMatrixPtr(), false, 1., *this->getMatrixPtr(), 1., false );
#elif defined HAVE_TRILINOS_EPETRAEXT // trilinos8
    EpetraExt::MatrixMatrix::Add( *matrix.getMatrixPtr(), false, 1., *this->getMatrixPtr(), 1. );
#else
#error error: do not have nor EpetraExt 6 nor 7 or 8
#endif

    return *this;
}

template<typename DataType>
EpetraMatrix<DataType>&  EpetraMatrix<DataType>::operator=( const EpetraMatrix& matrix )
{
    M_map        = matrix.M_map;
    M_domainMap  = matrix.M_domainMap;
    M_rangeMap   = matrix.M_rangeMap;
    *M_epetraCrs = *( matrix.M_epetraCrs );
    M_indexBase  = matrix.M_indexBase;

    return *this;
}

template<typename DataType>
typename EpetraMatrix<DataType>::vector_type
EpetraMatrix<DataType>::operator * ( const vector_type& vector ) const
{
    ASSERT_PRE( M_epetraCrs->Filled(),
                "EpetraMatrix::Operator*: globalAssemble(...) should be called first" );
    ASSERT_PRE( vec.getMap().MapsAreSimilar(*M_domainMap),
                "EpetraMatrix::Operator*: the map of vec is not the same of domainMap" );
    vector_type result(vector);

    M_epetraCrs->Apply( vector.getEpetraVector(), result.getEpetraVector() );

    return result;
}


template<typename DataType>
EpetraMatrix<DataType>&
EpetraMatrix<DataType>::operator *= ( const DataType value )
{
    M_epetraCrs->Scale( value );
    return *this;
}

template<typename DataType>
EpetraMatrix<DataType>
EpetraMatrix<DataType>::operator * ( const DataType scalar ) const
{
    EpetraMatrix<DataType> matrix( *this );
    return matrix *= scalar;
}

// ===================================================
// Methods
// ===================================================
template <typename DataType>
void EpetraMatrix<DataType>::openCrsMatrix()
{
    if ( M_epetraCrs->Filled() )
    {
        Int meanNumEntries = this->getMeanNumEntries();
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
void EpetraMatrix<DataType>::removeZeros()
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
            row = tmp->LRID( i+M_indexBase );
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
void EpetraMatrix<DataType>::swapCrsMatrix( matrix_ptrtype& matrixPtr )
{
    M_epetraCrs.swap( matrixPtr );
}


template <typename DataType>
void EpetraMatrix<DataType>::swapCrsMatrix( EpetraMatrix<DataType>& matrix )
{
    M_epetraCrs.swap( matrix.M_epetraCrs );
}

template <typename DataType>
Int EpetraMatrix<DataType>::Multiply( bool transposeCurrent,
                                      const EpetraMatrix<DataType> &matrix, bool transposeMatrix,
                                      EpetraMatrix<DataType> &result, bool callFillCompleteOnResult ) const
{
    Int errCode = EpetraExt::MatrixMatrix::Multiply( *M_epetraCrs, transposeCurrent,
                                                     *matrix.getMatrixPtr(), transposeMatrix,
                                                     *result.getMatrixPtr(), false );
    if (callFillCompleteOnResult)
        result.GlobalAssemble();

    return errCode;
}

template <typename DataType>
Int EpetraMatrix<DataType>::Multiply( bool transposeCurrent, const vector_type& vector1, vector_type &vector2 ) const
{
    ASSERT_PRE( M_epetraCrs->Filled(),
                "EpetraMatrix<DataType>::Multiply: GlobalAssemble(...) must be called first" );
    ASSERT_PRE( vector1.getMap().MapsAreSimilar(*M_domainMap),
                "EpetraMatrix<DataType>::Multiply: x map is different from M_domainMap" );
    ASSERT_PRE( vector2.getMap().MapsAreSimilar(*M_rangeMap),
                "EpetraMatrix<DataType>::Multiply: y map is different from M_rangeMap" );


    return M_epetraCrs->Multiply( transposeCurrent, vector1.getEpetraVector(), vector2.getEpetraVector() );
}

template <typename DataType>
void EpetraMatrix<DataType>::add ( const DataType scalar, const EpetraMatrix& matrix )
{
#if defined HAVE_TRILINOS_EPETRAEXT // trilinos8
    EpetraExt::MatrixMatrix::Add( *matrix.getMatrixPtr(), false, scalar, *this->getMatrixPtr(), 1. );
#else
#error error: do not have nor EpetraExt  8+
#endif
}

template <typename DataType>
void EpetraMatrix<DataType>::diagonalize ( std::vector<UInt> rVec,
                                           DataType const coefficient,
                                           UInt offset )
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
void EpetraMatrix<DataType>::diagonalize( UInt const row,
                                          DataType const coefficient,
                                          UInt offset )
{

    if ( !M_epetraCrs->Filled() )
    { // if not filled, I do not know how to diagonalize.
        ERROR_MSG( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap( M_epetraCrs->RowMap() );
    const Epetra_Map& colMap( M_epetraCrs->ColMap() );


    Int myCol = colMap.LID( row + M_indexBase + offset );

    // row: if r is mine, zero out values
    Int myRow = rowMap.LID( row + M_indexBase + offset );

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
void EpetraMatrix<DataType>::diagonalize( std::vector<UInt> rVec,
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
        Int lID = rowMap.LID(rVec[ii] + 1);
        if ( !( lID < 0 ) )
        {

            localIDs.push_back( rVec[ii] + 1 );
            localData.push_back( datumVec[ii] );
            localBC.insert( std::pair<Int, Real>( rVec[ii] + 1, datumVec[ii] ) );
        }
        else
        {
            remoteIDs.push_back( rVec[ii] + 1 );
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
    std::vector< std::vector<Int> > procToData( Comm.NumProc() );


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
                MPI_Send( &procToData[ii][0], length, MPI_INT, ii, 668, comm->Comm() );
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

                MPI_Recv( bufferData, length, MPI_INT, ii, 668, comm->Comm(), &status );
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
        diagonalize( im->first - 1, coefficient, rhs, im->second, offset );
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
void EpetraMatrix<DataType>::diagonalize( UInt const row,
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


    Int myCol = colMap.LID( row + M_indexBase + offset );

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
    Int myRow = rowMap.LID( row + M_indexBase + offset );

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
        rhs[row + M_indexBase + offset] = coefficient * datum; // correct right hand side for row r // BASEINDEX + M_indexBase

    }

}

template <typename DataType>
void EpetraMatrix<DataType>::matrixMarket( std::string const &fileName, const bool headers )
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
void EpetraMatrix<DataType>::spy( std::string const &fileName )
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
void EpetraMatrix<DataType>::showMe( std::ostream& output ) const
{
    output << "showMe must be implemented for the EpetraMatrix class" << std::endl;
}

template <typename DataType>
Int EpetraMatrix<DataType>::GlobalAssemble()
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
Int EpetraMatrix<DataType>::GlobalAssemble( const boost::shared_ptr<const EpetraMap> & domainMap,
                                            const boost::shared_ptr<const EpetraMap> & rangeMap )
{
    if ( M_epetraCrs->Filled() )
    {
        return -1;
    }

    M_domainMap = domainMap;
    M_rangeMap  = rangeMap;
    return  M_epetraCrs->GlobalAssemble( *domainMap->getMap(Unique), *rangeMap->getMap(Unique) );
}

template <typename DataType>
void
EpetraMatrix<DataType>::insertValueDiagonal( const DataType entry, const EpetraMap& Map, const UInt offset )
{
    for ( UInt i=0 ; i<Map.getMap(Unique)->NumMyElements(); ++i )//num from 1
    {
        set_mat_inc( offset + Map.getMap(Unique)->GID(i)-1 , offset + Map.getMap(Unique)->GID(i)-1, entry );
    }
}

template <typename DataType>
void EpetraMatrix<DataType>::insertValueDiagonal( const DataType& value, Int from, Int to )
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
void EpetraMatrix<DataType>::insertOneDiagonal( Int from, Int to )
{
    insertValueDiagonal( 1.0, from, to );
}

template <typename DataType>
void EpetraMatrix<DataType>::insertZeroDiagonal( Int from, Int to )
{
    insertValueDiagonal( 0.0, from, to );
}

template <typename DataType>
Real EpetraMatrix<DataType>::NormOne() const
{
    return M_epetraCrs->NormOne();
}

template <typename DataType>
Real EpetraMatrix<DataType>::NormInf() const
{
    return M_epetraCrs->NormInf();
}


// ===================================================
// Set Methods
// ===================================================
template <typename DataType>
void EpetraMatrix<DataType>::
set_mat( UInt row, UInt column, DataType localValue )
{
    // incrementing row and cols by indexBase;
    Int irow(    row + M_indexBase );
    Int icol( column + M_indexBase );

    Int ierr=M_epetraCrs->ReplaceGlobalValues( 1, &irow, 1, &icol, &localValue );
    if (ierr!=0)
        { std::cout << " error in matrix replacement " << ierr << std::endl; }

}

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat_inc( UInt row, UInt column, DataType localValue )
{

    // incrementing row and cols by indexBase;
    Int irow(    row + M_indexBase );
    Int icol( column + M_indexBase );

    Int ierr = M_epetraCrs->InsertGlobalValues( 1, &irow, 1, &icol, &localValue );

    if ( ierr < 0 ) std::cout << " error in matrix insertion " << ierr << std::endl;
}

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat_inc( Int const numRows, Int const numColumns,
             std::vector<Int> const row, std::vector<Int> const column,
             DataType* const* const localValue,
             Int format )
{

    // incrementing row and cols by indexBase;
    std::vector<Int> irow( numRows );
    std::vector<Int> icol( numColumns );

    std::vector<Int>::const_iterator pt;

    pt = row.begin();
    for ( std::vector<Int>::iterator i( irow.begin() ); i !=  irow.end() && pt != row.end(); ++i, ++pt )
        *i = *pt + M_indexBase;

    pt = column.begin();
    for ( std::vector<Int>::iterator i( icol.begin() ); i !=  icol.end() && pt != column.end(); ++i, ++pt )
        *i = *pt + M_indexBase;


    Int ierr = M_epetraCrs->InsertGlobalValues( numRows, &irow[0], numColumns, &icol[0], localValue, format );

    if ( ierr < 0 ) std::cout << " error in matrix insertion " << ierr << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
template <typename DataType>
Int EpetraMatrix<DataType>::getMeanNumEntries() const
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
Int EpetraMatrix<DataType>::MyPID()
{
    return  M_epetraCrs->Comm().MyPID();
}

template <typename DataType>
const EpetraMap& EpetraMatrix<DataType>::getMap() const
{
    ASSERT( M_map.get() != 0, "EpetraMatrix::getMap: Error: M_map pointer is null" );
    return *M_map;
}

template <typename DataType>
const EpetraMap& EpetraMatrix<DataType>::getDomainMap() const
{
    ASSERT( M_domainMap.get() != 0, "EpetraMatrix::getdomainMap: Error: M_domainMap pointer is null" );
    return *M_domainMap;
}

template <typename DataType>
const EpetraMap& EpetraMatrix<DataType>::getRangeMap() const
{
    ASSERT( M_rangeMap.get() != 0, "EpetraMatrix::getRangeMap: Error: M_rangeMap pointer is null" );
    return *M_rangeMap;
}

} // end namespace LifeV
//@@
//#undef OFFSET

#endif
