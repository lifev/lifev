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

#include <lifeconfig.h>
#include <Epetra_MpiComm.h>
#include <Epetra_FECrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
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
    EpetraMatrix( const EpetraMap&    _map, Int numEntries = 50, Int indexBase = 1 );

    //! Copy Constructor
    EpetraMatrix( const EpetraMatrix& _matrix);
    //     EpetraMatrix( const EpetraMatrix<DataType> &_epetra );

    //! Copies _matrix to a matrix which resides only on the processor "reduceToProc"
    EpetraMatrix( const EpetraMatrix& _matrix, const UInt reduceToProc);

    /*! Constructs an EpetraMatrix view of an Epetra_FECrsMatrix. This constructor can be used
     *  when we need to modify an Epetra_FECrsMatrix using a method of the class EpetraMatrix.
     */
    EpetraMatrix( matrix_ptrtype CRSMatrixPtr );

    ~EpetraMatrix() {};

    //@}


    //! @name Operators
    //@{

    //! Unary operator+
    EpetraMatrix& operator += (const EpetraMatrix& _matrix);

    //! Assignment operator
    EpetraMatrix&   operator=   (const EpetraMatrix& _matrix);

    //! Matrix-Vector multiplication
    vector_type     operator *  (const vector_type& vec) const;

    //! Unary operator* (it scales the matrix by the factor val)
    EpetraMatrix&   operator *= (const DataType     val);

    //! const operator* (it returns a rescaled matrix this)
    EpetraMatrix    operator *  (const DataType     val) const;

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
    void swapCrsMatrix(matrix_ptrtype& p);

    /*! Swap the matrix with the one given as argument
     *  @param B matrix which is swapped
     */
    void swapCrsMatrix(EpetraMatrix<DataType>& B);

    /*! Multiply the EpetraMatrix by the first given matrix and put the result in the second given matrix
     *  @param transposeA if true, it transposes the EpetraMatrix
     *  @param B matrix that multiply the EpetraMatrix
     *  @param transposeB if true, it transposes the matrix B
     *  @param C matrix to store the result
     *  @param call_FillComplete_on_result if true, the matrix C will be filled (i.e. closed) after the multiplication
     */
    Int Multiply(bool transposeA,
                 const EpetraMatrix<DataType> &B, bool transposeB,
                 EpetraMatrix<DataType> &C, bool call_FillComplete_on_result=true) const;

    /*! Multiply the first EpetraVector given as a parameter by the EpetraMatrix and put the result into the second given EpetraVector
     *  @param transposeA if true, it transposes the EpetraMatrix
     *  @param x vector that will be multiply by the EpetraMatrix
     *  @param y vector that will store the result
     */
    Int Multiply(bool transposeA, const vector_type& x, vector_type &y) const;

    //! Add a multiple of a given matrix:  *this += val*_matrix
    void add(const DataType val, const EpetraMatrix& _matrix);

    //! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
    void diagonalize ( std::vector<UInt> rVec, DataType const coeff, UInt offset=0 );

    //! set entry (r,r) to coeff and rest of row r to zero
    void diagonalize ( UInt const r, DataType const coeff, UInt offset=0 );

    /*! apply constraint on all rows rVec
     *  @param rVec vector of rows
     *  @param coeff value to set entry (r,r) at
     *  @param b right hand side Vector of the system to be adapted accordingly
     *  @param datumVec vector of values to constrain entry r of the solution at
     */
    void diagonalize( std::vector<UInt> rVec,
                      DataType const coeff,
                      vector_type &b,
                      std::vector<DataType> datumVec,
                      UInt offset=0 );

    /*! apply constraint on row r
     *  @param r row number
     *  @param coeff value to set entry (r,r) at
     *  @param b right hand side vector of the system to be adapted accordingly
     *  @param datum value to constrain entry r of the solution at
     */
    void diagonalize( UInt const r, DataType const coeff, vector_type &b,
                      DataType datum,
                      UInt offset=0 );

    /*! Save the matrix into a MatrixMarket (.mtx) file
     *  @param filename file where the matrix will be saved
     *  @param headers boolean to write the MM headers or not
     */

    void matrixMarket( std::string const &filename, const bool headers = true );

    /*! Save the matrix into a Matlab (.m) file
     *  @param filename file where the matrix will be saved
     */
    void spy    ( std::string const &filename );

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
    Int GlobalAssemble(const boost::shared_ptr<const EpetraMap> & domainMap,
                       const boost::shared_ptr<const EpetraMap> & rangeMap);

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
    void insertValueDiagonal(const DataType& value, Int from = -1, Int to = -2);

    //! insert ones into the diagonal to ensure the matrix' graph has a entry there
    //! Pay intention that this will add ones to the diagonal,
    //! so for later added values with set_mat_inc, the one
    //! will be added
    /** Inserts Zero on the local diagonal for diagonal elements >= from and < to;
    If from > to, process all diagonal entries entries
    If from = to, do nothing
    This methods works only if matrix is not closed.
    */
    void insertOneDiagonal(Int from = -1, Int to = -2);

    //! insert zeros into the diagonal to ensure the matrix' graph has a entry there
    //! This method does not remove non zero entries in the diagonal.
    /** Inserts Zero on the local diagonal for diagonal elements >= from and < to;
    If from > to, process all diagonal entries entries
    If from = to, do nothing
    This methods works only if matrix is not closed.
    */
    void insertZeroDiagonal(Int from = -1, Int to = -2);

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

    void set_mat( UInt row, UInt col, DataType loc_val );

    void set_mat_inc( Int const numRows, Int const numCols,
                      std::vector<Int> const row, std::vector<Int> const col,
                      DataType* const* const loc_val,
                      Int format = Epetra_FECrsMatrix::COLUMN_MAJOR);

    void set_mat_inc( UInt row, UInt col, DataType loc_val );

    //@}


    //! @name Get Methods
    //@{

    //! Return the shared_pointer of the Epetra_FECrsMatrix
    matrix_ptrtype& getMatrixPtr() {return M_epetraCrs;}

    //! Return the const shared_pointer of the Epetra_FECrsMatrix
    const matrix_ptrtype& getMatrixPtr() const {return M_epetraCrs;}

    Int getMeanNumEntries() const ;

    Int MyPID() { return  M_epetraCrs->Comm().MyPID(); }

    //! Return the row EpetraMap of the EpetraMatrix used in the assembling
    /*!
     * This method should be call when EpetraMap is still open.
     */
    const EpetraMap& getMap() const
    {
        ASSERT( M_map.get() != 0, "EpetraMatrix::getMap: Error: M_map pointer is null" );
        return *M_map;
    }

    //! Return the domain EpetraMap of the EpetraMatrix
    /*!
     * This function should be called only after EpetraMatrix<DataType>::GlobalAssemble(...) has been called.
     * If this is an open matrix that M_domainMap is an invalid pointer
     */
    const EpetraMap& getDomainMap() const
    {
        ASSERT( M_domainMap.get() != 0, "EpetraMatrix::getdomainMap: Error: M_domainMap pointer is null" );
        return *M_domainMap;
    }

    //! Return the range EpetraMap of the EpetraMatrix
    /*!
     * This function should be called only after EpetraMatrix<DataType>::GlobalAssemble(...) has been called.
     * If this is an open matrix that M_domainMap is an invalid pointer
     */
    const EpetraMap& getRangeMap() const
    {
        ASSERT( M_rangeMap.get() != 0, "EpetraMatrix::getRangeMap: Error: M_rangeMap pointer is null" );
        return *M_rangeMap;
    }

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
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMap& _map, Int numEntries, Int indexBase ):
        M_map      ( new EpetraMap  (_map)),
        M_epetraCrs   ( new matrix_type( Copy, *M_map->getMap(Unique), numEntries, false)),
        M_indexBase   ( indexBase )
{
}

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& _matrix):
        M_map(_matrix.M_map),
        M_domainMap(_matrix.M_domainMap),
        M_rangeMap(_matrix.M_rangeMap),
        M_epetraCrs(new matrix_type(*_matrix.M_epetraCrs)),
        M_indexBase(_matrix.M_indexBase )
{
}

// Copies _matrix to a matrix which resides only on the processor "reduceToProc"
template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& _matrix, const UInt reduceToProc):
        M_map   (_matrix.getMap().createRootMap(reduceToProc)),
        M_epetraCrs   ( new matrix_type( Copy, *M_map->getMap(Unique),
                                         numEntries(_matrix.M_epetraCrs->Map().Comm().MyPID() == reduceToProc) * 20,
                                         false) ),
        M_indexBase   (_matrix.M_indexBase)
{
    Int  me    = M_epetraCrs->Comm().MyPID();
    if (!me)
        std::cout << "_matrix.M_epetraCrs->Map().IndexBase() = "
                  << _matrix.M_epetraCrs->Map().IndexBase()
                  << std::endl
                  << "M_epetraCrs->Map().IndexBase() = "
                  << M_epetraCrs->Map().IndexBase()
                  << std::endl;
    //     std::cout << "reduced Map : " << M_epetraCrs->Map()
    //               <<    endl;

    //     M_epetraCrs->Comm().Barrier();

    //     std::cout << "original Map: " << _matrix.M_epetraCrs->Map()
    //               << std::endl;

    //     M_epetraCrs->Comm().Barrier();

    Epetra_Export reducedExport(M_epetraCrs->Map(), _matrix.M_epetraCrs->Map());
    M_epetraCrs->Import(*_matrix.M_epetraCrs, reducedExport, Add);

    if (M_epetraCrs->Filled())
    {
        M_domainMap = _matrix.getDomainMap().createRootMap(reduceToProc);
        M_rangeMap = _matrix.getRangeMap().createRootMap(reduceToProc);
    }
}

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( matrix_ptrtype CRSMatrixPtr ):
        M_map(),
        M_domainMap(),
        M_rangeMap(),
        M_indexBase(CRSMatrixPtr->Map().IndexBase())
{
    M_epetraCrs=CRSMatrixPtr;
}


// ===================================================
// Operators
// ===================================================
template <typename DataType>
EpetraMatrix<DataType>&
EpetraMatrix<DataType>::operator += ( const EpetraMatrix& _matrix)
{
    //    EpetraMatrix matrix(Copy, _matrix.RowMap(), _matrix.GlobalMaxNumEntries());

#ifdef HAVE_TRILINOS_EPETRAEXT_31 // trilinos6
    EpetraExt::MatrixMatrix::Add(*_matrix.getMatrixPtr(), false, 1., *this->getMatrixPtr(), 1., false);
#elif defined HAVE_TRILINOS_EPETRAEXT // trilinos8
    EpetraExt::MatrixMatrix::Add(*_matrix.getMatrixPtr(), false, 1., *this->getMatrixPtr(), 1.);
#else
#error error: do not have nor EpetraExt 6 nor 7 or 8
#endif

    return *this;
    //     return matrix;
}

template<typename DataType>
EpetraMatrix<DataType>&  EpetraMatrix<DataType>::operator=(const EpetraMatrix& _matrix)
{
    M_map       = _matrix.M_map;
    M_domainMap = _matrix.M_domainMap;
    M_rangeMap  = _matrix.M_rangeMap;
    *M_epetraCrs = *(_matrix.M_epetraCrs);
    M_indexBase = _matrix.M_indexBase;

    return *this;
}

template<typename DataType>
typename EpetraMatrix<DataType>::vector_type
EpetraMatrix<DataType>::operator * (const vector_type& vec) const
{
    ASSERT_PRE(M_epetraCrs->Filled(),
               "EpetraMatrix::Operator*: globalAssemble(...) should be called first");
    ASSERT_PRE(vec.getMap().MapsAreSimilar(*M_domainMap),
               "EpetraMatrix::Operator*: the map of vec is not the same of domainMap");
    vector_type result(vec);

    M_epetraCrs->Apply(vec.getEpetraVector(), result.getEpetraVector());

    return result;
}


template<typename DataType>
EpetraMatrix<DataType>&
EpetraMatrix<DataType>::operator *= (const DataType val)
{
    M_epetraCrs->Scale(val);
    return *this;
}

template<typename DataType>
EpetraMatrix<DataType>
EpetraMatrix<DataType>::operator * (const DataType val) const
{
    EpetraMatrix<DataType> matr(*this);
    return matr *= val;
}

// ===================================================
// Methods
// ===================================================
//Method to open again a matrix
template <typename DataType>
void EpetraMatrix<DataType>::openCrsMatrix()
{
    if (M_epetraCrs->Filled())
    {
        Int meanNumEntries = this->getMeanNumEntries();
        matrix_ptrtype tmp(M_epetraCrs);
        M_epetraCrs.reset(new matrix_type(Copy,M_epetraCrs->RowMap(), meanNumEntries ));

#ifdef HAVE_TRILINOS_EPETRAEXT_31 // trilinos6
        EpetraExt::MatrixMatrix::Add(*tmp, false, 1., *M_epetraCrs, 1., false);
#elif defined HAVE_TRILINOS_EPETRAEXT // trilinos8
        EpetraExt::MatrixMatrix::Add(*tmp, false, 1., *M_epetraCrs, 1.);
#else
#error error: do not have nor EpetraExt 6 nor 7 or 8
#endif
        M_domainMap.reset();
        M_rangeMap.reset();

    }
}


//Method to remove all the zeros contain in the matrix
template <typename DataType>
void EpetraMatrix<DataType>::removeZeros()
{
    if (M_epetraCrs->Filled())
    {
        Int meanNumEntries = this->getMeanNumEntries();
        matrix_ptrtype tmp(M_epetraCrs);
        M_epetraCrs.reset(new matrix_type(Copy,M_epetraCrs->RowMap(), meanNumEntries ));

        //Variables to store the informations
        Int NumEntries;
        Real* Values;
        Int* Indices;
        Int row(0);

        for (Int i(0); i<tmp->NumGlobalRows(); ++i)
        {
            row = tmp->LRID(i+M_indexBase);
            tmp->ExtractMyRowView(row, NumEntries, Values, Indices);

            Int Indices2[NumEntries];
            Real Values2[NumEntries];
            Int NumEntries2(0);

            for (Int j(0); j<NumEntries; ++j)
            {
                if (Values[j] != 0.0)
                {
                    Indices2[NumEntries2] = tmp->GCID(Indices[j]);
                    Values2[NumEntries2] = Values[j];
                    NumEntries2++;
                }
            }
            M_epetraCrs->InsertGlobalValues(row,NumEntries2,Values2,Indices2);
        }
        insertZeroDiagonal();
        M_epetraCrs->GlobalAssemble();
    }
}

//Swap the matrix with a new one

template <typename DataType>
void EpetraMatrix<DataType>::swapCrsMatrix(matrix_ptrtype& p)
{
    M_epetraCrs.swap(p);
}


template <typename DataType>
void EpetraMatrix<DataType>::swapCrsMatrix(EpetraMatrix<DataType>& B)
{
    M_epetraCrs.swap(B.M_epetraCrs);
}

template <typename DataType>
Int EpetraMatrix<DataType>::Multiply(bool transposeA,
                                     const EpetraMatrix<DataType> &B, bool transposeB,
                                     EpetraMatrix<DataType> &C, bool call_FillComplete_on_result) const
{
    //return EpetraExt::MatrixMatrix::Multiply(*M_epetraCrs,transposeA,*B.getMatrixPtr(),transposeB,*C.getMatrixPtr(),call_FillComplete_on_result);
    Int errCode = EpetraExt::MatrixMatrix::Multiply(*M_epetraCrs,transposeA,*B.getMatrixPtr(),transposeB,*C.getMatrixPtr(),false);
    if (call_FillComplete_on_result)
        C.GlobalAssemble();

    return errCode;
}


template <typename DataType>
Int EpetraMatrix<DataType>::Multiply(bool transposeA, const vector_type& x, vector_type &y) const
{
    ASSERT_PRE(M_epetraCrs->Filled(),
               "EpetraMatrix<DataType>::Multiply: GlobalAssemble(...) must be called first");
    ASSERT_PRE(x.getMap().MapsAreSimilar(*M_domainMap),
               "EpetraMatrix<DataType>::Multiply: x map is different from M_domainMap");
    ASSERT_PRE(y.getMap().MapsAreSimilar(*M_rangeMap),
               "EpetraMatrix<DataType>::Multiply: y map is different from M_rangeMap");


    return M_epetraCrs->Multiply(transposeA,x.getEpetraVector(),y.getEpetraVector());
}

template <typename DataType>
void EpetraMatrix<DataType>::add (const DataType val, const EpetraMatrix& _matrix)
{
    //    EpetraMatrix matrix(Copy, _matrix.RowMap(), _matrix.GlobalMaxNumEntries());
#if defined HAVE_TRILINOS_EPETRAEXT // trilinos8
    EpetraExt::MatrixMatrix::Add(*_matrix.getMatrixPtr(), false, val, *this->getMatrixPtr(), 1.);
#else
#error error: do not have nor EpetraExt  8+
#endif
}

//! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
template <typename DataType>
void EpetraMatrix<DataType>::diagonalize ( std::vector<UInt> rVec,
                                           DataType const coeff,
                                           UInt offset)
{

    const Epetra_Comm&  Comm(M_epetraCrs->Comm());
    Int numProcs(Comm.NumProc());
    Int MyPID   (Comm.MyPID()   );
    Int i;


    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    Int*     r;
    UInt*    Ur;


    // loop on all proc
    for ( Int p(0); p < numProcs; p++)
    {
        Int sizeVec(rVec.size());

        Comm.Broadcast(&sizeVec, 1, p);

        if ( p == MyPID )
        {
            Ur    =  &rVec    .front();
        }
        else
        {
            Ur    = new UInt    [sizeVec];
        }

        r     = (Int*) Ur;

        Comm.Broadcast(r,    sizeVec, p);

        for (i=0; i < sizeVec; i++)
            diagonalize( Ur[i], coeff, offset);

        if ( p != MyPID )
        {
            delete[] Ur;
        }

    }

}

//! set entry (r,r) to coeff and rest of row r to zero. NB: index r starts from 0
template <typename DataType>
void EpetraMatrix<DataType>::diagonalize( UInt const r,
                                          DataType const coeff,
                                          UInt offset)
{

    if ( !M_epetraCrs->Filled() )
    { //! if not filled, I do not know how to diagonalize.
        ERROR_MSG( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap(M_epetraCrs->RowMap());
    const Epetra_Map& colMap(M_epetraCrs->ColMap());


    Int myCol = colMap.LID(r + M_indexBase + offset);

    // row: if r is mine, zero out values
    Int myRow = rowMap.LID(r + M_indexBase + offset);

    if (myRow >= 0)  // I have this row
    {
        Int    NumEntries;
        Real* Values;
        Int* Indices;
        // Int globCol;

        M_epetraCrs->ExtractMyRowView(myRow, NumEntries, Values, Indices);

        for (Int i(0); i <  NumEntries; i++)
        {
            Values[i] = 0;
        }

        DataType coeff_(coeff);
        M_epetraCrs->ReplaceMyValues(myRow, 1, &coeff_, &myCol); // A(r,r) = coeff
    }

}

/*! Diagonalization of rows r_0 to r_n of the system. Done by merging and distributing
  the vector of rows w.r.t. all processors and calling diagonalize(r_i, ...)
  *  @param r vector of rows to diagonalize
  *  @param vector of coeff values to set the diagonal entry A(r,r) to
  *  @param vector of b right hand sides vectors to be corrected
  *  @param vector of datum value to set the fix the solution entry x(r) at
  */
template <typename DataType>
void EpetraMatrix<DataType>::diagonalize( std::vector<UInt> rVec,
                                          DataType const coeff,
                                          vector_type &b,
                                          std::vector<DataType> datumVec,
                                          UInt offset)
{


    const Epetra_Comm&  Comm(M_epetraCrs->Comm());
    Int numProcs(Comm.NumProc());
    Int MyPID   (Comm.MyPID()   );
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

    const Epetra_Map& rowMap(M_epetraCrs->RowMap());


    //Comm.Barrier();
    // we want to know which IDs are our or not

    for (Int ii = 0; ii < (Int)rVec.size(); ++ii)
    {
        Int lID = rowMap.LID(rVec[ii] + 1);
        if (!(lID < 0))
        {

            localIDs.push_back(rVec[ii] + 1);
            localData.push_back(datumVec[ii]);
            localBC.insert(std::pair<Int, Real>(rVec[ii] + 1, datumVec[ii]));
        }
        else
        {
            remoteIDs.push_back(rVec[ii] + 1);
            remoteData.push_back(datumVec[ii]);
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
                         LIDList);

    std::vector< std::vector<Int> > procToID  (Comm.NumProc());
    std::vector< std::vector<Int> > procToData(Comm.NumProc());


    for (Int ii = 0; ii < numIDs; ++ii)
    {
        Int pi = PIDList[ii];
        procToID[pi].push_back(remoteIDs[ii]);
        procToData[pi].push_back(remoteData[ii]);
    }

    // then, we send all the nodes where they belong


    const Epetra_MpiComm* comm = dynamic_cast<Epetra_MpiComm const*>(&Comm);

    assert(comm != 0);

    for (Int ii = 0; ii < (Int)procToID.size(); ++ii)
    {
        if (ii != me)
        {
            Int length;
            length = procToID[ii].size();
            //                    std::cout << me << " is sending " << *length << " to " << ii << std::endl;
            MPI_Send( &length, 1, MPI_INT, ii, 666, comm->Comm() );
            if (length > 0)
            {
                MPI_Send( &procToID[ii][0], length, MPI_INT, ii, 667, comm->Comm() );
                MPI_Send( &procToData[ii][0], length, MPI_INT, ii, 668, comm->Comm() );
                //std::cout << me << " has sent to " << ii << " : ";

                //for (Int jj = 0; jj < procToData[ii].size(); ++jj)
                //std::cout << procToID[ii][jj] << " ";
                //std::cout << " end sent" << std::endl;
            }
        }

    }

    for (Int ii = 0; ii < (Int)procToID.size(); ++ii)
    {
        if (ii != me)
        {
            Int length;
            MPI_Status status;
            MPI_Recv( &length, 1, MPI_INT, ii, 666, comm->Comm(), &status );
            //std::cout << me << " received " << *length << " from " << ii << std::endl;


            if (length > 0)
            {
                Int* bufferID = new Int[length];
                Int* ptrID(0);//    = new Int[length];

                MPI_Recv( bufferID, length, MPI_INT, ii, 667, comm->Comm(), &status );

                //std::cout << me << " has received ";
                ptrID = bufferID;

                //                             for (Int ii = 0; ii < *length; ++ii, ++ptrID)
                //                                 {
                //                                     std::cout << *ptrID << " ";
                //                                     localIDs.push_back(*ptrID);
                //                                 }

                //std::cout << me << " has received ";

                //                             for (Int ii = 0; ii < *length; ++ii, ++ptr)
                //                                 {
                //                                     //std::cout << *ptr << " ";
                //                                     localIDs.push_back(*ptr);
                //                                 }
                //std::cout << std::endl;

                Real* bufferData = new Real[length];
                Real* ptrData(0);

                MPI_Recv( bufferData, length, MPI_INT, ii, 668, comm->Comm(), &status );

                //                             std::cout << me << " has received ";
                ptrData = bufferData;

                for (Int ii = 0; ii < length; ++ii, ++ptrID, ++ptrData)
                {
                    localBC.insert(std::pair<Int, Real>
                                   (*ptrID, *ptrData));

                    //std::cout << *ptrID << " <-> " << *ptrData << std::endl;
                    //                                    std::cout << *ptr << " ";
                    //localData.push_back(*ptr);
                }
                //std::cout << std::endl;

                delete[] bufferID;
                delete[] bufferData;

            }

        }
    }

    delete[] PIDList;
    delete[] LIDList;

    for (im = localBC.begin(); im != localBC.end(); ++im)
    {
        //std::cout << me << " is filling " << im->first << " with ";
        //std::cout << im->second << std::endl;
        diagonalize( im->first - 1, coeff, b, im->second, offset);
    }
    //    std::cout << std::endl;


    return;
#endif


    // loop on all proc
    for ( Int p(0); p < numProcs; p++)
    {
        Int sizeVec(rVec.size());
        if ( sizeVec != Int(datumVec.size()))
        { //! vectors must be of the same size
            ERROR_MSG( "diagonalize: vectors must be of the same size\n" );
        }

        Comm.Broadcast(&sizeVec, 1, p);

        if ( p == MyPID )
        {
            Ur    =  &rVec    .front();
            datum = &datumVec.front();
        }
        else
        {
            Ur    = new UInt    [sizeVec];
            datum = new DataType[sizeVec];
        }

        r     = (Int*) Ur;

        Comm.Broadcast(r,    sizeVec, p);
        Comm.Broadcast(datum,sizeVec, p);

        for (i=0; i < sizeVec; i++)
            diagonalize( Ur[i], coeff, b, datum[i], offset);

        if ( p != MyPID )
        {
            delete[] Ur;
            delete[] datum;
        }

    }

}

/*! Diagonalization of row r of the system. Done by setting A(r,r) = coeff,
 *  A(r,j) = 0 and A(j,r) = 0 for j!=r, and suitably correcting the right hand
 *  side of the system.
 *  @param r row to diagonalize
 *  @param coeff value to set the diagonal entry A(r,r) to
 *  @param b right hand side vector to be corrected
 *  @param datum value to set the fix the solution entry x(r) at
 */
template <typename DataType>
void EpetraMatrix<DataType>::diagonalize( UInt const r,
                                          DataType const coeff,
                                          vector_type &b,
                                          DataType datum,
                                          UInt offset)
{

    if ( !M_epetraCrs->Filled() )
    { //! if not filled, I do not know how to diagonalize.
        ERROR_MSG( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap(M_epetraCrs->RowMap());
    const Epetra_Map& colMap(M_epetraCrs->ColMap());


    Int myCol = colMap.LID(r + M_indexBase + offset);

#ifdef EPETRAMATRIX_SYMMETRIC_DIAGONALIZE
    if (myCol >= 0)  // I have this column
    {
        Real zero(0);
        for (Int i(0); i < rowMap.NumMyElements(); i++)
            // Note that if a value is not already present for the specified location in the matrix,
            // the input value will be ignored and a positive warning code will be returned.
            M_epetraCrs->ReplaceMyValues(i,1, &zero, &myCol);
        //            b[ globCol ] -= Values[i] * datum; //@@ correct rhs : this is false, to be corrected
    }
#endif

    // row: if r is mine, zero out values
    Int myRow = rowMap.LID(r + M_indexBase + offset);

    if (myRow >= 0)  // I have this row
    {
        Int    NumEntries;
        Real* Values;
        Int* Indices;

        M_epetraCrs->ExtractMyRowView(myRow, NumEntries, Values, Indices);

        for (Int i(0); i <  NumEntries; i++)
        {
            Values[i] = 0;
        }

        DataType coeff_(coeff);

        M_epetraCrs->ReplaceMyValues(myRow, 1, &coeff_, &myCol); // A(r,r) = coeff
        b[ r + M_indexBase + offset] = coeff * datum; // correct right hand side for row r // BASEINDEX + M_indexBase

    }

}

template <typename DataType>
void EpetraMatrix<DataType>::matrixMarket( std::string const &filename, const bool headers)
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename;
    std::string desc = "Created by LifeV";

    //Int  me    = M_epetraCrs->Comm().MyPID();

    //
    // check on the file name
    //

    // std::ostringstream myStream;
    // myStream << me;

    nome = filename + ".mtx";

    EpetraExt::RowMatrixToMatrixMarketFile( nome.c_str(),
                                            *M_epetraCrs,
                                            nome.c_str(),
                                            desc.c_str(),
                                            headers
                                          );
}

template <typename DataType>
void EpetraMatrix<DataType>::spy( std::string const &filename)
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";

    Int  me    = M_epetraCrs->Comm().MyPID();

    //
    // check on the file name
    //

    std::ostringstream myStream;
    myStream << me;
    nome = filename + ".m";

    EpetraExt::RowMatrixToMatlabFile( nome.c_str(), *M_epetraCrs);

}

template <typename DataType>
Int EpetraMatrix<DataType>::GlobalAssemble()
{
    if ( M_epetraCrs->Filled ())
    {
        //         if (M_epetraCrs->Comm().MyPID() == 0)
        //             std::cout << "Matrix is already filled" << std::endl;
        return -1;
    }

    insertZeroDiagonal();
    M_domainMap = M_map;
    M_rangeMap  = M_map;
    return  M_epetraCrs->GlobalAssemble();
}

template <typename DataType>
Int EpetraMatrix<DataType>::GlobalAssemble(const boost::shared_ptr<const EpetraMap> & domainMap,
                                           const boost::shared_ptr<const EpetraMap> & rangeMap)
{
    if ( M_epetraCrs->Filled ())
    {
        //         if (M_epetraCrs->Comm().MyPID() == 0)
        //             std::cout << "Matrix is already filled" << std::endl;
        return -1;
    }

    M_domainMap = domainMap;
    M_rangeMap  = rangeMap;
    return  M_epetraCrs->GlobalAssemble(*domainMap->getMap(Unique), *rangeMap->getMap(Unique));
}


//! insert the given value into the diagonal according to a specified EpetraMap
//! Pay intention that this will add values to the diagonal,
//! so for later added values with set_mat_inc, the value
//! will be added
template <typename DataType>
void
EpetraMatrix<DataType>::insertValueDiagonal( const DataType entry, const EpetraMap& Map, const UInt offset )
{
    for (UInt i=0 ; i<Map.getMap(Unique)->NumMyElements(); ++i)//num from 1
    {
        set_mat_inc(  offset + Map.getMap(Unique)->GID(i)-1 ,   offset + Map.getMap(Unique)->GID(i)-1, entry);
    }
}


//! insert the given value into the diagonal
//! Pay intention that this will add values to the diagonal,
//! so for later added values with set_mat_inc, the value
//! will be added
template <typename DataType>
void EpetraMatrix<DataType>::insertValueDiagonal(const DataType& value, Int from, Int to)
{
    if ( M_epetraCrs->Filled ())
    {
        if (M_epetraCrs->Comm().MyPID() == 0)
            std::cout << "Matrix is already filled, it is impossible to insert the diagonal now" << std::endl;
        return;
    }

    if (to == from) return;

    if (to < from) // do all entries
    {
        from = M_epetraCrs->RowMap().MinMyGID ();
        to = M_epetraCrs->RowMap().MaxMyGID () + 1;
    }

    Int* p =  M_epetraCrs->RowMap().MyGlobalElements();
    Int ierr;

    for (Int i(0); i <  M_epetraCrs->RowMap().NumMyElements(); ++i, ++p)
    {
        if (*p < from || *p >= to) continue;

        ierr = M_epetraCrs->InsertGlobalValues (1, p, 1, p, &value);

        if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
    }
}

//! insert ones into the diagonal to ensure the matrix' graph has a entry there
//! Pay intention that this will add ones to the diagonal,
//! so for later added values with set_mat_inc, the one
//! will be added
template <typename DataType>
void EpetraMatrix<DataType>::insertOneDiagonal(Int from, Int to)
{
    insertValueDiagonal(1.0, from, to);
}

// Adds zeros into the diagonal to ensure the matrix' graph has a entry there
// This method does not remove non zero entries in the diagonal.
template <typename DataType>
void EpetraMatrix<DataType>::insertZeroDiagonal( Int from, Int to)
{
    insertValueDiagonal(0.0, from, to);
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
set_mat( UInt row, UInt col, DataType loc_val )
{
    //     if (M_epetraCrs->RowMap()

    // incrementing row and cols by indexBase;
    Int irow(row + M_indexBase);
    Int icol(col + M_indexBase);

    //    if (M_epetraCrs->RowMap().MyGID (row))
    Int ierr=M_epetraCrs->ReplaceGlobalValues (1, &irow, 1, &icol, &loc_val);
    if (ierr!=0)
        { std::cout << " error in matrix replacement " << ierr << std::endl;}

}

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat_inc( UInt row, UInt col, DataType loc_val )
{

    // incrementing row and cols by indexBase;
    Int irow(row + M_indexBase);
    Int icol(col + M_indexBase);

    //    Int  me    = M_epetraCrs->Comm().MyPID();

    //       std::cout << " -> ";
    //       std::cout << irow << " " << icol << " " << loc_val << std::endl;
    Int ierr = M_epetraCrs->InsertGlobalValues (1, &irow, 1, &icol, &loc_val);

    if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
    //    std::cout << ierr << std::endl;
}

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat_inc( Int const numRows, Int const numCols,
             std::vector<Int> const row, std::vector<Int> const col,
             DataType* const* const loc_val,
             Int format)
{

    // incrementing row and cols by indexBase;
    std::vector<Int> irow(numRows);
    std::vector<Int> icol(numCols);

    std::vector<Int>::const_iterator pt;

    pt = row.begin();
    for (std::vector<Int>::iterator i(irow.begin()); i !=  irow.end() && pt != row.end(); ++i, ++pt)
        *i = *pt + M_indexBase;

    pt = col.begin();
    for (std::vector<Int>::iterator i(icol.begin()); i !=  icol.end() && pt != col.end(); ++i, ++pt)
        *i = *pt + M_indexBase;


    Int ierr = M_epetraCrs->InsertGlobalValues (numRows, &irow[0], numCols, &icol[0], loc_val, format);

    if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
    //    std::cout << ierr << std::endl;
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

} // end namespace LifeV
//@@
//#undef OFFSET

#endif
