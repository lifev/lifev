/* -*- mode: c++ -*-

   This file is part of the LifeV library

   Author(s): Gilles Fourestey <gilles.fourestey@epfl.ch>
   Simone Deparis <simone.deparis@epfl.ch>
   Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   Date: 2006-10-04, 2009-11-03

   Copyright (C) 2006,2009 EPFL

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file EpetraMatrix.hpp
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
   \author Simone Deparis <simone.deparis@epfl.ch>
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date   2009-11-03

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

    typedef Epetra_FECrsMatrix             matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef EpetraVector                   vector_type;

    EpetraMatrix( const EpetraMap&    _map, int numEntries = 50, int indexBase = 1 );

    EpetraMatrix( const EpetraMatrix& _matrix);
    //     EpetraMatrix( const EpetraMatrix<DataType> &_epetra );

    //! Copies _matrix to a matrix which resides only on the processor "reduceToProc"
    EpetraMatrix( const EpetraMatrix& _matrix, const UInt reduceToProc);


    ~EpetraMatrix() {};

    //! Return the shared_pointer of the Epetra_FECrsMatrix
    matrix_ptrtype& getMatrixPtr(){return M_epetraCrs;}
    //! Return the const shared_pointer of the Epetra_FECrsMatrix
    const matrix_ptrtype& getMatrixPtr() const {return M_epetraCrs;}
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
    int Multiply(bool transposeA,
                 const EpetraMatrix<DataType> &B, bool transposeB,
                 EpetraMatrix<DataType> &C, bool call_FillComplete_on_result=true) const;
    /*! Multiply the first EpetraVector given as a parameter by the EpetraMatrix and put the result into the second given EpetraVector
     *  @param transposeA if true, it transposes the EpetraMatrix
     *  @param x vector that will be multiply by the EpetraMatrix
     *  @param y vector that will store the result
     */
    int Multiply(bool transposeA, const vector_type& x, vector_type &y) const;

    EpetraMatrix& operator += (const EpetraMatrix& _matrix);

    vector_type     operator *  (const vector_type& vec) const;
    EpetraMatrix&   operator *= (const DataType     val);
    EpetraMatrix    operator *  (const DataType     val) const;

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


    int getMeanNumEntries() const ;
    void set_mat( UInt row, UInt col, DataType loc_val );

    void set_mat_inc( int const numRows, int const numCols,
                      std::vector<int> const row, std::vector<int> const col,
                      DataType* const* const loc_val,
                      int format = Epetra_FECrsMatrix::COLUMN_MAJOR);

    void set_mat_inc( UInt row, UInt col, DataType loc_val );

    /*! Save the matrix into a Matlab (.m) file
     *  @param filename file where the matrix will be saved
     */
    void spy    ( std::string const &filename );

    // Calls insertZeroDiagonal and then epetra.globalAssemble;
    int GlobalAssemble();

    int MyPID() { return  M_epetraCrs->Comm().MyPID(); }

    //! Return the internal EpetraMap of the EpetraMatrix
    const EpetraMap& getMap() const
    {
        return *M_epetraMap;
    }

    //! insert ones into the diagonal to ensure the matrix' graph has a entry there
    //! Pay intention that this will add ones to the diagonal,
    //! so for later added values with set_mat_inc, the one
    //! will be added
    void insertOneDiagonal();

private:
    //! insert the given value into the diagonal
    //! Pay intention that this will add values to the diagonal,
    //! so for later added values with set_mat_inc, the one
    //! will be added
    void insertValueDiagonal(const DataType& value);

    //! insert zeros into the diagonal to ensure the matrix' graph has a entry there
    //! This method does not remove non zero entries in the diagonal.
    void insertZeroDiagonal();

    //! Shared pointer on an EpetraMap
    boost::shared_ptr< EpetraMap > M_epetraMap;

    //!Pointer on a Epetra_FECrsMatrix
    matrix_ptrtype  M_epetraCrs;

    //! The index base for the matrix
    int             M_indexBase;

};

//-------------------------------------------------------------------------------------------------------
// CSR - VALUES
//------------------------------------------------------------------------------------------------------

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& _matrix):
    M_epetraMap(_matrix.M_epetraMap),
    M_epetraCrs(new matrix_type(*_matrix.M_epetraCrs)),
    M_indexBase(_matrix.M_indexBase )
{
}

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMap& _map, int numEntries, int indexBase ):
    M_epetraMap   ( new EpetraMap  (_map)),
    M_epetraCrs   ( new matrix_type( Copy, *M_epetraMap->getMap(Unique), numEntries, false)),
    M_indexBase   ( indexBase )
{
}


// Copies _matrix to a matrix which resides only on the processor "reduceToProc"
template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& _matrix, const UInt reduceToProc):
    M_epetraMap   (_matrix.M_epetraMap->createRootMap(reduceToProc)),
    M_epetraCrs   ( new matrix_type( Copy, *M_epetraMap->getMap(Unique), numEntries
                                     (_matrix.M_epetraCrs->Map().Comm().MyPID() == reduceToProc) * 20,
                                     false) ),
    M_indexBase   (_matrix.M_indexBase)
{
    int  me    = M_epetraCrs->Comm().MyPID();
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
}

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
typename EpetraMatrix<DataType>::vector_type
EpetraMatrix<DataType>::operator * (const vector_type& vec) const
{
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

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat( UInt row, UInt col, DataType loc_val )
{
    //     if (M_epetraCrs->RowMap()

    // incrementing row and cols by indexBase;
    int irow(row + M_indexBase);
    int icol(col + M_indexBase);

    //    if (M_epetraCrs->RowMap().MyGID (row))
    int ierr=M_epetraCrs->ReplaceGlobalValues (1, &irow, 1, &icol, &loc_val);
    if(ierr!=0)
    { std::cout << " error in matrix replacement " << ierr << std::endl;}

}

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat_inc( UInt row, UInt col, DataType loc_val )
{

    // incrementing row and cols by indexBase;
    int irow(row + M_indexBase);
    int icol(col + M_indexBase);

    //    int  me    = M_epetraCrs->Comm().MyPID();

    //       std::cout << " -> ";
    //       std::cout << irow << " " << icol << " " << loc_val << std::endl;
    int ierr = M_epetraCrs->InsertGlobalValues (1, &irow, 1, &icol, &loc_val);

    if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
    //    std::cout << ierr << std::endl;
}

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat_inc( int const numRows, int const numCols,
             std::vector<int> const row, std::vector<int> const col,
             DataType* const* const loc_val,
             int format)
{

    // incrementing row and cols by indexBase;
    std::vector<int> irow(numRows);
    std::vector<int> icol(numCols);

    std::vector<int>::const_iterator pt;

    pt = row.begin();
    for (std::vector<int>::iterator i(irow.begin()); i !=  irow.end() && pt != row.end(); ++i, ++pt)
        *i = *pt + M_indexBase;

    pt = col.begin();
    for (std::vector<int>::iterator i(icol.begin()); i !=  icol.end() && pt != col.end(); ++i, ++pt)
        *i = *pt + M_indexBase;


    int ierr = M_epetraCrs->InsertGlobalValues (numRows, &irow[0], numCols, &icol[0], loc_val, format);

    if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
    //    std::cout << ierr << std::endl;
}

template <typename DataType>
int EpetraMatrix<DataType>::GlobalAssemble()
{
    if ( M_epetraCrs->Filled ())
    {
        //         if (M_epetraCrs->Comm().MyPID() == 0)
        //             std::cout << "Matrix is already filled" << std::endl;
        return -1;
    }

    insertZeroDiagonal();
    return  M_epetraCrs->GlobalAssemble();
}

//! insert the given value into the diagonal
//! Pay intention that this will add values to the diagonal,
//! so for later added values with set_mat_inc, the one
//! will be added
template <typename DataType>
void EpetraMatrix<DataType>::insertValueDiagonal(const DataType& value)
{
    if ( M_epetraCrs->Filled ())
    {
        if (M_epetraCrs->Comm().MyPID() == 0)
            std::cout << "Matrix is already filled, it is impossible to insert the diagonal now" << std::endl;
        return;
    }

    int* p =  M_epetraCrs->RowMap().MyGlobalElements();
    int ierr;

    for (int i(0); i <  M_epetraCrs->RowMap().NumMyElements(); ++i, ++p)
    {
        ierr = M_epetraCrs->InsertGlobalValues (1, p, 1, p, &value);

        if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
    }
}

//! insert ones into the diagonal to ensure the matrix' graph has a entry there
//! Pay intention that this will add ones to the diagonal,
//! so for later added values with set_mat_inc, the one
//! will be added
template <typename DataType>
void EpetraMatrix<DataType>::insertOneDiagonal()
{
    insertValueDiagonal(1.0);
}

// Adds zeros into the diagonal to ensure the matrix' graph has a entry there
// This method does not remove non zero entries in the diagonal.
template <typename DataType>
void EpetraMatrix<DataType>::insertZeroDiagonal()
{
    insertValueDiagonal(0.0);
}


//! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
template <typename DataType>
void EpetraMatrix<DataType>::diagonalize ( std::vector<UInt> rVec,
                                           DataType const coeff,
                                           UInt offset)
{

    const Epetra_Comm&  Comm(M_epetraCrs->Comm());
    int numProcs(Comm.NumProc());
    int MyPID   (Comm.MyPID()   );
    int i;


    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    int*     r;
    UInt*    Ur;


    // loop on all proc
    for ( int p(0); p < numProcs; p++)
    {
        int sizeVec(rVec.size());

        Comm.Broadcast(&sizeVec, 1, p);

        if ( p == MyPID )
        {
            Ur    =  &rVec    .front();
        }
        else
        {
            Ur    = new UInt    [sizeVec];
        }

        r     = (int*) Ur;

        Comm.Broadcast(r,    sizeVec, p);

        for (i=0; i < sizeVec; i++)
            diagonalize( Ur[i], coeff, offset);

        if ( p != MyPID )
        {
            delete[] Ur;
        }

    }

}

//! set entry (r,r) to coeff and rest of row r to zero
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


    int myCol = colMap.LID(r + M_indexBase + offset);

    // row: if r is mine, zero out values
    int myRow = rowMap.LID(r + M_indexBase + offset);

    if (myRow >= 0)  // I have this row
    {
        int    NumEntries;
        double* Values;
        int* Indices;
        // int globCol;

        M_epetraCrs->ExtractMyRowView(myRow, NumEntries, Values, Indices);

        for (int i(0); i <  NumEntries; i++)
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
    int numProcs(Comm.NumProc());
    int MyPID   (Comm.MyPID()   );
    int i;

    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    int*     r;
    UInt*    Ur;
    DataType* datum;

#if 1
    // 2 arrays to store le local and remote IDs

    int me = Comm.MyPID();

    std::vector<int> localIDs;
    std::vector<int> remoteIDs;

    std::vector<double> localData;
    std::vector<double> remoteData;

    std::map<int, double>           localBC;
    std::map<int, double>::iterator im;

    const Epetra_Map& rowMap(M_epetraCrs->RowMap());


    //Comm.Barrier();
    // we want to know which IDs are our or not

    for (int ii = 0; ii < (int)rVec.size(); ++ii)
    {
        int lID = rowMap.LID(rVec[ii] + 1);
        if (!(lID < 0))
        {

            localIDs.push_back(rVec[ii] + 1);
            localData.push_back(datumVec[ii]);
            localBC.insert(pair<int, double>(rVec[ii] + 1, datumVec[ii]));
        }
        else
        {
            remoteIDs.push_back(rVec[ii] + 1);
            remoteData.push_back(datumVec[ii]);
        }
    }

    // now, we have to fill our localIDs with IDs from other processors
    // first, we have to build the map of all the remoteIDs and their processor owner


    int numIDs = remoteIDs.size();

    int* PIDList = new int[numIDs];
    int* LIDList = new int[numIDs];

    rowMap.RemoteIDList( numIDs,
                         &remoteIDs[0],
                         PIDList,
                         LIDList);


    std::vector< std::vector<int> > procToID  (Comm.NumProc());
    std::vector< std::vector<int> > procToData(Comm.NumProc());


    for (int ii = 0; ii < numIDs; ++ii)
    {
        int pi = PIDList[ii];
        procToID[pi].push_back(remoteIDs[ii]);
        procToData[pi].push_back(remoteData[ii]);
    }

    // then, we send all the nodes where they belong


    const Epetra_MpiComm* comm = dynamic_cast<Epetra_MpiComm const*>(&Comm);

    assert(comm != 0);

    for (int ii = 0; ii < (int)procToID.size(); ++ii)
    {
        if (ii != me)
        {
            int length;
            length = procToID[ii].size();
            //                    std::cout << me << " is sending " << *length << " to " << ii << std::endl;
            MPI_Send( &length, 1, MPI_INT, ii, 666, comm->Comm() );
            if (length > 0)
            {
                MPI_Send( &procToID[ii][0], length, MPI_INT, ii, 667, comm->Comm() );
                MPI_Send( &procToData[ii][0], length, MPI_INT, ii, 668, comm->Comm() );
                //std::cout << me << " has sent to " << ii << " : ";

                //for (int jj = 0; jj < procToData[ii].size(); ++jj)
                //std::cout << procToID[ii][jj] << " ";
                //std::cout << " end sent" << std::endl;
            }
        }

    }

    for (int ii = 0; ii < (int)procToID.size(); ++ii)
    {
        if (ii != me)
        {
            int length;
            MPI_Status status;
            MPI_Recv( &length, 1, MPI_INT, ii, 666, comm->Comm(), &status );
            //std::cout << me << " received " << *length << " from " << ii << std::endl;


            if (length > 0)
            {
                int* bufferID = new int[length];
                int* ptrID(0);//    = new int[length];

                MPI_Recv( bufferID, length, MPI_INT, ii, 667, comm->Comm(), &status );

                //std::cout << me << " has received ";
                ptrID = bufferID;

                //                             for (int ii = 0; ii < *length; ++ii, ++ptrID)
                //                                 {
                //                                     std::cout << *ptrID << " ";
                //                                     localIDs.push_back(*ptrID);
                //                                 }

                //std::cout << me << " has received ";

                //                             for (int ii = 0; ii < *length; ++ii, ++ptr)
                //                                 {
                //                                     //std::cout << *ptr << " ";
                //                                     localIDs.push_back(*ptr);
                //                                 }
                //std::cout << std::endl;

                double* bufferData = new double[length];
                double* ptrData(0);

                MPI_Recv( bufferData, length, MPI_INT, ii, 668, comm->Comm(), &status );

                //                             std::cout << me << " has received ";
                ptrData = bufferData;

                for (int ii = 0; ii < length; ++ii, ++ptrID, ++ptrData)
                {
                    localBC.insert(pair<int, double>
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
    for ( int p(0); p < numProcs; p++)
    {
        int sizeVec(rVec.size());
        if ( sizeVec != int(datumVec.size()))
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

        r     = (int*) Ur;

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


    int myCol = colMap.LID(r + M_indexBase + offset);

#ifdef EPETRAMATRIX_SYMMETRIC_DIAGONALIZE
    if (myCol >= 0)  // I have this column
    {
        double zero(0);
        for (int i(0); i < rowMap.NumMyElements(); i++)
            // Note that if a value is not already present for the specified location in the matrix,
            // the input value will be ignored and a positive warning code will be returned.
            M_epetraCrs->ReplaceMyValues(i,1, &zero, &myCol);
        //            b[ globCol ] -= Values[i] * datum; //@@ correct rhs : this is false, to be corrected
    }
#endif

    // row: if r is mine, zero out values
    int myRow = rowMap.LID(r + M_indexBase + offset);

    if (myRow >= 0)  // I have this row
    {
        int    NumEntries;
        double* Values;
        int* Indices;

        M_epetraCrs->ExtractMyRowView(myRow, NumEntries, Values, Indices);

        for (int i(0); i <  NumEntries; i++)
        {
            Values[i] = 0;
        }

        DataType coeff_(coeff);

        M_epetraCrs->ReplaceMyValues(myRow, 1, &coeff_, &myCol); // A(r,r) = coeff
        b[ r + M_indexBase + offset] = coeff * datum; // correct right hand side for row r // BASEINDEX + M_indexBase

    }

}


template <typename DataType>
int EpetraMatrix<DataType>::getMeanNumEntries() const
{
    const int minEntries = M_epetraCrs->MaxNumEntries ()/2;
    if ( M_epetraCrs->NumMyRows() )
        return minEntries;

    int meanNumEntries = M_epetraCrs->NumMyNonzeros()/M_epetraCrs->NumMyRows();
    if ( meanNumEntries < minEntries || meanNumEntries > 2*minEntries )
        return minEntries;
    return meanNumEntries;
}

template <typename DataType>
void EpetraMatrix<DataType>::spy( std::string const &filename)
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";

    int  me    = M_epetraCrs->Comm().MyPID();

    //
    // check on the file name
    //

    std::ostringstream myStream;
    myStream << me;
    nome = filename + ".m";

    EpetraExt::RowMatrixToMatlabFile( nome.c_str(), *M_epetraCrs);

}

//Method to open again a matrix
template <typename DataType>
void EpetraMatrix<DataType>::openCrsMatrix()
{
    if(M_epetraCrs->Filled())
    {
        int meanNumEntries = this->getMeanNumEntries();
        matrix_ptrtype tmp(M_epetraCrs);
        M_epetraCrs.reset(new matrix_type(Copy,M_epetraCrs->RowMap(), meanNumEntries ));
        *M_epetraCrs += *tmp;
	}
}


//Method to remove all the zeros contain in the matrix
template <typename DataType>
void EpetraMatrix<DataType>::removeZeros()
{
    if(M_epetraCrs->Filled())
    {
        int meanNumEntries = this->getMeanNumEntries();
        matrix_ptrtype tmp(M_epetraCrs);
        M_epetraCrs.reset(new matrix_type(Copy,M_epetraCrs->RowMap(), meanNumEntries ));

        //Variables to store the informations
        int NumEntries;
        double* Values;
        int* Indices;
        int row(0);

        for(int i(0);i<tmp->NumGlobalRows();++i)
        {
            row = tmp->LRID(i+M_indexBase);
            tmp->ExtractMyRowView(row, NumEntries, Values, Indices);

            int Indices2[NumEntries];
            double Values2[NumEntries];
            int NumEntries2(0);

            for(int j(0);j<NumEntries;++j)
            {
                if(Values[j] != 0.0)
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
int EpetraMatrix<DataType>::Multiply(bool transposeA,
                                     const EpetraMatrix<DataType> &B, bool transposeB,
                                     EpetraMatrix<DataType> &C, bool call_FillComplete_on_result) const
{
    //return EpetraExt::MatrixMatrix::Multiply(*M_epetraCrs,transposeA,*B.getMatrixPtr(),transposeB,*C.getMatrixPtr(),call_FillComplete_on_result);
    int errCode = EpetraExt::MatrixMatrix::Multiply(*M_epetraCrs,transposeA,*B.getMatrixPtr(),transposeB,*C.getMatrixPtr(),false);
    if(call_FillComplete_on_result)
        C.GlobalAssemble();

    return errCode;
}


template <typename DataType>
int EpetraMatrix<DataType>::Multiply(bool transposeA, const vector_type& x, vector_type &y) const
{
    return M_epetraCrs->Multiply(transposeA,x.getEpetraVector(),y.getEpetraVector());
}





}
//@@
//#undef OFFSET

#endif
