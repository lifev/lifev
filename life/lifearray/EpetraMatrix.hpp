/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gilles Fourestey <gilles.fourestey@epfl.ch>
             Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-10-04

  Copyright (C) 2006 EPFL

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
             Simone Deparis <simone.deparis@epfl.ch>
   \date 2004-10-26
 */

#ifndef _EPETRAMATRIX_HPP_
#define _EPETRAMATRIX_HPP_


#include <lifeconfig.h>
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
////////////////////////////////////////////////////////////////
//
//  Epetra Matrix format Wrapper
//
///////////////////////////////////////////////////////////////

template <typename DataType>
class EpetraMatrix
{
public:

    typedef Epetra_FECrsMatrix     matrix_type;
    typedef EpetraVector vector_type;

//    EpetraMatrix(); //!< default constructor : NULL pattern
    //
    // Note that the constructors MUST be based on an existing pattern
    //
//     EpetraMatrix( const CSRPatt &_pattern,
//                   Epetra_Comm&   comm );
//     EpetraMatrix( const CSRPatt&    _pattern,
//                   const Epetra_Map& _mapRow );//,
//                  const Epetra_Map& _mapCol );

    EpetraMatrix( const EpetraMap&    _map, int numEntries = 50 );

    EpetraMatrix( const EpetraMatrix& _matrix);
//     EpetraMatrix( const EpetraMatrix<DataType> &_epetra );

//! Copies _matrix to a matrix which resides only on the processor "reduceToProc"
    EpetraMatrix( const EpetraMatrix& _matrix, const UInt reduceToProc);


    ~EpetraMatrix() {};

    matrix_type& getEpetraMatrix()       {return M_epetraCrs;}
    const matrix_type& getEpetraMatrix()  const {return M_epetraCrs;}

    EpetraMatrix& operator += (const EpetraMatrix& _matrix);

    vector_type     operator *  (const vector_type& vec) const;
    EpetraMatrix&   operator *= (const DataType     val);
    EpetraMatrix    operator *  (const DataType     val) const;

//    Epetra_FECrsMatrix const & getEpetraMatrix(){return M_epetraCrs;}

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

    void spy    ( std::string const &filename );

    // Calls insertZeroDiagonal and then epetra.globalAssemble;
    int GlobalAssemble();

    int MyPID() { return  M_epetraCrs.Comm().MyPID(); }

private:

    //! insert zeros into the diagonal to ensure the matrix' graph has a entry there
    //! This method does not remove non zero entries in the diagonal.
    void insertZeroDiagonal();

    matrix_type      M_epetraCrs;
};

//-------------------------------------------------------------------------------------------------------
// CSR - VALUES
//------------------------------------------------------------------------------------------------------

template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& _matrix):
    M_epetraCrs(_matrix.M_epetraCrs)
{
}


template <typename DataType>
EpetraMatrix<DataType>&
EpetraMatrix<DataType>::operator += ( const EpetraMatrix& _matrix)
{
//    EpetraMatrix matrix(Copy, _matrix.RowMap(), _matrix.GlobalMaxNumEntries());

#ifdef HAVE_TRILINOS_EPETRAEXT_31 // trilinos6
    EpetraExt::MatrixMatrix::Add(_matrix.getEpetraMatrix(), false, 1., (*this).getEpetraMatrix(), 1., false);
#elif defined HAVE_TRILINOS_EPETRAEXT // trilinos8
    EpetraExt::MatrixMatrix::Add(_matrix.getEpetraMatrix(), false, 1., (*this).getEpetraMatrix(), 1.);
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

    M_epetraCrs.Apply(vec.getEpetraVector(), result.getEpetraVector());

    return result;
}


template<typename DataType>
EpetraMatrix<DataType>&
EpetraMatrix<DataType>::operator *= (const DataType val)
{
    M_epetraCrs.Scale(val);
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
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMap& _map, int numEntries ):
    M_epetraCrs( Copy, *_map.getMap(Unique), numEntries, false)
//    M_epetraCrs( Copy, *_map.getEpetra_Map(), 0, false)
{
}

// template <typename DataType>
// EpetraMatrix<DataType>::
// EpetraMatrix( const EpetraMatrix<DataType> &_csr ):
//     M_epetraCrs( _csr.getEpetraMatrix() )
// {
// }


// Copies _matrix to a matrix which resides only on the processor "reduceToProc"
template <typename DataType>
EpetraMatrix<DataType>::EpetraMatrix( const EpetraMatrix& _matrix, const UInt reduceToProc):
    M_epetraCrs( Copy, Epetra_Map( _matrix.M_epetraCrs.Map().NumGlobalElements(),
                                   (_matrix.M_epetraCrs.Map().Comm().MyPID() == reduceToProc) * _matrix.M_epetraCrs.Map().NumGlobalElements(),
                                   _matrix.M_epetraCrs.Map().IndexBase(),
                                   _matrix.M_epetraCrs.Map().Comm()  ),
                 (_matrix.M_epetraCrs.Map().Comm().MyPID() == reduceToProc) * 20)
{
    int  me    = M_epetraCrs.Comm().MyPID();
    if (!me)
        std::cout << "_matrix.M_epetraCrs.Map().IndexBase() = "
                  << _matrix.M_epetraCrs.Map().IndexBase()
                  << std::endl
                  << "M_epetraCrs.Map().IndexBase() = "
                  << M_epetraCrs.Map().IndexBase()
                  << std::endl;
//     std::cout << "reduced Map : " << M_epetraCrs.Map()
//               <<    endl;

//     M_epetraCrs.Comm().Barrier();

//     std::cout << "original Map: " << _matrix.M_epetraCrs.Map()
//               << std::endl;

//     M_epetraCrs.Comm().Barrier();

    Epetra_Export reducedExport(M_epetraCrs.Map(), _matrix.M_epetraCrs.Map());
    M_epetraCrs.Import(_matrix.M_epetraCrs, reducedExport, Add);
}



template <typename DataType>
void EpetraMatrix<DataType>::
set_mat( UInt row, UInt col, DataType loc_val )
{
//     if (M_epetraCrs.RowMap()

    // incrementing row and cols by one;
    int irow(row + 1);
    int icol(col + 1);

//    if (M_epetraCrs.RowMap().MyGID (row))
    M_epetraCrs.ReplaceGlobalValues (1, &irow, 1, &icol, &loc_val);

}

template <typename DataType>
void EpetraMatrix<DataType>::
set_mat_inc( UInt row, UInt col, DataType loc_val )
{

    // incrementing row and cols by one;
    int irow(row + 1);
    int icol(col + 1);

//    int  me    = M_epetraCrs.Comm().MyPID();

//       std::cout << " -> ";
//       std::cout << irow << " " << icol << " " << loc_val << std::endl;
    int ierr = M_epetraCrs.InsertGlobalValues (1, &irow, 1, &icol, &loc_val);

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

    // incrementing row and cols by one;
    std::vector<int> irow(numRows);
    std::vector<int> icol(numCols);

    std::vector<int>::const_iterator pt;

    pt = row.begin();
    for (std::vector<int>::iterator i(irow.begin()); i !=  irow.end() && pt != row.end(); ++i, ++pt)
        *i = *pt + 1;

    pt = col.begin();
    for (std::vector<int>::iterator i(icol.begin()); i !=  icol.end() && pt != col.end(); ++i, ++pt)
        *i = *pt + 1;


    int ierr = M_epetraCrs.InsertGlobalValues (numRows, &irow[0], numCols, &icol[0], loc_val, format);

    if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
//    std::cout << ierr << std::endl;
}


template <typename DataType>
int EpetraMatrix<DataType>::GlobalAssemble()
{
    if ( M_epetraCrs.Filled ())
    {
//         if (M_epetraCrs.Comm().MyPID() == 0)
//             std::cout << "Matrix is already filled" << std::endl;
        return -1;
    }

    insertZeroDiagonal();
    return  M_epetraCrs.GlobalAssemble();
}

// Adds zeros into the diagonal to ensure the matrix' graph has a entry there
// This method does not remove non zero entries in the diagonal.
template <typename DataType>
void EpetraMatrix<DataType>::insertZeroDiagonal()
{

    if ( M_epetraCrs.Filled ())
    {
        if (M_epetraCrs.Comm().MyPID() == 0)
	    std::cout << "Matrix is already filled, it is impossible to insert the diagonal now" << std::endl;
        return;
    }

    int* p =  M_epetraCrs.RowMap().MyGlobalElements();
    int ierr;
    DataType const zero(0);

    for (int i(0); i <  M_epetraCrs.RowMap().NumMyElements(); ++i, ++p)
    {
        ierr = M_epetraCrs.InsertGlobalValues (1, p, 1, p, &zero);

        if (ierr < 0) std::cout << " error in matrix insertion " << ierr << std::endl;
    }

}


//! set entries (rVec(i),rVec(i)) to coeff and rest of row r(i) to zero
template <typename DataType>
void EpetraMatrix<DataType>::diagonalize ( std::vector<UInt> rVec,
                                           DataType const coeff,
                                           UInt offset)
{

    const Epetra_Comm&  Comm(M_epetraCrs.Comm());
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

    if ( !M_epetraCrs.Filled() )
    { //! if not filled, I do not know how to diagonalize.
      ERROR_MSG( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap(M_epetraCrs.RowMap());
    const Epetra_Map& colMap(M_epetraCrs.ColMap());


    int myCol = colMap.LID(r + 1 + offset);

    // row: if r is mine, zero out values
    int myRow = rowMap.LID(r + 1 + offset);

    if (myRow >= 0)  // I have this row
    {
        int    NumEntries;
        double* Values;
        int* Indices;
        // int globCol;

        M_epetraCrs.ExtractMyRowView(myRow, NumEntries, Values, Indices);

        for (int i(0); i <  NumEntries; i++)
        {
            Values[i] = 0;
        }

	DataType coeff_(coeff);
	M_epetraCrs.ReplaceMyValues(myRow, 1, &coeff_, &myCol); // A(r,r) = coeff


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

    const Epetra_Comm&  Comm(M_epetraCrs.Comm());
    int numProcs(Comm.NumProc());
    int MyPID   (Comm.MyPID()   );
    int i;

    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    int*     r;
    UInt*    Ur;
    DataType* datum;


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

    if ( !M_epetraCrs.Filled() )
    { //! if not filled, I do not know how to diagonalize.
      ERROR_MSG( "if not filled, I do not know how to diagonalize\n" );
    }

    const Epetra_Map& rowMap(M_epetraCrs.RowMap());
    const Epetra_Map& colMap(M_epetraCrs.ColMap());


    int myCol = colMap.LID(r + 1 + offset);

#ifdef EPETRAMATRIX_SYMMETRIC_DIAGONALIZE
    if (myCol >= 0)  // I have this column
    {
        double zero(0);
        for (int i(0); i < rowMap.NumMyElements(); i++)
            // Note that if a value is not already present for the specified location in the matrix,
            // the input value will be ignored and a positive warning code will be returned.
            M_epetraCrs.ReplaceMyValues(i,1, &zero, &myCol);
//            b[ globCol ] -= Values[i] * datum; //@@ correct rhs : this is false, to be corrected
    }
#endif

    // row: if r is mine, zero out values
    int myRow = rowMap.LID(r + 1 + offset);

    if (myRow >= 0)  // I have this row
    {
        int    NumEntries;
        double* Values;
        int* Indices;

        M_epetraCrs.ExtractMyRowView(myRow, NumEntries, Values, Indices);

        for (int i(0); i <  NumEntries; i++)
        {
            Values[i] = 0;
        }

        DataType coeff_(coeff);

        M_epetraCrs.ReplaceMyValues(myRow, 1, &coeff_, &myCol); // A(r,r) = coeff
        b[ r + 1 + offset] = coeff * datum; // correct right hand side for row r // BASEINDEX + 1

    }

}


template <typename DataType>
int EpetraMatrix<DataType>::getMeanNumEntries() const
{
    const int minEntries = M_epetraCrs.MaxNumEntries ()/2;
    if ( M_epetraCrs.NumMyRows() )
        return minEntries;

    int meanNumEntries = M_epetraCrs.NumMyNonzeros()/M_epetraCrs.NumMyRows();
    if ( meanNumEntries < minEntries || meanNumEntries > 2*minEntries )
        return minEntries;
    return meanNumEntries;
}

template <typename DataType>
void EpetraMatrix<DataType>::spy( std::string const &filename)
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";

     int  me    = M_epetraCrs.Comm().MyPID();

    //
    // check on the file name
    //

    std::ostringstream myStream;
    myStream << me;
    nome = filename + ".m";

    EpetraExt::RowMatrixToMatlabFile( nome.c_str(), M_epetraCrs);

}


}



//@@
//#undef OFFSET

#endif
