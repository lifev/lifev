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
   \file EpetraVector.hpp
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
             Simone Deparis <simone.deparis@epfl.ch>
   \date 2004-10-26
 */

#ifndef _EPETRAVECTOR_HPP_
#define _EPETRAVECTOR_HPP_

#include <life/lifealg/SolverAztec.hpp>

#include "Epetra_FEVector.h"
#include "life/lifecore/life.hpp"
//#include <life/lifearray/pattern.hpp>
#include "life/lifealg/EpetraMap.hpp"

#include <Epetra_Export.h>


//@@
//#define OFFSET 0

namespace LifeV
{
////////////////////////////////////////////////////////////////
//
//  Epetra Vector format Wrapper
//
///////////////////////////////////////////////////////////////





template <typename DataType = Real>
class EpetraVector
{
public:

    typedef Epetra_FEVector vector_type;

//    EpetraVector(); //!< default constructor : NULL pattern
    //
    // Note that the constructors MUST be based on an existing pattern
    //
//     EpetraVector( const CSRPatt &_pattern,
//                   Epetra_Comm&   comm );
//     EpetraVector( const CSRPatt&    _pattern,
//                   const Epetra_Map& _mapRow );//,
//                  const Epetra_Map& _mapCol );

    EpetraVector( const EpetraMap&          _map );
    EpetraVector( const Epetra_BlockMap&    _map );

    EpetraVector( const EpetraVector& _vector);

//! Copies _vector to a vector which resides only on the processor "reduceToProc"
    EpetraVector( const EpetraVector& _vector, const int reduceToProc);

 //! Copies _vector to a vector whith different map (usually the repeated map from EpetraMatrix)
    EpetraVector( const EpetraVector&     _vector,
                  const Epetra_BlockMap&  _map    );

    DataType& operator[] (const UInt row);
    const DataType& operator[] (const UInt row) const;
    DataType& operator() (const UInt row);
    const DataType& operator() (const UInt row) const;

    //! insert a global value. After insertion, you will have to call global assemble.
    int sumIntoGlobalValues (const int GID, const double value);

    const Epetra_BlockMap& Map() const;


    double Norm2()   const;
    double NormInf() const;

    void Norm2  (double* res) const;
    void NormInf(double* res) const;

    ~EpetraVector() {};
    vector_type& getEpetraVector()             {return M_epetraVector;}
    const vector_type& getEpetraVector() const {return M_epetraVector;}

    const int  size() const { return M_epetraVector.GlobalLength(); }

    void spy ( std::string const &filename );

    //! copies the value of a vector u. If the map is not the same,
    //! try to import the values. Calls Import with Add.
    EpetraVector& operator=(const EpetraVector& _vector);
    //! copies the value of a vector u. If the map is not the same,
    //! try to import the values. Let you decide wether to add or replace shared nodes:
    /*
        CombineMode Valuse:
        Add         Components on the receiving processor will be added together.
        Insert 	    Off-processor components will be inserted into locations on receiving processor replacing existing values.
        InsertAdd 	Off-processor components will be inserted into locations on receiving processor and added to existing values.
        Average 	Off-processor components will be averaged with existing components on the receiving processor.
        (+ Zero and AbsMax, probably never useful)
    */
    EpetraVector& Import   (const EpetraVector& _vector,
                            Epetra_CombineMode combineMode);

    /* adds to this the vector _vector with an offset.
       typically to do: (u,p) += p or (u,p) += u.
       Note: the nodes to add are taken by the map of _vector, hence:
       a) if this has a unique map: then _vector should also (otherwise run time error)
       b) if this has a repeated map: then _vector should also. (otherwise wrong)
     */
    EpetraVector&  add(const EpetraVector& _vector,
                       const int           offset = 0);


    /* set this to a subset of  _vector with an offset.
       typically to do: p = (u,p) or u = (u,p).
       Note: the nodes to add are taken by the map of this, hence:
       a) if _vector has a unique map: then this should also (otherwise run time error)
       b) if _vector has a repeated map: then this should also. (otherwise wrong)
     */
    EpetraVector<DataType>& subset(const EpetraVector& _vector,
                                   const int           offset = 0);

    //! if the map is not the same, try to import values
    EpetraVector& operator+=(const EpetraVector& _vector);
    EpetraVector& operator-=(const EpetraVector& _vector);

    EpetraVector& operator*=(DataType t);

    DataType operator*(EpetraVector const& a) const;


    int GlobalAssemble() { return  M_epetraVector.GlobalAssemble(); }

private:

    vector_type      M_epetraVector;
};

//-------------------------------------------------------------------------------------------------------
// CSR - VALUES
//------------------------------------------------------------------------------------------------------

template <typename DataType>
EpetraVector<DataType>::EpetraVector( const EpetraVector& _vector):
    M_epetraVector(_vector.M_epetraVector)
{
}

template <typename DataType>
EpetraVector<DataType>::EpetraVector( const EpetraMap& _map ):
    M_epetraVector( *_map.getUniqueEpetra_Map(), false)
{
}

template <typename DataType>
EpetraVector<DataType>::EpetraVector( const Epetra_BlockMap& _map ):
    M_epetraVector( _map, false)
{
}


// Copies _vector to a vector which resides only on the processor "reduceToProc"
template <typename DataType>
EpetraVector<DataType>::EpetraVector( const EpetraVector& _vector, const int reduceToProc):
    M_epetraVector( Epetra_BlockMap( _vector.Map().NumGlobalElements(),
                                     (_vector.Map().Comm().MyPID() == reduceToProc) * _vector.Map().NumGlobalElements(),
                                     _vector.Map().ElementSize(),
                                     _vector.Map().MinAllGID(),
                                     _vector.Map().Comm()  ),
                    false)
{
    operator = (_vector);

    /*
    Epetra_Export reducedExport(this->Map(), _vector.Map());
    M_epetraVector.Import(_vector.M_epetraVector, reducedExport, Add);
    */
}

// Copies _vector to a vector whitha different map (usually the repeated map from EpetraMatrix)
template <typename DataType>
EpetraVector<DataType>::EpetraVector( const EpetraVector&     _vector,
                                      const Epetra_BlockMap&  _map   ):
    M_epetraVector( _map, false )
{
    operator = (_vector);

    /*
    if (this->Map().SameAs(_vector.Map()) )
        M_epetraVector = _vector.M_epetraVector;
    else
        {
            Epetra_Export reducedExport(this->Map(), _vector.Map());
            M_epetraVector.Import(_vector.M_epetraVector, reducedExport, Add);
        }
    */
}


template <typename DataType>
DataType&
EpetraVector<DataType>::operator []( const UInt row )
{
    int lrow = Map().LID(row); // BASEINDEX + 1, row + 1

    // hint: with gdb: break LifeV::EpetraVector<double>::operator[](unsigned int)
    if (lrow < 0 )
    {
        std::cout << M_epetraVector.Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector<DataType>::operator [] ERROR : !! lrow < 0\n" );
    }
    return (M_epetraVector[0][lrow]);

}

template <typename DataType>
DataType&
EpetraVector<DataType>::operator ()( const UInt row )
{
    return operator[](row);
}


template <typename DataType>
const DataType&
EpetraVector<DataType>::operator ()( const UInt row ) const
{
    return operator[](row);
}


template <typename DataType>
const DataType&
EpetraVector<DataType>::operator [] ( const UInt row ) const
{
    int lrow = Map().LID(row); // BASEINDEX + 1 row+1

    if (lrow < 0 )
    {
        std::cout << M_epetraVector.Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector<DataType>::operator () ERROR : !! lrow < 0\n" );
    }

    return (M_epetraVector[0][lrow]);

}

template <typename DataType>
int
EpetraVector<DataType>::sumIntoGlobalValues (const int GID, const double value)
{
    return M_epetraVector.SumIntoGlobalValues(1, &GID, &value);
}


template <typename DataType>
const Epetra_BlockMap&
EpetraVector<DataType>::Map() const
{
    return M_epetraVector.Map();
}


template <typename DataType>
double
EpetraVector<DataType>::Norm2() const
{
    double res;
    M_epetraVector.Norm2(&res);
    return res;
}

template <typename DataType>
void
EpetraVector<DataType>::Norm2(double* res) const
{
    M_epetraVector.Norm2(res);
}

template <typename DataType>
double
EpetraVector<DataType>::NormInf() const
{
    double res;
    M_epetraVector.NormInf(&res);
    return res;
}

template <typename DataType>
void
EpetraVector<DataType>::NormInf(double* res) const
{
    M_epetraVector.NormInf(res);
}


template <typename DataType>
void EpetraVector<DataType>::spy( std::string const &filename )
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename;

    EpetraVector redVec(*this,0); // reduced vector (all at proc 0)

    int  me    = redVec.M_epetraVector.Comm().MyPID();

    if (me) return; // do not need other CPUs now

    //
    // check on the file name
    //

    std::ostringstream myStream;
    myStream << me;

    nome = filename + ".m";

    std::ofstream file_out( nome.c_str() );
    ASSERT( file_out, "Error: Output Vector (Values) file cannot be open" );

    file_out << "v = [";

    int           NumEntries = redVec.M_epetraVector.GlobalLength ();
    const double* Values     = redVec.M_epetraVector[0];
    const int*    Index      = redVec.Map().MyGlobalElements();

    for (int i = 0; i < NumEntries; i++)
    {
        file_out.width(20);
        file_out << Index[i] << " ";
        file_out << Values[i];
        file_out << "\n";
    }


    file_out << "];" << std::endl;

}



// copies the value of a vector u. If the map is not the same,
// try to import the values.
template <typename DataType>
EpetraVector<DataType>& EpetraVector<DataType>::
operator = (const EpetraVector& _vector)
{

    return Import(_vector, Add);

}

// copies the value of a vector u. If the map is not the same,
// try to import the values.
template <typename DataType>
EpetraVector<DataType>& EpetraVector<DataType>::
Import (const EpetraVector& _vector, Epetra_CombineMode combineMode)
{
    if (&_vector.getEpetraVector() == &this->getEpetraVector())
        return *this;

    if (this->Map().SameAs(_vector.Map()) )
        M_epetraVector = _vector.M_epetraVector;
    else
        {
            *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.
            Epetra_Export reducedExport(this->Map(), _vector.Map());
            M_epetraVector.Import(_vector.M_epetraVector, reducedExport, combineMode);
        }
    return *this;
}

template <typename DataType>
EpetraVector<DataType>& EpetraVector<DataType>::
add(const EpetraVector& _vector,
    const int           offset )
{
    if ( offset == 0 )
        return operator+= (_vector);

    int numMyEntries = _vector.M_epetraVector.MyLength ();
    const int*    gids       = _vector.M_epetraVector.Map().MyGlobalElements();

    // eg: (u,p) += p or (u,p) += u
    for (int i = 0; i < numMyEntries; ++i)
        {
            //        std::cout << gids[i] + offset << " " << gids[i] << std::endl;
            (*this)[gids[i]+offset] += _vector(gids[i]);
        }

    return *this;
}

template <typename DataType>
EpetraVector<DataType>& EpetraVector<DataType>::
subset(const EpetraVector& _vector,
       const int           offset )
{
    (*this) *= 0;

    int numMyEntries = M_epetraVector.MyLength ();
    const int*    gids       = M_epetraVector.Map().MyGlobalElements();

    // eg:  p = (u,p) or u = (u,p)
    for (int i = 0; i < numMyEntries; ++i)
        {
            //        std::cout << gids[i] + offset << " " << gids[i] << std::endl;
            (*this)[gids[i]] += _vector(gids[i]+offset);
        }

    return *this;
}


// if the map is not the same, try to import values
template <typename DataType>
EpetraVector<DataType>& EpetraVector<DataType>::
operator+=(const EpetraVector& _vector)
{
    if (this->Map().SameAs(_vector.Map()) )
        M_epetraVector.Update(1., _vector.M_epetraVector, 1.);
    else
        {
            EpetraVector<DataType> vCopy(_vector, this->Map());
            M_epetraVector.Update(1., vCopy.M_epetraVector, 1.);
        }

    return *this;
}

template <typename DataType>
EpetraVector<DataType>& EpetraVector<DataType>::
operator-=(const EpetraVector& _vector)
{
    if (this->Map().SameAs(_vector.Map()) )
        M_epetraVector.Update(-1., _vector.M_epetraVector, 1.);
    else
        {
            EpetraVector<DataType> vCopy(_vector, this->Map());
            M_epetraVector.Update(-1., vCopy.M_epetraVector, 1.);
        }

    return *this;
}

template <typename DataType>
EpetraVector<DataType>& EpetraVector<DataType>::
operator *= (DataType t)
{
    M_epetraVector.Scale(t);
    return *this;
}

template <typename DataType>
DataType EpetraVector<DataType>::
operator * ( EpetraVector<DataType> const& a) const
{
    DataType result;
    M_epetraVector.Dot(a.getEpetraVector(), &result);
    return result;
}


//! multiply scalar * vector.
template <typename DataType>
EpetraVector<DataType>
operator * (DataType t, const EpetraVector<DataType>& _vector)
{
    EpetraVector<DataType> result(_vector);
    return result *= t;
}


}
//@@
//#undef OFFSET

#endif
