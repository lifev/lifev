/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2008-03-11

  Copyright (C) 2008 EPFL

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
   \file EpetraVector.cpp
   \author Simone Deparis <simone.deparis@epfl.ch> Gilles Fourestey
   \date 2008-03-11
 */

#include <life/lifearray/EpetraVector.hpp>

namespace LifeV
{
////////////////////////////////////////////////////////////////
//
//  Epetra Vector format Wrapper
//
///////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------------------------------
// CSR - VALUES
//------------------------------------------------------------------------------------------------------

EpetraVector::EpetraVector( const EpetraVector& _vector):
    M_epetraVector(_vector.M_epetraVector),
    M_epetraMap   (_vector.M_epetraMap),
    M_maptype     (_vector.M_maptype)
{
}


EpetraVector::EpetraVector( const EpetraMap& _map, EpetraMapType maptype ):
    M_epetraVector(*_map.getMap(maptype)),
    M_epetraMap   (new EpetraMap(_map)),
    M_maptype     (maptype)
{
}

EpetraVector::EpetraVector( const EpetraVector& _vector, EpetraMapType maptype):
    M_epetraVector(*_vector.M_epetraMap->getMap(maptype)),
    M_epetraMap   (_vector.M_epetraMap),
    M_maptype     (maptype)
{
    operator = (_vector);
}

EpetraVector::EpetraVector( const EpetraVector& _vector, EpetraMapType maptype,
                            Epetra_CombineMode combineMode):
    M_epetraVector(*_vector.M_epetraMap->getMap(maptype)),
    M_epetraMap   (_vector.M_epetraMap),
    M_maptype     (maptype)
{

    if (maptype == _vector.M_maptype)
    {
        M_epetraVector = _vector.getEpetraVector();
        return;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    switch (M_maptype) {
    case Unique:
        M_epetraVector.Export(_vector.M_epetraVector, M_epetraMap->getImporter(), combineMode);
        return ;
    case Repeated:
        M_epetraVector.Import(_vector.M_epetraVector, M_epetraMap->getExporter(), combineMode);
        return ;
    }
}

/*
EpetraVector::EpetraVector( const Epetra_BlockMap& _map ):
    M_epetraVector( _map, false)
{
}
*/


// Copies _vector to FEvector that comes as Multivector
EpetraVector::EpetraVector( const Epetra_MultiVector&    _vector,
                            boost::shared_ptr<EpetraMap> _map,
                            EpetraMapType                maptype )
:
    M_epetraVector( *_map->getMap(maptype)),
    M_epetraMap   ( _map),
    M_maptype     (maptype)
{

    assert(this->BlockMap().SameAs(_vector.Map()) );

    M_epetraVector.Update(1., _vector, 0.);

}


// Copies _vector to a vector which resides only on the processor "reduceToProc"
EpetraVector::EpetraVector( const EpetraVector& _vector, const int reduceToProc):
    M_epetraVector(_vector.M_epetraMap->getRootMap(reduceToProc)),
    M_epetraMap   (),
    M_maptype     (Unique)
{
    operator = (_vector);
}

/*
// Copies _vector to a vector whitha different map (usually the repeated map from EpetraMatrix)
EpetraVector::EpetraVector( const EpetraVector&     _vector,
                            const Epetra_BlockMap&  _map   ):
    M_epetraVector( _map, false )
{
    operator = (_vector);
}
*/


EpetraVector::data_type&
EpetraVector::operator []( const UInt row )
{
    int lrow = BlockMap().LID(row); // BASEINDEX + 1, row + 1

    // hint: with gdb: break LifeV::EpetraVector<double>::operator[](unsigned int)
    if (lrow < 0 )
    {
        std::cout << M_epetraVector.Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::operator [] ERROR : !! lrow < 0\n" );
    }
    return (M_epetraVector[0][lrow]);

}

EpetraVector::data_type&
EpetraVector::operator ()( const UInt row )
{
    return operator[](row);
}


const EpetraVector::data_type&
EpetraVector::operator ()( const UInt row ) const
{
    return operator[](row);
}


const EpetraVector::data_type&
EpetraVector::operator [] ( const UInt row ) const
{
    int lrow = BlockMap().LID(row); // BASEINDEX + 1 row+1

    if (lrow < 0 )
    {
        std::cout << M_epetraVector.Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::operator () ERROR : !! lrow < 0\n" );
    }

    return (M_epetraVector[0][lrow]);

}

//! if row is mine returns the LID
//! if row is not mine and if the numCpus > 1, returns -1
//! if row is not mine and if the numCpus == 1, asserts
int EpetraVector::checkLID(const UInt row) const
{
    int lrow = BlockMap().LID(row); // BASEINDEX + 1, row + 1

    if (lrow < 0 && BlockMap().Comm().NumProc() == 1)
    {
        std::cout << M_epetraVector.Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::checkLID ERROR : !! lrow < 0\n" );
    }

    return lrow;

}

//! if row is mine sets this[row] = value and return true
//! if row is not mine and if the numCpus > 1, returns false
//! if row is not mine and if the numCpus == 1, asserts
bool EpetraVector::checkAndSet(const UInt row, const data_type& value, UInt offset)
{
    int lrow = checkLID(row + offset);
    if (lrow < 0)
        return false;

    M_epetraVector[0][lrow] = value;
    return true;
}

//! Set the row row of the vector to value. If it isn't on this processor,
//! store it and send it and send it at next GlobalAssemble
int EpetraVector::replaceGlobalValues(std::vector<int>& rVec, std::vector<double>& datumVec)
{
    ASSERT( rVec.size() == datumVec.size(), "Error: rVec and datumVec should have the same size" );
    ASSERT( M_maptype == Unique, "Error: Vector must have a unique map" );

    // Coding this part by hands, in fact I do not trust the following line (Simone, June 2008)
    // return M_epetraVector.ReplaceGlobalValues(rVec.size(), &rVec.front(), &datumVec.front());

    const Epetra_Comm&  Comm(M_epetraVector.Comm());
    int numProcs(Comm.NumProc());
    int MyPID   (Comm.MyPID()   );
    int i;

    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    int*       r;
    data_type* datum;


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
                    r    =  &rVec    .front();
                    datum = &datumVec.front();
                }
            else
                {
                    r    = new int     [sizeVec];
                    datum = new data_type[sizeVec];
                }

            Comm.Broadcast(r,    sizeVec, p);
            Comm.Broadcast(datum,sizeVec, p);

            // row: if r is mine, Assign values
            for (i=0; i < sizeVec; ++i)
                checkAndSet(r[i], datum[i]);

            if ( p != MyPID )
                {
                    delete[] r;
                    delete[] datum;
                }

        }

    return 0;

}

int
EpetraVector::sumIntoGlobalValues (const int GID, const double value)
{
    return M_epetraVector.SumIntoGlobalValues(1, &GID, &value);
}


const Epetra_BlockMap&
EpetraVector::BlockMap() const
{
    return M_epetraVector.Map();
}

void
EpetraVector::MeanValue(double* res) const
{
    M_epetraVector.MeanValue(res);
}

double
EpetraVector::Norm1() const
{
    double res;
    M_epetraVector.Norm1(&res);
    return res;
}

void
EpetraVector::Norm1(double* res) const
{
    M_epetraVector.Norm1(res);
}

double
EpetraVector::Norm2() const
{
    double res;
    M_epetraVector.Norm2(&res);
    return res;
}

void
EpetraVector::Norm2(double* res) const
{
    M_epetraVector.Norm2(res);
}

double
EpetraVector::NormInf() const
{
    double res;
    M_epetraVector.NormInf(&res);
    return res;
}

void
EpetraVector::NormInf(double* res) const
{
    M_epetraVector.NormInf(res);
}


void EpetraVector::spy( std::string const &filename ) const
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

    nome = "load" + filename + ".m";

    std::ofstream file_out( nome.c_str() );
    ASSERT( file_out, "Error: Output Vector (Values) file cannot be open" );

    file_out << filename << " =  [";

    int           NumEntries = redVec.M_epetraVector.GlobalLength ();
    const double* Values     = redVec.M_epetraVector[0];
    const int*    Index      = redVec.BlockMap().MyGlobalElements();

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
EpetraVector& EpetraVector::
operator = (const EpetraVector& _vector)
{
    if (&_vector.getEpetraVector() == &this->getEpetraVector())
        return *this;

    if (this->BlockMap().SameAs(_vector.BlockMap()) )
    {
        M_epetraVector = _vector.getEpetraVector();
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    // vector have the same underlying EpetraMap, we then use the existing importer/exporter
    if ( M_epetraMap.get()  && _vector.M_epetraMap.get() &&
         M_epetraMap->MapsAreSimilar( *_vector.M_epetraMap ) )
    {
        // we shouldn't get here if we have the same maptype!
        assert(M_maptype != _vector.M_maptype);

        switch (M_maptype) {
        case Unique:
            M_epetraVector.Export(_vector.M_epetraVector, M_epetraMap->getImporter(), Add);
            return *this;
        case Repeated:
             M_epetraVector.Import(_vector.M_epetraVector, M_epetraMap->getExporter(), Add);
            return *this;
        }
    }
    switch (_vector.M_maptype) {
    case Repeated:
        //
        if (M_maptype != Repeated)
            return Export(_vector.M_epetraVector, Add);
    case Unique:
            return Import(_vector.M_epetraVector, Add);
    }

    // if we get here, it means that we have two different repeated maps.
    // To hande this case, we have to create a unique copy first:

    std::cout << "Tentative of export import from two repeated vectors based on different maps."
              << std::endl;

    EpetraVector vectorUnique(*M_epetraMap, Unique);
    vectorUnique.Export(_vector.M_epetraVector, Add);
    M_epetraVector.Import(vectorUnique.M_epetraVector, M_epetraMap->getExporter(), Add);


    return *this;

}

//! copies the value of a Epetra_MultiVector u (assumed of width 1). If the map is not the same,
//! try to import the values. Calls Import with Add.
EpetraVector& EpetraVector::
operator = (const Epetra_MultiVector& _vector)
{
    Epetra_FEVector const* feVec (dynamic_cast<Epetra_FEVector const*>(&_vector));
    assert( feVec );

    // We hope we are guessing right
    switch (M_maptype) {
    case Unique:
        return Export(*feVec, Add);
    case Repeated:
        return Import(*feVec, Add);
    }

    return *this;

}


// copies the value of a vector u. If the map is not the same,
// try to import the values.
EpetraVector& EpetraVector::
Import (const Epetra_FEVector& _vector, Epetra_CombineMode combineMode)
{
    if (&_vector == &this->getEpetraVector())
        return *this;

    if (this->BlockMap().SameAs(_vector.Map()) )
    {
        M_epetraVector = _vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Export reducedExport(this->BlockMap(), _vector.Map());
    M_epetraVector.Import(_vector, reducedExport, combineMode);

    return *this;
}

// copies the value of a vector u. If the map is not the same,
// try to import the values.
EpetraVector& EpetraVector::
Export (const Epetra_FEVector& _vector, Epetra_CombineMode combineMode)
{
    if (&_vector == &this->getEpetraVector())
        return *this;

    if (this->BlockMap().SameAs(_vector.Map()) )
    {
        M_epetraVector = _vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Import reducedImport( _vector.Map(), this->BlockMap());
    M_epetraVector.Export(_vector, reducedImport, combineMode);

    return *this;
}

EpetraVector& EpetraVector::
add(const EpetraVector& _vector,
    const int           offset )
{
    if ( offset == 0 )
        return operator+= (_vector);

    int numMyEntries = _vector.M_epetraVector.MyLength ();
    const int*    gids       = _vector.BlockMap().MyGlobalElements();

    // eg: (u,p) += p or (u,p) += u
    for (int i = 0; i < numMyEntries; ++i)
        {
            //        std::cout << gids[i] + offset << " " << gids[i] << std::endl;
            (*this)[gids[i]+offset] += _vector(gids[i]);
        }

    return *this;
}

EpetraVector& EpetraVector::
subset(const EpetraVector& _vector,
       const int           offset )
{
    (*this) *= 0;

    int numMyEntries = M_epetraVector.MyLength ();
    const int*    gids       = BlockMap().MyGlobalElements();

    // eg:  p = (u,p) or u = (u,p)
    for (int i = 0; i < numMyEntries; ++i)
        {
            //        std::cout << gids[i] + offset << " " << gids[i] << std::endl;
            (*this)[gids[i]] += _vector(gids[i]+offset);
        }

    return *this;
}


// if the map is not the same, try to import values
EpetraVector& EpetraVector::
operator+=(const EpetraVector& _vector)
{
    if (this->BlockMap().SameAs(_vector.BlockMap()) )
        M_epetraVector.Update(1., _vector.M_epetraVector, 1.);
    else
        {
            EpetraVector vCopy(_vector, M_maptype);
            M_epetraVector.Update(1., vCopy.M_epetraVector, 1.);
        }

    return *this;
}

EpetraVector& EpetraVector::
operator-=(const EpetraVector& _vector)
{
    if (this->BlockMap().SameAs(_vector.BlockMap()) )
        M_epetraVector.Update(-1., _vector.M_epetraVector, 1.);
    else
        {
            EpetraVector vCopy(_vector, M_maptype);
            M_epetraVector.Update(-1., vCopy.M_epetraVector, 1.);
        }

    return *this;
}

EpetraVector& EpetraVector::
operator *= (data_type t)
{
    M_epetraVector.Scale(t);
    return *this;
}

EpetraVector& EpetraVector::
operator = (data_type t)
{
    M_epetraVector.PutScalar(t);
    return *this;
}

EpetraVector::data_type
EpetraVector::
operator * ( EpetraVector const& a) const
{
    data_type result;
    M_epetraVector.Dot(a.getEpetraVector(), &result);
    return result;
}


//! multiply scalar * vector.
EpetraVector
operator * (EpetraVector::data_type t, const EpetraVector& _vector)
{
    EpetraVector result(_vector);
    return result *= t;
}



}  // end namespace LifeV


