//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Model 1D
 *
 *  @version 1.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @author Simone Deparis <simone.deparis@epfl.ch>
 *  @date 04-10-2006
 *
 *  @version 1.10
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 15-10-2009
 */

#include <life/lifearray/EpetraVector.hpp>
#include <EpetraExt_MultiVectorOut.h>

namespace LifeV {

// ===================================================
// Constructors
// ===================================================
EpetraVector::EpetraVector( const EpetraMapType& maptype ):
    M_epetraMap   (),
    M_maptype     ( maptype ),
    M_epetraVector(),
    M_combineMode ( Add )
{
}

EpetraVector::EpetraVector( const EpetraMap& map, const EpetraMapType& maptype ):
    M_epetraMap   ( new EpetraMap( map ) ),
    M_maptype     ( maptype ),
    M_epetraVector( new vector_type( *M_epetraMap->getMap(M_maptype) )),
    M_combineMode ( Add )
{
}

EpetraVector::EpetraVector( const boost::shared_ptr<EpetraMap>& map, const EpetraMapType& maptype ):
    M_epetraMap   ( map ),
    M_maptype     ( maptype ),
    M_epetraVector( new vector_type( *M_epetraMap->getMap(M_maptype) ) ),
    M_combineMode ( Add )
{
}

EpetraVector::EpetraVector( const EpetraVector& vector):
    M_epetraMap   ( vector.M_epetraMap ),
    M_maptype     ( vector.M_maptype ),
    M_epetraVector( new vector_type( vector.getEpetraVector() ) ), //This make a true copy!
    M_combineMode ( vector.M_combineMode )
{
}

EpetraVector::EpetraVector( const EpetraVector& vector, const EpetraMapType& maptype):
    M_epetraMap   ( vector.M_epetraMap ),
    M_maptype     ( maptype ),
    M_epetraVector( new vector_type( *M_epetraMap->getMap( M_maptype ) ) ),
    M_combineMode ( Add )
{
    operator = (vector);
}

EpetraVector::EpetraVector( const EpetraVector& vector, const EpetraMapType& maptype,
                            const Epetra_CombineMode& combineMode ):
    M_epetraMap   ( vector.M_epetraMap ),
    M_maptype     ( maptype ),
    M_epetraVector( new vector_type( *M_epetraMap->getMap( M_maptype ) ) ),
    M_combineMode ( Add )
{
    if (maptype == vector.M_maptype)
    {
        *M_epetraVector = vector.getEpetraVector();
        return;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    switch (M_maptype)
    {
    case Unique:
        M_epetraVector->Export( vector.getEpetraVector(), M_epetraMap->getImporter(), combineMode );
        return ;
    case Repeated:
        M_epetraVector->Import( vector.getEpetraVector(), M_epetraMap->getExporter(), combineMode );
        return ;
    }
}

EpetraVector::EpetraVector( const Epetra_MultiVector&           vector,
                            const boost::shared_ptr<EpetraMap>  map,
                            const EpetraMapType&                maptype ):
    M_epetraMap   ( map ),
    M_maptype     ( maptype ),
    M_epetraVector( new vector_type( *M_epetraMap->getMap( M_maptype ) ) ),
    M_combineMode ( Add )
{
    assert( this->BlockMap().SameAs(vector.Map()) );
    M_epetraVector->Update(1., vector, 0.);
}

EpetraVector::EpetraVector( const EpetraVector& vector, const int& reduceToProc):
    M_epetraMap   ( vector.M_epetraMap->createRootMap( reduceToProc ) ),
    M_maptype     ( Unique ),
    M_epetraVector( new vector_type( *M_epetraMap->getMap( M_maptype ) ) ),
    M_combineMode ( Add )
{
    operator = ( vector );
}

// ===================================================
// Methods
// ===================================================
int EpetraVector::checkLID( const UInt row ) const
{
    int lrow = BlockMap().LID(row); // BASEINDEX + 1, row + 1

    if (lrow < 0 && BlockMap().Comm().NumProc() == 1)
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::checkLID ERROR : !! lrow < 0\n" );
    }

    return lrow;
}

bool EpetraVector::checkAndSet(const UInt row, const data_type& value, UInt offset)
{
    int lrow = checkLID(row + offset);
    if (lrow < 0)
        return false;

    (*M_epetraVector)[0][lrow] = value;
    return true;
}

int EpetraVector::replaceGlobalValues(std::vector<int>& rVec, std::vector<double>& datumVec)
{
    ASSERT( rVec.size() == datumVec.size(), "Error: rVec and datumVec should have the same size" );
    ASSERT( M_maptype == Unique, "Error: Vector must have a unique map" );

    // Coding this part by hands, in fact I do not trust the following line (Simone, June 2008)
    // return M_epetraVector->ReplaceGlobalValues(rVec.size(), &rVec.front(), &datumVec.front());

    const Epetra_Comm&  Comm(M_epetraVector->Comm());
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
        int sizeVec( static_cast<int>( rVec.size() ) );
        if ( sizeVec != static_cast<int>( datumVec.size() ) )
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
EpetraVector::sumIntoGlobalValues ( const int GID, const double value )
{
    return M_epetraVector->SumIntoGlobalValues(1, &GID, &value);
}

EpetraVector&
EpetraVector::add( const EpetraVector& vector, const int offset )
{
    if ( offset == 0 )
        return operator+= (vector);

    int numMyEntries = vector.M_epetraVector->MyLength ();
    const int*    gids       = vector.BlockMap().MyGlobalElements();

    // eg: (u,p) += p or (u,p) += u
    for (int i = 0; i < numMyEntries; ++i)
    {
        //        std::cout << gids[i] + offset << " " << gids[i] << std::endl;
        (*this)[gids[i]+offset] += vector(gids[i]);
    }

    return *this;
}

EpetraVector&
EpetraVector::subset( const EpetraVector& vector,
                      const UInt          offset )
{
    return this->subset( vector, getMap(), offset, static_cast<UInt> (0) );
}

EpetraVector&
EpetraVector::subset( const EpetraVector& vector,
                      const EpetraMap&    map,
                      const UInt          offset1,
                      const UInt          offset2 )
{
    if (M_maptype==Repeated && vector.getMaptype()==Unique )
    {
        return subset(EpetraVector(vector, Repeated), map, offset1, offset2);
    }
    return subset(vector.getEpetraVector(), map, offset1, offset2);
}

EpetraVector&
EpetraVector::subset( const Epetra_MultiVector& vector,
                      const EpetraMap&    map,
                      const UInt          offset1,
                      const UInt          offset2,
                      const UInt          column )
{
    const int*    gids        = map.getMap(M_maptype)->MyGlobalElements();
    const UInt    numMyEntries = map.getMap(M_maptype)->NumMyElements();

    int lid1 ;
    int lid2 ;

    // eg:  p = (u,p) or u = (u,p)
    for (UInt i = 0; i < numMyEntries; ++i)
    {
        lid1 = vector.Map().LID(gids[i]+offset1);
        lid2 = BlockMap().LID(gids[i]+offset2);
        ASSERT( ( lid2 >= 0 ) && ( lid1 >= 0 ), "EpetraVector::subset ERROR : !! lid < 0\n" );
        //        std::cout << gids[i] + offset << " " << gids[i] << std::endl;
        (*M_epetraVector)[0][lid2] = vector[column][lid1];
    }

    return *this;
}

void
EpetraVector::MeanValue(double* res) const
{
    M_epetraVector->MeanValue(res);
}

double
EpetraVector::Norm1() const
{
    double res;
    M_epetraVector->Norm1(&res);
    return res;
}

void
EpetraVector::Norm1(double* res) const
{
    M_epetraVector->Norm1(res);
}

void
EpetraVector::Norm1(double& res) const
{
    M_epetraVector->Norm1(&res);
}

double
EpetraVector::Norm2() const
{
    double res;
    M_epetraVector->Norm2(&res);
    return res;
}

void
EpetraVector::Norm2(double* res) const
{
    M_epetraVector->Norm2(res);
}

void
EpetraVector::Norm2( double& res ) const
{
    M_epetraVector->Norm2( &res );
}

double
EpetraVector::NormInf() const
{
    double res;
    M_epetraVector->NormInf(&res);
    return res;
}

void
EpetraVector::NormInf(double* res) const
{
    M_epetraVector->NormInf(res);
}

void
EpetraVector::NormInf(double& res) const
{
    M_epetraVector->NormInf(&res);
}

double
EpetraVector::MinValue() const
{
    double res;
    M_epetraVector->MinValue(&res);
    return res;
}

double
EpetraVector::MaxValue() const
{
    double res;
    M_epetraVector->MaxValue(&res);
    return res;
}

void
EpetraVector::MinValue(double* res) const
{
    M_epetraVector->MinValue(res);
}

void
EpetraVector::MaxValue(double* res) const
{
    M_epetraVector->MaxValue(res);
}

void
EpetraVector::MinValue(double& res) const
{
    M_epetraVector->MinValue(&res);
}

void
EpetraVector::MaxValue(double& res) const
{
    M_epetraVector->MaxValue(&res);
}

void
EpetraVector::Abs( void )
{
    M_epetraVector->Abs( *M_epetraVector );
}

void
EpetraVector::Abs( EpetraVector& vector )
{
    vector.M_epetraVector->Abs( *M_epetraVector );
}

// Scalar Products
EpetraVector::data_type
EpetraVector::Dot( const EpetraVector& vector ) const
{
    data_type scalarProduct;
    M_epetraVector->Dot( vector.getEpetraVector(), &scalarProduct );

    return scalarProduct;
}

// Scalar Products
void
EpetraVector::Dot( const EpetraVector& vector, data_type& scalarProduct )
{
    M_epetraVector->Dot( vector.getEpetraVector(), &scalarProduct );
}



void EpetraVector::matrixMarket( std::string const &filename, const bool headers ) const
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename;
    std::string desc = "Created by LifeV";

    // int  me    = M_epetraVector->Comm().MyPID();

    if (M_maptype == Repeated)
    {
        EpetraVector unique(*this, Unique, Zero);
        unique.spy(filename);
        return;
    }

    // check on the file name
    // std::ostringstream myStream;
    // myStream << me;

    nome = filename + ".mtx";

    EpetraExt::MultiVectorToMatrixMarketFile(nome.c_str(),
                                             *M_epetraVector,
                                             nome.c_str(),
                                             desc.c_str(),
                                             headers
                                             );
}


void EpetraVector::spy( std::string const &filename ) const
{
    // Purpose: Matlab dumping and spy
    std::string nome = filename, uti = " , ";

    int  me    = M_epetraVector->Comm().MyPID();

    if (M_maptype == Repeated)
    {
        EpetraVector unique(*this, Unique, Zero);
        unique.spy(filename);
        return;
    }

    // check on the file name
    std::ostringstream myStream;
    myStream << me;
    nome = filename + ".m";

    EpetraExt::MultiVectorToMatlabFile(nome.c_str(), *M_epetraVector);
}


void EpetraVector::ShowMe( std::ostream& output ) const
{
    EpetraVector redVec( *this, 0 ); // reduced vector (all at proc 0)

    if ( redVec.M_epetraVector->Comm().MyPID() )
        return; // do not need other CPUs now

    const double* Values = redVec.getEpetraVector()[0];
    for ( int i = 0; i < redVec.M_epetraVector->GlobalLength () ; ++i )
        output << Values[i] << std::endl;
}

// ===================================================
// Operators
// ===================================================
EpetraVector::data_type&
EpetraVector::operator[]( const UInt row )
{
    int lrow = BlockMap().LID(row); // BASEINDEX + 1, row + 1

    // hint: with gdb: break LifeV::EpetraVector<double>::operator[](unsigned int)
    if (lrow < 0 )
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::operator [] ERROR : !! lrow < 0\n" );
    }

    return (*M_epetraVector)[0][lrow];
}

const EpetraVector::data_type&
EpetraVector::operator[]( const UInt row ) const
{
    int lrow = BlockMap().LID(row); // BASEINDEX + 1 row+1

    if (lrow < 0 )
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::operator () ERROR : !! lrow < 0\n" );

    }

    return ((*M_epetraVector)[0][lrow]);
}

EpetraVector::data_type&
EpetraVector::operator()( const UInt row )
{
    return operator[](row);
}


const EpetraVector::data_type&
EpetraVector::operator()( const UInt row ) const
{
    return operator[](row);
}

// copies the value of a vector u. If the map is not the same,
// try to import the values.
EpetraVector&
EpetraVector::operator=( const EpetraVector& vector )
{
    if (&vector.getEpetraVector() == &this->getEpetraVector())
        return *this;

    if (BlockMap().SameAs(vector.BlockMap()) )
    {
        *M_epetraVector = vector.getEpetraVector();
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    // vector have the same underlying EpetraMap, we then use the existing importer/exporter
    if ( M_epetraMap.get()  && vector.M_epetraMap.get() &&
         M_epetraMap->MapsAreSimilar( *vector.M_epetraMap ) )
    {
        // we shouldn't get here if we have the same maptype!
        assert(M_maptype != vector.M_maptype);

        switch (M_maptype)
        {
        case Unique:
            M_epetraVector->Export(vector.getEpetraVector(), M_epetraMap->getImporter(), M_combineMode);
            return *this;
        case Repeated:
             M_epetraVector->Import(vector.getEpetraVector(), M_epetraMap->getExporter(), M_combineMode);
            return *this;
        }
    }
    switch (vector.M_maptype)
    {
    case Repeated:
        //
        if (M_maptype != Repeated)
            return Export(vector.getEpetraVector(), M_combineMode);
    case Unique:
            return Import(vector.getEpetraVector(), M_combineMode);
    }

    // if we get here, it means that we have two different repeated maps.
    // To hande this case, we have to create a unique copy first:

    std::cout << "Tentative of export import from two repeated vectors based on different maps."
              << std::endl;

    EpetraVector vectorUnique(*M_epetraMap, Unique);
    vectorUnique.Export(vector.getEpetraVector(), M_combineMode);
    M_epetraVector->Import(vectorUnique.getEpetraVector(), M_epetraMap->getExporter(), M_combineMode);

    return *this;
}

//! copies the value of a Epetra_MultiVector u (assumed of width 1). If the map is not the same,
//! try to import the values. Calls Import with Add.
EpetraVector&
EpetraVector::operator=( const Epetra_MultiVector& vector )
{
    Epetra_FEVector const* feVec (dynamic_cast<Epetra_FEVector const*>(&vector));
    assert( feVec );

    // We hope we are guessing right
    switch (M_maptype)
    {
    case Unique:
        return Export(*feVec, M_combineMode);
    case Repeated:
        return Import(*feVec, M_combineMode);
    }

    return *this;
}

EpetraVector&
EpetraVector::operator=( data_type t )
{
    M_epetraVector->PutScalar(t);

    return *this;
}

EpetraVector&
EpetraVector::operator+=( const EpetraVector& vector )
{
    if ( this->BlockMap().SameAs( vector.BlockMap() ) )
        M_epetraVector->Update( 1., vector.getEpetraVector(), 1. );
    else
    {
        EpetraVector vCopy( vector, M_maptype );
        M_epetraVector->Update( 1., vCopy.getEpetraVector(), 1. );
    }

    return *this;
}

EpetraVector&
EpetraVector::operator-=( const EpetraVector& vector )
{
    if ( this->BlockMap().SameAs( vector.BlockMap() ) )
        M_epetraVector->Update( -1., vector.getEpetraVector(), 1. );
    else
    {
        EpetraVector vCopy(vector, M_maptype);
        M_epetraVector->Update( -1., vCopy.getEpetraVector(), 1. );
    }

    return *this;
}

EpetraVector&
EpetraVector::operator*=( const EpetraVector& vector )
{
    if ( this->BlockMap().SameAs( vector.BlockMap() ) )
        M_epetraVector->Multiply( 1.0, vector.getEpetraVector(), *M_epetraVector, 0.0 );
    else
    {
        EpetraVector vectorCopy( vector, M_maptype );
        M_epetraVector->Multiply( 1.0, vectorCopy.getEpetraVector(), *M_epetraVector, 0.0 );
    }

    return *this;
    /*
  int numMyEntries;
  const int*  gids;

  if (M_maptype == Unique)
    {
      numMyEntries = this->M_epetraVector->MyLength();
      gids  = this->BlockMap().MyGlobalElements();

      if (!this->BlockMap().SameAs(vector.BlockMap()))
    {

      EpetraVector vCopy(vector, M_maptype, Insert);
      vector = vCopy;
    }
    }

  if (M_maptype == Repeated)
    {
      numMyEntries = vector.M_epetraVector->MyLength();
      gids   = vector.BlockMap().MyGlobalElements();

      if (!this->BlockMap().SameAs(vector.BlockMap()))
    {
      EpetraVector vCopy(vector, M_maptype, Zero);
      vector = vCopy;
    }
    }


    for (int i = 0; i < numMyEntries; ++i)
    {
      (*this)[gids[i]] *= vector(gids[i]);
    }


  return *this;
  */
}

// Element by element division
EpetraVector&
EpetraVector::operator/=( const EpetraVector& vector )
{
    if ( this->BlockMap().SameAs( vector.BlockMap() ) )
        M_epetraVector->ReciprocalMultiply( 1.0, vector.getEpetraVector(), *M_epetraVector, 0.0 );
    else
    {
        EpetraVector vectorCopy( vector, M_maptype );
        M_epetraVector->ReciprocalMultiply( 1.0, vectorCopy.getEpetraVector(), *M_epetraVector, 0.0 );
    }

    return *this;
}

// Element by element sum
const EpetraVector
EpetraVector::operator+( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy += vector;

    return MyVectorCopy;
}

// Element by element minus
const EpetraVector
EpetraVector::operator-( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy -= vector;

    return MyVectorCopy;
}

// Element by element multiplication
const EpetraVector
EpetraVector::operator*( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy *= vector;

    return MyVectorCopy;
}

// Element by element division
const EpetraVector
EpetraVector::operator/( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy /= vector;

    return MyVectorCopy;
}

// Add a scalar quantity
EpetraVector&
EpetraVector::operator+=( const data_type& scalar )
{
    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            (*M_epetraVector)[i][j] += scalar;

    return *this;
}

// Remove a scalar quantity
EpetraVector&
EpetraVector::operator-=( const data_type& scalar )
{
    this->operator+=( -scalar );

    return *this;
}

// Multiply by a scalar quantity
EpetraVector&
EpetraVector::operator*=( const data_type& scalar )
{
    M_epetraVector->Scale( scalar );

    return *this;
}

// Divide by a scalar quantity
EpetraVector&
EpetraVector::operator/=( const data_type& scalar )
{
    this->operator*=( 1. / scalar );

    return *this;
}

// Add a scalar quantity
const EpetraVector
EpetraVector::operator+( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy += scalar;

    return MyVectorCopy;
}

// Remove a scalar quantity
const EpetraVector
EpetraVector::operator-( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy -= scalar;

    return MyVectorCopy;
}

// Multiply by a scalar quantity
const EpetraVector
EpetraVector::operator*( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy *= scalar;

    return MyVectorCopy;
}

// Divide by a scalar quantity
const EpetraVector
EpetraVector::operator/( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy /= scalar;

    return MyVectorCopy;
}

//! - vector.
EpetraVector
operator-( const EpetraVector& vector )
{
    EpetraVector VectorCopy( vector );

    return VectorCopy *= static_cast<EpetraVector::data_type> ( -1.0 );
}

//! scalar + vector.
EpetraVector
operator+( const EpetraVector::data_type& scalar, const EpetraVector& vector )
{
    EpetraVector VectorCopy( vector );

    return VectorCopy += scalar;
}

//! scalar - vector.
EpetraVector
operator-( const EpetraVector::data_type& scalar, const EpetraVector& vector )
{
    EpetraVector VectorCopy( -vector );

    return VectorCopy += scalar;
}

//! scalar * vector.
EpetraVector
operator*( const EpetraVector::data_type& scalar, const EpetraVector& vector )
{
    EpetraVector VectorCopy( vector );

    return VectorCopy *= scalar;
}

// Comparison with a scalar value
EpetraVector
EpetraVector::operator==( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] == scalar ? true : false;

    return comparisonVector;
}

// Comparison with a scalar value
EpetraVector
EpetraVector::operator!=( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] != scalar ? true : false;

    return comparisonVector;
}

// Comparison with a scalar value
EpetraVector
EpetraVector::operator<( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] < scalar ? true : false;

    return comparisonVector;
}

// Comparison with a scalar value
EpetraVector
EpetraVector::operator>( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] > scalar ? true : false;

    return comparisonVector;
}

// Comparison with a scalar value
EpetraVector
EpetraVector::operator<=( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] <= scalar ? true : false;

    return comparisonVector;
}

// Comparison with a scalar value
EpetraVector
EpetraVector::operator>=( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] >= scalar ? true : false;

    return comparisonVector;
}

// Logical comparison
EpetraVector
EpetraVector::operator&&( const EpetraVector& vector )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] && vector.getEpetraVector()[i][j];

    return comparisonVector;
}

// Logical comparison
EpetraVector
EpetraVector::operator||( const EpetraVector& vector )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = (*M_epetraVector)[i][j] || vector.getEpetraVector()[i][j];

    return comparisonVector;
}

// Logical comparison
EpetraVector
EpetraVector::operator!( void )
{
    EpetraVector comparisonVector( *M_epetraMap, M_maptype );

    int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.getEpetraVector()[i][j] = !(*M_epetraVector)[i][j];

    return comparisonVector;
}


// ===================================================
// Set Methods
// ===================================================
void
EpetraVector::setCombineMode( Epetra_CombineMode combineMode )
{
    M_combineMode = combineMode;
}

void
EpetraVector::setDefaultCombineMode( )
{
    setCombineMode(Add);
}

void
EpetraVector::setMap( const EpetraMap& map )
{
    M_epetraMap.reset( new EpetraMap( map ) );
    M_epetraVector.reset( new vector_type( *M_epetraMap->getMap(M_maptype) ) );
}

// ===================================================
// Private Methods
// ===================================================
EpetraVector&
EpetraVector::Import (const Epetra_FEVector& vector, Epetra_CombineMode combineMode)
{
    if (&vector == &this->getEpetraVector())
        return *this;

    if (BlockMap().SameAs(vector.Map()) )
    {
        *M_epetraVector = vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Export reducedExport(BlockMap(), vector.Map());
    M_epetraVector->Import(vector, reducedExport, combineMode);

    return *this;
}

EpetraVector&
EpetraVector::Export (const Epetra_FEVector& vector, Epetra_CombineMode combineMode)
{
    if (&vector == &this->getEpetraVector())
        return *this;

    if (BlockMap().SameAs(vector.Map()) )
    {
        *M_epetraVector = vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Import reducedImport( vector.Map(), BlockMap());
    M_epetraVector->Export(vector, reducedImport, combineMode);

    return *this;
}

}  // end namespace LifeV
