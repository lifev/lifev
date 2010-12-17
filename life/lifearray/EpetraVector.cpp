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
    @brief EpetraVector

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 04-10-2006
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <EpetraExt_MultiVectorOut.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/life.hpp>
#include <life/lifearray/EpetraVector.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
EpetraVector::EpetraVector( const EpetraMapType& mapType ):
        M_epetraMap   (),
        M_mapType     ( mapType ),
        M_epetraVector(),
        M_combineMode ( Add )
{
}

EpetraVector::EpetraVector( const EpetraMap& map, const EpetraMapType& mapType ):
        M_epetraMap   ( new EpetraMap( map ) ),
        M_mapType     ( mapType ),
        M_epetraVector( new vector_type( *M_epetraMap->map(M_mapType) ) ),
        M_combineMode ( Add )
{
}

EpetraVector::EpetraVector( const boost::shared_ptr<EpetraMap>& map, const EpetraMapType& mapType ):
        M_epetraMap   ( map ),
        M_mapType     ( mapType ),
        M_epetraVector( new vector_type( *M_epetraMap->map(M_mapType) ) ),
        M_combineMode ( Add )
{
}

EpetraVector::EpetraVector( const EpetraVector& vector):
        M_epetraMap   ( vector.M_epetraMap ),
        M_mapType     ( vector.M_mapType ),
        M_epetraVector( new vector_type( vector.epetraVector() ) ), //This make a true copy!
        M_combineMode ( vector.M_combineMode )
{
}

EpetraVector::EpetraVector( const EpetraVector& vector, const EpetraMapType& mapType):
        M_epetraMap   ( vector.M_epetraMap ),
        M_mapType     ( mapType ),
        M_epetraVector( new vector_type( *M_epetraMap->map( M_mapType ) ) ),
        M_combineMode ( Add )
{
    operator = (vector);
}

EpetraVector::EpetraVector( const EpetraVector& vector, const EpetraMapType& mapType,
                            const Epetra_CombineMode& combineMode ):
        M_epetraMap   ( vector.M_epetraMap ),
        M_mapType     ( mapType ),
        M_epetraVector( new vector_type( *M_epetraMap->map( M_mapType ) ) ),
        M_combineMode ( Add )
{
    if (mapType == vector.M_mapType)
    {
        *M_epetraVector = vector.epetraVector();
        return;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    switch (M_mapType)
    {
    case Unique:
        M_epetraVector->Export( vector.epetraVector(), M_epetraMap->importer(), combineMode );
        return ;
    case Repeated:
        M_epetraVector->Import( vector.epetraVector(), M_epetraMap->exporter(), combineMode );
        return ;
    }
}

EpetraVector::EpetraVector( const Epetra_MultiVector&          vector,
                            const boost::shared_ptr<EpetraMap> map,
                            const EpetraMapType&               mapType ):
        M_epetraMap   ( map ),
        M_mapType     ( mapType ),
        M_epetraVector( new vector_type( *map->map( mapType ) ) ),
        M_combineMode ( Add )
{
    assert( this->blockMap().SameAs(vector.Map()) );
    M_epetraVector->Update(1., vector, 0.);
}

EpetraVector::EpetraVector( const EpetraVector& vector, const Int& reduceToProc):
        M_epetraMap   ( vector.M_epetraMap->createRootMap( reduceToProc ) ),
        M_mapType     ( Unique ),
        M_epetraVector( new vector_type( *M_epetraMap->map( M_mapType ) ) ),
        M_combineMode ( Add )
{
    operator = ( vector );
}


// ===================================================
// Operators
// ===================================================
EpetraVector::data_type&
EpetraVector::operator[]( const UInt row )
{
    Int lrow = blockMap().LID(row); // BASEINDEX + 1, row + 1

    // hint: with gdb: break LifeV::EpetraVector<Real>::operator[](unsigned Int)
    if ( lrow < 0 )
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::operator [] ERROR : !! lrow < 0\n" );
    }

    return (*M_epetraVector)[0][lrow];
}

const EpetraVector::data_type&
EpetraVector::operator[]( const UInt row ) const
{
    Int lrow = blockMap().LID(row); // BASEINDEX + 1 row+1

    if ( lrow < 0 )
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
    if ( &vector.epetraVector() == &this->epetraVector() )
        return *this;

    if ( blockMap().SameAs( vector.blockMap() ) )
    {
        *M_epetraVector = vector.epetraVector();
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    // vector have the same underlying EpetraMap, we then use the existing importer/exporter
    if ( M_epetraMap.get()  && vector.M_epetraMap.get() &&
            M_epetraMap->mapsAreSimilar( *vector.M_epetraMap ) )
    {
        // we shouldn't get here if we have the same maptype!
        assert( M_mapType != vector.M_mapType );

        switch ( M_mapType )
        {
        case Unique:
            M_epetraVector->Export(vector.epetraVector(), M_epetraMap->importer(), M_combineMode);
            return *this;
        case Repeated:
            M_epetraVector->Import(vector.epetraVector(), M_epetraMap->exporter(), M_combineMode);
            return *this;
        }
    }
    switch ( vector.M_mapType )
    {
    case Repeated:
        if ( M_mapType != Repeated )
            return Export(vector.epetraVector(), M_combineMode);
    case Unique:
        return Import(vector.epetraVector(), M_combineMode);
    }

    // if we get here, it means that we have two different repeated maps.
    // To hande this case, we have to create a unique copy first:

    std::cout << "Tentative of export import from two repeated vectors based on different maps."
              << std::endl;

    EpetraVector vectorUnique( *M_epetraMap, Unique );
    vectorUnique.Export( vector.epetraVector(), M_combineMode );
    M_epetraVector->Import( vectorUnique.epetraVector(), M_epetraMap->exporter(), M_combineMode );

    return *this;
}

// copies the value of a Epetra_MultiVector u (assumed of width 1). If the map is not the same,
// try to import the values. Calls Import with Add.
EpetraVector&
EpetraVector::operator=( const Epetra_MultiVector& vector )
{
    Epetra_FEVector const* feVec ( dynamic_cast<Epetra_FEVector const*>(&vector) );
    assert( feVec );

    // We hope we are guessing right
    switch ( M_mapType )
    {
    case Unique:
        return Export( *feVec, M_combineMode );
    case Repeated:
        return Import( *feVec, M_combineMode );
    }

    return *this;
}

EpetraVector&
EpetraVector::operator=( data_type scalar )
{
    M_epetraVector->PutScalar(scalar);

    return *this;
}

EpetraVector&
EpetraVector::operator+=( const EpetraVector& vector )
{
    if ( this->blockMap().SameAs( vector.blockMap() ) )
        M_epetraVector->Update( 1., vector.epetraVector(), 1. );
    else
    {
        EpetraVector vCopy( vector, M_mapType );
        M_epetraVector->Update( 1., vCopy.epetraVector(), 1. );
    }

    return *this;
}

EpetraVector&
EpetraVector::operator-=( const EpetraVector& vector )
{
    if ( this->blockMap().SameAs( vector.blockMap() ) )
        M_epetraVector->Update( -1., vector.epetraVector(), 1. );
    else
    {
        EpetraVector vCopy(vector, M_mapType);
        M_epetraVector->Update( -1., vCopy.epetraVector(), 1. );
    }

    return *this;
}

EpetraVector&
EpetraVector::operator*=( const EpetraVector& vector )
{
    if ( this->blockMap().SameAs( vector.blockMap() ) )
        M_epetraVector->Multiply( 1.0, vector.epetraVector(), *M_epetraVector, 0.0 );
    else
    {
        EpetraVector vectorCopy( vector, M_mapType );
        M_epetraVector->Multiply( 1.0, vectorCopy.epetraVector(), *M_epetraVector, 0.0 );
    }

    return *this;
}

EpetraVector&
EpetraVector::operator/=( const EpetraVector& vector )
{
    if ( this->blockMap().SameAs( vector.blockMap() ) )
        M_epetraVector->ReciprocalMultiply( 1.0, vector.epetraVector(), *M_epetraVector, 0.0 );
    else
    {
        EpetraVector vectorCopy( vector, M_mapType );
        M_epetraVector->ReciprocalMultiply( 1.0, vectorCopy.epetraVector(), *M_epetraVector, 0.0 );
    }

    return *this;
}

const EpetraVector
EpetraVector::operator+( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy += vector;

    return MyVectorCopy;
}

const EpetraVector
EpetraVector::operator-( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy -= vector;

    return MyVectorCopy;
}

const EpetraVector
EpetraVector::operator*( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy *= vector;

    return MyVectorCopy;
}

const EpetraVector
EpetraVector::operator/( const EpetraVector& vector ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy /= vector;

    return MyVectorCopy;
}

EpetraVector&
EpetraVector::operator+=( const data_type& scalar )
{
    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            (*M_epetraVector)[i][j] += scalar;

    return *this;
}

EpetraVector&
EpetraVector::operator-=( const data_type& scalar )
{
    this->operator+=( -scalar );

    return *this;
}

EpetraVector&
EpetraVector::operator*=( const data_type& scalar )
{
    M_epetraVector->Scale( scalar );

    return *this;
}

EpetraVector&
EpetraVector::operator/=( const data_type& scalar )
{
    this->operator*=( 1. / scalar );

    return *this;
}

const EpetraVector
EpetraVector::operator+( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy += scalar;

    return MyVectorCopy;
}

const EpetraVector
EpetraVector::operator-( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy -= scalar;

    return MyVectorCopy;
}

const EpetraVector
EpetraVector::operator*( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy *= scalar;

    return MyVectorCopy;
}

const EpetraVector
EpetraVector::operator/( const data_type& scalar ) const
{
    EpetraVector MyVectorCopy( *this );

    MyVectorCopy /= scalar;

    return MyVectorCopy;
}

EpetraVector
EpetraVector::operator==( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] == scalar ? true : false;

    return comparisonVector;
}

EpetraVector
EpetraVector::operator!=( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] != scalar ? true : false;

    return comparisonVector;
}

EpetraVector
EpetraVector::operator<( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] < scalar ? true : false;

    return comparisonVector;
}

EpetraVector
EpetraVector::operator>( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] > scalar ? true : false;

    return comparisonVector;
}

EpetraVector
EpetraVector::operator<=( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] <= scalar ? true : false;

    return comparisonVector;
}

EpetraVector
EpetraVector::operator>=( const Real& scalar )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] >= scalar ? true : false;

    return comparisonVector;
}

EpetraVector
EpetraVector::operator&&( const EpetraVector& vector )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] && vector.epetraVector()[i][j];

    return comparisonVector;
}

EpetraVector
EpetraVector::operator||( const EpetraVector& vector )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = (*M_epetraVector)[i][j] || vector.epetraVector()[i][j];

    return comparisonVector;
}

EpetraVector
EpetraVector::operator!( void )
{
    EpetraVector comparisonVector( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i=0; i < M_epetraVector->NumVectors(); ++i )
        for ( j=0; j < M_epetraVector->MyLength(); ++j )
            comparisonVector.epetraVector()[i][j] = !(*M_epetraVector)[i][j];

    return comparisonVector;
}


// ===================================================
// Methods
// ===================================================
Int EpetraVector::globalToLocalRowId( const UInt row ) const
{
    Int lrow = blockMap().LID(row); // BASEINDEX + 1, row + 1

    if ( lrow < 0 && blockMap().Comm().NumProc() == 1 )
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG( "EpetraVector::globalToLocalRowId ERROR : !! lrow < 0\n" );
    }

    return lrow;
}

bool EpetraVector::setCoefficient( const UInt row, const data_type& value, UInt offset )
{
    Int lrow = globalToLocalRowId(row + offset);
    if ( lrow < 0 )
        return false;

    (*M_epetraVector)[0][lrow] = value;
    return true;
}

Int EpetraVector::setCoefficients( std::vector<Int>& rowsVector, std::vector<Real>& valuesVector )
{
    ASSERT( rowsVector.size() == valuesVector.size(), "Error: rowsVector and valuesVector should have the same size" );
    ASSERT( M_mapType == Unique, "Error: Vector must have a unique map" );

    // Coding this part by hands, in fact I do not trust the following line (Simone, June 2008)
    // return M_epetraVector->ReplaceGlobalValues(rowsVector.size(), &rowsVector.front(), &valuesVector.front());

    const Epetra_Comm&  Comm( M_epetraVector->Comm() );
    Int numProcs( Comm.NumProc() );
    Int MyPID   ( Comm.MyPID() );
    Int i;

    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    Int*       r;
    data_type* datum;

    // loop on all proc
    for ( Int p(0); p < numProcs; p++ )
    {
        Int sizeVec( static_cast<Int>( rowsVector.size() ) );
        if ( sizeVec != static_cast<Int>( valuesVector.size() ) )
        { //! vectors must be of the same size
            ERROR_MSG( "diagonalize: vectors must be of the same size\n" );
        }

        Comm.Broadcast(&sizeVec, 1, p);

        if ( p == MyPID )
        {
            r     =  &rowsVector.front();
            datum = &valuesVector.front();
        }
        else
        {
            r     = new Int      [sizeVec];
            datum = new data_type[sizeVec];
        }

        Comm.Broadcast(r,    sizeVec, p);
        Comm.Broadcast(datum,sizeVec, p);

        // row: if r is mine, Assign values
        for (i=0; i < sizeVec; ++i)
            setCoefficient(r[i], datum[i]);

        if ( p != MyPID )
        {
            delete[] r;
            delete[] datum;
        }
    }

    return 0;

}

Int
EpetraVector::sumIntoGlobalValues ( const Int GID, const Real value )
{
    return M_epetraVector->SumIntoGlobalValues( 1, &GID, &value );
}

EpetraVector&
EpetraVector::add( const EpetraVector& vector, const Int offset )
{

    Int numMyEntries = vector.M_epetraVector->MyLength ();
    const Int* gids  = vector.blockMap().MyGlobalElements();

    // eg: (u,p) += p or (u,p) += u
    for ( Int i = 0; i < numMyEntries; ++i )
    {
        (*this)[gids[i]+offset] += vector(gids[i]);
    }

    return *this;
}

EpetraVector&
EpetraVector::subset( const EpetraVector& vector,
                      const UInt          offset )
{
    return this->subset( vector, map(), offset, static_cast<UInt> (0) );
}

EpetraVector&
EpetraVector::subset( const EpetraVector& vector,
                      const EpetraMap&    map,
                      const UInt          offset1,
                      const UInt          offset2 )
{
    if ( M_mapType == Repeated && vector.mapType() == Unique )
    {
        return subset(EpetraVector( vector, Repeated), map, offset1, offset2 );
    }
    return subset( vector.epetraVector(), map, offset1, offset2 );
}

EpetraVector&
EpetraVector::subset( const Epetra_MultiVector& vector,
                      const EpetraMap&    map,
                      const UInt          offset1,
                      const UInt          offset2,
                      const UInt          column )
{
    const Int*    gids         = map.map(M_mapType)->MyGlobalElements();
    const UInt    numMyEntries = map.map(M_mapType)->NumMyElements();

    Int lid1 ;
    Int lid2 ;

    // eg:  p = (u,p) or u = (u,p)
    for ( UInt i = 0; i < numMyEntries; ++i )
    {
        lid1 = vector.Map().LID(gids[i]+offset1);
        lid2 = blockMap().LID(gids[i]+offset2);
        ASSERT( ( lid2 >= 0 ) && ( lid1 >= 0 ), "EpetraVector::subset ERROR : !! lid < 0\n" );
        (*M_epetraVector)[0][lid2] = vector[column][lid1];
    }

    return *this;
}

void
EpetraVector::meanValue( Real* result ) const
{
    M_epetraVector->MeanValue(result);
}

Real
EpetraVector::norm1() const
{
    Real result;
    M_epetraVector->Norm1(&result);
    return result;
}

void
EpetraVector::norm1( Real* result ) const
{
    M_epetraVector->Norm1(result);
}

void
EpetraVector::norm1( Real& result ) const
{
    M_epetraVector->Norm1(&result);
}

Real
EpetraVector::norm2() const
{
    Real result;
    M_epetraVector->Norm2(&result);
    return result;
}

void
EpetraVector::norm2( Real* result ) const
{
    M_epetraVector->Norm2(result);
}

void
EpetraVector::norm2( Real& result ) const
{
    M_epetraVector->Norm2( &result );
}

Real
EpetraVector::normInf() const
{
    Real result;
    M_epetraVector->NormInf(&result);
    return result;
}

void
EpetraVector::normInf( Real* result ) const
{
    M_epetraVector->NormInf(result);
}

void
EpetraVector::normInf( Real& result ) const
{
    M_epetraVector->NormInf(&result);
}

Real
EpetraVector::minValue() const
{
    Real result;
    M_epetraVector->MinValue(&result);
    return result;
}

Real
EpetraVector::maxValue() const
{
    Real result;
    M_epetraVector->MaxValue(&result);
    return result;
}

void
EpetraVector::minValue( Real* result ) const
{
    M_epetraVector->MinValue(result);
}

void
EpetraVector::maxValue( Real* result ) const
{
    M_epetraVector->MaxValue(result);
}

void
EpetraVector::minValue( Real& result ) const
{
    M_epetraVector->MinValue(&result);
}

void
EpetraVector::maxValue( Real& result ) const
{
    M_epetraVector->MaxValue(&result);
}

void
EpetraVector::abs( void )
{
    M_epetraVector->Abs( *M_epetraVector );
}

void
EpetraVector::abs( EpetraVector& vector )
{
    vector.M_epetraVector->Abs( *M_epetraVector );
}

// Scalar Products
EpetraVector::data_type
EpetraVector::dot( const EpetraVector& vector ) const
{
    data_type scalarProduct;
    M_epetraVector->Dot( vector.epetraVector(), &scalarProduct );

    return scalarProduct;
}

void
EpetraVector::dot( const EpetraVector& vector, data_type& scalarProduct )
{
    M_epetraVector->Dot( vector.epetraVector(), &scalarProduct );
}

void EpetraVector::matrixMarket( std::string const &fileName, const bool headers ) const
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName;
    std::string desc = "Created by LifeV";

    if (M_mapType == Repeated)
    {
        EpetraVector unique(*this, Unique, Zero);
        unique.spy(fileName);
        return;
    }

    name = fileName + ".mtx";

    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(),
                                              *M_epetraVector,
                                              name.c_str(),
                                              desc.c_str(),
                                              headers );
}

void EpetraVector::spy( std::string const &fileName ) const
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName, uti = " , ";

    Int  me    = M_epetraVector->Comm().MyPID();

    if (M_mapType == Repeated)
    {
        EpetraVector unique( *this, Unique, Zero );
        unique.spy( fileName );
        return;
    }

    // check on the file name
    std::ostringstream myStream;
    myStream << me;
    name = fileName + ".m";

    EpetraExt::MultiVectorToMatlabFile( name.c_str(), *M_epetraVector );
}


void EpetraVector::showMe( std::ostream& output ) const
{
    EpetraVector redVec( *this, 0 ); // reduced vector (all at proc 0)

    if ( redVec.M_epetraVector->Comm().MyPID() )
        return; // do not need other CPUs now

    const Real* Values = redVec.epetraVector()[0];
    for ( Int i = 0; i < redVec.M_epetraVector->GlobalLength () ; ++i )
        output << Values[i] << std::endl;
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
EpetraVector::setDefaultCombineMode()
{
    setCombineMode(Add);
}

void
EpetraVector::setMap( const EpetraMap& map )
{
    M_epetraMap.reset( new EpetraMap( map ) );
    M_epetraVector.reset( new vector_type( *M_epetraMap->map(M_mapType) ) );
}

// ===================================================
// Get Methods
// ===================================================
Int
EpetraVector::size() const
{
    if ( M_epetraVector.get() )
        return M_epetraVector->GlobalLength();
    return 0;
}

// ===================================================
// Private Methods
// ===================================================
EpetraVector&
EpetraVector::Import (const Epetra_FEVector& vector, Epetra_CombineMode combineMode )
{
    if ( &vector == &this->epetraVector() )
        return *this;

    if ( blockMap().SameAs(vector.Map()) )
    {
        *M_epetraVector = vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Export reducedExport( blockMap(), vector.Map() );
    M_epetraVector->Import( vector, reducedExport, combineMode );

    return *this;
}

EpetraVector&
EpetraVector::Export ( const Epetra_FEVector& vector, Epetra_CombineMode combineMode )
{
    if ( &vector == &this->epetraVector() )
        return *this;

    if ( blockMap().SameAs(vector.Map()) )
    {
        *M_epetraVector = vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Import reducedImport( vector.Map(), blockMap());
    M_epetraVector->Export(vector, reducedImport, combineMode);

    return *this;
}

EpetraVector
operator-( const EpetraVector& vector )
{
    EpetraVector VectorCopy( vector );

    return VectorCopy *= static_cast<EpetraVector::data_type> ( -1.0 );
}

EpetraVector
operator+( const EpetraVector::data_type& scalar, const EpetraVector& vector )
{
    EpetraVector VectorCopy( vector );

    return VectorCopy += scalar;
}

EpetraVector
operator-( const EpetraVector::data_type& scalar, const EpetraVector& vector )
{
    EpetraVector VectorCopy( -vector );

    return VectorCopy += scalar;
}

EpetraVector
operator*( const EpetraVector::data_type& scalar, const EpetraVector& vector )
{
    EpetraVector VectorCopy( vector );

    return VectorCopy *= scalar;
}

}  // end namespace LifeV
