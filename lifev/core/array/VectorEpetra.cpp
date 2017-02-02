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
    @brief VectorEpetra

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 04-10-2006
 */


#include <EpetraExt_MultiVectorOut.h>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
VectorEpetra::VectorEpetra ( const MapEpetraType& mapType, const combineMode_Type combineMode ) :
    M_epetraMap   (),
    M_mapType     ( mapType ),
    M_epetraVector(),
    M_combineMode ( combineMode )
{
}

VectorEpetra::VectorEpetra ( const MapEpetra& map,
                             const MapEpetraType& mapType,
                             const combineMode_Type combineMode ) :
    M_epetraMap   ( new MapEpetra ( map ) ),
    M_mapType     ( mapType ),
    M_combineMode ( combineMode )
{
    ASSERT (M_epetraMap->map(M_mapType).get()!=0, "Error! The stored MapEpetra does not have valid map pointer.\n");

    M_epetraVector.reset( new vector_type ( *M_epetraMap->map (M_mapType) ) );
}

VectorEpetra::VectorEpetra ( const std::shared_ptr<MapEpetra>& map,
                             const MapEpetraType& mapType,
                             const combineMode_Type combineMode ) :
    M_epetraMap   ( map ),
    M_mapType     ( mapType ),
    M_combineMode ( combineMode )
{
    ASSERT (M_epetraMap->map(M_mapType).get()!=0, "Error! The stored MapEpetra does not have valid map pointer.\n");

    M_epetraVector.reset( new vector_type ( *M_epetraMap->map (M_mapType) ) );
}

VectorEpetra::VectorEpetra ( const VectorEpetra& vector) :
    M_epetraMap   ( vector.M_epetraMap ),
    M_mapType     ( vector.M_mapType ),
    M_combineMode ( vector.M_combineMode )
{
    if (vector.epetraVectorPtr().get()!=0)
        M_epetraVector.reset( new vector_type ( vector.epetraVector() ) );   //This make a true copy!
}

VectorEpetra::VectorEpetra ( const VectorEpetra& vector, const MapEpetraType& mapType) :
    M_epetraMap   ( vector.M_epetraMap ),
    M_mapType     ( mapType ),
    M_combineMode ( vector.M_combineMode )
{
    ASSERT (M_epetraMap->map(M_mapType).get()!=0, "Error! The stored MapEpetra does not have valid map pointer.\n");

    M_epetraVector.reset( new vector_type ( *M_epetraMap->map ( M_mapType ) ) );

    operator = (vector);
}

VectorEpetra::VectorEpetra ( const VectorEpetra& vector, const MapEpetraType& mapType,
                             const combineMode_Type& combineMode ) :
    M_epetraMap   ( vector.M_epetraMap ),
    M_mapType     ( mapType ),
    M_combineMode ( vector.M_combineMode )
{
    ASSERT (M_epetraMap->map(M_mapType).get()!=0, "Error! The stored MapEpetra does not have valid map pointer.\n");

    M_epetraVector.reset( new vector_type ( *M_epetraMap->map ( M_mapType ) ) );

    if (mapType == vector.M_mapType)
    {
        *M_epetraVector = vector.epetraVector();
        return;
    }

    *this = 0.; // because of a buggy behaviour in case of multidefined indeces.

    switch (M_mapType)
    {
        case Unique:
            M_epetraVector->Export ( vector.epetraVector(), M_epetraMap->importer(), combineMode );
            return ;
        case Repeated:
            M_epetraVector->Import ( vector.epetraVector(), M_epetraMap->exporter(), combineMode );
            return ;
    }
}

VectorEpetra::VectorEpetra ( const Epetra_MultiVector&          vector,
                             const std::shared_ptr<MapEpetra> map,
                             const MapEpetraType&               mapType,
                             const combineMode_Type             combineMode ) :
    M_epetraMap   ( map ),
    M_mapType     ( mapType ),
    M_combineMode ( combineMode )
{
    ASSERT (M_epetraMap->map(M_mapType).get()!=0, "Error! The stored MapEpetra does not have valid map pointer.\n");

    assert ( M_epetraMap->map ( mapType )->SameAs (vector.Map() ) );

    M_epetraVector.reset( new vector_type ( *M_epetraMap->map ( mapType ) ) );
    M_epetraVector->Update (1., vector, 0.);
}

VectorEpetra::VectorEpetra ( const VectorEpetra& vector, const Int& reduceToProc) :
    M_mapType     ( Unique ),
    M_combineMode ( vector.M_combineMode )
{
    ASSERT (vector.map().map(M_mapType).get()!=0, "Error! The stored MapEpetra does not have valid map pointer.\n");

    M_epetraMap = vector.M_epetraMap->createRootMap ( reduceToProc );

    M_epetraVector.reset( new vector_type ( *M_epetraMap->map (M_mapType) ) );

    operator = ( vector );
}

// ===================================================
// Operators
// ===================================================
VectorEpetra::data_type&
VectorEpetra::operator[] ( const UInt row )
{
    Int lrow = blockMap().LID (static_cast<EpetraInt_Type> (row));

#ifdef HAVE_LIFEV_DEBUG
    if ( lrow < 0 )
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG ( "VectorEpetra::operator [] ERROR : !! lrow < 0\n" );
    }
#endif

    return (*M_epetraVector) [0][lrow];
}

const VectorEpetra::data_type&
VectorEpetra::operator[] ( const UInt row ) const
{
    Int lrow = blockMap().LID (static_cast<EpetraInt_Type> (row));

#ifdef HAVE_LIFEV_DEBUG
    if ( lrow < 0 )
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG ( "VectorEpetra::operator () ERROR : !! lrow < 0\n" );

    }
#endif

    return ( (*M_epetraVector) [0][lrow]);
}

VectorEpetra::data_type&
VectorEpetra::operator() ( const UInt row )
{
    return operator[] (static_cast<EpetraInt_Type> (row));
}


const VectorEpetra::data_type&
VectorEpetra::operator() ( const UInt row ) const
{
    return operator[] (static_cast<EpetraInt_Type> (row));
}

// copies the value of a vector u. If the map is not the same,
// try to import the values.
VectorEpetra&
VectorEpetra::operator= ( const VectorEpetra& vector )
{
    if ( &vector.epetraVector() == &this->epetraVector() )
    {
        return *this;
    }

    if ( blockMap().SameAs ( vector.blockMap() ) )
    {
        *M_epetraVector = vector.epetraVector();
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    // vector have the same underlying MapEpetra, we then use the existing importer/exporter
    if ( M_epetraMap.get()  && vector.M_epetraMap.get() &&
            M_epetraMap->mapsAreSimilar ( *vector.M_epetraMap ) )
    {
        // we shouldn't get here if we have the same maptype!
        assert ( M_mapType != vector.M_mapType );

        switch ( M_mapType )
        {
            case Unique:
                M_epetraVector->Export (vector.epetraVector(), M_epetraMap->importer(), M_combineMode);
                return *this;
            case Repeated:
                M_epetraVector->Import (vector.epetraVector(), M_epetraMap->exporter(), M_combineMode);
                return *this;
        }
    }
    switch ( vector.M_mapType )
    {
        case Repeated:
            if ( M_mapType != Repeated )
            {
                return Export (vector.epetraVector(), M_combineMode);
            }
        case Unique:
            return Import (vector.epetraVector(), M_combineMode);
    }

    // if we get here, it means that we have two different repeated maps.
    // To handle this case, we have to create a unique copy first:

    std::cout << "Tentative of export import from two repeated vectors based on different maps."
              << std::endl;

    VectorEpetra vectorUnique ( *M_epetraMap, Unique );
    vectorUnique.Export ( vector.epetraVector(), M_combineMode );
    M_epetraVector->Import ( vectorUnique.epetraVector(), M_epetraMap->exporter(), M_combineMode );

    return *this;
}

// copies the value of a Epetra_MultiVector u (assumed of width 1). If the map is not the same,
// try to import the values. Calls Import with Add.
VectorEpetra&
VectorEpetra::operator= ( const Epetra_MultiVector& vector )
{
    Epetra_FEVector const* feVec ( dynamic_cast<Epetra_FEVector const*> (&vector) );
    assert ( feVec );

    // We hope we are guessing right
    switch ( M_mapType )
    {
        case Unique:
            return Export ( *feVec, M_combineMode );
        case Repeated:
            return Import ( *feVec, M_combineMode );
    }

    return *this;
}

VectorEpetra&
VectorEpetra::operator= ( data_type scalar )
{
    M_epetraVector->PutScalar (scalar);

    return *this;
}

VectorEpetra&
VectorEpetra::operator+= ( const VectorEpetra& vector )
{
    if ( this->blockMap().SameAs ( vector.blockMap() ) )
    {
        M_epetraVector->Update ( 1., vector.epetraVector(), 1. );
    }
    else
    {
        VectorEpetra vCopy ( vector, M_mapType, M_combineMode );
        M_epetraVector->Update ( 1., vCopy.epetraVector(), 1. );
    }

    return *this;
}

VectorEpetra&
VectorEpetra::operator-= ( const VectorEpetra& vector )
{
    if ( this->blockMap().SameAs ( vector.blockMap() ) )
    {
        M_epetraVector->Update ( -1., vector.epetraVector(), 1. );
    }
    else
    {
        VectorEpetra vCopy (vector, M_mapType);
        M_epetraVector->Update ( -1., vCopy.epetraVector(), 1. );
    }

    return *this;
}

VectorEpetra&
VectorEpetra::operator*= ( const VectorEpetra& vector )
{
    if ( this->blockMap().SameAs ( vector.blockMap() ) )
    {
        M_epetraVector->Multiply ( 1.0, vector.epetraVector(), *M_epetraVector, 0.0 );
    }
    else
    {
        VectorEpetra vectorCopy ( vector, M_mapType );
        M_epetraVector->Multiply ( 1.0, vectorCopy.epetraVector(), *M_epetraVector, 0.0 );
    }

    return *this;
}

VectorEpetra&
VectorEpetra::operator/= ( const VectorEpetra& vector )
{
    if ( this->blockMap().SameAs ( vector.blockMap() ) )
    {
        M_epetraVector->ReciprocalMultiply ( 1.0, vector.epetraVector(), *M_epetraVector, 0.0 );
    }
    else
    {
        VectorEpetra vectorCopy ( vector, M_mapType );
        M_epetraVector->ReciprocalMultiply ( 1.0, vectorCopy.epetraVector(), *M_epetraVector, 0.0 );
    }

    return *this;
}

const VectorEpetra
VectorEpetra::operator+ ( const VectorEpetra& vector ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy += vector;

    return MyVectorCopy;
}

const VectorEpetra
VectorEpetra::operator- ( const VectorEpetra& vector ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy -= vector;

    return MyVectorCopy;
}

const VectorEpetra
VectorEpetra::operator* ( const VectorEpetra& vector ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy *= vector;

    return MyVectorCopy;
}

const VectorEpetra
VectorEpetra::operator/ ( const VectorEpetra& vector ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy /= vector;

    return MyVectorCopy;
}

VectorEpetra&
VectorEpetra::operator+= ( const data_type& scalar )
{
    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            (*M_epetraVector) [i][j] += scalar;
        }

    return *this;
}

VectorEpetra&
VectorEpetra::operator-= ( const data_type& scalar )
{
    this->operator+= ( -scalar );

    return *this;
}

VectorEpetra&
VectorEpetra::operator*= ( const data_type& scalar )
{
    M_epetraVector->Scale ( scalar );

    return *this;
}

VectorEpetra&
VectorEpetra::operator/= ( const data_type& scalar )
{
    this->operator*= ( 1. / scalar );

    return *this;
}

const VectorEpetra
VectorEpetra::operator+ ( const data_type& scalar ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy += scalar;

    return MyVectorCopy;
}

const VectorEpetra
VectorEpetra::operator- ( const data_type& scalar ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy -= scalar;

    return MyVectorCopy;
}

const VectorEpetra
VectorEpetra::operator* ( const data_type& scalar ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy *= scalar;

    return MyVectorCopy;
}

const VectorEpetra
VectorEpetra::operator/ ( const data_type& scalar ) const
{
    VectorEpetra MyVectorCopy ( *this );

    MyVectorCopy /= scalar;

    return MyVectorCopy;
}

VectorEpetra
VectorEpetra::operator== ( const Real& scalar ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] == scalar ? true : false;
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator!= ( const Real& scalar ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] != scalar ? true : false;
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator< ( const Real& scalar ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] < scalar ? true : false;
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator> ( const Real& scalar ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] > scalar ? true : false;
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator<= ( const Real& scalar ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] <= scalar ? true : false;
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator>= ( const Real& scalar ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] >= scalar ? true : false;
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator&& ( const VectorEpetra& vector ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] && vector.epetraVector() [i][j];
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator|| ( const VectorEpetra& vector ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = (*M_epetraVector) [i][j] || vector.epetraVector() [i][j];
        }

    return comparisonVector;
}

VectorEpetra
VectorEpetra::operator! ( void ) const
{
    VectorEpetra comparisonVector ( *M_epetraMap, M_mapType );

    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            comparisonVector.epetraVector() [i][j] = ! (*M_epetraVector) [i][j];
        }

    return comparisonVector;
}


// ===================================================
// Methods
// ===================================================
bool VectorEpetra::isGlobalIDPresent (const UInt row) const
{
    return blockMap().LID ( static_cast<EpetraInt_Type> (row) ) >= 0;
}

Int VectorEpetra::globalToLocalRowId ( const UInt row ) const
{
    Int lrow = blockMap().LID ( static_cast<EpetraInt_Type> (row) );

    if ( lrow < 0 && blockMap().Comm().NumProc() == 1 )
    {
        std::cout << M_epetraVector->Comm().MyPID() << " " << row << " " << lrow << std::endl;
        ERROR_MSG ( "VectorEpetra::globalToLocalRowId ERROR : !! lrow < 0\n" );
    }

    return lrow;
}

bool VectorEpetra::setCoefficient ( const UInt row, const data_type& value, UInt offset )
{
    Int lrow = globalToLocalRowId (row + offset);
    if ( lrow < 0 )
    {
        return false;
    }

    (*M_epetraVector) [0][lrow] = value;
    return true;
}

Int VectorEpetra::setCoefficients ( std::vector<Int>& rowsVector, std::vector<Real>& valuesVector )
{
    ASSERT ( rowsVector.size() == valuesVector.size(), "Error: rowsVector and valuesVector should have the same size" );
    ASSERT ( M_mapType == Unique, "Error: Vector must have a unique map" );

    // Coding this part by hands, in fact I do not trust the following line (Simone, June 2008)
    // return M_epetraVector->ReplaceGlobalValues(rowsVector.size(), &rowsVector.front(), &valuesVector.front());

    const Epetra_Comm&  Comm ( M_epetraVector->Comm() );
    Int numProcs ( Comm.NumProc() );
    Int MyPID   ( Comm.MyPID() );
    Int i;

    // Note: Epetra_Comm::broadcast does not support passing of uint, hence
    //       I define an int pointer to make the broadcast but then come back to an
    //       UInt pointer to insert the data
    Int*       r;
    data_type* datum;

    // loop on all proc
    for ( Int p (0); p < numProcs; p++ )
    {
        Int sizeVec ( static_cast<Int> ( rowsVector.size() ) );
        if ( sizeVec != static_cast<Int> ( valuesVector.size() ) )
        {
            //! vectors must be of the same size
            ERROR_MSG ( "diagonalize: vectors must be of the same size\n" );
        }

        Comm.Broadcast (&sizeVec, 1, p);

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

        Comm.Broadcast (r,    sizeVec, p);
        Comm.Broadcast (datum, sizeVec, p);

        // row: if r is mine, Assign values
        for (i = 0; i < sizeVec; ++i)
        {
            setCoefficient (r[i], datum[i]);
        }

        if ( p != MyPID )
        {
            delete[] r;
            delete[] datum;
        }
    }

    return 0;

}

Int
VectorEpetra::sumIntoGlobalValues ( const Int GID, const Real value )
{
    return M_epetraVector->SumIntoGlobalValues ( 1, &GID, &value );
}

VectorEpetra&
VectorEpetra::add ( const VectorEpetra& vector, const Int offset )
{

    Int numMyEntries = vector.M_epetraVector->MyLength ();
    const Int* gids  = vector.blockMap().MyGlobalElements();

    // eg: (u,p) += p or (u,p) += u
    for ( Int i = 0; i < numMyEntries; ++i )
    {
        (*this) [gids[i] + offset] += vector (gids[i]);
    }

    return *this;
}

VectorEpetra&
VectorEpetra::replace ( const VectorEpetra& vector, const Int& offset )
{
    // Definitions
    Int numMyEntries = vector.M_epetraVector->MyLength ();
    const Int* globalIDs = vector.blockMap().MyGlobalElements();

    // Replace part of the vector
    for ( Int i (0); i < numMyEntries; ++i )
    {
        ( *this ) [globalIDs[i] + offset] = vector ( globalIDs[i] );
    }

    return *this;
}

VectorEpetra&
VectorEpetra::subset ( const VectorEpetra& vector,
                       const UInt          offset )
{
    return this->subset ( vector, map(), offset, static_cast<UInt> (0) );
}

VectorEpetra&
VectorEpetra::subset ( const VectorEpetra& vector,
                       const MapEpetra&    map,
                       const UInt          offset1,
                       const UInt          offset2 )
{
    if ( M_mapType == Repeated && vector.mapType() == Unique )
    {
        return subset (VectorEpetra ( vector, Repeated), map, offset1, offset2 );
    }
    return subset ( vector.epetraVector(), map, offset1, offset2 );
}

VectorEpetra&
VectorEpetra::subset ( const Epetra_MultiVector& vector,
                       const MapEpetra&    map,
                       const UInt          offset1,
                       const UInt          offset2,
                       const UInt          column )
{
    const Int*    gids         = map.map (M_mapType)->MyGlobalElements();
    const UInt    numMyEntries = map.map (M_mapType)->NumMyElements();

    Int lid1 ;
    Int lid2 ;

    // eg:  p = (u,p) or u = (u,p)
    for ( UInt i = 0; i < numMyEntries; ++i )
    {
        lid1 = vector.Map().LID ( static_cast<EpetraInt_Type> (gids[i] + offset1) );
        lid2 = blockMap().LID ( static_cast<EpetraInt_Type> (gids[i] + offset2) );
        ASSERT ( ( lid2 >= 0 ) && ( lid1 >= 0 ), "VectorEpetra::subset ERROR : !! lid < 0\n" );
        (*M_epetraVector) [0][lid2] = vector[column][lid1];
    }

    return *this;
}

void
VectorEpetra::meanValue ( Real* result ) const
{
    M_epetraVector->MeanValue (result);
}

Real
VectorEpetra::norm1() const
{
    Real result;
    M_epetraVector->Norm1 (&result);
    return result;
}

void
VectorEpetra::norm1 ( Real* result ) const
{
    this->norm1 (*result);
}

void
VectorEpetra::norm1 ( Real& result ) const
{

    if (this->mapType() == Repeated)
    {
        VectorEpetra vUnique (*this, Unique, M_combineMode);
        vUnique.norm1 ( &result );
        return;
    }
    M_epetraVector->Norm1 (&result);

}

Real
VectorEpetra::norm2() const
{
    Real result;
    M_epetraVector->Norm2 (&result);
    return result;
}

void
VectorEpetra::norm2 ( Real* result ) const
{
    this->norm2 (*result);
}

void
VectorEpetra::norm2 ( Real& result ) const
{
    if (this->mapType() == Repeated)
    {
        VectorEpetra vUnique (*this, Unique, M_combineMode);
        vUnique.norm2 ( &result );
        return;
    }

    M_epetraVector->Norm2 ( &result );
}

Real
VectorEpetra::normInf() const
{
    Real result;
    M_epetraVector->NormInf (&result);
    return result;
}

void
VectorEpetra::normInf ( Real* result ) const
{
    M_epetraVector->NormInf (result);
}

void
VectorEpetra::normInf ( Real& result ) const
{
    M_epetraVector->NormInf (&result);
}

Real
VectorEpetra::minValue() const
{
    Real result;
    M_epetraVector->MinValue (&result);
    return result;
}

Real
VectorEpetra::maxValue() const
{
    Real result;
    M_epetraVector->MaxValue (&result);
    return result;
}

void
VectorEpetra::minValue ( Real* result ) const
{
    M_epetraVector->MinValue (result);
}

void
VectorEpetra::maxValue ( Real* result ) const
{
    M_epetraVector->MaxValue (result);
}

void
VectorEpetra::minValue ( Real& result ) const
{
    M_epetraVector->MinValue (&result);
}

void
VectorEpetra::maxValue ( Real& result ) const
{
    M_epetraVector->MaxValue (&result);
}

void
VectorEpetra::abs ( void )
{
    M_epetraVector->Abs ( *M_epetraVector );
}

void
VectorEpetra::abs ( VectorEpetra& vector )
{
    vector.M_epetraVector->Abs ( *M_epetraVector );
}

void
VectorEpetra::sqrt ()
{
    for ( Int i (0); i < M_epetraVector->NumVectors(); ++i )
    {
        for ( Int j (0); j < M_epetraVector->MyLength(); ++j )
        {
            (*M_epetraVector) [i][j] = std::sqrt ( (*M_epetraVector) [i][j] );
        }
    }
}

// Scalar Products
VectorEpetra::data_type
VectorEpetra::dot ( const VectorEpetra& vector ) const
{
    data_type scalarProduct;
    M_epetraVector->Dot ( vector.epetraVector(), &scalarProduct );

    return scalarProduct;
}

void
VectorEpetra::dot ( const VectorEpetra& vector, data_type& scalarProduct )
{
    M_epetraVector->Dot ( vector.epetraVector(), &scalarProduct );
}

void VectorEpetra::matrixMarket ( std::string const& fileName, const bool headers ) const
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName;
    std::string desc = "Created by LifeV";

    if (M_mapType == Repeated)
    {
        VectorEpetra unique (*this, Unique, Zero);
        unique.spy (fileName);
        return;
    }

    name = fileName + ".mtx";

    EpetraExt::MultiVectorToMatrixMarketFile ( name.c_str(),
                                               *M_epetraVector,
                                               name.c_str(),
                                               desc.c_str(),
                                               headers );
}

void VectorEpetra::spy ( std::string const& fileName ) const
{
    // Purpose: Matlab dumping and spy
    std::string name = fileName;

    Int  me    = M_epetraVector->Comm().MyPID();

    if (M_mapType == Repeated)
    {
        VectorEpetra unique ( *this, Unique, Zero );
        unique.spy ( fileName );
        return;
    }

    // check on the file name
    std::ostringstream myStream;
    myStream << me;
    name = fileName + ".m";

    EpetraExt::MultiVectorToMatlabFile ( name.c_str(), *M_epetraVector );
}


void VectorEpetra::showMe ( std::ostream& output ) const
{
    VectorEpetra redVec ( *this, 0 ); // reduced vector (all at proc 0)

    if ( redVec.M_epetraVector->Comm().MyPID() )
    {
        return;    // do not need other CPUs now
    }

    const Real* Values = redVec.epetraVector() [0];
    for ( Int i = 0; i < redVec.M_epetraVector->GlobalLength () ; ++i )
    {
        output << Values[i] << std::endl;
    }
}

void VectorEpetra::apply (const std::function<Real (Real)>& f)
{
    Int i, j;
    for ( i = 0; i < M_epetraVector->NumVectors(); ++i )
        for ( j = 0; j < M_epetraVector->MyLength(); ++j )
        {
            (*M_epetraVector) [i][j] = f ( (*M_epetraVector) [i][j]);
        }
}


// ===================================================
// Set Methods
// ===================================================
void
VectorEpetra::setCombineMode ( combineMode_Type combineMode )
{
    M_combineMode = combineMode;
}

void
VectorEpetra::setDefaultCombineMode()
{
    setCombineMode (Add);
}

void
VectorEpetra::setMap ( const MapEpetra& map )
{
    M_epetraMap.reset ( new MapEpetra ( map ) );
    M_epetraVector.reset ( new vector_type ( *M_epetraMap->map (M_mapType) ) );
}

// ===================================================
// Get Methods
// ===================================================
Int
VectorEpetra::size() const
{
    if ( M_epetraVector.get() )
    {
        return M_epetraVector->GlobalLength();
    }
    return 0;
}

// ===================================================
// Private Methods
// ===================================================
VectorEpetra&
VectorEpetra::Import (const Epetra_FEVector& vector, combineMode_Type combineMode )
{
    if ( &vector == &this->epetraVector() )
    {
        return *this;
    }

    if ( blockMap().SameAs (vector.Map() ) )
    {
        *M_epetraVector = vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Export reducedExport ( blockMap(), vector.Map() );
    M_epetraVector->Import ( vector, reducedExport, combineMode );

    return *this;
}

VectorEpetra&
VectorEpetra::Export ( const Epetra_FEVector& vector, combineMode_Type combineMode )
{
    if ( &vector == &this->epetraVector() )
    {
        return *this;
    }

    if ( blockMap().SameAs (vector.Map() ) )
    {
        *M_epetraVector = vector;
        return *this;
    }

    *this *= 0.; // because of a buggy behaviour in case of multidefined indeces.

    Epetra_Import reducedImport ( vector.Map(), blockMap() );
    M_epetraVector->Export (vector, reducedImport, combineMode);

    return *this;
}

VectorEpetra
operator- ( const VectorEpetra& vector )
{
    VectorEpetra VectorCopy ( vector );

    return VectorCopy *= static_cast<VectorEpetra::data_type> ( -1.0 );
}

VectorEpetra
operator+ ( const VectorEpetra::data_type& scalar, const VectorEpetra& vector )
{
    VectorEpetra VectorCopy ( vector );

    return VectorCopy += scalar;
}

VectorEpetra
operator- ( const VectorEpetra::data_type& scalar, const VectorEpetra& vector )
{
    VectorEpetra VectorCopy ( -vector );

    return VectorCopy += scalar;
}

VectorEpetra
operator* ( const VectorEpetra::data_type& scalar, const VectorEpetra& vector )
{
    VectorEpetra VectorCopy ( vector );

    return VectorCopy *= scalar;
}

}  // end namespace LifeV

