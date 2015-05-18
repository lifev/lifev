/*
 * MeshVolumeSubdivision.hpp
 *
 *  Created on: Apr 21, 2015
 *      Author: dalsanto
 *        Mail: niccolo.dalsanto@epfl.ch
 */

#ifndef _MeshVolumeSubdivision_HPP_
#define _MeshVolumeSubdivision_HPP_

namespace LifeV
{
//! MeshVolumeSubdivision.
/*!
    @author
    Constructs local arrays containing the volumes ID corresponding to the Physical Entities specified by different regionflags

    @warning
    it supposes the regionFlags are > 0

 */
template<typename MeshType>
class MeshVolumeSubdivision
{
public:
    //! @name Public Types
    //@{
    typedef MeshType                             mesh_Type;
    typedef boost::shared_ptr<MeshType>          meshPtr_Type;

    MeshVolumeSubdivision();
    MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm );
    MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, bool _verbose = false );
    MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, UInt _numSubregions = 1, bool _verbose = false );
    MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, meshPtr_Type _mesh,
                     UInt _numSubregions = 1, bool _verbose = false );
    MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, meshPtr_Type _mesh,
                           Epetra_IntSerialDenseVector _regionFlags, UInt _numSubregions = 1, bool _verbose = false );

    ~MeshVolumeSubdivision();

private:

    UInt countElementPerFlag();
    UInt assignRegionFlags( Epetra_IntSerialDenseVector _regionFlags );
    UInt allocatePerElements();
    UInt fillElementPerFlag();

public:

    void printFlags();
    void printNumElementPerFlag();
    void printElementPerFlag();
    UInt makeSubDivision();
    const UInt getNumElements( UInt flag ) const;
    const UInt * getSubmesh( UInt flag ) const;

private:

    bool M_verbose;
    boost::shared_ptr< Epetra_Comm > M_comm;

    UInt                             M_numSubregions;
    UInt **                          M_elements;
    Epetra_IntSerialDenseVector      M_regionFlags;
    Epetra_IntSerialDenseVector      M_numElementPerFlag;

    meshPtr_Type                     M_mesh;
    bool                             M_readRegionFlags;

};

template<typename MeshType>
MeshVolumeSubdivision<MeshType>::
MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm )
:
M_comm( _comm )
{}

template<typename MeshType>
MeshVolumeSubdivision<MeshType>::
MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, bool _verbose )
:
M_comm(_comm),
M_verbose( _verbose )
{}

template<typename MeshType>
MeshVolumeSubdivision<MeshType>::
MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, UInt _numSubregions, bool _verbose )
:
M_comm(_comm),
M_verbose( _verbose ),
M_numSubregions( _numSubregions )
{
    M_elements = new UInt*[M_numSubregions];

    for( UInt iRegion; iRegion < M_numSubregions; iRegion++ )
    {
        M_numElementPerFlag( iRegion ) = 0;
    }


}

template<typename MeshType>
MeshVolumeSubdivision<MeshType>::
MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, meshPtr_Type _mesh,
                                            UInt _numSubregions, bool _verbose )
:
M_comm(_comm),
M_verbose( _verbose ),
M_mesh( _mesh ),
M_numSubregions( _numSubregions ),
M_regionFlags( _numSubregions ),
M_numElementPerFlag( _numSubregions ),
M_readRegionFlags( false )
{
    M_elements = new UInt*[M_numSubregions];

    for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
    {
        M_numElementPerFlag( iRegion ) = 0;
    }

}


template<typename MeshType>
MeshVolumeSubdivision<MeshType>::
MeshVolumeSubdivision( boost::shared_ptr< Epetra_Comm > _comm, meshPtr_Type _mesh,
                                            Epetra_IntSerialDenseVector _regionFlags, UInt _numSubregions, bool _verbose )
:
M_comm(_comm),
M_verbose( _verbose ),
M_mesh( _mesh ),
M_numSubregions( _numSubregions ),
M_regionFlags( _regionFlags ),
M_numElementPerFlag( _numSubregions ),
M_readRegionFlags( true )
{
    M_elements = new UInt*[M_numSubregions];
}

template<typename MeshType>
MeshVolumeSubdivision<MeshType>::
~MeshVolumeSubdivision()
{
    for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
    {
        delete M_elements[iRegion];
    }
    delete M_elements;
}

template<typename MeshType>
UInt
MeshVolumeSubdivision<MeshType>::
assignRegionFlags( Epetra_IntSerialDenseVector _regionFlags )
{
    M_regionFlags = _regionFlags;
    M_readRegionFlags = true;
}




template<typename MeshType>
UInt
MeshVolumeSubdivision<MeshType>::
countElementPerFlag()
{

    if ( !M_readRegionFlags )
    {
        std::cout << "I HAVE NOT READ THE FLAGS" << std::endl;
        assert(false);
    }

    UInt nbElements ( M_mesh->numElements() );
    UInt oldMarkerID = 0;
    UInt numRegion = 0;

    std::cout << std::endl << "TOTAL NUMBER OF ELEMENTS proc " << M_comm->MyPID()  << ": " << nbElements << std::endl << std::endl;

    for (UInt iElement(0); iElement < nbElements; iElement++)
    {
        // Extracting the marker
        UInt markerID = M_mesh->element( iElement ).markerID( );

        if( oldMarkerID != markerID )
        {
                for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
                {
                    if( M_regionFlags( iRegion ) == markerID )
                    {
                        numRegion = iRegion;
                        iRegion = M_numSubregions;
                    }
                }
        }

        M_numElementPerFlag( numRegion ) = M_numElementPerFlag( numRegion ) + 1;
        oldMarkerID = markerID;

    }

    printNumElementPerFlag();

}


template<typename MeshType>
UInt
MeshVolumeSubdivision<MeshType>::
fillElementPerFlag()
{

    if ( !M_readRegionFlags )
    {
        std::cout << "I HAVE NOT READ THE FLAGS" << std::endl;
        assert(false);
    }

    std::vector<int> counters(M_numSubregions, 0);

    UInt nbElements( M_mesh->numElements( ) );
    UInt oldMarkerID = 0;
    UInt numRegion = 0;

    for (UInt iElement(0); iElement < nbElements; iElement++)
    {
        // Extracting the marker
        UInt markerID = M_mesh->element( iElement ).markerID( );

        if( oldMarkerID != markerID )
        {
                for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
                {
                    if( M_regionFlags( iRegion ) == markerID )
                    {
                        numRegion = iRegion;
                        iRegion = M_numSubregions;
                    }
                }
        }

            (M_elements[numRegion])[ counters[numRegion] ] = iElement;
            counters[numRegion]++;
            oldMarkerID = markerID;
    }

}

template<typename MeshType>
void
MeshVolumeSubdivision<MeshType>::
printFlags()
{

    for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
    {
        std::cout << "ID " << M_comm->MyPID() << " Region: " << iRegion
                  << " flag " << M_regionFlags( iRegion ) << std::endl;
    }


}

template<typename MeshType>
void
MeshVolumeSubdivision<MeshType>::
printNumElementPerFlag()
{

    std::cout << "Printing number of element per each flag " << std::endl;

    for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
    {
        std::cout << "ID " << M_comm->MyPID()
                  << " Region: " << iRegion
                  << " flag " << M_regionFlags( iRegion )
                  << " numElements: " << M_numElementPerFlag( iRegion ) << std::endl;
    }

}

template<typename MeshType>
UInt
MeshVolumeSubdivision<MeshType>::
allocatePerElements()
{

    for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
    {
        if( M_numElementPerFlag( iRegion ) > 0 )
        {
            M_elements[iRegion] = new UInt[M_numElementPerFlag( iRegion ) ];
        }
        else
        {
            M_elements[iRegion] = nullptr;
        }
    }


}

template<typename MeshType>
void
MeshVolumeSubdivision<MeshType>::
printElementPerFlag()
{

    std::cout << "Printing number of element per each flag " << std::endl;

    for( UInt iRegion(0); iRegion < M_numSubregions; iRegion++ )
    {
        for( UInt iElement(0); iElement < M_numElementPerFlag( iRegion ); iElement++ )
        {
            std::cout << "ID " << M_comm->MyPID()
                      << " Region: " << iRegion
                      << " flag " << M_regionFlags( iRegion )
                      << " numElements: " << M_numElementPerFlag( iRegion )
                      << " element " << iElement
                      << " globally in the mpi local mesh " << (M_elements[iRegion])[ iElement ]
                      << std::endl;

        }
    }

}


template<typename MeshType>
UInt
MeshVolumeSubdivision<MeshType>::
makeSubDivision()
{
    countElementPerFlag();
    allocatePerElements();
    fillElementPerFlag();

}

template<typename MeshType>
const UInt
MeshVolumeSubdivision<MeshType>::
getNumElements( UInt flag ) const
{
    return M_numElementPerFlag[flag];
}

template<typename MeshType>
const UInt *
MeshVolumeSubdivision<MeshType>::
getSubmesh( UInt flag ) const
{
    return M_elements[flag];
}













} // end LfeV namespace

#endif /* LIFEV_REDUCED_BASIS_SOLVER_MeshVolumeSubdivision_HPP_ */
