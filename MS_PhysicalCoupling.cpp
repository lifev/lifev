/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-09-02

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_PhysicalCoupling.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-09-02
 */

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>

namespace LifeV {

std::map< std::string, couplingsTypes > couplingsMap;

UInt MS_PhysicalCoupling::M_couplingsNumber = 0;

// ===================================================
//! Constructors
// ===================================================
MS_PhysicalCoupling::MS_PhysicalCoupling() :
    M_ID                     (),
    M_type                   (),
    M_dataFile               (),
    M_models                 (),
    M_flags                  (),
    M_couplingName           (),
    M_comm                   (),
    M_displayer              (),
    M_couplingIndex          (),
    M_LocalCouplingVariables (),
    M_LocalCouplingResiduals (),
    M_LocalJacobianProduct   (),
    M_LocalDeltaCouplingVariables ()
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::MS_PhysicalCoupling() \n";
#endif

    M_ID = M_couplingsNumber++;
}

MS_PhysicalCoupling::MS_PhysicalCoupling( const MS_PhysicalCoupling& coupling ) :
    M_ID                     ( coupling.M_ID ),
    M_type                   ( coupling.M_type ),
    M_dataFile               ( coupling.M_dataFile ),
    M_models                 ( coupling.M_models ),
    M_flags                  ( coupling.M_flags ),
    M_couplingName           ( coupling.M_couplingName ),
    M_comm                   ( coupling.M_comm ),
    M_displayer              ( coupling.M_displayer ),
    M_couplingIndex          ( coupling.M_couplingIndex ),
    M_LocalCouplingVariables ( coupling.M_LocalCouplingVariables ),
    M_LocalCouplingResiduals ( coupling.M_LocalCouplingResiduals ),
    M_LocalJacobianProduct   ( coupling.M_LocalJacobianProduct ),
    M_LocalDeltaCouplingVariables ( coupling.M_LocalDeltaCouplingVariables )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::MS_PhysicalCoupling( coupling ) \n";
#endif

    M_ID = M_couplingsNumber++;
}

// ===================================================
//! Methods
// ===================================================
MS_PhysicalCoupling&
MS_PhysicalCoupling::operator=( const MS_PhysicalCoupling& coupling )
{
    if ( this != &coupling )
    {
        M_ID                     = coupling.M_ID;
        M_type                   = coupling.M_type;
        M_dataFile               = coupling.M_dataFile;
        M_models                 = coupling.M_models;
        M_flags                  = coupling.M_flags;
        M_couplingName           = coupling.M_couplingName;
        M_comm                   = coupling.M_comm;
        M_displayer              = coupling.M_displayer;
        M_couplingIndex          = coupling.M_couplingIndex;
        M_LocalCouplingVariables = coupling.M_LocalCouplingVariables;
        M_LocalCouplingResiduals = coupling.M_LocalCouplingResiduals;
        M_LocalJacobianProduct   = coupling.M_LocalJacobianProduct;
        M_LocalDeltaCouplingVariables = coupling.M_LocalDeltaCouplingVariables;
    }
    return *this;
}

void
MS_PhysicalCoupling::BuildCouplingVectorsMap( UInt&             couplingVariablesGlobalNumber,
                                              std::vector<int>& MyGlobalElements )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::BuildCouplingVectorsMap( couplingVariablesGlobalNumber, MyGlobalElements ) \n";
#endif

    M_couplingIndex.second = couplingVariablesGlobalNumber;

    for ( UInt i = 0 ; i < M_couplingIndex.first ; ++i, ++couplingVariablesGlobalNumber )
        MyGlobalElements.push_back( couplingVariablesGlobalNumber );
}

void
MS_PhysicalCoupling::ImportCouplingVariables( const VectorType& CouplingVariables )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ImportCouplingVariables( CouplingVariables ) \n";
#endif

    ImportCouplingVector( CouplingVariables, *M_LocalCouplingVariables );
}

void
MS_PhysicalCoupling::ExportCouplingVariables( VectorType& CouplingVariables )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ExportCouplingVariables( CouplingVariables ) \n";
#endif

    ExportCouplingVector( *M_LocalCouplingVariables, CouplingVariables );
}

void
MS_PhysicalCoupling::ExportJacobianProduct( const VectorType& deltaCouplingVariables, VectorType& JacobianProduct )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_Model_MultiScale::ExportJacobianProduct() \n";
#endif

    VectorType LocalDeltaCouplingVariables( *M_LocalCouplingVariables );
    ImportCouplingVector( deltaCouplingVariables, LocalDeltaCouplingVariables );

    this->ComputeJacobianProduct( LocalDeltaCouplingVariables );

    ExportCouplingVector( *M_LocalJacobianProduct, JacobianProduct );
}

// ===================================================
//! Set Methods
// ===================================================
void
MS_PhysicalCoupling::SetDataFile( const std::string& dataFile )
{

#ifdef DEBUG
    Debug( 8100 ) << "MultiScale_PhysicalMCoupling::SetDataFile( dataFile ) \n";
#endif

    M_dataFile = GetPot( dataFile );

    // Read multiscale parameters
    M_couplingName = M_dataFile( "MultiScale/couplingName", "couplingName" );

    UInt componentSize = M_dataFile.vector_variable_size( "MultiScale/couplingFlags" );

    M_flags.reserve( componentSize );
    for ( UInt j( 0 ); j < componentSize; ++j )
        M_flags.push_back( M_dataFile( "MultiScale/couplingFlags", 0, j ) ); // flags
}

void
MS_PhysicalCoupling::SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm.get() ) );
}

// ===================================================
//! MultiScale PhysicalCoupling Virtual Functions
// ===================================================
void
MS_PhysicalCoupling::ShowMe( void )
{
    std::cout << "Coupling id         = " << M_ID << std::endl
              << "Coupling name       = " << M_couplingName << std::endl
              << "Coupling type       = " << Enum2String( M_type, couplingsMap ) << std::endl << std::endl;

    std::cout << "Models number       = " << modelsNumber() << std::endl;
    std::cout << "Models type(s)      = ";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        std::cout << Enum2String( M_models[i]->GetType(), modelsMap ) << " ";
    std::cout << std::endl << std::endl;

    std::cout << "Flags number        = " << modelsNumber() << std::endl;
    std::cout << "Flags list          = ";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        std::cout << M_flags[i] << " ";
    std::cout << std::endl << std::endl;
}


// ===================================================
//! Protected Methods
// ===================================================
void
MS_PhysicalCoupling::switchErrorMessage( const PhysicalModel_ptr& model )
{
    std::cout << "ERROR: Invalid model type [" << Enum2String( model->GetType(), modelsMap )
              << "] for coupling type [" << Enum2String( M_type, couplingsMap ) << "]" << std::endl;
}

void
MS_PhysicalCoupling::ImportCouplingVector( const VectorType& globalVector, VectorType& localVector )
{
    for ( UInt i(0) ; i < M_couplingIndex.first ; ++i )
        localVector[i] = globalVector[ M_couplingIndex.second + i ];
}

void
MS_PhysicalCoupling::ExportCouplingVector( const VectorType& localVector, VectorType& globalVector )
{
    for ( UInt i(0) ; i < M_couplingIndex.first ; ++i )
        globalVector[ M_couplingIndex.second + i ] = localVector[i];
}

} // Namespace LifeV
