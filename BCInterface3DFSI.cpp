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
 *  @file
 *  @brief File containing the BCInterface3DFSI class
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/BCInterface3DFSI.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterface3DFSI< FSIOperator >::BCInterface3DFSI() :
        M_FSIFunction           (),
        M_physicalSolver        (),
        M_name                  (),
        M_flag                  (),
        M_type                  (),
        M_mode                  (),
        M_comV                  (),
        M_vectorFunctionRobin   (),
        M_robinRHS              (),
        M_robinAlphaCoefficient (),
        M_robinBetaCoefficient  ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFSI::BCInterface3DFSI()" << "\n";
#endif

}

BCInterface3DFSI< FSIOperator >::BCInterface3DFSI( const data_Type& data ) :
        M_FSIFunction           (),
        M_physicalSolver        (),
        M_name                  (),
        M_flag                  (),
        M_type                  (),
        M_mode                  (),
        M_comV                  (),
        M_vectorFunctionRobin   (),
        M_robinRHS              (),
        M_robinAlphaCoefficient (),
        M_robinBetaCoefficient  ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFSI::BCInterface3DFSI( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Methods
// ===================================================
void
BCInterface3DFSI< FSIOperator >::exportData( data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFSI::exportData" << "\n";
#endif

    data.setName( M_name );
    data.setFlag( M_flag );
    data.setType( M_type );
    data.setMode( M_mode );
    data.setComV( M_comV );
}

void
BCInterface3DFSI< FSIOperator >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFSI::updatePhysicalSolverVariables" << "\n";
#endif

    switch ( M_FSIFunction )
    {
    case RobinWall:
    {
        if ( !M_physicalSolver->isSolid() )
            return;

        // Update the physical solver variables
        for ( UInt i( 0 ); i < M_vectorFunctionRobin.size(); ++i )
        {
            boost::shared_ptr< BCInterfaceFunctionSolver< physicalSolver_Type > > castedFunctionSolver =
                boost::dynamic_pointer_cast< BCInterfaceFunctionSolver< physicalSolver_Type > > ( M_vectorFunctionRobin[i] );

            if ( castedFunctionSolver != 0 )
                castedFunctionSolver->updatePhysicalSolverVariables();
        }

        // Set coefficients
        Int gid;
        Real x, y, z;
        Real alpha, beta;
        Real t( M_physicalSolver->dataSolid()->getdataTime()->time() );
        Real timeStep( M_physicalSolver->dataSolid()->getdataTime()->timeStep() );

        Int verticesGlobalNumber( M_physicalSolver->solidMeshPart().meshPartition()->numGlobalVertices() );
        for ( UInt i(0) ; i < M_physicalSolver->solidMeshPart().meshPartition()->numVertices() ; ++i )
        {
            gid = M_physicalSolver->solidMeshPart().meshPartition()->pointInitial( i ).id();

            x   = M_physicalSolver->solidMeshPart().meshPartition()->pointInitial( i ).x();
            y   = M_physicalSolver->solidMeshPart().meshPartition()->pointInitial( i ).y();
            z   = M_physicalSolver->solidMeshPart().meshPartition()->pointInitial( i ).z();

            alpha = M_vectorFunctionRobin[0]->functionTimeSpace( t, x, y, z, 0 );
            beta  = M_vectorFunctionRobin[1]->functionTimeSpace( t, x, y, z, 0 );

            alpha += 2 / timeStep * beta;

            (*M_robinAlphaCoefficient)[gid] = alpha;
            (*M_robinBetaCoefficient)[gid]  = beta;

            (*M_robinAlphaCoefficient)[gid + verticesGlobalNumber] = alpha;
            (*M_robinBetaCoefficient)[gid + verticesGlobalNumber]  = beta;

            (*M_robinAlphaCoefficient)[gid + verticesGlobalNumber * 2] = alpha;
            (*M_robinBetaCoefficient)[gid + verticesGlobalNumber * 2]  = beta;
        }

        // Set displacement and velocity at time tn (mid-point scheme for the solid)
        FSIOperator::vector_Type displacementTn( M_physicalSolver->dFESpace().map(), Repeated, Zero );
        FSIOperator::vector_Type velocityTn( M_physicalSolver->dFESpace().map(), Repeated, Zero );

        M_physicalSolver->exportFluidDisplacement( displacementTn );
        M_physicalSolver->exportSolidVelocity( velocityTn );

        *M_robinRHS = 2 / timeStep * displacementTn + velocityTn;
    }
    default:
        ;// Do nothing
    }
}


// ===================================================
// Set Methods
// ===================================================
void
BCInterface3DFSI< FSIOperator >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFSI::setData" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, FSIFunction > mapFunction;
    mapFunction["DerFluidLoadToFluid"]              = DerFluidLoadToFluid;
    mapFunction["DerFluidLoadToStructure"]          = DerFluidLoadToStructure;
    mapFunction["DerHarmonicExtensionVelToFluid"]   = DerHarmonicExtensionVelToFluid;
    mapFunction["DerStructureDispToSolid"]          = DerStructureDispToSolid;
    mapFunction["FluidInterfaceDisp"]               = FluidInterfaceDisp;
    mapFunction["FluidLoadToStructure"]             = FluidLoadToStructure;
    mapFunction["HarmonicExtensionVelToFluid"]      = HarmonicExtensionVelToFluid;
    mapFunction["SolidLoadToStructure"]             = SolidLoadToStructure;
    mapFunction["StructureDispToHarmonicExtension"] = StructureDispToHarmonicExtension;
    mapFunction["StructureDispToSolid"]             = StructureDispToSolid;
    mapFunction["StructureToFluid"]                 = StructureToFluid;
    mapFunction["RobinWall"]                        = RobinWall;

    // Retrieving the strings
    M_FSIFunction = mapFunction[ data.baseString() ];

    M_name = data.name();
    M_flag = data.flag();
    M_type = data.type();
    M_mode = data.mode();
    M_comV = data.comV();

    if ( M_FSIFunction == RobinWall )
    {
        factory_Type factory;
        M_vectorFunctionRobin.reserve(2);
        data_Type temporaryData ( data );

        // Create the mass term function
        temporaryData.setRobinBaseAlpha();
        M_vectorFunctionRobin.push_back( factory.createFunction( temporaryData ) );

        // Create the RHS
        temporaryData.setRobinBaseBeta();
        M_vectorFunctionRobin.push_back( factory.createFunction( temporaryData ) );
    }
}

} // Namespace LifeV
