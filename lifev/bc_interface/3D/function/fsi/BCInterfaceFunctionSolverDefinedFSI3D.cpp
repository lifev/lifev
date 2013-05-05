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
 *  @brief File containing the BCInterfaceFunctionSolverDefined class and specializations
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/bc_interface/3D/function/fsi/BCInterfaceFunctionSolverDefinedFSI3D.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceFunctionSolverDefined< BCHandler, FSIOperator >::BCInterfaceFunctionSolverDefined() :
    M_FSIFunction           (),
    M_physicalSolver        (),
    M_name                  (),
    M_flag                  (),
    M_type                  (),
    M_mode                  (),
    M_componentsVector      (),
    M_vectorFunctionRobin   (),
    M_robinRHS              (),
    M_robinAlphaCoefficient (),
    M_robinBetaCoefficient  ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::BCInterfaceFunctionSolverDefined()" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< BCHandler, FSIOperator >::exportData ( dataPtr_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::exportData" << "\n";
#endif

    data->setName ( M_name );
    data->setFlag ( M_flag );
    data->setType ( M_type );
    data->setMode ( M_mode );
    data->setComponentsVector ( M_componentsVector );
}

void
BCInterfaceFunctionSolverDefined< BCHandler, FSIOperator >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::updatePhysicalSolverVariables" << "\n";
#endif

    switch ( M_FSIFunction )
    {
        case RobinWall:
        {
            if ( !M_physicalSolver->isSolid() )
            {
                return;
            }

            // Update the physical solver variables
            for ( UInt i ( 0 ); i < M_vectorFunctionRobin.size(); ++i )
            {
                functionParserSolverPtr_Type castedFunctionSolver = boost::dynamic_pointer_cast< functionParserSolver_Type > ( M_vectorFunctionRobin[i] );

                if ( castedFunctionSolver != 0 )
                {
                    castedFunctionSolver->updatePhysicalSolverVariables();
                }
            }

            // Set coefficients
            Int gid;
            Real x, y, z;
            Real alpha, beta;
            Real t ( M_physicalSolver->dataSolid()->dataTime()->time() );
            Real timeStep ( M_physicalSolver->dataSolid()->dataTime()->timeStep() );

            Int verticesGlobalNumber ( M_physicalSolver->solidLocalMesh().numGlobalVertices() );
            for ( UInt i (0) ; i < M_physicalSolver->solidLocalMesh().numVertices() ; ++i )
            {
                gid = M_physicalSolver->solidLocalMesh().meshTransformer().pointInitial ( i ).id();

                x   = M_physicalSolver->solidLocalMesh().meshTransformer().pointInitial ( i ).x();
                y   = M_physicalSolver->solidLocalMesh().meshTransformer().pointInitial ( i ).y();
                z   = M_physicalSolver->solidLocalMesh().meshTransformer().pointInitial ( i ).z();

                alpha = M_vectorFunctionRobin[0]->functionTimeSpace ( t, x, y, z, 0 );
                beta  = M_vectorFunctionRobin[1]->functionTimeSpace ( t, x, y, z, 0 );

                alpha += M_physicalSolver->solidTimeAdvance()->coefficientFirstDerivative ( 0 ) / timeStep * beta;

                (*M_robinAlphaCoefficient) [gid] = alpha;
                (*M_robinBetaCoefficient) [gid]  = beta;

                (*M_robinAlphaCoefficient) [gid + verticesGlobalNumber] = alpha;
                (*M_robinBetaCoefficient) [gid + verticesGlobalNumber]  = beta;

                (*M_robinAlphaCoefficient) [gid + verticesGlobalNumber * 2] = alpha;
                (*M_robinBetaCoefficient) [gid + verticesGlobalNumber * 2]  = beta;
            }

            M_physicalSolver->solidTimeAdvance()->updateRHSFirstDerivative ( timeStep );
            if ( M_physicalSolver->data().method().compare ("monolithicGE") == 0 || M_physicalSolver->data().method().compare ("monolithicGI") == 0 )
            {
                M_robinRHS->subset ( M_physicalSolver->solidTimeAdvance()->rhsContributionFirstDerivative(),
                                     boost::dynamic_pointer_cast< FSIMonolithic > ( M_physicalSolver )->offset() );
            }
            else
            {
                *M_robinRHS = M_physicalSolver->solidTimeAdvance()->rhsContributionFirstDerivative();
            }
            break;
        }
        default:
            break;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< BCHandler, FSIOperator >::setData ( const dataPtr_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::setData" << "\n";
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
    M_FSIFunction = mapFunction[ data->baseString() ];

    M_name = data->name();
    M_flag = data->flag();
    M_type = data->type();
    M_mode = data->mode();
    M_componentsVector = data->componentsVector();

    if ( M_FSIFunction == RobinWall )
    {
        factory_Type factory;
        M_vectorFunctionRobin.reserve (2);
        dataPtr_Type temporaryData ( new data_Type ( *data ) );

        // Create the mass term function
        temporaryData->setRobinBaseAlpha();
        M_vectorFunctionRobin.push_back ( factory.createFunctionParser ( temporaryData ) );

        // Create the RHS
        temporaryData->setRobinBaseBeta();
        M_vectorFunctionRobin.push_back ( factory.createFunctionParser ( temporaryData ) );
    }
}


// ===================================================
// Get Methods
// ===================================================
baseContainer_Type
BCInterfaceFunctionSolverDefined< BCHandler, FSIOperator >::baseType() const
{
    switch ( M_FSIFunction )
    {
        case DerFluidLoadToFluid:
        case DerFluidLoadToStructure:
        case DerHarmonicExtensionVelToFluid:
        case DerStructureDispToSolid:
        case FluidInterfaceDisp:
        case FluidLoadToStructure:
        case HarmonicExtensionVelToFluid:
        case SolidLoadToStructure:
        case StructureDispToHarmonicExtension:
        case StructureDispToSolid:
        case StructureToFluid:

            return BASEVectorInterface3D;

        case RobinWall:

            return BASEVector3D;

        default:

            std::cout << " !!! Error: " << M_FSIFunction << " is not available as a FSIFunction !!!" << std::endl;
            return BASEDefault;
    }
}

} // Namespace LifeV
