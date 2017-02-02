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

#include <lifev/bc_interface/3D/function/solid/BCInterfaceFunctionSolverDefinedSolid3D.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::BCInterfaceFunctionSolverDefined() :
    M_solid3DFunction       (),
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
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::exportData ( dataPtr_Type& data )
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
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::updatePhysicalSolverVariables" << "\n";
#endif

    switch ( M_solid3DFunction )
    {
        case RobinWall:
        {
            // Update the physical solver variables
            for ( UInt i ( 0 ); i < M_vectorFunctionRobin.size(); ++i )
            {
                functionParserSolverPtr_Type castedFunctionSolver = std::dynamic_pointer_cast< functionParserSolver_Type > ( M_vectorFunctionRobin[i] );

                if ( castedFunctionSolver != 0 )
                {
                    castedFunctionSolver->updatePhysicalSolverVariables();
                }
            }

            // Set coefficients
            Int gid;
            Real x, y, z;
            Real alpha, beta;
            Real t ( M_physicalSolver->data()->dataTime()->time() );
            Real timeStep ( M_physicalSolver->data()->dataTime()->timeStep() );

            Int verticesGlobalNumber ( M_physicalSolver->dispFESpace().mesh()->numGlobalVertices() );
            for ( UInt i (0) ; i < M_physicalSolver->dispFESpace().mesh()->numVertices() ; ++i )
            {
                gid = M_physicalSolver->dispFESpace().mesh()->meshTransformer().pointInitial ( i ).id();

                x   = M_physicalSolver->dispFESpace().mesh()->meshTransformer().pointInitial ( i ).x();
                y   = M_physicalSolver->dispFESpace().mesh()->meshTransformer().pointInitial ( i ).y();
                z   = M_physicalSolver->dispFESpace().mesh()->meshTransformer().pointInitial ( i ).z();

                alpha = M_vectorFunctionRobin[0]->functionTimeSpace ( t, x, y, z, 0 );
                beta  = M_vectorFunctionRobin[1]->functionTimeSpace ( t, x, y, z, 0 );

                alpha += M_physicalSolver->timeAdvancePtr()->coefficientFirstDerivative ( 0 ) / timeStep * beta;

                (*M_robinAlphaCoefficient) [gid] = alpha;
                (*M_robinBetaCoefficient) [gid]  = beta;

                (*M_robinAlphaCoefficient) [gid + verticesGlobalNumber] = alpha;
                (*M_robinBetaCoefficient) [gid + verticesGlobalNumber]  = beta;

                (*M_robinAlphaCoefficient) [gid + verticesGlobalNumber * 2] = alpha;
                (*M_robinBetaCoefficient) [gid + verticesGlobalNumber * 2]  = beta;
            }

            M_physicalSolver->timeAdvancePtr()->updateRHSFirstDerivative ( timeStep );
            *M_robinRHS = M_physicalSolver->timeAdvancePtr()->rhsContributionFirstDerivative();
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
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::setData ( const dataPtr_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::setData" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, Solid3DFunction > mapFunction;
    mapFunction["RobinWall"] = RobinWall;

    // Retrieving the strings
    M_solid3DFunction = mapFunction[ data->baseString() ];

    M_name = data->name();
    M_flag = data->flag();
    M_type = data->type();
    M_mode = data->mode();
    M_componentsVector = data->componentsVector();

    if ( M_solid3DFunction == RobinWall )
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
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::baseType() const
{
    switch ( M_solid3DFunction )
    {
        case RobinWall:

            return BASEVector3D;

        default:

            std::cout << " !!! Error: " << M_solid3DFunction << " is not available as a Solid3DFunction !!!" << std::endl;
            return BASEDefault;
    }
}



// ===================================================
// Private Methods
// ===================================================
void
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::checkFunction ( BCVectorInterface& /*base*/ )
{
    switch ( M_solid3DFunction )
    {
        default:

            std::cout << " !!! Error: " << M_solid3DFunction << " is not available as a BCVectorInterface !!!" << std::endl;

            break;
    }
}

void
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::checkFunction ( BCVector& base )
{
    switch ( M_solid3DFunction )
    {
        case RobinWall:

#ifdef HAVE_LIFEV_DEBUG
            debugStream ( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          RobinWall" << "\n";
#endif

            // Define the vectors
            M_robinRHS.reset ( new physicalSolver_Type::vector_Type ( M_physicalSolver->dispFESpace().map(), Repeated, Zero ) );
            M_robinAlphaCoefficient.reset ( new physicalSolver_Type::vector_Type ( M_physicalSolver->dispFESpace().map(), Repeated, Zero ) );
            M_robinBetaCoefficient.reset ( new physicalSolver_Type::vector_Type ( M_physicalSolver->dispFESpace().map(), Repeated, Zero ) );

            // Set the vectors (still empty)
            base.setRhsVector ( *M_robinRHS, M_physicalSolver->dispFESpace().dof().numTotalDof(), 0 );
            base.setRobinCoeffVector ( *M_robinAlphaCoefficient );
            base.setBetaCoeffVector ( *M_robinBetaCoefficient );

            // Set the physical solver in the Robin functions for alpha and beta
            for ( UInt i ( 0 ); i < M_vectorFunctionRobin.size(); ++i )
            {
                functionParserSolverPtr_Type castedFunctionSolver = std::dynamic_pointer_cast< functionParserSolver_Type > ( M_vectorFunctionRobin[i] );

                if ( castedFunctionSolver != 0 )
                {
                    castedFunctionSolver->setPhysicalSolver ( M_physicalSolver );
                }
            }

            break;

        default:

            std::cout << " !!! Error: " << M_solid3DFunction << " is not available as a BCVector !!!" << std::endl;

            break;
    }
}

void
BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >::checkFunction ( BCFunctionBase& /*base*/ )
{
    switch ( M_solid3DFunction )
    {
        default:

            std::cout << " !!! Error: " << M_solid3DFunction << " is not available as a BCFunction !!!" << std::endl;

            return;
    }
}

} // Namespace LifeV
