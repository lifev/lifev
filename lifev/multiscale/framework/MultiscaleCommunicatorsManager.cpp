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
 *  @brief File containing the Multiscale Communicators Manager
 *
 *  @date 13-04-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/framework/MultiscaleCommunicatorsManager.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors
// ===================================================
MultiscaleCommunicatorsManager::MultiscaleCommunicatorsManager() :
    M_comm               (),
    M_commContainer      (),
    M_serialModelsID     (),
    M_serialProcesses    (),
    M_parallelModelsID   (),
    M_parallelModelsLoad (),
    M_parallelProcesses  ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8005 ) << "MultiscaleCommunicatorsManager::MultiscaleCommunicatorsManager() \n";
#endif

}

// ===================================================
// Methods
// ===================================================
void
MultiscaleCommunicatorsManager::splitCommunicator()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8005 ) << "MultiscaleCommunicatorsManager::splitCommunicator() \n";
#endif

    // Preliminaries
    Int myPID = M_comm->MyPID();
    Int numberOfProcesses = M_comm->NumProc();
    MPI_Comm comm = ( boost::dynamic_pointer_cast< Epetra_MpiComm > ( M_comm ) )->Comm();

    // Group initialization
    MPI_Group commGroup;
    MPI_Comm_group ( comm, &commGroup );

    // Serial models: assign processes number to the models
    M_serialProcesses.resize ( M_serialModelsID.size(), std::vector< Int > ( 1, 0 ) );
    for ( Int i (0) ; i < static_cast <Int> ( M_serialProcesses.size() ) ; ++i )
    {
        M_serialProcesses[i][0] = i % numberOfProcesses;
    }

    // Serial models: create communicators
    Int serialMembers[1] = { myPID };
    MPI_Group serialCommGroup;
    MPI_Group_incl ( commGroup, 1, serialMembers, &serialCommGroup );

    MPI_Comm serialComm;
    MPI_Comm_create ( comm, serialCommGroup, &serialComm );

    for ( Int i (0) ; i < static_cast <Int> ( M_serialProcesses.size() ) ; ++i )
        if ( M_serialProcesses[i][0] == myPID )
        {
            M_commContainer[M_serialModelsID[i]].reset ( new Epetra_MpiComm ( serialComm ) );
        }

    // Parallel models: identify number of processes per model
    std::vector<Real> localNumberOfProcesses ( M_parallelModelsID.size(), 0 );
    parallelProcessesDistribution ( localNumberOfProcesses, numberOfProcesses );

    // Parallel models: assign processes number to the models
    parallelProcessesAssignment ( M_parallelProcesses, localNumberOfProcesses, numberOfProcesses );

    // Parallel models: create communicators
    for ( UInt i (0) ; i < M_parallelModelsID.size() ; ++i )
    {
        // Definitions
        bool myComm ( false );
        Int parallelMembers[M_parallelProcesses[i].size()];

        // Fill parallel members
        for ( UInt j (0) ; j < M_parallelProcesses[i].size() ; ++j )
        {
            parallelMembers[j] = M_parallelProcesses[i][j];
            if ( parallelMembers[j] == myPID )
            {
                myComm = true;
            }
        }

        // Create parallel group
        MPI_Group parallelCommGroup;
        MPI_Group_incl ( commGroup, M_parallelProcesses[i].size(), parallelMembers, &parallelCommGroup );

        // Create parallel comm
        MPI_Comm localParallelComm;
        MPI_Comm_create ( comm, parallelCommGroup, &localParallelComm );

        // Assign parallel comm
        if ( myComm )
        {
            M_commContainer[M_parallelModelsID[i]].reset ( new Epetra_MpiComm ( localParallelComm ) );
        }
    }
}

bool
MultiscaleCommunicatorsManager::myModel ( const UInt& modelID ) const
{
    if ( M_commContainer.find ( modelID ) != M_commContainer.end() )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void
MultiscaleCommunicatorsManager::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        std::cout << "Serial models number    = " << M_serialModelsID.size() << std::endl;
        std::cout << "Serial models list      = ";
        for ( UInt i ( 0 ) ; i < M_serialModelsID.size() ; ++i )
        {
            std::cout << M_serialModelsID[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Serial models processes = ";
        for ( UInt i ( 0 ) ; i < M_serialModelsID.size() ; ++i )
        {
            std::cout << M_serialProcesses[i][0] << " ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Parallel models number  = " << M_parallelModelsID.size() << std::endl;
        for ( UInt i ( 0 ) ; i < M_parallelModelsID.size() ; ++i )
        {
            std::cout << "Model " << M_parallelModelsID[i]
                      << ", load " << M_parallelModelsLoad[i]
                      << "%, processes: ";
            for ( UInt j (0) ; j < M_parallelProcesses[i].size() ; ++j )
            {
                std::cout << M_parallelProcesses[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
MultiscaleCommunicatorsManager::addGroup ( const Real& load, const modelsID_Type& modelsID )
{
    if ( load < 0 )
    {
        M_serialModelsID.insert ( M_serialModelsID.end(), modelsID.begin(), modelsID.end() );
    }
    else
    {
        // We sort the models from the cheapest to the most expensive
        modelsIDIterator_Type modelsIDIterator     = M_parallelModelsID.begin();
        modelsLoadIterator_Type modelsLoadIterator = M_parallelModelsLoad.begin();
        for ( ; modelsIDIterator != M_parallelModelsID.end() ; ++modelsIDIterator, ++modelsLoadIterator )
            if ( load < *modelsLoadIterator )
            {
                break;
            }

        M_parallelModelsID.insert ( modelsIDIterator, modelsID.begin(), modelsID.end() );

        modelsLoad_Type loadVector ( modelsID.size(), load );
        M_parallelModelsLoad.insert ( modelsLoadIterator, loadVector.begin(), loadVector.end() );
    }
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleCommunicatorsManager::parallelProcessesDistribution ( std::vector<Real>& localNumberOfProcesses, const Int& numberOfProcesses )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8005 ) << "MultiscaleCommunicatorsManager::parallelProcessesDistribution() \n";
#endif

    // Preliminary distribution
    for ( UInt i (0) ; i < M_parallelModelsID.size() ; ++i )
    {
        localNumberOfProcesses[i] = numberOfProcesses * ( M_parallelModelsLoad[i] / 100 );
    }

    //    if ( M_comm->MyPID() == 0 )
    //    {
    //        std::cout << "Preliminary distribution" << std::endl;
    //        for ( UInt i(0) ; i < M_parallelModelsID.size() ; ++i )
    //            std::cout << M_parallelModelsID[i] << " " << localNumberOfProcesses[i] << std::endl;
    //        std::cout << std::endl;
    //    }

    // Definitions
    Real availableResource ( 0 );
    bool optimize ( false );

    // Set all the model as unoptimized
    std::vector<bool> unoptimized ( M_parallelModelsID.size(), true );

    // Check that all the models have at least one processor
    for ( UInt i (0) ; i < M_parallelModelsID.size() ; ++i )
    {
        if ( localNumberOfProcesses[i] < 1 )
        {
            availableResource -= 1 - localNumberOfProcesses[i];
            localNumberOfProcesses[i] = 1;

            unoptimized[i] = false;
        }
        optimize |= unoptimized[i];
    }

    //    if ( M_comm->MyPID() == 0 )
    //    {
    //        std::cout << "After one CPU check" << std::endl;
    //        std::cout << "Sum:        " << availableResource << std::endl;
    //        for ( UInt i(0) ; i < M_parallelModelsID.size() ; ++i )
    //            std::cout << M_parallelModelsID[i] << " " << localNumberOfProcesses[i] << std::endl;
    //        std::cout << std::endl;
    //    }

    // Optimize the other models (using the cores per node information)
    Int  totalResources (0);
    Real optimizationResources (0);
    Real resourcesToBelow (0);
    Real resourcesToAbove (0);
    std::vector<Real> delta ( M_parallelModelsID.size(), 0 );

    while ( optimize )
    {

        //        if ( M_comm->MyPID() == 0 )
        //        {
        //            std::cout << "Optimization progress" << std::endl;
        //            std::cout << "Sum:        " << availableResource << std::endl;
        //        }

        for ( UInt i (0) ; i < M_parallelModelsID.size() ; ++i )
            if ( unoptimized[i] )
            {
                if ( std::floor ( localNumberOfProcesses[i] / multiscaleCoresPerNode ) > 0 )
                {
                    resourcesToBelow = std::fmod ( localNumberOfProcesses[i], multiscaleCoresPerNode );
                    resourcesToAbove = resourcesToBelow - multiscaleCoresPerNode;
                }
                else
                {
                    resourcesToBelow  = std::fmod ( localNumberOfProcesses[i], 1 );
                    resourcesToAbove = resourcesToBelow - 1;
                }
                //                if ( M_comm->MyPID() == 0 )
                //                {
                //                    std::cout << "to below:   " << resourcesToBelow << std::endl;
                //                    std::cout << "to above:   " << resourcesToAbove << std::endl;
                //                }
                // Compute the proposed optimization
                if ( std::fabs ( resourcesToBelow ) >= std::fabs ( resourcesToAbove ) )
                {
                    delta[i] += resourcesToAbove;
                }
                else
                {
                    delta[i] += resourcesToBelow;
                }
                optimizationResources += delta[i];
            }


        //        if ( M_comm->MyPID() == 0 )
        //            std::cout << "currentSum: " << optimizationResources << std::endl;
        totalResources = roundToInteger ( availableResource + optimizationResources );
        if ( totalResources >= 0 )
        {
            //            if ( M_comm->MyPID() == 0 )
            //                std::cout << "TRUE" << std::endl;

            for ( UInt i (0) ; i < M_parallelModelsID.size() ; ++i )
                if ( unoptimized[i] )
                {
                    localNumberOfProcesses[i] -= delta[i];

                    if ( totalResources > 0 )
                    {
                        Int availableSlot = multiscaleCoresPerNode - std::fmod ( localNumberOfProcesses[i], multiscaleCoresPerNode );
                        for ( ; ( availableSlot > 0 && totalResources > 0 ) ; ++localNumberOfProcesses[i], --totalResources );
                    }

                    unoptimized[i] = false;
                }
        }
        else  // Proposed optimization refused
        {
            //            if ( M_comm->MyPID() == 0 )
            //                std::cout << "FALSE" << std::endl;
            UInt ID (0);
            optimizationResources = multiscaleCoresPerNode + 1;
            for ( UInt i (0) ; i < M_parallelModelsID.size() ; ++i )
                if ( unoptimized[i] && ( delta[i] < optimizationResources ) )
                {
                    ID = i;
                    optimizationResources = delta[ID];
                }

            if ( optimizationResources >= 0 )
            {
                availableResource += optimizationResources;
                localNumberOfProcesses[ID] -= optimizationResources;
            }
            else
            {
                if ( std::floor ( localNumberOfProcesses[ID] / multiscaleCoresPerNode ) > 0 )
                {
                    availableResource += multiscaleCoresPerNode + optimizationResources;
                    localNumberOfProcesses[ID] -= multiscaleCoresPerNode + optimizationResources;
                }
                else
                {
                    availableResource += 1 + optimizationResources;
                    localNumberOfProcesses[ID] -= 1 + optimizationResources;
                }
            }
            unoptimized[ID] = false;
        }

        // Reset variables
        optimizationResources = 0;
        delta.assign ( M_parallelModelsID.size(), 0 );

        // Update optimize
        optimize =  false;
        for ( UInt i (0) ; i < M_parallelModelsID.size() ; ++i )
        {
            optimize |= unoptimized[i];
        }

        //        if ( M_comm->MyPID() == 0 )
        //        {
        //            for ( UInt i(0) ; i < M_parallelModelsID.size() ; ++i )
        //                std::cout << M_parallelModelsID[i] << " " << localNumberOfProcesses[i] << std::endl;
        //            std::cout << std::endl;
        //        }
    }
}

void
MultiscaleCommunicatorsManager::parallelProcessesAssignment ( std::vector< std::vector< Int > >& parallelProcesses, const std::vector<Real>& localNumberOfProcesses, const Int& numberOfProcesses )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8005 ) << "MultiscaleCommunicatorsManager::parallelProcessesAssignment() \n";
#endif

    // Initialize the vector
    parallelProcesses.resize ( M_parallelModelsID.size() );
    for ( UInt i (0) ; i < parallelProcesses.size() ; ++i )
    {
        parallelProcesses[i] = std::vector< Int > ( localNumberOfProcesses[i], 0 );
    }

    // We place the models on the nodes starting from the most expensive
    UInt processID (0);
    for ( Int i ( parallelProcesses.size() - 1 ) ; i > -1; --i )
    {
        for ( UInt j (0) ; j < localNumberOfProcesses[i] ; ++j, ++processID )
        {
            parallelProcesses[i][j] = processID % numberOfProcesses;
        }
    }

    //    if ( M_comm->MyPID() == 0 )
    //    {
    //        std::cout << "parallelProcessesAssignment " << std::endl;
    //        for ( UInt i(0) ; i < parallelProcesses.size() ; ++i )
    //        {
    //            std::cout << M_parallelModelsID[i] << ": ";
    //            for ( UInt j(0) ; j < parallelProcesses[i].size() ; ++j )
    //                std::cout << parallelProcesses[i][j] << " ";
    //            std::cout << std::endl;
    //        }
    //    }
}

} // Namespace multiscale
} // Namespace LifeV
