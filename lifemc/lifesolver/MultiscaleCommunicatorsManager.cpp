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

#include <lifemc/lifesolver/MultiscaleCommunicatorsManager.hpp>

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
        M_parallelModelsID   (),
        M_parallelModelsLoad ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleCommunicatorsManager::MultiscaleCommunicatorsManager() \n";
#endif

}

// ===================================================
// Methods
// ===================================================
void
MultiscaleCommunicatorsManager::splitCommunicators()
{
//    MPI_Comm localComm;
//    MPI_Comm_split( ( dynamic_cast<Epetra_MpiComm*> ( &(*M_comm) ) )->Comm(), M_comm->MyPID(), M_comm->MyPID(), &localComm );
//    M_comm.reset( new Epetra_MpiComm( localComm ) );
}

void
MultiscaleCommunicatorsManager::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        std::cout << std::endl;
        std::cout << "================= Communicators Information =================" << std::endl << std::endl;

        std::cout << "Serial models number = " << M_serialModelsID.size() << std::endl;
        std::cout << "Serial models list   = ";
        for ( UInt i( 0 ) ; i < M_serialModelsID.size() ; ++i )
            std::cout << M_serialModelsID[i] << " ";
        std::cout << std::endl << std::endl;

        std::cout << "Parallel models number = " << M_parallelModelsID.size() << std::endl;
        std::cout << "Parallel models list   = ";
        for ( UInt i( 0 ) ; i < M_parallelModelsID.size() ; ++i )
            std::cout << M_parallelModelsID[i] << " ";
        std::cout << std::endl;
        std::cout << "Parallel models load   = ";
        for ( UInt i( 0 ) ; i < M_parallelModelsID.size() ; ++i )
            std::cout << M_parallelModelsLoad[i] << " ";
        std::cout << std::endl << std::endl;

        std::cout << "=============================================================" << std::endl << std::endl;
    }
}

} // Namespace multiscale
} // Namespace LifeV
