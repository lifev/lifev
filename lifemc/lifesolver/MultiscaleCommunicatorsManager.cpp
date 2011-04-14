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
//@HEADERR

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
        M_comm              (),
        M_commContainer     (),
        M_loadContainer     (),
        M_modelsIDContainer ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleCommunicatorsManager::MultiscaleCommunicatorsManager() \n";
#endif

}

// ===================================================
// Methods
// ===================================================
void
MultiscaleCommunicatorsManager::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        std::cout << std::endl;
        std::cout << "================= Communicators Information =================" << std::endl << std::endl;

        std::cout << "Groups number       = " << M_loadContainer.size() << std::endl;
        std::cout << "Load / Models       = ";
        for ( UInt i( 0 ); i < M_loadContainer.size(); ++i )
        {
            std::cout << M_loadContainer[i] << " / ";
            for ( UInt j( 0 ); j < M_modelsIDContainer[i].size(); ++j )
                std::cout << M_modelsIDContainer[i][j] << " ";
            std::cout << std::endl << "                      ";
        }

        std::cout << std::endl;

        std::cout << "=============================================================" << std::endl << std::endl;
    }
}

} // Namespace multiscale
} // Namespace LifeV
