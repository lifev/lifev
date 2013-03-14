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
 *  @brief File containing the Multiscale Model 0D
 *
 *  @version 1.0
 *  @date 30-09-2011
 *  @author Mahmoud Jafargholi <mahmoud.jafargholi@epfl.ch>
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/models/MultiscaleModel0D.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModel0D::MultiscaleModel0D() :
    M_data          ( new data_Type() ),
    M_solver        (),
    M_bc            ( new bcInterface_Type() )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::MultiscaleModel0D() \n";
#endif

    M_type = ZeroDimensional;
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModel0D::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::setupData( fileName ) \n";
#endif

    GetPot dataFile ( fileName );

    std::string circuitDataFile = dataFile ( "0D_Model/CircuitDataFile", "./inputFile.dat" );
    M_bc->createHandler();
    M_bc->fillHandler ( circuitDataFile, "Files" );

    M_data->setup (dataFile, M_bc->handler() );
    if ( M_globalData.get() )
    {
        setupGlobalData ( fileName );
    }

    // The 0D solver requires Rythmos/NOX/Thyra for now
#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
    M_solver.reset ( new solver_Type ( M_data->unknownCounter(), M_comm, M_data->circuitData()  ) );
    M_solver->setup ( M_data->solverData() );
#endif
}

void
MultiscaleModel0D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::setupModel() \n";
#endif

    M_data->initializeSolution();
}

void
MultiscaleModel0D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::buildModel() \n";
#endif

}

void
MultiscaleModel0D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::updateModel() \n";
#endif

}

void
MultiscaleModel0D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::solveModel() \n";
#endif

    // The 0D solver requires Rythmos/NOX/Thyra for now
#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
    M_solver->takeStep ( M_data->dataTime()->previousTime(), M_data->dataTime()->time() );
#endif

}

void
MultiscaleModel0D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::updateSolution() \n";
#endif

}

void
MultiscaleModel0D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::saveSolution() \n";
#endif

    M_data->saveSolution();
}

void
MultiscaleModel0D::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();

        M_data->showMe();
    }
}

Real
MultiscaleModel0D::checkSolution() const
{
    return M_data->circuitData()->Nodes()->nodeListAt ( 1 )->voltage() + M_data->circuitData()->Elements()->elementListAt ( 1 )->current();
}



// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModel0D::setupGlobalData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8160 ) << "MultiscaleModel0D::setupGlobalData( fileName ) \n";
#endif

    GetPot dataFile ( fileName );

    //Global data time
    M_data->setTimeData ( M_globalData->dataTime() );
}

} // Namespace Multiscale
} // Namespace LifeV
