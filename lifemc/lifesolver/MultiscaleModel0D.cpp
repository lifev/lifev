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
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <MultiscaleModel0D.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModel0D::MultiscaleModel0D():
    M_data				(new ZeroDimensionalData),
    M_bc                                ( new bcInterface_Type() )
{
}

MultiscaleModel0D::~MultiscaleModel0D()
{
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModel0D::setupData( const std::string& fileName )
{
    GetPot dataFile( fileName );
    //Set Input Output files
    std::string circuitDataFile = dataFile( "0D_Model/CircuitDataFile", "./inputFile.dat" );
    M_bc->createHandler();
    M_bc->fillHandler( circuitDataFile, "Files" );
    std::cout<<"\n Im at BC handler \n";
    M_data->setup(dataFile, M_bc );
    M_solver.reset(new ZeroDimensionalSolver(M_data->unknownCounter(),
                                             M_comm,
                                             M_data->circuitData()  ) );
    M_solver->setup(M_data->solverData());
}

void
MultiscaleModel0D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8151 ) << "MultiscaleModel0D::setupModel() \n";
#endif
    M_data->initializeSolution();
    M_Tn = M_data->dataTime() ->time();
    M_TnPlus = M_Tn;
    updateModel();
}

void
MultiscaleModel0D::buildModel()
{   
    M_Tn = M_TnPlus;
    M_TnPlus = M_data->dataTime() ->time();

}

void
MultiscaleModel0D::updateModel()
{
    M_Tn = M_TnPlus;
    M_TnPlus = M_data->dataTime() ->time();
}

void
MultiscaleModel0D::solveModel()
{
    M_solver->takeStep(M_Tn,M_TnPlus);
}

void
MultiscaleModel0D::saveSolution()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8151 ) << "MultiscaleModel0D::saveSolution() \n";
#endif
    M_data->saveSolution();
}

void
MultiscaleModel0D::showMe()
{
    std::cout << "This MultiscaleModel0D----------------------------"<< std::endl;
    M_data->showMe();
}

Real
MultiscaleModel0D::checkSolution() const
{
    return (0.0);
}

} // Namespace Multiscale
} // Namespace LifeV
