//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Model 1D
 *
 *  @author gilles fourestey <gilles.fourestey@epfl.ch>
 *  @date 02-25-2010
 */

#include "MS_Model_1D.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_1D::MS_Model_1D() :
    super                          (),
    M_output                       (),
    M_FESpace                      (),
    M_solver                       (),
    M_linearSolver                 (),
    M_BC                           (),
    M_data                         (),
    M_solution                     ()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_1D::MS_Model_1D() \n";
#endif

    //M_type = oneD;
}



// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================

void
MS_Model_1D::SetupData()
{
    return;
}

void
MS_Model_1D::SetupData(const GetPot& dataFile, std::string section)
{
    Debug( 8120 ) << "MS_Model_1D::SetupData( ) \n";

    //Fluid data
    //this->SetDataFile(dataFile);

    M_data.reset  ( new DataType( dataFile, section ) );
    M_params.reset(new Params(dataFile, section));


    M_linearSolver.reset(new solver_type(*M_comm));
    M_linearSolver->setUpPrec        (dataFile, section + "/prec");
    M_linearSolver->setDataFromGetPot(dataFile, section + "/solver");

}

void
MS_Model_1D::SetupModel()
{
    Debug( 8120 ) << "MS_Model_1D::SetupProblem() \n";

    M_flux.reset  (new Flux(*M_params));
    M_source.reset(new Source(*M_params));

    setUpFESpace();

    M_solver.reset (new OneDModelSolver<Params, Flux, Source> (*M_data, *M_params, *M_flux, *M_source, *M_FESpace, *(this->M_comm)));
    M_solver->setup();

    M_solver->setLinearSolver(M_linearSolver);

    M_BC.reset(new BCType(M_solver->U_thistime(), *M_flux, M_FESpace->dim()));

}

void
MS_Model_1D::BuildSystem()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_1D::BuildSystem() \n";
#endif

    M_solver->initialize();
}


void
MS_Model_1D::UpdateSystem()
{

    Debug( 8120 ) << "MS_Model_1D::UpdateSystem() \n";
    M_data->timeAdvance();
    M_solver->timeAdvance( M_data->time() );

}

void
MS_Model_1D::SolveSystem()
{

    Debug( 8120 ) << "MS_Model_1D::SolveSystem() \n";

    int    count = 0;

    M_solver->iterate( *M_BC, M_data->time() , count );

}

void
MS_Model_1D::SaveSolution()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_1D::SaveSolution() \n";
#endif

    //Post-processing

//     *M_FluidSolution = M_Fluid->solution();
//     M_output->postProcess( M_FluidData->getTime() );

    //MPI Barrier
    //M_comm->Barrier();
}

void
MS_Model_1D::ShowMe()
{
}

// ===================================================
// Methods
// ===================================================
void
MS_Model_1D::SetupLinearData()
{

}

void
MS_Model_1D::SetupLinearModel()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_1D::SetupLinearModel( ) \n";
#endif
}

void
MS_Model_1D::UpdateLinearModel()
{
}

void
MS_Model_1D::SolveLinearModel( bool& SolveLinearSystem )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_1D::SolveLinearModel() \n";
#endif
}

// ===================================================
// Get Methods
// ===================================================
MS_Model_1D::BCType&
MS_Model_1D::GetBC()
{
    return *M_BC;
}

MS_Model_1D::FESpaceType&
MS_Model_1D::GetFESpace()
{
    return *M_FESpace;
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_1D::setUpFESpace()
{
    Debug( 8120 ) << "MS_Model_1D::setupFEspace() \n";

    const RefFE*    refFE = &feSegP1;
    const QuadRule* qR    = &quadRuleSeg3pt;
    const QuadRule* bdQr  = &quadRuleSeg1pt;

    M_FESpace.reset(new FESpace<MeshType, EpetraMap>(M_data->mesh(), *refFE, *qR, *bdQr, 1, *(this->M_comm)));
}


} // Namespace LifeV
