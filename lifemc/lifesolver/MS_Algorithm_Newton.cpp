/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-26

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
 \file MS_Algorithm_Newton.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-26
 */

#include <lifemc/lifesolver/MS_Algorithm_Newton.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
MS_Algorithm_Newton::MS_Algorithm_Newton() :
    super               (),
    M_solver            (),
    M_operator          ()
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Newton::MS_Algorithm_Newton() \n";
#endif

    M_type = Newton;
}

MS_Algorithm_Newton::MS_Algorithm_Newton( const MS_Algorithm_Newton& algorithm ) :
    super               ( algorithm ),
    M_solver            ( algorithm.M_solver ),
    M_operator          ( algorithm.M_operator )
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Newton::MS_Algorithm_Newton( algorithm ) \n";
#endif

}

// ===================================================
//! Methods
// ===================================================
MS_Algorithm_Newton&
MS_Algorithm_Newton::operator=( const MS_Algorithm_Newton& algorithm )
{
    if ( this != &algorithm )
    {
        super::operator=( algorithm );
        M_solver          = algorithm.M_solver;
        M_operator        = algorithm.M_operator;
    }
    return *this;
}

// ===================================================
//! Virtual Methods
// ===================================================
void
MS_Algorithm_Newton::SetupData( const GetPot& DataFile )
{

#ifdef DEBUG
    Debug( 8012 ) << "MS_Algorithm_Newton::SetupData( DataFile ) \n";
#endif

    super::SetupData( DataFile );

    M_operator.SetAlgorithm( this );

    M_solver.SetCommunicator( *M_comm );
    M_solver.setDataFromGetPot( DataFile, "Solver/Algorithm/Newton_method/AztecOO" );
    //M_solver.setUpPrec( DataFile, "Solver/Algorithm/Newton_method/Preconditioner" );

    M_solver.setOperator( M_operator );

    M_solver.SetParameters( true );
}

void
MS_Algorithm_Newton::SubIterate( void )
{

#ifdef DEBUG
    Debug( 8012 ) << "MS_Algorithm_Newton::SubIterate( tolerance, subITMax ) \n";
#endif

    if ( M_displayer->isLeader() )
        std::cout << " MS-  Newton Algorithm \n" << std::endl;

    M_multiscale->ExportCouplingVariables( *M_couplingVariables );

    Real residual = 0;
    VectorType delta( *M_couplingResiduals ); delta = 0.0;
    VectorType minusCouplingResidual( *M_couplingResiduals ); minusCouplingResidual = 0.0;
    for ( UInt subIT = 0; subIT < M_SubiterationsMaximumNumber; ++subIT )
    {
        M_multiscale->ExportCouplingResiduals( *M_couplingResiduals );
        residual = M_couplingResiduals->Norm2();

        if ( M_displayer->isLeader() )
        {
            std::cout << " MS-  Sub-iteration n.:                        " << subIT << std::endl;
            std::cout << " MS-  Residual:                                " << residual << std::endl;
        }

        // Verify tolerance
        if ( residual <= M_Tolerance )
            break;

        // To be moved in a post-processing class
        std::cout << " MS-  CouplingVariables:\n" << std::endl;
        M_couplingVariables->ShowMe();
        std::cout << " MS-  CouplingResiduals:\n" << std::endl;
        M_couplingResiduals->ShowMe();

        //Compute delta using -R
        minusCouplingResidual = -( *M_couplingResiduals );
        M_solver.solve( delta, minusCouplingResidual );

        // Update Coupling Variables using the Newton Method
        *M_couplingVariables += delta;

        std::cout << " MS-  New CouplingVariables:\n" << std::endl;
        M_couplingVariables->ShowMe();

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->ImportCouplingVariables( *M_couplingVariables );

        // SolveSystem
        M_multiscale->SolveSystem();
    }
}

void
MS_Algorithm_Newton::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();
    }

    //MPI Barrier
    M_comm->Barrier();
}

} // Namespace LifeV
