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
 *  @brief MultiScale Newton Algorithm
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 26-10-2009
 */

#include <lifemc/lifesolver/MS_Algorithm_Newton.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Algorithm_Newton::MS_Algorithm_Newton() :
    super               (),
    M_solver            (),
    M_Jacobian          ()
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Newton::MS_Algorithm_Newton() \n";
#endif

    M_type = Newton;
}

MS_Algorithm_Newton::MS_Algorithm_Newton( const MS_Algorithm_Newton& algorithm ) :
    super               ( algorithm ),
    M_solver            ( algorithm.M_solver ),
    M_Jacobian          ( algorithm.M_Jacobian )
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Newton::MS_Algorithm_Newton( algorithm ) \n";
#endif

}

// ===================================================
// Operators
// ===================================================
MS_Algorithm_Newton&
MS_Algorithm_Newton::operator=( const MS_Algorithm_Newton& algorithm )
{
    if ( this != &algorithm )
    {
        super::operator=( algorithm );
        M_solver          = algorithm.M_solver;
        M_Jacobian        = algorithm.M_Jacobian;
    }
    return *this;
}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MS_Algorithm_Newton::SetupData( const GetPot& DataFile )
{

#ifdef DEBUG
    Debug( 8012 ) << "MS_Algorithm_Newton::SetupData( DataFile ) \n";
#endif

    super::SetupData( DataFile );

    M_solver.SetCommunicator( *M_comm );
    M_solver.setDataFromGetPot( DataFile, "Solver/Algorithm/Newton_method/AztecOO" );
    //M_solver.setUpPrec(         DataFile, "Solver/Algorithm/Newton_method/Preconditioner" );
}

void
MS_Algorithm_Newton::SubIterate()
{

#ifdef DEBUG
    Debug( 8012 ) << "MS_Algorithm_Newton::SubIterate( tolerance, subITMax ) \n";
#endif

    // Check if it is necessary to performs subiterations
    M_multiscale->ExportCouplingResiduals( *M_couplingResiduals );
    Real residual = M_couplingResiduals->Norm2();

    if ( M_displayer->isLeader() )
        std::cout << " MS-  Residual:                                " << residual << std::endl;

    if ( residual <= M_Tolerance )
        return;

    // Starting Newton Algorithm
    if ( M_displayer->isLeader() )
        std::cout << " MS-  Newton Algorithm" << std::endl;

    M_multiscale->ExportCouplingVariables( *M_couplingVariables );

    VectorType delta( *M_couplingResiduals ); delta = 0.0;
    VectorType minusCouplingResidual( *M_couplingResiduals ); minusCouplingResidual = 0.0;

    for ( UInt subIT = 1; subIT <= M_SubiterationsMaximumNumber; ++subIT )
    {
        // Compute the Jacobian
        if ( subIT == 1 )
        {
            M_Jacobian.reset( new MatrixType( M_couplingVariables->getMap(), 50, 0 ) );
            M_multiscale->ExportJacobian( *M_Jacobian );
            M_Jacobian->GlobalAssemble();
            M_solver.setMatrix( *M_Jacobian );

            //M_Jacobian->spy( "Jacobian" );
        }

        // To be moved in a post-processing class
        //std::cout << " MS-  CouplingVariables:\n" << std::endl;
        //M_couplingVariables->ShowMe();
        //std::cout << " MS-  CouplingResiduals:\n" << std::endl;
        //M_couplingResiduals->ShowMe();

        //Compute delta using -R
        minusCouplingResidual = -( *M_couplingResiduals );

        M_solver.solve( delta, minusCouplingResidual );
        //M_solver.solveSystem( minusCouplingResidual, delta, M_Jacobian, false );

        // Update Coupling Variables using the Newton Method
        *M_couplingVariables += delta;

        //std::cout << " MS-  New CouplingVariables:\n" << std::endl;
        //M_couplingVariables->ShowMe();

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->ImportCouplingVariables( *M_couplingVariables );

        // SolveSystem
        M_multiscale->SolveSystem();

        // Check the new residual
        M_multiscale->ExportCouplingResiduals( *M_couplingResiduals );
        Real residual = M_couplingResiduals->Norm2();

        // Display subiteration information
        if ( M_displayer->isLeader() )
        {
            std::cout << " MS-  Sub-iteration n.:                        " << subIT << std::endl;
            std::cout << " MS-  Residual:                                " << residual << std::endl;
        }

        // Verify tolerance
        if ( residual <= M_Tolerance )
            return;
    }

    MS_ErrorCheck( MS_Tolerance, "Newton algorithm residual: " + number2string( residual ) +
                                 " (required: " + number2string( M_Tolerance ) + ")\n" );
}

void
MS_Algorithm_Newton::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();
    }

    //MPI Barrier
    M_comm->Barrier();
}

} // Namespace LifeV
