/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-23

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
 \file MS_Algorithm_Aitken.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-23
 */

#include <lifemc/lifesolver/MS_Algorithm_Aitken.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
MS_Algorithm_Aitken::MS_Algorithm_Aitken() :
    super               (),
    M_methodMap         (),
    M_method            (),
    M_inverseOmega      (),
    M_generalizedAitken ()
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Aitken::MS_Algorithm_Aitken() \n";
#endif

    M_type = Aitken;

    //Set Map
    M_methodMap["Scalar"]  = Scalar;
    M_methodMap["Vectorial"] = Vectorial;
    //M_methodMap["VectorialBlock"] = VectorialBlock;
}

MS_Algorithm_Aitken::MS_Algorithm_Aitken( const MS_Algorithm_Aitken& algorithm ) :
    super               ( algorithm ),
    M_methodMap         ( algorithm.M_methodMap ),
    M_method            ( algorithm.M_method ),
    M_inverseOmega      ( algorithm.M_inverseOmega ),
    M_generalizedAitken ( algorithm.M_generalizedAitken )
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Aitken::MS_Algorithm_Aitken( algorithm ) \n";
#endif

}

// ===================================================
//! Methods
// ===================================================
MS_Algorithm_Aitken&
MS_Algorithm_Aitken::operator=( const MS_Algorithm_Aitken& algorithm )
{
    if ( this != &algorithm )
    {
        super::operator=( algorithm );
        M_methodMap         = algorithm.M_methodMap;
        M_method            = algorithm.M_method;
        M_inverseOmega      = algorithm.M_inverseOmega;
        M_generalizedAitken = algorithm.M_generalizedAitken;
    }
    return *this;
}

// ===================================================
//! Virtual Methods
// ===================================================
void
MS_Algorithm_Aitken::SetupData( const GetPot& DataFile )
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Aitken::SetupData( algorithm ) \n";
#endif

    super::SetupData( DataFile );

    M_generalizedAitken.setDefault( DataFile( "Solver/Algorithm/Aitken_method/Omega", 1.e-3 ) );
    M_generalizedAitken.UseDefaultOmega( DataFile( "Solver/Algorithm/Aitken_method/fixedOmega",   false ) );
    M_method       = M_methodMap[ DataFile( "Solver/Algorithm/Aitken_method/method", "Vectorial" ) ];
    M_inverseOmega = DataFile( "Solver/Algorithm/Aitken_method/inverseOmega", true );
}

void
MS_Algorithm_Aitken::SubIterate( void )
{

#ifdef DEBUG
    Debug( 8011 ) << "MS_Algorithm_Aitken::SubIterate( tolerance, subITMax ) \n";
#endif

    if ( M_displayer->isLeader() )
        std::cout << " MS-  Aitken Algorithm \n" << std::endl;

    M_multiscale->ExportCouplingVariables( *M_couplingVariables );

    Real residual = 0;
    M_generalizedAitken.restart( true );

    // Temporary Computation of a Block Vector - Testing purpose
    //VectorType blocksVector( M_couplingVariables ); blocksVector = 0.0;
    //for ( UInt i = 1 ; i < blocksVector.size() ; i = i+2)
    //    blocksVector[i] = 1.0;
    //std::cout << "blocksVector: " << std::endl;
    //blocksVector.ShowMe();

    for ( UInt subIT = 0; subIT < M_SubiterationsMaximumNumber; ++subIT )
    {
        M_multiscale->ExportCouplingResiduals( *M_couplingResiduals );
        residual = M_couplingResiduals->Norm2();

        if ( M_displayer->isLeader() )
        {
            std::cout << " MS-  Sub-iteration n.:                        " << subIT << std::endl;
            std::cout << " MS-  Residual:                                " << residual << std::endl;
        }

        // To be moved in a post-processing class
        std::cout << " MS-  CouplingVariables:\n" << std::endl;
        M_couplingVariables->ShowMe();
        std::cout << " MS-  CouplingResiduals:\n" << std::endl;
        M_couplingResiduals->ShowMe();

        // Verify tolerance
        if ( residual <= M_Tolerance )
            break;

        // Update Coupling Variables
        switch ( M_method )
        {
            case Scalar:

                *M_couplingVariables += M_generalizedAitken.computeDeltaLambdaScalar( *M_couplingVariables, *M_couplingResiduals, M_inverseOmega );

                break;

            case Vectorial:

                *M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVector( *M_couplingVariables, *M_couplingResiduals, M_inverseOmega, true );

                break;

            case VectorialBlock:

                //*M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVectorBlock( *M_couplingVariables, *M_couplingResiduals, blocksVector, 2, true );

                break;
        }

        std::cout << " MS-  New CouplingVariables:\n" << std::endl;
        M_couplingVariables->ShowMe();

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->ImportCouplingVariables( *M_couplingVariables );

        // SolveSystem
        M_multiscale->SolveSystem();
    }
}

void
MS_Algorithm_Aitken::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Aitken Method       = " << Enum2String( M_method, M_methodMap ) << std::endl
                  << "Minimize inv. Omega = " << ( M_inverseOmega ? "true" : "false" ) << std::endl;
        std::cout << std::endl << std::endl;
    }

    //MPI Barrier
    M_comm->Barrier();
}

} // Namespace LifeV
