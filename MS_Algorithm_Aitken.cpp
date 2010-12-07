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
 *  @brief MultiScale Aitken Algorithm
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 23-10-2009
 */

#include <lifemc/lifesolver/MS_Algorithm_Aitken.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Algorithm_Aitken::MS_Algorithm_Aitken() :
        super               (),
        M_methodMap         (),
        M_method            (),
        M_generalizedAitken ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MS_Algorithm_Aitken::MS_Algorithm_Aitken() \n";
#endif

    M_type = Aitken;

    //Set Map
    M_methodMap["Scalar"]  = Scalar;
    M_methodMap["Vectorial"] = Vectorial;
    //M_methodMap["VectorialBlock"] = VectorialBlock;
}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MS_Algorithm_Aitken::SetupData( const std::string& FileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MS_Algorithm_Aitken::SetupData( algorithm ) \n";
#endif

    super::SetupData( FileName );

    GetPot DataFile( FileName );

    M_generalizedAitken.setDefaultOmega( DataFile( "Solver/Algorithm/Aitken_method/Omega", 1.e-3 ) );
    M_generalizedAitken.useDefaultOmega( DataFile( "Solver/Algorithm/Aitken_method/fixedOmega",   false ) );
    M_generalizedAitken.setOmegaMin( DataFile( "Solver/Algorithm/Aitken_method/range", M_generalizedAitken.defaultOmegaFluid()/1024, 0 ) );
    M_generalizedAitken.setOmegaMax( DataFile( "Solver/Algorithm/Aitken_method/range", M_generalizedAitken.defaultOmegaFluid()*1024, 1 ) );
    M_generalizedAitken.setMinimizationType( DataFile( "Solver/Algorithm/Aitken_method/inverseOmega", true ) );
    M_method = M_methodMap[ DataFile( "Solver/Algorithm/Aitken_method/method", "Vectorial" ) ];
}

void
MS_Algorithm_Aitken::SubIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MS_Algorithm_Aitken::SubIterate( tolerance, subITMax ) \n";
#endif

    super::SubIterate();

    // Verify tolerance
    if ( ToleranceSatisfied() )
    {
        Save( 0, M_couplingResiduals->Norm2() );
        return;
    }

    M_multiscale->ExportCouplingVariables( *M_couplingVariables );

    M_generalizedAitken.restart();

    // Temporary Computation of a Block Vector - Testing purpose
    //VectorType blocksVector( M_couplingVariables ); blocksVector = 0.0;
    //for ( UInt i = 1 ; i < blocksVector.size() ; i = i+2)
    //    blocksVector[i] = 1.0;
    //std::cout << "blocksVector: " << std::endl;
    //blocksVector.ShowMe();

    for ( UInt subIT = 1; subIT <= M_SubiterationsMaximumNumber; ++subIT )
    {
        // To be moved in a post-processing class
        //std::cout << " MS-  CouplingVariables:\n" << std::endl;
        //M_couplingVariables->ShowMe();
        //std::cout << " MS-  CouplingResiduals:\n" << std::endl;
        //M_couplingResiduals->ShowMe();

        // Update Coupling Variables
        switch ( M_method )
        {
        case Scalar:

            *M_couplingVariables += M_generalizedAitken.computeDeltaLambdaScalar( *M_couplingVariables, *M_couplingResiduals );

            break;

        case Vectorial:

            *M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVector( *M_couplingVariables, *M_couplingResiduals, true );

            break;

        case VectorialBlock:

            //*M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVectorBlock( *M_couplingVariables, *M_couplingResiduals, blocksVector, 2 );

            break;
        }

        //std::cout << " MS-  New CouplingVariables:\n" << std::endl;
        //M_couplingVariables->ShowMe();

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->ImportCouplingVariables( *M_couplingVariables );

        // SolveSystem
        M_multiscale->SolveSystem();

        // Display subiteration information
        if ( M_displayer->isLeader() )
            std::cout << " MS-  Sub-iteration n.:                        " << subIT << std::endl;

        // Verify tolerance
        if ( ToleranceSatisfied() )
        {
            Save( subIT, M_couplingResiduals->Norm2() );
            return;
        }
    }

    Save( M_SubiterationsMaximumNumber, M_couplingResiduals->Norm2() );

    MS_ErrorCheck( MS_Tolerance, "Aitken algorithm residual: " + number2string( M_couplingResiduals->Norm2() ) +
                   " (required: " + number2string( M_Tolerance ) + ")\n" );
}

void
MS_Algorithm_Aitken::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Aitken Method       = " << Enum2String( M_method, M_methodMap ) << std::endl;
        std::cout << std::endl << std::endl;
    }
}

} // Namespace LifeV
