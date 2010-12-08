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
 *  @brief MultiScale Algorithm
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 23-10-2009
 */

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV
{

std::map< std::string, algorithmsTypes > MS_algorithmsMap;

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Algorithm::MS_Algorithm() :
        M_type                       (),
        M_multiscale                 (),
        M_couplingVariables          (),
        M_couplingResiduals          (),
        M_comm                       (),
        M_displayer                  (),
        M_SubiterationsMaximumNumber (),
        M_Tolerance                  ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::MS_Algorithm() \n";
#endif

}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MS_Algorithm::SetupData( const std::string& FileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::SetupData( DataFile ) \n";
#endif

    GetPot DataFile( FileName );

    M_SubiterationsMaximumNumber = DataFile( "Solver/Algorithm/subITMax", 100 );
    M_Tolerance                  = DataFile( "Solver/Algorithm/tolerance", 1.e-10 );
}

void
MS_Algorithm::SubIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::SubIterate( tolerance, subITMax ) \n";
#endif

    // Algorithm Type
    if ( M_displayer->isLeader() )
        std::cout << " MS-  " << Enum2String( M_type, MS_algorithmsMap ) << " Algorithm" << std::endl;
}

void
MS_Algorithm::UpdateCouplingVariables()
{
    // The default approach is to extrapolate the next coupling variables
    M_multiscale->ExtrapolateCouplingVariables();
}

void
MS_Algorithm::ShowMe()
{
    std::cout << "=================== Algorithm Information ===================" << std::endl << std::endl;

    std::cout << "Algorithm type      = " << Enum2String( M_type, MS_algorithmsMap ) << std::endl
              << "Max Sub-iterations  = " << M_SubiterationsMaximumNumber << std::endl
              << "Tolerance           = " << M_Tolerance << std::endl << std::endl;
    std::cout << std::endl << std::endl;
}

// ===================================================
// Methods
// ===================================================
void
MS_Algorithm::InitializeCouplingVariables()
{
    M_multiscale->InitializeCouplingVariables();
}

Real
MS_Algorithm::Residual() const
{
    // Compute residual
    M_multiscale->ExportCouplingResiduals( *M_couplingResiduals );
    return M_couplingResiduals->Norm2();
}

// ===================================================
// Set Methods
// ===================================================
void
MS_Algorithm::SetCommunicator( const MS_Comm_PtrType& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MS_Algorithm::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

void
MS_Algorithm::SetModel( const MS_Model_PtrType model )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::SetMultiScaleProblem( multiscale ) \n";
#endif

    M_multiscale = boost::dynamic_pointer_cast< MS_Model_MultiScale >( model );

    // Build coupling variables and residuals vectors
    std::vector<int> MyGlobalElements(0);
    EpetraMap couplingMap( -1, static_cast<int> ( MyGlobalElements.size() ), &MyGlobalElements[0], 0, M_comm );
    M_multiscale->CreateCouplingMap( couplingMap );

    M_couplingVariables.reset( new EpetraVector( couplingMap, Unique ) );
    M_couplingResiduals.reset( new EpetraVector( couplingMap, Unique ) );
}

// ===================================================
// Get Methods
// ===================================================
const algorithmsTypes&
MS_Algorithm::GetType() const
{
    return M_type;
}

const boost::shared_ptr< MS_Model_MultiScale >
MS_Algorithm::GetMultiScaleProblem() const
{
    return M_multiscale;
}

const MS_Vector_PtrType
MS_Algorithm::GetCouplingVariables() const
{
    return M_couplingVariables;
}

const MS_Vector_PtrType
MS_Algorithm::GetCouplingResiduals() const
{
    return M_couplingResiduals;
}

const MS_Comm_PtrType
MS_Algorithm::GetCommunicator() const
{
    return M_comm;
}

const UInt&
MS_Algorithm::GetSubiterationsMaximumNumber() const
{
    return M_SubiterationsMaximumNumber;
}

const Real&
MS_Algorithm::GetTolerance() const
{
    return M_Tolerance;
}

// ===================================================
// Protected Methods
// ===================================================
void
MS_Algorithm::Save( const UInt& SubiterationsNumber, const Real& residual )
{
    std::ofstream output;
    output << std::scientific << std::setprecision( 15 );

    if ( M_comm->MyPID() == 0 )
    {
        std::string filename = MS_ProblemFolder + "Step_" + number2string( MS_ProblemStep ) + "_Algorithm.mfile";

        if ( M_multiscale->GetGlobalData()->GetDataTime()->isFirstTimeStep() )
        {
            output.open( filename.c_str(), std::ios::trunc );
            output << "% Algorithm Type: " << Enum2String( M_type, MS_algorithmsMap ) << std::endl;
            output << "% Subiteration maximum number: " << M_SubiterationsMaximumNumber << std::endl;
            output << "% Tolerance: " << M_Tolerance << std::endl << std::endl;
            output << "% TIME                     Subiterations      Residual" << std::endl;
        }
        else
            output.open( filename.c_str(), std::ios::app );

        output << M_multiscale->GetGlobalData()->GetDataTime()->getTime()
        << "      "<< SubiterationsNumber
        << "                  "<< residual << std::endl;

        output.close();
    }
}

bool
MS_Algorithm::ToleranceSatisfied()
{
    // Compute residual
    Real residual ( Residual() );

    // Display residual value
    if ( M_displayer->isLeader() )
        std::cout << " MS-  Residual:                                " << residual << std::endl;

    // Is the tolerance satisfied?
    if ( residual <= M_Tolerance )
        return true;
    else
        return false;
}

} // Namespace LifeV
