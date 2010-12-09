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
 *  @brief File containing the MultiScale Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
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
        M_subiterationsMaximumNumber (),
        M_tolerance                  ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::MS_Algorithm() \n";
#endif

}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MS_Algorithm::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::SetupData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    M_subiterationsMaximumNumber = dataFile( "Solver/Algorithm/subITMax", 100 );
    M_tolerance                  = dataFile( "Solver/Algorithm/tolerance", 1.e-10 );
}

void
MS_Algorithm::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::SubIterate() \n";
#endif

    // Algorithm Type
    if ( M_displayer->isLeader() )
        std::cout << " MS-  " << Enum2String( M_type, MS_algorithmsMap ) << " Algorithm" << std::endl;
}

void
MS_Algorithm::updateCouplingVariables()
{
    // The default approach is to extrapolate the next coupling variables
    M_multiscale->ExtrapolateCouplingVariables();
}

void
MS_Algorithm::showMe()
{
    std::cout << "=================== Algorithm Information ===================" << std::endl << std::endl;

    std::cout << "Algorithm type      = " << Enum2String( M_type, MS_algorithmsMap ) << std::endl
              << "Max Sub-iterations  = " << M_subiterationsMaximumNumber << std::endl
              << "Tolerance           = " << M_tolerance << std::endl << std::endl;
    std::cout << std::endl << std::endl;
}

// ===================================================
// Methods
// ===================================================
void
MS_Algorithm::initializeCouplingVariables()
{
    M_multiscale->InitializeCouplingVariables();
}

Real
MS_Algorithm::computeResidual() const
{
    // Compute computeResidual
    M_multiscale->ExportCouplingResiduals( *M_couplingResiduals );
    return M_couplingResiduals->Norm2();
}

// ===================================================
// Set Methods
// ===================================================
void
MS_Algorithm::setCommunicator( const MS_Comm_PtrType& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MS_Algorithm::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

void
MS_Algorithm::setModel( const MS_Model_PtrType model )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MS_Algorithm::SetMultiScaleProblem( multiscale ) \n";
#endif

    M_multiscale = boost::dynamic_pointer_cast< MS_Model_MultiScale >( model );

    // Build coupling variables and residuals vectors
    std::vector<Int> myGlobalElements(0);
    EpetraMap couplingMap( -1, static_cast<Int> ( myGlobalElements.size() ), &myGlobalElements[0], 0, M_comm );
    M_multiscale->CreateCouplingMap( couplingMap );

    M_couplingVariables.reset( new EpetraVector( couplingMap, Unique ) );
    M_couplingResiduals.reset( new EpetraVector( couplingMap, Unique ) );
}

// ===================================================
// Get Methods
// ===================================================
const algorithmsTypes&
MS_Algorithm::type() const
{
    return M_type;
}

const MS_Algorithm::multiscaleModelPtr_Type
MS_Algorithm::multiScaleProblem() const
{
    return M_multiscale;
}

const MS_Vector_PtrType
MS_Algorithm::couplingVariables() const
{
    return M_couplingVariables;
}

const MS_Vector_PtrType
MS_Algorithm::couplingResiduals() const
{
    return M_couplingResiduals;
}

const MS_Comm_PtrType
MS_Algorithm::communicator() const
{
    return M_comm;
}

const UInt&
MS_Algorithm::subiterationsMaximumNumber() const
{
    return M_subiterationsMaximumNumber;
}

const Real&
MS_Algorithm::tolerance() const
{
    return M_tolerance;
}

// ===================================================
// Protected Methods
// ===================================================
void
MS_Algorithm::save( const UInt& subiterationsNumber, const Real& residual )
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
            output << "% Subiteration maximum number: " << M_subiterationsMaximumNumber << std::endl;
            output << "% Tolerance: " << M_tolerance << std::endl << std::endl;
            output << "% TIME                     Subiterations      Residual" << std::endl;
        }
        else
            output.open( filename.c_str(), std::ios::app );

        output << M_multiscale->GetGlobalData()->GetDataTime()->getTime() << "      "
               << subiterationsNumber << "                  " << residual << std::endl;

        output.close();
    }
}

bool
MS_Algorithm::toleranceSatisfied()
{
    // Compute computeResidual
    Real residual ( computeResidual() );

    // Display computeResidual value
    if ( M_displayer->isLeader() )
        std::cout << " MS-  Residual:                                " << residual << std::endl;

    // Is the tolerance satisfied?
    if ( residual <= M_tolerance )
        return true;
    else
        return false;
}

} // Namespace LifeV
