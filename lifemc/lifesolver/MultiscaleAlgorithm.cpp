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

#include <lifemc/lifesolver/MultiscaleAlgorithm.hpp>

namespace LifeV
{
namespace multiscale
{

std::map< std::string, algorithms_Type > multiscaleAlgorithmsMap;

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithm::MultiscaleAlgorithm() :
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
    Debug( 8010 ) << "MultiscaleAlgorithm::MultiscaleAlgorithm() \n";
#endif

}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithm::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MultiscaleAlgorithm::SetupData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    M_subiterationsMaximumNumber = dataFile( "Solver/Algorithm/subITMax", 100 );
    M_tolerance                  = dataFile( "Solver/Algorithm/tolerance", 1.e-10 );
}

void
MultiscaleAlgorithm::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MultiscaleAlgorithm::SubIterate() \n";
#endif

    // Algorithm Type
    if ( M_displayer->isLeader() )
        std::cout << " MS-  " << Enum2String( M_type, multiscaleAlgorithmsMap ) << " Algorithm" << std::endl;
}

void
MultiscaleAlgorithm::showMe()
{
    std::cout << "=================== Algorithm Information ===================" << std::endl << std::endl;

    std::cout << "Algorithm type      = " << Enum2String( M_type, multiscaleAlgorithmsMap ) << std::endl
              << "Max Sub-iterations  = " << M_subiterationsMaximumNumber << std::endl
              << "Tolerance           = " << M_tolerance << std::endl << std::endl;
    std::cout << std::endl << std::endl;
}

// ===================================================
// Methods
// ===================================================
Real
MultiscaleAlgorithm::computeResidual() const
{
    // Compute computeResidual
    M_multiscale->exportCouplingResiduals( *M_couplingResiduals );
    return M_couplingResiduals->norm2();
}

// ===================================================
// Set Methods
// ===================================================
void
MultiscaleAlgorithm::setCommunicator( const multiscaleCommPtr_Type& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MultiscaleAlgorithm::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

void
MultiscaleAlgorithm::setModel( const multiscaleModelPtr_Type model )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MultiscaleAlgorithm::SetMultiScaleProblem( multiscale ) \n";
#endif

    M_multiscale = boost::dynamic_pointer_cast< MultiscaleModelMultiscale >( model );

    // Build coupling variables and residuals vectors
    std::vector<Int> myGlobalElements(0);
    EpetraMap couplingMap( -1, static_cast<Int> ( myGlobalElements.size() ), &myGlobalElements[0], 0, M_comm );
    M_multiscale->createCouplingMap( couplingMap );

    M_couplingVariables.reset( new EpetraVector( couplingMap, Unique ) );
    M_couplingResiduals.reset( new EpetraVector( couplingMap, Unique ) );
}

// ===================================================
// Protected Methods
// ===================================================
void
MultiscaleAlgorithm::save( const UInt& subiterationsNumber, const Real& residual )
{
    std::ofstream output;
    output << std::scientific << std::setprecision( 15 );

    if ( M_comm->MyPID() == 0 )
    {
        std::string filename = multiscaleProblemFolder + "Step_" + number2string( multiscaleProblemStep ) + "_Algorithm.mfile";

        if ( M_multiscale->globalData()->dataTime()->isFirstTimeStep() )
        {
            output.open( filename.c_str(), std::ios::trunc );
            output << "% Algorithm Type: " << Enum2String( M_type, multiscaleAlgorithmsMap ) << std::endl;
            output << "% Subiteration maximum number: " << M_subiterationsMaximumNumber << std::endl;
            output << "% Tolerance: " << M_tolerance << std::endl << std::endl;
            output << "% TIME                     Subiterations      Residual" << std::endl;
        }
        else
            output.open( filename.c_str(), std::ios::app );

        output << M_multiscale->globalData()->dataTime()->getTime() << "      "
               << subiterationsNumber << "                  " << residual << std::endl;

        output.close();
    }
}

bool
MultiscaleAlgorithm::toleranceSatisfied()
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

} // Namespace multiscale
} // Namespace LifeV
