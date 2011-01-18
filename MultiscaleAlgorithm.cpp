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
 *  @brief File containing the Multiscale Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleAlgorithm.hpp>

namespace LifeV
{
namespace Multiscale
{

std::map< std::string, algorithms_Type > multiscaleAlgorithmsMap;

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithm::MultiscaleAlgorithm() :
        M_type                       (),
        M_name                       (),
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
// Multiscale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithm::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MultiscaleAlgorithm::setupData( fileName ) \n";
#endif

    GetPot dataFile( fileName );
    M_name = dataFile( "Multiscale/algorithmName", "algorithmName" );
}

void
MultiscaleAlgorithm::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MultiscaleAlgorithm::subIterate() \n";
#endif

    // Algorithm Type
    if ( M_displayer->isLeader() )
        std::cout << " MS-  " << enum2String( M_type, multiscaleAlgorithmsMap ) << " Algorithm" << std::endl;
}

void
MultiscaleAlgorithm::showMe()
{
    std::cout << "=================== Algorithm Information ===================" << std::endl << std::endl;

    std::cout << "Algorithm type      = " << enum2String( M_type, multiscaleAlgorithmsMap ) << std::endl
              << "Algorithm name      = " << M_name << std::endl
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
    Debug( 8100 ) << "MultiscaleAlgorithm::setCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

void
MultiscaleAlgorithm::setModel( const multiscaleModelPtr_Type model )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8010 ) << "MultiscaleAlgorithm::setMultiscaleProblem( multiscale ) \n";
#endif

    M_multiscale = boost::dynamic_pointer_cast< MultiscaleModelMultiscale >( model );

    // Build coupling variables and residuals vectors
    std::vector<Int> myGlobalElements(0);
    MapEpetra couplingMap( -1, static_cast<Int> ( myGlobalElements.size() ), &myGlobalElements[0],  M_comm );
    M_multiscale->createCouplingMap( couplingMap );

    M_couplingVariables.reset( new VectorEpetra( couplingMap, Unique ) );
    M_couplingResiduals.reset( new VectorEpetra( couplingMap, Unique ) );
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
            output << "% Algorithm Type: " << enum2String( M_type, multiscaleAlgorithmsMap ) << std::endl;
            output << "% Subiteration maximum number: " << M_subiterationsMaximumNumber << std::endl;
            output << "% Tolerance: " << M_tolerance << std::endl << std::endl;
            output << "% TIME                     Subiterations      Residual" << std::endl;
        }
        else
            output.open( filename.c_str(), std::ios::app );

        output << M_multiscale->globalData()->dataTime()->time() << "      "
               << subiterationsNumber << "                  " << residual << std::endl;

        output.close();
    }
}

bool
MultiscaleAlgorithm::checkResidual( const UInt& subIT )
{
    // Compute computeResidual
    Real residual ( computeResidual() );

    // Display subIT and residual values
    if ( M_displayer->isLeader() )
    {
        if ( subIT > 0 )
            std::cout << " MS-  Sub-iteration n.:                        " << subIT << std::endl;
        std::cout << " MS-  Residual:                                " << residual << std::endl;
    }
    // Is the tolerance satisfied?
    if ( residual <= M_tolerance )
    {
        save( subIT, M_couplingResiduals->norm2() );
        return true;
    }
    else
        return false;
}

} // Namespace Multiscale
} // Namespace LifeV
