/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright ( C ) 2001, 2002, 2003, 2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or ( at your option ) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/*!

\author M.A. Fernandez
\date 01/01/2004

\brief Main program for solving he Oseen's equation with reaction using
P1/P1 or P2/P2 FEM with Interior Penalty (IP) stabilization
(See Burman-Fernandez-Hansbo 2004).

*/

#include <life.hpp>
#include <NavierStokesSolverIP.hpp>
#include <ud_functions.hpp>
#include <ethierSteinman.hpp>
// #include <simple.hpp>
#include <picard.hpp>

#include <iostream>
#include <set>

using namespace LifeV;

std::set<UInt> parseList( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    UInt commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );
        setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
        stringList = stringList.substr( commaPos+1 );
    }
    setList.insert( atoi( stringList.c_str() ) );
    return setList;
}

int main( int argc, char** argv )
{
    // Reading from data file
    //
    GetPot commandLine( argc, argv );
    const char* dataFileName = commandLine.follow( "data", 2, "-f", "--file" );
    GetPot dataFile( dataFileName );

    // Problem definition
    typedef EthierSteinmanSteady Problem;
    Problem::setParamsFromGetPot( dataFile );

    // Boundary conditions
    std::string dirichletList = dataFile( "fluid/problem/dirichletList", "" );
    std::set<UInt> dirichletMarkers = parseList( dirichletList );
    std::string neumannList = dataFile( "fluid/problem/neumannList", "" );
    std::set<UInt> neumannMarkers = parseList( neumannList );

    BCHandler::BCHints hint = neumannMarkers.size() != 0 ?
        BCHandler::HINT_BC_NONE : BCHandler::HINT_BC_ONLY_ESSENTIAL;
    BCHandler bcH( 0, hint );
    BCFunctionBase uWall( Problem::uexact );
    BCFunctionBase uNeumann( Problem::fNeumann );

    for (std::set<UInt>::const_iterator it = dirichletMarkers.begin();
         it != dirichletMarkers.end(); ++it)
    {
        bcH.addBC( "Wall", *it, Essential, Full, uWall, 3 );
    }
    for (std::set<UInt>::const_iterator it = neumannMarkers.begin();
         it != neumannMarkers.end(); ++it)
    {
        bcH.addBC( "Flux", *it, Natural, Full, uNeumann, 3 );
    }

    //bcH.addBC( "Wall", 2, Essential, Component, bcf, icomp );
    //BCFunctionBase bcf( afZero );
    //BCFunctionBase inFlow( u2 );
    //bcH.addBC( "Wall", 1, Essential, Full, bcf, 3 );
    //bcH.addBC( "InFlow", 2, Natural, Full, inFlow, 3 );
    //bcH.addBC( "Edges", 20, Essential, Full, bcf, 3 );

    //bcH.showMe();

    // fluid solver
    const RefFE* refFE;
    const QuadRule* qR;
    const QuadRule* bdQr;

    if ( dataFile( "fluid/discretization/p_order", 1 ) == 2 )
    {
        // P2 only if specified 2
        refFE = &feTetraP2;
        qR = &quadRuleTetra15pt; // DoE 5
        bdQr = &quadRuleTria3pt; // DoE 1
    }
    else
    {
        // P1 by default
        refFE = &feTetraP1;
        qR = &quadRuleTetra4pt; // DoE 2
        bdQr = &quadRuleTria3pt; // DoE 2
    }

    NavierStokesSolverIP< RegionMesh3D<LinearTetra> >
        fluid( dataFile, *refFE, *qR, *bdQr, bcH );
    //fluid.showMe();

    // linearization of fluid
    if ( dataFile( "fluid/discretization/linearized", 0 ) )
    {
        fluid.linearize( Problem::uexact );
    }

    // steady or unsteady solve
    if ( dataFile( "fluid/miscellaneous/steady", 1 ) != 0 )
    {
        // Picard-Aitken iterations: steady version
        //

        Real abstol = dataFile( "fluid/picard/abstol", 1.e6 );
        Real reltol = dataFile( "fluid/picard/reltol", 0.0 );
        int maxiter = dataFile( "fluid/picard/maxiter", 300 );
        int method = dataFile( "fluid/picard/method", 1 ); // 1=Aitken
        Real omega  = dataFile( "fluid/picard/omega", 0 );

        UInt lVec = nDimensions * fluid.uDof().numTotalDof();

        Vector fx1( lVec );
        Vector fx0( lVec );
        Vector gx1( lVec );
        Vector gx0( lVec );
        Vector x1( lVec );
        Vector x0( lVec );

        fx1 = ZeroVector( lVec );
        fx0 = ZeroVector( lVec );
        gx1 = ZeroVector( lVec );
        gx0 = ZeroVector( lVec );
        x1  = ZeroVector( lVec );
        x0  = ZeroVector( lVec );

        // Compute right hand side
        fluid.initialize( Problem::xexact );
        fluid.timeAdvance( Problem::f, 0.0 );

        int  status = picard( &fluid, norm_inf_adaptor(), fx1, fx0, gx1, gx0,
                              x1, x0, abstol, reltol, maxiter, method, omega );

        if( status == 1 )
        {
            std::cout << "Inner iterations failed" << std::endl;
            exit( 1 );
        }
        else
        {
            // std::cout << "End of time "<< time << std::endl;
            std::cout << "Number of inner iterations       : " << maxiter
                      << std::endl;
            std::cout << "      - L2 pressure error = "
                      << fluid.pErrorL2( Problem::pexact, 0. )
                      << std::endl;
            std::cout << "      - L2 velocity error = "
                      << fluid.uErrorL2( Problem::uexact, 0. )
                      << std::endl;
            fluid.postProcess();
        }
    }
    else
    {
        // bdf: unsteady version

        // Initialization

        Real dt = fluid.timestep();
        Real t0 = fluid.inittime();
        Real tFinal = fluid.endtime();
        fluid.initialize( Problem::xexact, t0, dt );

        // Temporal loop

        for ( Real time = t0+dt ; time <= tFinal+dt/2; time+=dt )
        {
            fluid.timeAdvance( Problem::f, time );
            fluid.iterate( time );

            Real epr;
            Real epL2 = fluid.pErrorL2( Problem::pexact, time, &epr );
            Real eur;
            Real euL2 = fluid.uErrorL2( Problem::uexact, time, &eur );

            std::cout << "      - L2 pressure error = "
                      << epL2 << std::endl;
            std::cout << "      - L2 velocity error = "
                      << euL2 << std::endl;
            std::cout << "      - L2 p error (rel.) = "
                      << epr << std::endl;
            std::cout << "      - L2 u error (rel.) = "
                      << eur << std::endl;

            // save result on file
            std::ostringstream indexout;
            indexout << ( (int)(time/dt+0.5) );
            std::string voutname;
            voutname = "fluid.res"+indexout.str();
            std::fstream resFile( voutname.c_str(),
                                  std::ios::out | std::ios::binary );
            resFile.write( ( char* )&fluid.u()( 1 ),
                           fluid.u().size()*sizeof( double ) );
            resFile.write( ( char* )&fluid.p()( 1 ),
                           fluid.p().size()*sizeof( double ) );
            resFile.close();

            fluid.postProcess();

            //// save exact solution and error
            Vector uExactVec( fluid.u().size() );

            fluid.uInterpolate( Problem::pexact, uExactVec, time );
            Vector pExactVec( fluid.p().size() );
            for( UInt i=0; i<pExactVec.size(); ++i )
            {
                pExactVec[ i ] = uExactVec[ i ];
            }
            fluid.removeMean( pExactVec, 1 );

            fluid.uInterpolate( Problem::uexact, uExactVec, time );

            wr_gmv_ascii( "test_exact_" + indexout.str() + ".inp",
                          fluid.mesh(),
                          fluid.uDof().numTotalDof(),
                          &(*uExactVec.begin()), &(*pExactVec.begin()) );

            uExactVec -= fluid.u();
            uExactVec *= -1;
            pExactVec -= fluid.p();
            pExactVec *= -1;
            wr_gmv_ascii( "test_error_" + indexout.str() + ".inp",
                          fluid.mesh(),
                          fluid.uDof().numTotalDof(),
                          &(*uExactVec.begin()), &(*pExactVec.begin()) );
        } // temporal loop

    } // steady or unsteady

    return 0;

}
