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

\brief Main program for solving he Ossen's equation with reaction using
P1/P1 or P2/P2 FEM with Interior Penalty ( IP ) stabilization
( See Burman-Fernandez-Hansbo 2004 ).

*/

#include <lifeV.hpp>
#include <NavierStokesSolverIP.hpp>
#include <ud_functions.hpp>
//#include <ethierSteinman.hpp>
//#include <simple.hpp>
#include <picard.hpp>

int main( int argc, char** argv )
{
    using namespace LifeV;
    //using namespace std;

    // Reading from data file
    //
    GetPot commandLine( argc, argv );
    const char* dataFileName = commandLine.follow( "data", 2, "-f", "--file" );
    GetPot dataFile( dataFileName );

    // Number of boundary conditions for the fluid velocity,
    // solid displacement, and fluid mesh motion
    //
    BCHandler bcH( 1 );

    NavierStokesSolverIP< RegionMesh3D<LinearTetra> >
        fluid( dataFile, feTetraP1, quadRuleTetra4pt, quadRuleTria3pt, bcH );
    //fluid.showMe();

    // Boundary conditions for the fluid velocity
    //vector<ID> icomp( 1 );
    //icomp[0] = 1;
    BCFunctionBase uWall( uexact );
    // BCFunctionBase bcf( fZero );
    bcH.addBC( "Wall", 1, Essential, Full, uWall, 3 );
    //bcH.addBC( "Wall", 2, Essential, Component, bcf, icomp );
    //BCFunctionBase bcf( afZero );
    //BCFunctionBase inFlow( u2 );
    //bcH.addBC( "Wall", 1, Essential, Full, bcf, 3 );
    //bcH.addBC( "InFlow", 2, Natural, Full, inFlow, 3 );
    //bcH.addBC( "Edges", 20, Essential, Full, bcf, 3 );

    //bcH.showMe();

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
        fluid.initialize( xexact );
        fluid.timeAdvance( f, 0.0 );

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
                      << fluid.pErrorL2(pexact, 0.) << std::endl;
            std::cout << "      - L2 velocity error = "
                      << fluid.uErrorL2(uexact, 0.) << std::endl;
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
        fluid.initialize( xexact, t0, dt );

        // Temporal loop

        for ( Real time = t0+dt ; time <= tFinal; time+=dt )
        {
            fluid.timeAdvance( f, time );
            fluid.iterate( time );

            std::cout << "      - L2 pressure error = "
                      << fluid.pErrorL2(pexact, time) << std::endl;
            std::cout << "      - L2 velocity error = "
                      << fluid.uErrorL2(uexact, time) << std::endl;

            // save result on file
            std::ostringstream indexout;
            indexout << ( time*100 );
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
        }

    }

    return 0;

}
