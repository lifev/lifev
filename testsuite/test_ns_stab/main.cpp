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
\date 01/05/2005

\brief Flow around a cylinder using Navier-Stokes.
       The result file "cl" contains the 
       lift coefficient at each time, so 
       that the Strouhal number can be computed

*/

#include <life/lifecore/life.hpp>
#include <NavierStokesSolver.hpp>
#include <ud_functions.hpp>
#include <life/lifefilters/medit.hpp>

#include <iostream>
#include <set>

using namespace LifeV;



int main( int argc, char** argv )
{
    // Reading from data file
    //
    GetPot commandLine( argc, argv );
    const char* dataFileName = commandLine.follow( "data", 2, "-f", "--file" );
    GetPot dataFile( dataFileName );


    // boundary conditions
    BCHandler bcH;    
    BCFunctionBase bcf( fZero );
    BCFunctionBase inFlow( u_in );
    std::vector<ID> icomp(1);
    icomp[0]=3;
    bcH.addBC( "Sym", 20, Essential, Component, bcf, icomp );
    bcH.addBC( "Wall", 2, Essential, Full, bcf, 3 );
    bcH.addBC( "Wall", 1, Essential, Full, inFlow, 3 );
    

    // fluid solver
    NavierStokesSolver< RegionMesh3D<LinearTetra> >
      fluid( dataFile, feTetraP1, quadRuleTetra4pt, quadRuleTria3pt, bcH);
    
    fluid.showMe();

   // post-processing
    Medit medit(fluid);
    
    
    // Initialization
    
    Real dt = fluid.timestep();
    Real t0 = fluid.inittime();
    Real tFinal = fluid.endtime();
    fluid.initialize( u0, t0, dt );
    //fluid.initialize("sol0");
    
    std::ofstream out_cl("cl");
    
    // Temporal loop
    
    for ( Real time = t0+dt ; time <= tFinal+dt/2; time+=dt )
      {
	fluid.timeAdvance( fZero, time );
	fluid.iterate( time );
	medit.postProcess(fluid);
	
	out_cl << time << " " << fluid.liftCoeff(2) << std::endl;

      } // temporal loop



    return 0;

}
