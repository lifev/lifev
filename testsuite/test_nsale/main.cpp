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

\brief Fluid code to be coupled with an external solid solver using Jean-Fred 
pvm coupler (masterFSI). PVM required. Tested with SD (P1/P1 and Q1/Q1) or IP stabilization (P1/P1).

*/

#include <life/lifecore/life.hpp>
#include <life/lifefem/regionMesh3D_ALE.hpp>
#include <life/lifesolver/NavierStokesAleSolver.hpp>
#include <ud_functions.hpp>
#include <life/lifefilters/medit.hpp>
#include <iostream>
#include <set>

#include <life/lifesolver/fluidToMaster.hpp>
#include <pvm3.h>

using namespace LifeV;

int main( int argc, char** argv )
{
    // Reading from data file
    //
    GetPot commandLine( argc, argv );
    const char* dataFileName = commandLine.follow( "data", 2, "-f", "--file" );
    GetPot dataFile( dataFileName );

    // recover FSI test case
    UInt FSItest= dataFile( "FSItest", 0);
    UInt FSIalgo= dataFile( "FSIalgo", 0);

    BCHandler bcHu, bcHmesh; 
    BCHandler bcHdu, bcHdp, bcHdz; 


    //bcH.showMe();

    // fluid solver

    //
    //     NavierStokesAleSolver< RegionMesh3D_ALE<LinearHexa> >
    //fluid( dataFile, feHexaQ1, quadRuleHexa8pt, quadRuleQuad4pt, bcHu, bcHmesh );
   
    NavierStokesAleSolver< RegionMesh3D_ALE<LinearTetra> >
      fluid( dataFile, feTetraP1, quadRuleTetra4pt, quadRuleTria3pt, bcHu, bcHmesh );



    fluid.showMe();

    int masterId = pvm_parent();
    
    if( masterId == PvmNoParent )
      {
	cout << " I am a poor lonesome job\n ";
	exit(1);
      } 
    else if( masterId == 0 )
      {
	cout << "PVM is not ok\n";
	exit(1);
      }
    
    
    FluidToMaster< NavierStokesAleSolver< RegionMesh3D_ALE<LinearTetra> > >  fluidToMaster(fluid,masterId,-1);
		      
    std::cout <<"%%%%%%%%%%% \n";
    std::cout << " I am the slave " << pvm_mytid()<< endl;
    std::cout << " My master is " << masterId << endl;


    Medit medit(fluid);

    
    BCFunctionBase bcf0( fZero );

    UInt dim_fluid = fluid.uDof().numTotalDof();
    
    Vector disp_interf(3*dim_fluid); 
    Vector minus_rho_da_interf(3*dim_fluid); 
    disp_interf = ZeroVector( disp_interf.size() ); 
    minus_rho_da_interf = ZeroVector( minus_rho_da_interf.size() ); 

    BCVector dispbd(  disp_interf      , dim_fluid); 
    BCVector u_wall(  fluid.w()        , dim_fluid);
    BCVector du_wall( fluid.dw()       , dim_fluid);
    BCVector da_wall( minus_rho_da_interf, dim_fluid, 2); 
    BCVector dg_wall( fluid.dp(),          dim_fluid, 1);
		
    BCFunctionBase inFlow1(g_aneurism); 
    BCFunctionBase inFlow2(g_caronew); 
    BCFunctionBase inFlow3(g_branch);
    BCFunctionBase inFlow(g);
    BCFunctionBase inFlow7(u_in);
    vector<ID> comp(1);
    comp[0]=3;
 
    switch( FSItest ) 
      {
      case 1:
	//
	// Aneurism test case
	//
	
	std::cout << "\n\n";
	std::cout << " ================--------------================\n";
	std::cout << " ================ Test Case 1  ================\n";
   	std::cout << " ================--------------================\n\n";

	bcHu.addBC("Wall",    4,  Essential, Full, u_wall,   3);
	bcHu.addBC("Wall",    5,  Essential, Full, u_wall,   3); 
	bcHu.addBC("Wall",    6,  Essential, Full, u_wall,   3);   
	bcHu.addBC("Wall",    7,  Essential, Full, u_wall,   3);
	bcHu.addBC("OutFlow", 1,  Natural,   Full, bcf0,     3);
	bcHu.addBC("OutFlow", 2,  Natural,   Full, bcf0,     3);
	bcHu.addBC("InFlow",  3,  Natural,   Full, inFlow1,  3);
	
	bcHdu.addBC("Wall",    4,  Essential, Full, du_wall, 3);
	bcHdu.addBC("Wall",    5,  Essential, Full, du_wall, 3); 
	bcHdu.addBC("Wall",    6,  Essential, Full, du_wall, 3);   
	bcHdu.addBC("Wall",    7,  Essential, Full, du_wall, 3);
	bcHdu.addBC("OutFlow", 1,  Natural,   Full, bcf0,    3);  
	bcHdu.addBC("OutFlow", 2,  Natural,   Full, bcf0,    3);
	bcHdu.addBC("InFlow",  3,  Natural,   Full, bcf0,    3);

	bcHdp.addBC("Wall",    4,  Natural,   Scalar, da_wall);
	bcHdp.addBC("Wall",    5,  Essential, Scalar, bcf0); 
	bcHdp.addBC("Wall",    6,  Essential, Scalar, bcf0);   
	bcHdp.addBC("Wall",    7,  Essential, Scalar, bcf0);
	bcHdp.addBC("OutFlow", 1,  Essential, Scalar, bcf0);
	bcHdp.addBC("OutFlow", 2,  Essential, Scalar, bcf0);
	bcHdp.addBC("InFlow",  3,  Essential, Scalar, bcf0);

	bcHmesh.addBC("Wall",   4, Essential, Full, dispbd,  3);  //  displacement given on the lateral points
	bcHmesh.addBC("Wall",   5, Essential, Full, dispbd,  3); 
	bcHmesh.addBC("Wall",   6, Essential, Full, dispbd,  3); 
	bcHmesh.addBC("Wall",   7, Essential, Full, dispbd,  3); 
	bcHmesh.addBC("InLet",  3, Essential, Full, bcf0,    3);
	bcHmesh.addBC("OutLet", 1, Essential, Full, bcf0,    3);
	bcHmesh.addBC("OutLet", 2, Essential, Full, bcf0,    3);
	
	bcHdz.addBC("Interface", 4, Natural, Full, dg_wall, 3); 

	break;

      case 2:
	//
	// Carotid test case
	//

	std::cout << "\n\n";
	std::cout << " ================--------------================\n";
	std::cout << " ================ Test Case 2  ================\n";
   	std::cout << " ================--------------================\n\n";
	
	bcHu.addBC("Wall",      4,  Essential, Full, u_wall, 3);
	bcHu.addBC("Wall",      5,  Essential, Full, u_wall, 3); 
	bcHu.addBC("Wall",      6,  Essential, Full, u_wall, 3);   
	bcHu.addBC("Wall",      7,  Essential, Full, u_wall, 3);
	bcHu.addBC("OutFlow",   2,  Natural,   Full, bcf0,   3);
	bcHu.addBC("OutFlow",   3,  Natural,   Full, bcf0,   3);
	bcHu.addBC("InFlow",    1,  Natural,   Full, inFlow, 3);
	
	bcHdu.addBC("Wall",     4,  Essential, Full, du_wall, 3);
	bcHdu.addBC("Wall",     5,  Essential, Full, du_wall, 3); 
	bcHdu.addBC("Wall",     6,  Essential, Full, du_wall, 3);   
	bcHdu.addBC("Wall",     7,  Essential, Full, du_wall, 3);
	bcHdu.addBC("OutFlow",  2,  Natural,   Full, bcf0,    3);  
	bcHdu.addBC("OutFlow",  3,  Natural,   Full, bcf0,    3);
	bcHdu.addBC("InFlow",   1,  Natural,   Full, bcf0,    3);

	bcHdp.addBC("Wall",     4,  Natural, Scalar, da_wall);
	bcHdp.addBC("Wall",     5,  Essential, Scalar, bcf0);  
	bcHdp.addBC("Wall",     6,  Essential, Scalar, bcf0);  
	bcHdp.addBC("Wall",     7,  Essential, Scalar, bcf0);  
	bcHdp.addBC("OutFlow",  2,  Essential, Scalar, bcf0); 
	bcHdp.addBC("OutFlow",  3,  Essential, Scalar, bcf0); 
	bcHdp.addBC("InFlow",   1,  Essential, Scalar, bcf0); 

	bcHmesh.addBC("Wall",   4, Essential, Full, dispbd, 3);  //  displacement given on the lateral points
	bcHmesh.addBC("Wall",   5, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("Wall",   6, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("Wall",   7, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("InLet",  1, Essential, Full, bcf0,   3);
	bcHmesh.addBC("OutLet", 2, Essential, Full, bcf0,   3);
	bcHmesh.addBC("OutLet", 3, Essential, Full, bcf0,   3);
	
	bcHdz.addBC("Interface", 4, Natural, Full, dg_wall, 3); 

	break;

      case 3:
	//
	// Straight cylinder test case
	//	
	
	std::cout << "\n\n";
	std::cout << " ================--------------================\n";
	std::cout << " ================ Test Case 3  ================\n";
   	std::cout << " ================--------------================\n\n";


	bcHu.addBC("Wall",      2,  Essential, Full, u_wall, 3);
	bcHu.addBC("Wall",      4,  Essential, Full, u_wall, 3); 
	bcHu.addBC("Wall",      5,  Essential, Full, u_wall, 3);   
	bcHu.addBC("OutFlow",   3,  Natural,   Full, bcf0,   3);
	bcHu.addBC("InFlow",    1,  Natural,   Full, inFlow, 3);
	
	bcHdu.addBC("Wall",     2,  Essential, Full, du_wall, 3);
	bcHdu.addBC("Wall",     4,  Essential, Full, du_wall, 3); 
	bcHdu.addBC("Wall",     5,  Essential, Full, du_wall, 3);   
	bcHdu.addBC("OutFlow",  3,  Natural,   Full, bcf0,    3);
	bcHdu.addBC("InFlow",   1,  Natural,   Full, bcf0,    3);

	bcHdp.addBC("Wall",     2,  Natural,   Scalar, da_wall);
	bcHdp.addBC("Wall",     4,  Essential, Scalar, bcf0); 
	bcHdp.addBC("Wall",     5,  Essential, Scalar, bcf0); 
	bcHdp.addBC("OutFlow",  3,  Essential, Scalar, bcf0); 
	bcHdp.addBC("InFlow",   1,  Essential, Scalar, bcf0); 
		
	bcHmesh.addBC("Wall",   2, Essential, Full, dispbd, 3);  //  displacement given on the lateral points
	bcHmesh.addBC("Wall",   4, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("Wall",   5, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("InLet",  1, Essential, Full, bcf0,   3);
	bcHmesh.addBC("OutLet", 3, Essential, Full, bcf0,   3);
	
	bcHdz.addBC("Interface", 2, Natural, Full, dg_wall, 3); 
	
	break;

      case 4:

	bcHu.addBC("Wall",            4,  Essential, Full, u_wall,   3);
	bcHu.addBC("Wall_Edges",    15,  Essential, Full, u_wall,   3); 
	bcHu.addBC("Wall_Edges",     9,  Essential, Full, u_wall,   3); 
	bcHu.addBC("Wall_Edges",    16,  Essential, Full, u_wall,   3);
 	bcHu.addBC("Wall_Edges",     6,  Essential, Full, u_wall,   3);
	bcHu.addBC("OutFlow",        11,  Natural,   Full, bcf0,     3);
	bcHu.addBC("OutFlow",        14,  Natural,   Full, bcf0,     3);
	bcHu.addBC("OutFlow",         8,  Natural,   Full, bcf0,     3);
	bcHu.addBC("InFlow",         17,  Natural,   Full, inFlow2,  3);
      
	bcHdu.addBC("Wall",            4,  Essential, Full, du_wall,   3);
	bcHdu.addBC("Wall_Edges",    15,  Essential, Full, du_wall,   3); 
	bcHdu.addBC("Wall_Edges",     9,  Essential, Full, du_wall,   3); 
	bcHdu.addBC("Wall_Edges",    16,  Essential, Full, du_wall,   3);
 	bcHdu.addBC("Wall_Edges",     6,  Essential, Full, du_wall,   3);
	bcHdu.addBC("OutFlow",        11,  Natural,   Full, bcf0,     3);
	bcHdu.addBC("OutFlow",        14,  Natural,   Full, bcf0,     3);
	bcHdu.addBC("OutFlow",         8,  Natural,   Full, bcf0,     3);
	bcHdu.addBC("InFlow",         17,  Natural,   Full, bcf0,     3);

	bcHdp.addBC("Wall",        4,  Natural,   Scalar,  da_wall);
	bcHdp.addBC("InFlow",      17, Essential, Scalar,  bcf0);
	bcHdp.addBC("Wall_Edges",  15, Essential, Scalar,  bcf0);
	bcHdp.addBC("Wall_Edges",  9,  Essential, Scalar,  bcf0);
	bcHdp.addBC("Wall_Edges",  16, Essential, Scalar,  bcf0);
	bcHdp.addBC("Wall_Edges",  6,  Essential, Scalar,  bcf0);
	bcHdp.addBC("OutFlow",     11, Essential, Scalar,  bcf0);
	bcHdp.addBC("OutFlow",     14, Essential, Scalar,  bcf0);
	bcHdp.addBC("OutFlow",     8,  Essential, Scalar,  bcf0);


	bcHmesh.addBC("Interface", 4, Essential, Full, dispbd, 3);
	bcHmesh.addBC("Edges",    15, Essential, Full, dispbd, 3);
	bcHmesh.addBC("Edges",     9, Essential, Full, dispbd, 3);
	bcHmesh.addBC("Edges",    16, Essential, Full, dispbd, 3);
	bcHmesh.addBC("Edges",     6, Essential, Full, dispbd, 3);
	bcHmesh.addBC("Base",     17, Essential, Full, bcf0,  3);
	bcHmesh.addBC("Base",      8, Essential, Full, bcf0,  3);
	bcHmesh.addBC("Base",     11, Essential, Full, bcf0,  3);
	bcHmesh.addBC("Base",     14, Essential, Full, bcf0,  3);
	
	
	bcHdz.addBC("Interface", 4, Natural, Full, dg_wall, 3);
	
      
        break;	


      case 5:
	//
	// Straight cylinder test case
	//	
	
	std::cout << "\n\n";
	std::cout << " ================--------------================\n";
	std::cout << " ================ Test Case 4  ================\n";
   	std::cout << " ================--------------================\n\n";


	bcHu.addBC("Wall",      1,  Essential, Full, u_wall, 3);
	bcHu.addBC("Wall",      5,  Essential, Full, u_wall, 3); 
	bcHu.addBC("Wall",     10,  Essential, Full, u_wall, 3);   
	bcHu.addBC("OutFlow",   3,  Natural,   Full, bcf0,   3);
	bcHu.addBC("InFlow",    2,  Natural,   Full, inFlow, 3);
	
	bcHdu.addBC("Wall",     1,  Essential, Full, du_wall, 3);
	bcHdu.addBC("Wall",     5,  Essential, Full, du_wall, 3); 
	bcHdu.addBC("Wall",    10,  Essential, Full, du_wall, 3);   
	bcHdu.addBC("OutFlow",  3,  Natural,   Full, bcf0,    3);
	bcHdu.addBC("InFlow",   2,  Natural,   Full, bcf0,    3);

	bcHdp.addBC("Wall",     1,  Natural,   Scalar, da_wall);
	bcHdp.addBC("Wall",     5,  Essential, Scalar, bcf0); 
	bcHdp.addBC("Wall",    10,  Essential, Scalar, bcf0); 
	bcHdp.addBC("OutFlow",  3,  Essential, Scalar, bcf0); 
	bcHdp.addBC("InFlow",   2,  Essential, Scalar, bcf0); 
		
	bcHmesh.addBC("Wall",   1, Essential, Full, dispbd, 3);  //  displacement given on the lateral points
	bcHmesh.addBC("Wall",   5, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("Wall",  10, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("InLet",  2, Essential, Full, bcf0,   3);
	bcHmesh.addBC("OutLet", 3, Essential, Full, bcf0,   3);
	
	bcHdz.addBC("Interface", 1, Natural, Full, dg_wall, 3); 
	
	break;


      case 6:
	//
	// Aneurism test case
	//
	
	std::cout << "\n\n";
	std::cout << " ================--------------================\n";
	std::cout << " ================ Test Case 6  ================\n";
   	std::cout << " ================--------------================\n\n";

	bcHu.addBC("Wall",    9,  Essential, Full, u_wall,   3);
	bcHu.addBC("Wall",    8,  Essential, Full, u_wall,   3); 
	bcHu.addBC("Wall",   10,  Essential, Full, u_wall,   3); 
	bcHu.addBC("Wall",    4,  Essential, Full, bcf0,     3);	
	bcHu.addBC("Sym",     5,  Essential, Component, bcf0, comp);
	bcHu.addBC("OutFlow", 3,  Natural,   Full, bcf0,     3);
	bcHu.addBC("OutFlow", 2,  Natural,   Full, bcf0,     3);
	bcHu.addBC("InFlow",  1,  Natural,   Full, inFlow3,  3);
	
	bcHdu.addBC("Wall",    9,  Essential, Full, du_wall,   3);
	bcHdu.addBC("Wall",    8,  Essential, Full, du_wall,   3);
	bcHdu.addBC("Wall",   10,  Essential, Full, du_wall,   3); 
	bcHdu.addBC("Wall",    4,  Essential, Full, bcf0,     3);	
	bcHdu.addBC("Sym",     5,  Essential, Component, bcf0, comp);
	bcHdu.addBC("OutFlow", 3,  Natural,   Full, bcf0,     3);
	bcHdu.addBC("OutFlow", 2,  Natural,   Full, bcf0,     3);
	bcHdu.addBC("InFlow",  1,  Natural,   Full, bcf0,  3);

	bcHdp.addBC("Wall",    8,  Natural,   Scalar, da_wall);
	bcHdp.addBC("Wall",    9,  Essential, Scalar, bcf0); 
	bcHdp.addBC("Wall",   10,  Essential, Scalar, bcf0); 
	bcHdp.addBC("Wall",    4,  Essential, Scalar, bcf0);   
	bcHdp.addBC("OutFlow", 3,  Essential, Scalar, bcf0);
	bcHdp.addBC("OutFlow", 2,  Essential, Scalar, bcf0);
	bcHdp.addBC("InFlow",  1,  Essential, Scalar, bcf0);

	bcHmesh.addBC("Wall",   9, Essential, Full, dispbd,  3);  //  displacement given on the lateral points
	bcHmesh.addBC("Wall",  10, Essential, Full, dispbd,  3); 	
	bcHmesh.addBC("Wall",   8, Essential, Full, dispbd,  3); 
	bcHmesh.addBC("Wall",   5, Essential, Component, bcf0,  comp); 
	bcHmesh.addBC("Wall",   4, Essential, Full, bcf0,  3); 
	bcHmesh.addBC("InLet",  1, Essential, Full, bcf0,    3);
	bcHmesh.addBC("OutLet", 2, Essential, Full, bcf0,    3);
	bcHmesh.addBC("OutLet", 3, Essential, Full, bcf0,    3);
	
	bcHdz.addBC("Interface", 8, Natural, Full, dg_wall, 3); 

	break;

     case 7:
	//
	// Straight cylinder test case
	//	
	
	std::cout << "\n\n";
	std::cout << " ================--------------================\n";
	std::cout << " ================ Test Case 7  ================\n";
   	std::cout << " ================--------------================\n\n";


	bcHu.addBC("Wall",      4,  Essential, Full, u_wall, 3);
	bcHu.addBC("Wall",      5,  Essential, Full, u_wall, 3); 
	bcHu.addBC("Wall",      6,  Essential, Full, u_wall, 3);
   	bcHu.addBC("Rigid-Wall",      3,  Essential, Full, bcf0, 3);
	bcHu.addBC("OutFlow",   2,  Natural,   Full, bcf0,   3);
	bcHu.addBC("InFlow",    1,  Essential, Full, inFlow7,   3);
	
	bcHdu.addBC("Wall",     4,  Essential, Full, du_wall, 3);
	bcHdu.addBC("Wall",     5,  Essential, Full, du_wall, 3); 
	bcHdu.addBC("Wall",     6,  Essential, Full, du_wall, 3);	
	bcHdu.addBC("Rigid-Wall",3,  Essential, Full, bcf0, 3);   
	bcHdu.addBC("OutFlow",  2,  Natural,   Full, bcf0,    3);
	bcHdu.addBC("InFlow",   1,  Essential,   Full, bcf0,    3);

	bcHdp.addBC("Wall",     4,  Natural,   Scalar, da_wall);
	bcHdp.addBC("Wall",     5,  Essential, Scalar, bcf0); 
	bcHdp.addBC("Wall",     6,  Essential, Scalar, bcf0); 	
	bcHdp.addBC("Rigid-Wall",3,  Essential, Scalar, bcf0);   
	bcHdp.addBC("OutFlow",  2,  Essential, Scalar, bcf0); 
	bcHdp.addBC("InFlow",   1,  Essential, Scalar, bcf0); 
		
	bcHmesh.addBC("Wall",   4, Essential, Full, dispbd, 3);  //  displacement given on the lateral points
	bcHmesh.addBC("Wall",   5, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("Wall",   6, Essential, Full, dispbd, 3); 
	bcHmesh.addBC("Rigid-Wall",3,  Essential, Full, bcf0, 3); 
	bcHmesh.addBC("InLet",  1, Essential, Full, bcf0,   3);
	bcHmesh.addBC("OutLet", 2, Essential, Full, bcf0,   3);
	
	bcHdz.addBC("Interface", 4, Natural, Full, dg_wall, 3); 
	
	break;


      default:
	ERROR_MSG("This FSI test is not yet implemented");
      }
    
 
    fluidToMaster.rvFromMasterFSI(140);
    fluidToMaster.sdCoorToMasterFSI();

     
    // Initialization
    fluid.initialize( u0 );
      
    // Temporal loop
    
    Real time = 0.0;
    
    do  {
      
      fluidToMaster.rvFromMasterFSI(110);
         
      switch ( fluidToMaster.status() ) {
	
      case -1: // Normal exit
	cout << "Normal fluid exit\n" << flush;  
	medit.postProcessALE(fluid);
	exit(1); 
	break;
      case 1: // New time step

	medit.postProcessALE(fluid);
	// new time
	time+=fluid.timestep();

	// advancing fluid
	fluid.timeAdvance( fZero, time );
	
	// Recovering displacement at the interface
	// from pvm
	fluidToMaster.dep_fluid_interf(disp_interf); 	

	std::cout << "norm( analytic )  = " << norm_inf( disp_interf ) << std::endl;
    	
	// moving mesh and solving fluid
	fluid.updateMesh( time );
	fluid.iterate( time );

	// sending fluid force to master
	fluidToMaster.sdToMasterFSI();
	break;

      case 2: // Solving Tangent fluid
	switch ( FSIalgo ) 
	  {
	  case 1: // Newton: with shape derivative terms
	    fluidToMaster.dep_fluid_interf(disp_interf); 
	    std::cout << "norm( disp_interf )  = " << norm_inf( disp_interf ) << std::endl;
	    fluid.updateDispVelo();
	    fluid.iterateLin( time, bcHdu );
	    break;
	  case 2: // Quasi-Newton: without shape derivative terms (Linearized NS_
	    fluidToMaster.dep_fluid_interf(disp_interf); 
	    std::cout << "norm( disp_interf )  = " << norm_inf( disp_interf ) << std::endl;
	    fluid.updateDispVelo();
	    fluid.iterateLin( time, bcHdu, 0 );
	    break;
	  case 3: // Quasi-Newton: simplified fluid 
	    fluidToMaster.minus_rho_da_interf(minus_rho_da_interf); 
	    std::cout << "norm( minus_rho_da_interf )  = " << norm_inf( minus_rho_da_interf ) << std::endl;
	    fluid.iterateReducedLin( time, bcHdp, bcHdz );
	    break;
	  default:
	    ERROR_MSG("This FSI algorithm is not yet implemented");
	  }

	fluidToMaster.sdToMasterFSI();
	break; 
      case 0: // Inner fluid solver evaluation
	fluidToMaster.dep_fluid_interf(disp_interf); 
	std::cout << "norm( analytic )  = " << norm_inf( disp_interf ) << std::endl;
	fluid.updateMesh( time );
	fluid.iterate( time );
	fluidToMaster.sdToMasterFSI();
	break;
      default:
	cout << " This status number is not allowed.\n" << endl << flush;
	exit(1);	
      }    // close switch 
    }
    while ( 1 );  // close infinite do 
    return 0;
}
