/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

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
\date 10/03/2005

\brief Fluid code to be coupled with an external solid solver using Jean-Fred 
pvm coupler (masterFSI), using the semi-implicit scheme introduced in  
INRIA Research Report 5700, "A projection semi-implicit scheme for the coupling of 
an elastic structure with an incompressible fluid", by M.A. Fernández, J.-F. Gebeau and C. Grandmont.
PVM required. Tested with SD (P1/P1 and Q1/Q1) or IP stabilization (P1/P1).

*/

#include <life/lifecore/life.hpp>
#include <life/lifefem/regionMesh3D_ALE.hpp>
#include <life/lifesolver/NavierStokesAleSolverCT.hpp>
#include <ud_functions.hpp>
#include <life/lifefilters/medit.hpp>
#include <iostream>

#include <life/lifesolver/fluidToMaster.hpp>
#include <pvm3.h>


 
using namespace LifeV;
using namespace std;

int main(int argc, char** argv) {
  
 
  
  // Reading from data file
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  UInt FSItest= data_file( "FSItest", 0);
  UInt FSIalgo= data_file( "FSIalgo", 0);
  UInt FSIsche= data_file( "FSIsche", 0);

  // Boundary conditions handlers
  //
  BCHandler bcHu; 
  BCHandler bcHd;   
  BCHandler bcHdz;
  BCHandler bcHp; 
  BCHandler bcHdu;     
  BCHandler bcHdp; 
  BCHandler bcHextOper; 
  
  
  // The NavierStokes ALE solver
  /*
  NavierStokesAleSolverCT< RegionMesh3D_ALE<LinearHexa> > fluid(data_file,
								 feHexaQ1,    
								 quadRuleHexa8pt, 
								 quadRuleQuad4pt, 
								 bcHu,
								 bcHp,
								 bcHextOper);
  */
  NavierStokesAleSolverCT< RegionMesh3D_ALE<LinearTetra> > fluid(data_file,
								 feTetraP1,
								 quadRuleTetra4pt,
								 quadRuleTria3pt, 
								 bcHu,
								 bcHp,
								 bcHextOper);
  // 
  // Outputs
  fluid.showMe();
  
  // postprocess
  Medit medit(fluid);
  
  // solid and fluid number of total dof
  UInt dim_fluid = fluid.uDof().numTotalDof();
  
  //======================================================================
  //  BOUNDARY CONDITIONS
  //======================================================================
  
  BCFunctionBase bcf0(fZero); 
  BCFunctionBase bcfg(g); 
  BCFunctionBase bcfg_anev(g_anev); 
  BCFunctionBase bcfg_caro(g_caro);  
  BCFunctionBase inFlow7(u_in);

  Vector disp_interf(3*dim_fluid); 
  Vector minus_rho_da_interf(3*dim_fluid);  
  Vector minus_rho_a_interf(3*dim_fluid); 

  disp_interf = ZeroVector( disp_interf.size() ); 
  minus_rho_da_interf = ZeroVector( minus_rho_da_interf.size() ); 
  minus_rho_a_interf = ZeroVector( minus_rho_a_interf.size() );


  

  BCVector a_wall(  minus_rho_a_interf,  dim_fluid, 2); 
  BCVector da_wall( minus_rho_da_interf, dim_fluid, 2); 

  BCVector displ(   disp_interf ,        dim_fluid); 
  BCVector g_wall(  fluid.p(),           dim_fluid, 1); 
  BCVector dg_wall( fluid.dp(),          dim_fluid, 1);
  BCVector u_wall(  fluid.w(),           dim_fluid); 
  BCVector du_wall( fluid.dw(),         dim_fluid); 

  switch( FSItest )
    {  
    case 1:
      
      bcHp.addBC("Wall",        4, Natural,   Scalar, a_wall);   
      bcHp.addBC("Wall_Edges",  5, Essential, Scalar, bcf0);    
      bcHp.addBC("Wall_Edges",  6, Essential, Scalar, bcf0);   
      bcHp.addBC("Wall_Edges",  7, Essential, Scalar, bcfg_anev);
      bcHp.addBC("OutFlow",     1, Essential, Scalar, bcf0);   
      bcHp.addBC("OutFlow",     2, Essential, Scalar, bcf0);
      bcHp.addBC("InFlow",      3, Essential, Scalar, bcfg_anev);
      
      bcHdp.addBC("Wall",        4, Natural,   Scalar, da_wall);   
      bcHdp.addBC("OutFlow",     2, Essential, Scalar, bcf0);   
      bcHdp.addBC("OutFlow",     3, Essential, Scalar, bcf0);
      bcHdp.addBC("Wall_Edges",  5, Essential, Scalar, bcf0);    
      bcHdp.addBC("Wall_Edges",  6, Essential, Scalar, bcf0);   
      bcHdp.addBC("Wall_Edges",  7, Essential, Scalar, bcf0);
      bcHdp.addBC("InFlow",      1, Essential, Scalar, bcf0);
      
      
      bcHu.addBC("Wall",        4, Essential, Full, u_wall,  3);    
      bcHu.addBC("Wall_Edges",  5, Essential, Full, u_wall,  3); 
      bcHu.addBC("Wall_Edges",  6, Essential, Full, u_wall,  3);
      bcHu.addBC("Wall_Edges",  7, Essential, Full, u_wall,  3);
      
      bcHdu.addBC("Wall",        4, Essential, Full, du_wall,  3);    
      bcHdu.addBC("Wall_Edges",  5, Essential, Full, du_wall,  3); 
      bcHdu.addBC("Wall_Edges",  6, Essential, Full, du_wall,  3);
      bcHdu.addBC("Wall_Edges",  7, Essential, Full, du_wall,  3);
      
      // Boundary conditions: computation of the pressure stress 
      bcHdz.addBC("Interface", 4, Natural, Full, dg_wall, 3); 
      bcHd.addBC("Interface", 4, Natural, Full, g_wall, 3);
      
      
      bcHextOper.addBC("Interface", 4, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     5, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     6, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     7, Essential, Full, displ, 3);
      bcHextOper.addBC("Axis",      1, Essential, Full, bcf0,  3); 
      bcHextOper.addBC("Axis",      2, Essential, Full, bcf0,  3); 
      bcHextOper.addBC("Axis",      3, Essential, Full, bcf0,  3); 
      break; 
 
    case 3:
      
      bcHp.addBC("Wall",            2, Natural,   Scalar, a_wall);   
      bcHp.addBC("Wall_Edges_Out",  5, Essential, Scalar, bcf0);    
      bcHp.addBC("Wall_Edges_In",   4, Essential, Scalar, bcfg);   
      bcHp.addBC("OutFlow",         3, Essential, Scalar, bcf0);
      bcHp.addBC("InFlow",          1, Essential, Scalar, bcfg);
      
      bcHdp.addBC("Wall",            2, Natural,   Scalar, da_wall);   
      bcHdp.addBC("Wall_Edges_Out",  5, Essential, Scalar, bcf0);    
      bcHdp.addBC("Wall_Edges_In",   4, Essential, Scalar, bcf0);   
      bcHdp.addBC("OutFlow",         3, Essential, Scalar, bcf0);
      bcHdp.addBC("InFlow",          1, Essential, Scalar, bcf0);

      bcHu.addBC("Wall",        2, Essential, Full, u_wall,  3);    
      bcHu.addBC("Wall_Edges",  4, Essential, Full, u_wall,  3); 
      bcHu.addBC("Wall_Edges",  5, Essential, Full, u_wall,  3);

      bcHdu.addBC("Wall",        2, Essential, Full, du_wall, 3);    
      bcHdu.addBC("Wall_Edges",  4, Essential, Full, du_wall, 3); 
      bcHdu.addBC("Wall_Edges",  5, Essential, Full, du_wall, 3);


      // Boundary conditions: computation of the pressure stress 
      bcHdz.addBC("Interface", 2, Natural, Full, dg_wall, 3); 
      bcHd.addBC("Interface", 2, Natural, Full, g_wall, 3);
      
      
      bcHextOper.addBC("Interface", 2, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     4, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     5, Essential, Full, displ, 3);
      bcHextOper.addBC("Axis",      1, Essential, Full, bcf0,  3); 
      bcHextOper.addBC("Axis",      3, Essential, Full, bcf0,  3); 
      break;

   case 4:
      
      
      bcHp.addBC("Wall",        4, Natural,   Scalar, a_wall);
      bcHp.addBC("InFlow",      17, Essential, Scalar, bcfg_caro);
      bcHp.addBC("Wall_Edges",  15, Essential, Scalar, bcfg_caro);
      bcHp.addBC("Wall_Edges",  9, Essential, Scalar, bcf0);
      bcHp.addBC("Wall_Edges",  16, Essential, Scalar, bcf0);
      bcHp.addBC("Wall_Edges",  6, Essential, Scalar, bcf0);
      bcHp.addBC("OutFlow",     11, Essential, Scalar, bcf0);
      bcHp.addBC("OutFlow",     14, Essential, Scalar, bcf0);
      bcHp.addBC("OutFlow",     8, Essential, Scalar, bcf0);
      
      bcHdp.addBC("Wall",        4,  Natural,   Scalar,  da_wall);
      bcHdp.addBC("InFlow",      17, Essential, Scalar,   bcf0);
      bcHdp.addBC("Wall_Edges",  15, Essential, Scalar,   bcf0);
      bcHdp.addBC("Wall_Edges",  9,  Essential, Scalar,  bcf0);
      bcHdp.addBC("Wall_Edges",  16, Essential, Scalar,   bcf0);
      bcHdp.addBC("Wall_Edges",  6,  Essential, Scalar,  bcf0);
      bcHdp.addBC("OutFlow",     11, Essential, Scalar,   bcf0);
      bcHdp.addBC("OutFlow",     14, Essential, Scalar,   bcf0);
      bcHdp.addBC("OutFlow",     8,  Essential, Scalar,  bcf0);

      
      bcHu.addBC("Wall",        4, Essential, Full, u_wall,  3);
      bcHu.addBC("Wall_Edges", 15, Essential, Full, u_wall,  3);
      bcHu.addBC("Wall_Edges",  9, Essential, Full, u_wall,  3);
      bcHu.addBC("Wall_Edges", 16, Essential, Full, u_wall,  3);
      bcHu.addBC("Wall_Edges",  6, Essential, Full, u_wall,  3);
      
      bcHdu.addBC("Wall",        4, Essential, Full, du_wall,  3);
      bcHdu.addBC("Wall_Edges", 15, Essential, Full, du_wall,  3);
      bcHdu.addBC("Wall_Edges",  9, Essential, Full, du_wall,  3);
      bcHdu.addBC("Wall_Edges", 16, Essential, Full, du_wall,  3);
      bcHdu.addBC("Wall_Edges",  6, Essential, Full, du_wall,  3);
      
      // Boundary conditions: computation of the pressure stress
      bcHdz.addBC("Interface", 4, Natural, Full, dg_wall, 3);
      bcHd.addBC("Interface", 4, Natural, Full, g_wall, 3);
      
      
      bcHextOper.addBC("Interface", 4, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",    15, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     9, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",    16, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     6, Essential, Full, displ, 3);
      bcHextOper.addBC("Base",     17, Essential, Full, bcf0,  3);
      bcHextOper.addBC("Base",      8, Essential, Full, bcf0,  3);
      bcHextOper.addBC("Base",     11, Essential, Full, bcf0,  3);
      bcHextOper.addBC("Base",     14, Essential, Full, bcf0,  3);
      break;



case 5:
      
      bcHp.addBC("Wall",            1, Natural,   Scalar, a_wall);   
      bcHp.addBC("Wall_Edges_Out",  5, Essential, Scalar, bcf0);    
      bcHp.addBC("Wall_Edges_In",  10, Essential, Scalar, bcfg);   
      bcHp.addBC("OutFlow",         3, Essential, Scalar, bcf0);
      bcHp.addBC("InFlow",          2, Essential, Scalar, bcfg);
      
      bcHdp.addBC("Wall",            1, Natural,   Scalar, da_wall);   
      bcHdp.addBC("Wall_Edges_Out",  5, Essential, Scalar, bcf0);    
      bcHdp.addBC("Wall_Edges_In",  10, Essential, Scalar, bcf0);   
      bcHdp.addBC("OutFlow",         3, Essential, Scalar, bcf0);
      bcHdp.addBC("InFlow",          2, Essential, Scalar, bcf0);

      bcHu.addBC("Wall",        1, Essential, Full, u_wall,  3);    
      bcHu.addBC("Wall_Edges",  5, Essential, Full, u_wall,  3); 
      bcHu.addBC("Wall_Edges", 10, Essential, Full, u_wall,  3);
      
      bcHdu.addBC("Wall",        1, Essential, Full, du_wall,  3);    
      bcHdu.addBC("Wall_Edges",  5, Essential, Full, du_wall,  3); 
      bcHdu.addBC("Wall_Edges", 10, Essential, Full, du_wall,  3);

      // Boundary conditions: computation of the pressure stress 
      bcHdz.addBC("Interface", 1, Natural, Full, dg_wall, 3); 
      bcHd.addBC("Interface", 1, Natural, Full, g_wall, 3);
      
      
      bcHextOper.addBC("Interface", 1, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",     5, Essential, Full, displ, 3);
      bcHextOper.addBC("Edges",    10, Essential, Full, displ, 3);
      bcHextOper.addBC("Axis",      2, Essential, Full, bcf0,  3); 
      bcHextOper.addBC("Axis",      3, Essential, Full, bcf0,  3); 
      break;


    case 7:
      
      bcHp.addBC("Wall",            4, Natural,   Scalar, a_wall);   
      bcHp.addBC("Wall_Edges_Out",  7, Essential, Scalar, bcf0);    
      bcHp.addBC("OutFlow",         2, Essential, Scalar, bcf0);
       
      bcHdp.addBC("Wall",            4, Natural,   Scalar, da_wall);   
      bcHdp.addBC("Wall_Edges_Out",  7, Essential, Scalar, bcf0);    
      bcHdp.addBC("OutFlow",         2, Essential, Scalar, bcf0);
   

      bcHu.addBC("Wall",      4,  Essential, Full, u_wall, 3);
      bcHu.addBC("Wall",      5,  Essential, Full, u_wall, 3); 
      bcHu.addBC("Wall",      6,  Essential, Full, u_wall, 3);
      bcHu.addBC("Rigid-Wall",      3,  Essential, Full, bcf0, 3);
      bcHu.addBC("Rigid-Wall",      7,  Essential, Full, bcf0, 3);
      bcHu.addBC("OutFlow",   2,  Natural,   Full, bcf0,   3);
      bcHu.addBC("InFlow",    1,  Essential, Full, inFlow7,   3);

      bcHdu.addBC("Wall",     4,  Essential, Full, du_wall, 3);
      bcHdu.addBC("Wall",     5,  Essential, Full, du_wall, 3); 
      bcHdu.addBC("Wall",     6,  Essential, Full, du_wall, 3);	
      bcHdu.addBC("Rigid-Wall",3,  Essential, Full, bcf0, 3);   
      bcHdu.addBC("Rigid-Wall",      7,  Essential, Full, bcf0, 3);
      bcHdu.addBC("OutFlow",  2,  Natural,   Full, bcf0,    3);
      bcHdu.addBC("InFlow",   1,  Essential, Full, bcf0,    3);


      // Boundary conditions: computation of the pressure stress 
      bcHdz.addBC("Interface", 4, Natural, Full, dg_wall, 3); 
      bcHd.addBC("Interface", 4, Natural, Full, g_wall, 3);
      	
      bcHextOper.addBC("Wall",   4, Essential, Full, displ, 3);  //  displacement given on the lateral points
      bcHextOper.addBC("Wall",   5, Essential, Full, displ, 3); 
      bcHextOper.addBC("Wall",   6, Essential, Full, displ, 3); 
      bcHextOper.addBC("Rigid-Wall",3,  Essential, Full, bcf0, 3); 
      bcHextOper.addBC("InLet",  1, Essential, Full, bcf0,   3);
      bcHextOper.addBC("OutLet", 2, Essential, Full, bcf0,   3);
      bcHextOper.addBC("Wall_OutLet", 7, Essential, Full, bcf0,   3);
      
      break;


    default:
      ERROR_MSG("This FSI test is not yet implemented");
    }
  
  int masterId = pvm_parent();
  if(masterId==PvmNoParent)
    {
      cout << "I am a poor lonesome job\n";
      exit(1);
    } 
  else if(masterId==0)
    {
      cout << "PVM is not ok\n";
      exit(1);
    }


  FluidToMaster< NavierStokesAleSolverCT< RegionMesh3D_ALE<LinearTetra> > > 
    fluidToMaster(fluid,masterId,-1);
  
  cout <<"%%%%%%%%%%% \n";
  cout << " I am the slave " << pvm_mytid()<< endl;
  cout << " My master is " << masterId << endl;
  fluidToMaster.rvFromMasterFSI(140);
  fluidToMaster.sdCoorToMasterFSI(); 
 
  // initial conditions
  fluid.initialize(u0);
      
  // Temporal loop
  Real time=0.0;
  
  do 
    {
      
      // receive from master
      fluidToMaster.rvFromMasterFSI(110); 
      
      //============================================================= 
      // switch on status
      switch ( fluidToMaster.status() ) {
	
	//===========================================================
	// Normal exit
	//
      case -1:
	// post-processing last converged
	medit.postProcessALE(fluid);
	cout << "Normal fluid exit\n" << flush;  
	exit(1); 
	break;

	//===========================================================
	// New time step (continue to case 0)
	//
      case 1: 
	
	time+=fluid.timestep();

	// post-processing last converged fluid
	medit.postProcessALE(fluid); 
	
	fluid.timeAdvance( f, time);

	if ( FSIsche == 2 ) // semi-implicit
	  {
	    fluidToMaster.dep_fluid_interf( disp_interf );
	    //fluidToMaster.dep_fluid_interf_1( disp_interf );
	    fluid.updateMesh(time);
	    fluid.iterateVelocity(time);
	  }
	
	//============================================================
	// inner evaluation (Picard and Newton)
      case 0:
	 
	fluidToMaster.minus_rho_a_interf( minus_rho_a_interf );
	
	if ( FSIsche == 2 ) // semi-implicit
	    fluid.iteratePressure(time);
	else // implicit
	  {
	    fluidToMaster.dep_fluid_interf( disp_interf );
	    fluid.updateMesh(time);
	    fluid.iterate(time);
	  }
	fluid.updateForce(bcHd);
	break;

	//============================================================
	// Solving tangent problem: for Newton of quas-Newton
      case 2:  
	if  ( FSIsche == 2 ) // semi-implicit
	  {
	    fluidToMaster.minus_rho_da_interf( minus_rho_da_interf );
	    fluid.iteratePressureLin(time,bcHdp);
	  }
	else // implicit
	  {
	    switch ( FSIalgo ) 
	      {
	      case 1: // newton: shape
		
		fluidToMaster.dep_fluid_interf( disp_interf );
		fluidToMaster.minus_rho_da_interf( minus_rho_da_interf );
		fluid.dwUpdate(time);
		fluid.iterateLin(time,bcHdu,bcHdp);
		break;

	      case 2: // quasi newton: linearized ns, no-shape
		
		fluidToMaster.dep_fluid_interf( disp_interf );
		fluidToMaster.minus_rho_da_interf( minus_rho_da_interf );
		fluid.dwUpdate(time);
		fluid.iterateLin(time,bcHdu,bcHdp,0);
		break;
	  
	      case 3: // quasi-newton: laplacian
		fluidToMaster.minus_rho_da_interf( minus_rho_da_interf );
		fluid.iteratePressureLin(time,bcHdp);
		break;
		
	      default:
		ERROR_MSG("This FSI algorithm is not yet implemented");
	      }
	  }
	// update linearized force
	fluid.updateLinForce(bcHdz);
	break; 
      default:
	cout << " This status number is not allowed.\n" << endl << flush;
	exit(1);
      }
      
      //===============================================================
      // at the end send to master
      fluidToMaster.sdToMasterFSI();
    }
    while ( 1 ); 
    return 0;
}
