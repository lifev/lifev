#ifndef __CT_CASE_HH
#define __CT_CASE_HH 1

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <ChorinTemam.hpp>
#include <iostream>

using namespace LifeV;

/*!
 * \struct CTcase
 * \brief Interface for Chorin-Temam / projection case study
 * \see Howto add a specific case study below.
 *
 */

struct CT::CTcase
{
    
    /*********************************************************************** 
     * Mandatory part:
     *   this part should be sufficiently general for not to be modified 
     *   by specific case studies.
     *
     ***********************************************************************/

    // typedefs
    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID&)> func_type;
    
    // mandatory members
    GetPot C_data;
    bool C_hasdata;
    int NUM_BCS;
    BCHandler C_bcHu;
    BCHandler C_bcHp;
    Epetra_Comm *C_comm;

    // mandatory methods : methods called by CT upstream class
    BCHandler *get_bcHu() {return &C_bcHu;}
    BCHandler *get_bcHp() {return &C_bcHp;}
    void set_data(const GetPot& _data, Epetra_Comm *_comm) 
    {
        C_comm = _comm;
	C_data = _data;
	C_hasdata = true;
	set_user_data();
    }
    // void set_user_data();
    // void set_bcs();

    // constructor
    CTcase() : C_hasdata(false),
               NUM_BCS(5),
    	       C_bcHu(NUM_BCS, BCHandler::HINT_BC_NONE),
	       C_bcHp(NUM_BCS, BCHandler::HINT_BC_NONE)
	       {}

    // handy die
    void die(std::string msg, Epetra_Comm *comm)
    { 
    	if (!comm->MyPID()) {
	  std::cout << "  --  ERROR: " << msg << std::endl;
	  std::cout << "  --  Exiting ..." << std::endl;
	}
	exit(1);
    }

    /**********************************************************************
     * Specific case part: 
     *   place where user implements specific case study stuff.
     *
     **********************************************************************/

    // user defined members (those read in the data file)
    Real C_density;
    Real C_viscosity;
    Real C_dimX;
    Real C_dimY;
    Real C_dimZ;
    Real C_diameter;
    Real C_ampl;
    Real C_amplstep;
    
    // user defined members (eg. reference used in the mesh)
    int INLET;
    int WALL;
    int SLIPWALL;
    int OUTLET;
    int CYLINDER;
    std::vector<ID> zComp;

    // user defined (boost) functions
    func_type u3D_zero;
    func_type u3D_in;

    // user defined data settings
    void set_user_data()
    {
    	if (!C_hasdata) die("Data case has not be loaded", C_comm);
	// now read whatever data we need
	C_density   = C_data ("fluid/physics/density", 1.);
	C_viscosity = C_data ("fluid/physics/vicosity", 1.);
	C_dimX = C_data ("fluid/problem/dimX", 40.);
	C_dimY = C_data ("fluid/problem/dimY", 20.);
	C_dimZ = C_data ("fluid/problem/dimZ", 4.);
	C_diameter = C_data ("fluid/problem/diameter", 1.);
	C_ampl = C_data ("fluid/problem/amplitude", 1.);
	C_amplstep = C_data ("fluid/problem/amplstep", 2.);
	// specify boundary vertices references
	INLET    = 40;
	WALL     = 60;
	SLIPWALL = 61;
	OUTLET   = 50;
	CYLINDER = 70;
	// bind some locally defined methods to (boost) functions for use in BC, etc.
	u3D_zero = boost::bind(&CT::CTcase::u3DZero, this, _1, _2, _3, _4, _5);
	u3D_in = boost::bind(&CT::CTcase::u3Dcyl_dyn, this, _1, _2, _3, _4, _5);
    }

    // user defined methods to be used like functions in BC, etc.
      
    // zero 
    Real u3DZero(const Real&, const Real&, const Real&, const Real&, const ID&)
    {return(0);}

    // time independant parabolic profile
    Real u3Dcyl(const Real&, const Real&, const Real& y, const Real&, const ID& i)
    {
      if (i == 1)
        return ( C_ampl / (C_dimY * C_dimY)*(y + C_dimY)*(C_dimY-y) );
      else
	return(0.);
    }
    
    // time dependant parabolic profile
    Real u3Dcyl_dyn(const Real& t, const Real&, const Real& y, const Real&, const ID& i)
    {
      if (i==1)
      {
        if (t < C_amplstep)
	  return ( (t/C_amplstep) * C_ampl / (C_dimY * C_dimY)*(y + C_dimY)*(C_dimY-y) );
	else 
	  return ( C_ampl / (C_dimY * C_dimY)*(y + C_dimY)*(C_dimY-y) );
      }
      else
	  return(0.);
    }
    
    // set boundary conditions
    void set_bcs()
    {
      zComp.resize(1);
      zComp[0] = 3;
      BCFunctionBase uIn(u3D_in);
      BCFunctionBase uZero(u3D_zero);

      // bc for the velocity
      C_bcHu.addBC("Inlet",    INLET, 	 Essential, Full,      uIn,   3);
      C_bcHu.addBC("Outlet",   OUTLET,	 Natural,   Full,      uZero, 3);
      C_bcHu.addBC("Wall",     WALL, 	 Essential, Full,      uZero, 3);
      C_bcHu.addBC("SlipWall", SLIPWALL, Essential, Component, uZero, zComp);
      C_bcHu.addBC("Cylinder", CYLINDER, Essential, Full,      uZero, 3);

      // bc for the pressure
      C_bcHp.addBC("Inlet",    INLET,    Natural,   Scalar,    uZero);
      C_bcHp.addBC("Outlet",   OUTLET,   Essential, Scalar,    uZero);
      C_bcHp.addBC("Wall",     WALL,     Natural,   Scalar,    uZero);
      C_bcHp.addBC("SlipWall", SLIPWALL, Natural,   Scalar,    uZero);
      C_bcHp.addBC("Cylinder", CYLINDER, Natural,   Scalar,    uZero);
    
    }
    
}; // struct CT::CTcase

#endif /* __CT_CASE_HH */
