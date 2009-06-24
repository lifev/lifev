/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  L. Iapichino  <laura.iapichino@epfl.ch>
              C. Malossi    <cristiano.malossi@epfl.ch>
              A. Manzoni    <andrea.manzoni@epfl.ch>
       Date: 2009-03-24

  Copyright (C) 2009 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file laplacian.cpp
   \author L. Iapichino <laura.iapichino@epfl.ch>, C. Malossi <cristiano.malossi@epfl.ch>, A. Manzoni <andrea.manzoni@epfl.ch>
   \date 2009-03-24
 */

/* ========================================================

Simple Laplacian test with Dirichlet Boundary condition

Solve the problem

               - \Delta u = f

               u = 0 on the boundary


 3D: with the source term f = 12 \pi^2 sin(2 \pi x) sin(2 \pi y) sin (2 \pi z) on a cube
 2D: with the source term f = 8 \pi^2 sin(2 \pi x) sin(2 \pi y) on a square

 the rhs is computed as rhs = Mass_Matrix * f_iterpolated


 More generally this test can solve the problem:

               - \nu \Delta u + \beta \nabla u + \sigma u = f

               u = g on the boundary

 being \nu and \sigma constants defined in the data file and \beta interpolated.

*/


// ===================================================
//! Includes
// ===================================================
#ifdef EPETRA_MPI
	#include "Epetra_MpiComm.h"
	#include <mpi.h>
#else
	#include "Epetra_SerialComm.h"
#endif

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataADR.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifesolver/AdvectionDiffusionReactionSolver.hpp>

#include "laplacian.hpp"

#define POSTPROCESS 0


// ===================================================
//! Namespaces & Define
// ===================================================
using namespace LifeV;

#ifdef TWODIM

	const int TOP    = 3;
	const int BOTTOM = 1;
	const int LEFT   = 4;
	const int RIGHT  = 2;

    typedef RegionMesh2D<LinearTriangle> RegionMesh;

    const std::string discretization_section="adr/space_discretization2D";

#elif defined THREEDIM

    const int TOP    = 6;
    const int BOTTOM = 5;
    const int LEFT   = 3;
    const int RIGHT  = 4;
    const int FRONT  = 2;
    const int BACK   = 1;

    typedef RegionMesh3D<LinearTetra> RegionMesh;

    const std::string discretization_section="adr/space_discretization3D";


#endif

// ===================================================
//! User functions
// ===================================================

#ifdef TWODIM

class AnalyticalSol
{
public:
	inline Real operator()(Real /*t*/, Real x,Real y,Real /*z*/, UInt /*ic*/=0) const {
		return sin(2*Pi*x)*sin(2*Pi*y);
	}
	inline Real grad(UInt icoor, Real /*t*/, Real x,Real y,Real /*z*/, UInt /*ic*/=0) const {
		switch(icoor)
		{
  	  	case 1: //der_x
  	  		return 2*Pi*cos(2*Pi*x)*sin(2*Pi*y);
  	  	case 2: //der_y
  	  		return 2*Pi*sin(2*Pi*x)*cos(2*Pi*y);
  	  	default:
  	  		return 0;
		}
	}
};


Real source_in( const Real& /* t */,
                const Real& x,
                const Real& y,
                const Real& /*z*/,
                const ID&   icomp)
{
    return 8*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y);
}

//solution on the boundary
Real g( const Real& /* t */,
           const Real& x,
           const Real& y,
           const Real& z,
           const ID&   )
{
    return sin(2*Pi*x)*sin(2*Pi*y);
}

Real beta( const Real& /* t */,
           const Real& ,
           const Real& ,
           const Real& ,
           const ID& icomp )
{
	switch(icomp)
	{
		case 1:
			return 0;
		case 2:
			return 0;
		default:
			return 0;
	}
}



#elif defined THREEDIM
class AnalyticalSol
{
public:
	inline Real operator()(Real /*t*/, Real x,Real y,Real z, UInt /*ic*/=0) const {
		return sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z);
	}
	inline Real grad(UInt icoor, Real /*t*/, Real x,Real y,Real z, UInt /*ic*/=0) const {
		switch(icoor)
		{
  	  	case 1: //der_x
  	  		return 2*Pi*cos(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z);
  	  	case 2: //der_y
  	  		return 2*Pi*sin(2*Pi*x)*cos(2*Pi*y)*sin(2*Pi*z);
  	  	case 3:
  	  		return 2*Pi*sin(2*Pi*x)*sin(2*Pi*y)*cos(2*Pi*z);
  	  	default:
  	  		return 0;
		}
	}
};

//solution on the boundary
Real g( const Real& /* t */,
           const Real& x,
           const Real& y,
           const Real& z,
           const ID&   icomp)
{
    return sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z);
}


Real source_in( const Real& /* t */,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  icomp)
{
    return 12*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z);
}

Real beta( const Real& /* t */,
           const Real& ,
           const Real& ,
           const Real& ,
           const ID& icomp )
{
	switch(icomp)
	{
		case 1:
			return 0;
		case 2:
			return 0;
		case 3:
			return 0;
		default:
			return 0;
	}
}
#endif


Real UOne( const Real& /* t */,
           const Real& ,
           const Real& ,
           const Real& ,
           const ID&   )
{
    return 1.;
}

Real UZero( const Real& /* t */,
           const Real& ,
           const Real& ,
           const Real& ,
           const ID&   )
{
    return 0.;
}




// ===================================================
//! Private Members
// ===================================================
struct laplacian::Private
{
    Private() {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    std::string    data_file_name;
    Epetra_Comm*   comm;

    /**
     *
     *
     * Function Types
     */

    fct_type getUOne()
    {
    	fct_type f;
    	f = boost::bind(&UOne, _1, _2, _3, _4, _5);
        return f;
        }

    fct_type getUZero()
    {
    	fct_type f;
    	f = boost::bind(&UZero, _1, _2, _3, _4, _5);
    	return f;
    }
};





// ===================================================
//! Constructors
// ===================================================
laplacian::laplacian( int argc,
                      char** argv,
                      LifeV::AboutData const& /*ad*/,
                      LifeV::po::options_description const& /*od*/ ): Members( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;

	#ifdef EPETRA_MPI
		std::cout << "Epetra Initialization" << std::endl;

		Members->comm = new Epetra_MpiComm( MPI_COMM_WORLD );
	#else
		Members->comm = new Epetra_SerialComm();
	#endif
}





// ===================================================
//! Methods
// ===================================================
void
laplacian::run()
{
	typedef ADRSolver< RegionMesh >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type>        vector_ptrtype;

    // Reading from data file
    //
    GetPot dataFile( Members->data_file_name.c_str() );

    bool verbose = (Members->comm->MyPID() == 0);

    //
    // The Laplacian Solver
    //

    if (verbose) std::cout << "The Laplacian Solver" << std::flush;

    // the boundary conditions
    BCFunctionBase g_Ess(g);
    BCFunctionBase uZero( Members->getUZero() );
    BCFunctionBase uOne ( Members->getUOne()  );

#ifdef TWODIM

	BCHandler bcADR( 4, BCHandler::HINT_BC_ONLY_ESSENTIAL );

    bcADR.addBC( "Top",     TOP,    Essential, Full,      g_Ess, 1 );
    bcADR.addBC( "Bottom",  BOTTOM, Essential, Full,      g_Ess, 1 );
    bcADR.addBC( "Left",    LEFT,   Essential, Full,      g_Ess, 1 );
    bcADR.addBC( "Right",   RIGHT,  Essential, Full,      g_Ess, 1 );

#elif defined THREEDIM

	BCHandler bcADR( 6, BCHandler::HINT_BC_ONLY_ESSENTIAL );

    bcADR.addBC( "Top",     TOP,    Essential, Full,      uZero, 1 );
    bcADR.addBC( "Bottom",  BOTTOM, Essential, Full,      uZero, 1 );
    bcADR.addBC( "Left",    LEFT,   Essential, Full,      uZero, 1 );
    bcADR.addBC( "Right",   RIGHT,  Essential, Full,      uZero, 1 );
    bcADR.addBC( "Front",   FRONT,  Essential, Full,      uZero, 1 );
    bcADR.addBC( "Back",    BACK,   Essential, Full,      uZero, 1 );
#endif

    DataADR<RegionMesh> dataADR( dataFile, discretization_section);

    partitionMesh< RegionMesh>   meshPart(* dataADR.mesh(), *Members->comm);
    dataADR.setMesh(meshPart.mesh());

    if(verbose) dataADR.showMe();

    //finite element space of the solution

    std::string adrOrder =  dataFile( (discretization_section+"/order").data(), "P1");
    FESpace< RegionMesh, EpetraMap > adrFESpace(meshPart,adrOrder,1,*Members->comm);

    //finite element space of the advection term
    std::string betaOrder = "P1_HIGH";
    FESpace< RegionMesh, EpetraMap > betaFESpace(meshPart,betaOrder,nDimensions,*Members->comm);

    //instantiation of the AdvectionDiffusionReactionSolver class

    ADRSolver< RegionMesh > adr (dataADR,
                                 adrFESpace,
                                 betaFESpace,
                                 *Members->comm);

    Chrono chrono;

    chrono.start();
    adr.setUp(dataFile);
    adr.buildSystem();

    // Laplacian Solver

#if POSTPROCESS
    boost::shared_ptr< Exporter<RegionMesh> > exporter;
    exporter.reset( new Ensight<RegionMesh> ( dataFile, meshPart.mesh(), "u", Members->comm->MyPID()) );
#endif

    dataADR.setTime(0);

    EpetraMap fullAdrMap(adr.getMap());

    //computing the iterpolation of the advection vector
    vector_type betaFluid( betaFESpace.map() );
    betaFESpace.interpolate(beta, betaFluid);

    //computing the rhs
    vector_type rhsADR ( fullAdrMap );
    adrFESpace.interpolate(source_in, rhsADR);
    rhsADR = adr.matrMass()*rhsADR;


    //updating the system with the reaction and advection terms and the rhs
    //(including the treatment of the boundary conditions)
    Real sigma = dataFile( "adr/physics/sigma", 0);
    adr.updateSystem(sigma, betaFluid, rhsADR);

    //solve the linear system
    adr.iterate(bcADR);

    chrono.stop();

    if (verbose) std::cout << "\n \n -- Total time = " << chrono.diff() << std::endl << std::endl;

    // post processing setup

#if POSTPROCESS
    vector_ptrtype temperature  ( new vector_type(adr.solution(), Repeated ) );
    exporter->addVariable( ExporterData::Scalar, "temperature", temperature,
                           UInt(0), UInt(adrFESpace.dof().numTotalDof()));
    exporter->postProcess( 0 );
#endif

    // Error L2 and H1 Norms

    AnalyticalSol uExact;
    vector_type uComputed(adr.solution(), Repeated );

    Real H1_Error, H1_RelError, L2_Error, L2_RelError;

    L2_Error = adrFESpace.L2Error(uExact, uComputed, 0 ,&L2_RelError);
    H1_Error = adrFESpace.H1Error(uExact, uComputed, 0 ,&H1_RelError);

    if (verbose)
    	std::cout << "Error Norm L2: " << L2_Error <<
		"\nRelative Error Norm L2: " << L2_RelError<<
        "\nError Norm H1: " << H1_Error <<
        "\nRelative Error Norm H1: " << H1_RelError<<std::endl;
}
