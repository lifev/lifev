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





// ===================================================
//! Includes
// ===================================================
#ifdef EPETRA_MPI
	#include "Epetra_MpiComm.h"
	#include <mpi.h>
#else
	#include "Epetra_SerialComm.h"
#endif

#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
//#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifesolver/dataADR.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifesolver/AdvectionDiffusionReactionSolver.hpp>

#include "laplacian.hpp"





// ===================================================
//! Namespaces & Define
// ===================================================
using namespace LifeV;

const int TOP    = 6;
const int BOTTOM = 5;
const int LEFT   = 3;
const int RIGHT  = 4;
const int FRONT  = 2;
const int BACK   = 1;





// ===================================================
//! User functions
// ===================================================

Real USinExact( const Real& /* t */,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID&  )
{
    Real pi = 3.141592653589793;
    return sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
    //return sin(pi*x)*sin(pi*y)*sin(pi*z);
}

Real USin( const Real& /* t */,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID&  )
{
    Real pi = 3.141592653589793;
    return 12*pi*pi*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
    //return pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z);
}

Real UOne( const Real& /* t */,
           const Real&,
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
            const ID&  )
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
    typedef ADRSolver< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type>                   vector_ptrtype;

    // Reading from data file
    //
    GetPot dataFile( Members->data_file_name.c_str() );

    bool verbose = (Members->comm->MyPID() == 0);

    //
    // The Laplacian Solver
    //
    if (verbose) std::cout << "The Laplacian Solver" << std::flush;


    // the boundary conditions
    BCHandler bcADR( 6, BCHandler::HINT_BC_ONLY_ESSENTIAL );
    BCFunctionBase uZero( Members->getUZero() );
    BCFunctionBase uOne ( Members->getUOne()  );

    // cube
    bcADR.addBC( "Top",     TOP,    Essential, Full,      uZero, 1 );
    bcADR.addBC( "Bottom",  BOTTOM, Essential, Full,      uZero, 1 );
    bcADR.addBC( "Left",    LEFT,   Essential, Full,      uZero, 1 );
    bcADR.addBC( "Right",   RIGHT,  Essential, Full,      uZero, 1 );
    bcADR.addBC( "Front",   FRONT,  Essential, Full,      uZero, 1 );
    bcADR.addBC( "Back",    BACK,   Essential, Full,      uZero, 1 );

    DataADR<RegionMesh3D<LinearTetra> > dataADR( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(* dataADR.mesh(), *Members->comm);
    dataADR.setMesh(meshPart.mesh());

    // Advection Velocity Space: u0 P1 interpolation

    const RefFE*    refFE_vel  (0);
    const QuadRule* qR_vel     (0);
    const QuadRule* bdQr_vel   (0);

    refFE_vel = &feTetraP1;
    qR_vel    = &quadRuleTetra4pt;  // DoE 2
    bdQr_vel  = &quadRuleTria3pt;   // DoE 2

    // Scalar Temperature Space: adr
    std::string adrOrder =  dataFile( "adr/discretization/order", "P1");
    const RefFE*    refFE_adr(0);
    const QuadRule* qR_adr(0);
    const QuadRule* bdQr_adr(0);

    if ( adrOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "  Space order : P1" << std::flush;
            refFE_adr = &feTetraP1;
            qR_adr    = &quadRuleTetra4pt; // DoE 5
            bdQr_adr  = &quadRuleTria3pt;   // DoE 2
        }
    else
        if ( adrOrder.compare("P2") == 0 )
            {
                if (verbose) std::cout << " Space order : P2 interpolation ";
                refFE_adr = &feTetraP2;
                qR_adr    = &quadRuleTetra15pt;  // DoE 2
                bdQr_adr  = &quadRuleTria3pt;   // DoE 2
            }

    if (verbose) std::cout << std::endl;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > adrFESpace(meshPart,
                                                               *refFE_adr,
                                                               *qR_adr,
                                                               *bdQr_adr,
                                                               1,
                                                               *Members->comm);

    // Laplacian u0=(0,0,0)
    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > u0FESpace(meshPart,
                                                             *refFE_vel,
                                                             *qR_vel,
                                                             *bdQr_vel,
                                                             3,
                                                             *Members->comm);


    UInt totalADRDof   = adrFESpace.map().getMap(Unique)->NumGlobalElements();

    ADRSolver< RegionMesh3D<LinearTetra> > adr (dataADR,
                                                adrFESpace,
                                                u0FESpace,
                                                *Members->comm);

    adr.setUp(dataFile);
    adr.buildSystem();

    // Laplacian Solver
    vector_type betaFluid( u0FESpace.map() );
    betaFluid *= 0.;

    MPI_Barrier(MPI_COMM_WORLD);

    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;
    exporter.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "temperature", Members->comm->MyPID()) );

    dataADR.setTime(0);

    EpetraMap fullAdrMap(adr.getMap());
    vector_type rhsADR ( fullAdrMap );

    rhsADR  *= 0.;
    adrFESpace.interpolate(USin, rhsADR);
    adrFESpace.interpolateBC(bcADR, rhsADR, 0.);

    rhsADR = adr.matrMass()*rhsADR;
    adr.updateSystem(1., betaFluid, rhsADR);
    adr.iterate(bcADR);

    adr.resetPrec();

    // post processing setup
    vector_type uComputed(adr.solution(), Repeated );

    vector_ptrtype temperature  ( new vector_type(adr.solution(), Repeated ) );
    exporter->addVariable( ExporterData::Scalar, "temperature", temperature,
                           UInt(0), UInt(adrFESpace.dof().numTotalDof()));
    exporter->postProcess( 0 );

    // Error H1 Norm
    vector_type uExact ( fullAdrMap, Repeated );
    uExact  *= 0.;
    adrFESpace.interpolate(USinExact, uExact);

    uExact -= uComputed;
    Real H1error = adrFESpace.H1Norm( uExact );
    std::cout << "Error Norm H1: " << H1error << std::endl;
}
