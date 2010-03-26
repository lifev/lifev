/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-19

  Copyright (C) 2005 EPFL

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
   \file cylinder.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-19
 */

// #include <lifeconfig.h>

// #include <life/lifecore/life.hpp>
// #include <life/lifecore/GetPot.hpp>
// #include <life/lifecore/debug.hpp>

// #include <life/lifefilters/importer.hpp>

//#include "NavierStokesSolverIP.hpp"
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataElasticStructure.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>

#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifefilters/noexport.hpp>


#include "structure.hpp"
#include <iostream>

using namespace LifeV;



Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    if (t == 0)
    {
        switch(i) {
            case 1:
	      return -z*(z - 5.)*x/50.;
               break;
            case 2:
	      return -z*(z - 5.)*y/50.;
	      break;
            case 3:
                return 0.0;
                break;
            default:
                ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
                break;
        }
    }

    return 0.0;
}

Real w0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{

  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break;
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }
  return 0.0;

}


Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

struct Structure::Private
{
    Private() :
        rho(1), poisson(1), young(1)
        {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;
    double rho, poisson, young;

    std::string data_file_name;

    Epetra_Comm*     comm;


};

Structure::Structure( int                                   argc,
                      char**                                argv,
                      Epetra_Comm &                         structComm,
                      LifeV::AboutData const&               ad,
                      LifeV::po::options_description const& od ):
    parameters( new Private() )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->rho = dataFile( "solid/physics/density", 1. );
    parameters->young = dataFile( "solid/physics/young", 1. );
    parameters->poisson  = dataFile( "solid/physics/poisson", 1. );

    std::cout << "density = " << parameters->rho << std::endl
              << "young   = " << parameters->young << std::endl
              << "poisson = " << parameters->poisson << std::endl;

    parameters->comm = &structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID()) std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;

}

void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}

void
Structure::run3d()
{

    typedef VenantKirchhofSolver< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    bool verbose = (parameters->comm->MyPID() == 0);
    // Number of boundary conditions for the velocity and mesh motion
    //
    BCHandler BCh(2);
    BCFunctionBase dZero( zero_scalar );
    //
    // dataElasticStructure
    //

    GetPot dataFile( parameters->data_file_name.c_str() );

    DataElasticStructure< RegionMesh3D<LinearTetra > >    dataStructure( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >            meshPart( *dataStructure.mesh(),
                                                                    *parameters->comm );

//    meshPart.rebuildMesh();


//    exit(0);

    dataStructure.setMesh(meshPart.mesh());

    std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");
    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > dFESpace(meshPart,dOrder,3,*parameters->comm);
    if (verbose) std::cout << std::endl;

    EpetraMap structMap(dFESpace.refFE(), meshPart, *parameters->comm);

    EpetraMap fullMap;

    for (UInt ii = 0; ii < nDimensions; ++ii)
    {
         fullMap += structMap;
    }

    VenantKirchhofSolver< RegionMesh3D<LinearTetra> > solid( dataStructure,
                                                             dFESpace,
                                                             *parameters->comm);
    solid.setUp(dataFile);
    solid.buildSystem();
    //
    // Boundary conditions for the displacement
    //
    BCFunctionBase fixed(dZero);

    BCh.addBC("Base2 ", 2 , Essential, Full, fixed, 3);
    BCh.addBC("Base3 ", 3 , Essential, Full, fixed, 3);
   
    //
    // Temporal data and initial conditions
    //
    Real dt = dataStructure.getTimeStep();
    Real T  = dataStructure.getEndTime();

    vector_ptrtype disp(new vector_type(solid.disp(), Unique));
    vector_ptrtype vel(new vector_type(solid.vel(), Unique));

    dFESpace.interpolate(d0, *disp, 0.0);
    dFESpace.interpolate(w0, *vel , 0.0);

    if (verbose) std::cout << "S- initialization ... ";

//    solid.initialize(d0, w0);
    solid.initialize(disp,vel); // displacement and velocity

    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose ) std::cout << "ok." << std::endl;
    //if (parameters->comm->NumProc() == 1 )  solid.postProcess();


    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
        exporter.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, "structure" ) );
    }
    else
#endif
    {
        if (exporterType.compare("none") == 0)
        {
            exporter.reset( new NoExport<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "structure", parameters->comm->MyPID()) );
        }
        else
        {
            exporter.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "structure", parameters->comm->MyPID()) );
        }
    }

    exporter->setDirectory( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( meshPart.mesh(), parameters->comm->MyPID() );

    vector_ptrtype solidDisp ( new vector_type(solid.disp(), exporter->mapType() ) );
    vector_ptrtype solidVel  ( new vector_type(solid.vel(),  exporter->mapType() ) );

    exporter->addVariable( ExporterData::Vector, "displacement", solidDisp,
                           UInt(0), dFESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::Vector, "velocity", solidVel,
                           UInt(0), dFESpace.dof().numTotalDof() );

    exporter->postProcess( 0 );


    //
    // Temporal loop
    //

    for (Real time = dt; time <= T; time += dt)
    {

      dataStructure.setTime(time);

        if (verbose)
            {
                std::cout << std::endl;
                std::cout << "S- Now we are at time " << dataStructure.getTime() << " s." << std::endl;
            }

        //solid.updateSystem(dZero);    // Computes the rigth hand side
        solid.updateSystem();    // Computes the rigth hand side
        solid.iterate( BCh );                  // Computes the matrices and solves the system
        //if (parameters->comm->NumProc() == 1 )  solid.postProcess(); // Post-presssing

        *solidDisp = solid.disp();
        *solidVel  = solid.vel();
        exporter->postProcess( time );

    }


}


//////////////////////


