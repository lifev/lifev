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

// #include <life/lifecore/life.hpp>

// #include <life/lifecore/life.hpp>
// #include <life/lifecore/GetPot.hpp>
// #include <life/lifecore/debug.hpp>

// #include <life/lifefilters/importer.hpp>

//#include "NavierStokesSolverIP.hpp"

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataElasticStructure.hpp>
#include <life/lifesolver/LinearVenantKirchhofSolver.hpp>

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
        switch (i)
        {
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

    switch (i)
    {
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

    boost::shared_ptr<Epetra_Comm>     comm;


};

Structure::Structure( int                                   argc,
                      char**                                argv,
                      boost::shared_ptr<Epetra_Comm>        structComm):
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

    parameters->comm = structComm;
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

    typedef VenantKirchhoffSolver< RegionMesh3D<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    bool verbose = (parameters->comm->MyPID() == 0);
    // Number of boundary conditions for the velocity and mesh motion
    //
    boost::shared_ptr<BCHandler> BCh( new BCHandler() );
    BCFunctionBase dZero( zero_scalar );
    //
    // dataElasticStructure
    //

    GetPot dataFile( parameters->data_file_name.c_str() );

    boost::shared_ptr<VenantKirchhoffElasticData> dataStructure(new VenantKirchhoffElasticData( ));
    dataStructure->setup(dataFile);

    DataMesh             dataMesh;
    dataMesh.setup(dataFile, "solid/space_discretization");

    boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);
    readMesh(*fullMeshPtr, dataMesh);


    partitionMesh< RegionMesh3D<LinearTetra> > meshPart( fullMeshPtr, parameters->comm );

//    meshPart.rebuildMesh();


//    exit(0);


    std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");

    typedef FESpace< RegionMesh3D<LinearTetra>, EpetraMap > solidFESpace_type;
    typedef boost::shared_ptr<solidFESpace_type> solidFESpace_ptrtype;
    solidFESpace_ptrtype dFESpace( new solidFESpace_type(meshPart,dOrder,3,parameters->comm) );
    if (verbose) std::cout << std::endl;

    EpetraMap structMap(dFESpace->refFE(), meshPart, parameters->comm);

    EpetraMap fullMap;

    for (UInt ii = 0; ii < nDimensions; ++ii)
    {
        fullMap += structMap;
    }

    VenantKirchhoffSolverLinear< RegionMesh3D<LinearTetra> > solid;
    solid.setup(dataStructure,
                dFESpace,
                parameters->comm);

    solid.setDataFromGetPot(dataFile);
    solid.buildSystem();
    //
    // Boundary conditions for the displacement
    //
    BCFunctionBase fixed(dZero);

    BCh->addBC("Base2 ", 2 , Essential, Full, fixed, 3);
    BCh->addBC("Base3 ", 3 , Essential, Full, fixed, 3);

    //
    // Temporal data and initial conditions
    //

    Real dt = dataStructure->getDataTime()->timeStep();
    Real T  = dataStructure->getDataTime()->endTime();

    vectorPtr_Type disp(new vector_Type(solid.getDisplacement(), Unique));
    vectorPtr_Type vel(new vector_Type(solid.getVelocity(), Unique));

    dFESpace->interpolate(d0, *disp, 0.0);
    dFESpace->interpolate(w0, *vel , 0.0);

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
            exporter.reset( new NoExport<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID()) );
        }
        else
        {
            exporter.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID()) );
        }
    }

    exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( meshPart.meshPartition(), parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type(solid.getDisplacement(), exporter->mapType() ) );
    vectorPtr_Type solidVel  ( new vector_Type(solid.getVelocity(),  exporter->mapType() ) );

    exporter->addVariable( ExporterData::Vector, "displacement", solidDisp,
                           UInt(0), dFESpace->dof().numTotalDof() );

    exporter->addVariable( ExporterData::Vector, "velocity", solidVel,
                           UInt(0), dFESpace->dof().numTotalDof() );

    exporter->postProcess( 0 );


    //
    // Temporal loop
    //

    for (Real time = dt; time <= T; time += dt)
    {

        dataStructure->getDataTime()->setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "S- Now we are at time " << dataStructure->getDataTime()->time() << " s." << std::endl;
        }

        //solid.updateSystem(dZero);    // Computes the rigth hand side
        solid.updateSystem();    // Computes the rigth hand side
        solid.iterate( BCh );                  // Computes the matrices and solves the system
        //if (parameters->comm->NumProc() == 1 )  solid.postProcess(); // Post-presssing

        *solidDisp = solid.getDisplacement();
        *solidVel  = solid.getVelocity();
            CheckResults(solid.getDisplacement().norm2(),time);

        exporter->postProcess( time );

    }


}


void Structure::CheckResults(const LifeV::Real& dispNorm,const LifeV::Real& time)
{
    if ( time == 0.005 && abs(dispNorm-1.55991)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
    else if ( time == 0.01  && abs(dispNorm-1.49237)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
    else if ( time == 0.015  && abs(dispNorm-1.34538)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
    else if ( time == 0.02  && abs(dispNorm-1.11341)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
}

LifeV::UInt Structure::RESULT_CHANGED_EXCEPTION(const LifeV::Real time)
  {
      LifeV::UInt value = EXIT_SUCCESS;
      std::cout << "Some modifications led to changes in the l2 norm of the solution at time" << time << std::endl;
      return value = EXIT_FAILURE;
  }


//////////////////////


