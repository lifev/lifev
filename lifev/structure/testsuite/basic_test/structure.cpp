//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER
/**
   \file cylinder.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-19
 */

// #include <life/lifecore/LifeV.hpp>
// #include <life/lifefilters/GetPot.hpp>
// #include <life/lifecore/LifeDebug.hpp>

// #include <life/lifefilters/importer.hpp>

//#include "NavierStokesSolverIP.hpp"

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/structure/solver/VenantKirchhoffElasticData.hpp>
#include <lifev/structure/solver/VenantKirchhoffSolverLinear.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>


#include "structure.hpp"
#include <iostream>

using namespace LifeV;



Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    if (t == 0)
    {
        switch (i)
        {
        case 0:
            return -z*(z - 5.)*x/50.;
            break;
        case 1:
            return -z*(z - 5.)*y/50.;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
            break;
        }
    }

    return 0.0;
}

Real w0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch (i)
    {
    case 0:
        return 0.0;
        break;
    case 1:
        return 0.0;
        break;
    case 2:
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
    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef VenantKirchhoffSolver< mesh_Type >::vector_Type  vector_Type;
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

    MeshData             meshData;
    meshData.setup(dataFile, "solid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr( new mesh_Type( *( parameters->comm ) ) );
    readMesh(*fullMeshPtr, meshData);

    boost::shared_ptr<mesh_Type > localMeshPtr;
    {
        MeshPartitioner< mesh_Type > meshPart( fullMeshPtr, parameters->comm );
        localMeshPtr = meshPart.meshPartition();
    }

//    meshPart.rebuildMesh();


//    exit(0);


    std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");

    typedef FESpace< mesh_Type, MapEpetra > solidFESpace_type;
    typedef boost::shared_ptr<solidFESpace_type> solidFESpace_ptrtype;
    solidFESpace_ptrtype dFESpace( new solidFESpace_type(localMeshPtr,dOrder,3,parameters->comm) );
    if (verbose) std::cout << std::endl;

    VenantKirchhoffSolverLinear< mesh_Type > solid;
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

    Real dt = dataStructure->getdataTime()->timeStep();
    Real T  = dataStructure->getdataTime()->endTime();

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


    boost::shared_ptr< Exporter<mesh_Type > > exporter;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
        exporter.reset( new ExporterHDF5<mesh_Type > ( dataFile, "structure" ) );
    }
    else
#endif
    {
        if (exporterType.compare("none") == 0)
        {
            exporter.reset( new ExporterEmpty<mesh_Type > ( dataFile, localMeshPtr, "structure", parameters->comm->MyPID()) );
        }
        else
        {
            exporter.reset( new ExporterEnsight<mesh_Type > ( dataFile, localMeshPtr, "structure", parameters->comm->MyPID()) );
        }
    }

    exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( localMeshPtr, parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type(solid.getDisplacement(), exporter->mapType() ) );
    vectorPtr_Type solidVel  ( new vector_Type(solid.getVelocity(),  exporter->mapType() ) );

    exporter->addVariable( ExporterData<mesh_Type>::VectorField, "displacement",
                           dFESpace, solidDisp, UInt(0) );

    exporter->addVariable( ExporterData<mesh_Type>::VectorField, "velocity",
                           dFESpace, solidVel, UInt(0) );

    exporter->postProcess( 0 );


    //
    // Temporal loop
    //

    for (Real time = dt; time <= T; time += dt)
    {

        dataStructure->getdataTime()->setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "S- Now we are at time " << dataStructure->getdataTime()->time() << " s." << std::endl;
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
    if ( time == 0.005 && std::fabs(dispNorm-1.55991)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
    else if ( time == 0.01  && std::fabs(dispNorm-1.49237)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
    else if ( time == 0.015  && std::fabs(dispNorm-1.34538)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
    else if ( time == 0.02  && std::fabs(dispNorm-1.11341)>1e-4 )
        RESULT_CHANGED_EXCEPTION(time);
}

LifeV::UInt Structure::RESULT_CHANGED_EXCEPTION(const LifeV::Real time)
{
    LifeV::UInt value = EXIT_SUCCESS;
    std::cout << "Some modifications led to changes in the l2 norm of the solution at time" << time << std::endl;
    return value = EXIT_FAILURE;
}


//////////////////////


