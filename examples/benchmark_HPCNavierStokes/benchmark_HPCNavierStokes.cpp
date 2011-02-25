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

/*!
    @file
    @brief Benchmark to test the performance of LifeV for HPC

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 2010-08-30
 */

#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifesolver/OseenSolver.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>
//#include <fstream> // To create an output for the flux

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

typedef LifeV::MatrixEpetra<LifeV::Real> matrix_type;
typedef LifeV::RegionMesh3D<LifeV::LinearTetra>       mesh_type;
typedef LifeV::MapEpetra                              map_type;
typedef LifeV::FESpace<mesh_type,map_type>            fespace_type;
typedef boost::shared_ptr< fespace_type >             fespace_ptr;
typedef LifeV::OseenSolver< mesh_type >               fluid_type;
typedef fluid_type::vector_Type                       vector_type;
typedef boost::shared_ptr<vector_type>                vector_ptrtype;
typedef LifeV::Preconditioner                         basePrec_Type;

// +-----------------------------------------------+
// | Data and functions for the boundary conditions|
// +-----------------------------------------------+
const int INFLOW   = 1;
const int OUTFLOW1 = 2;
const int OUTFLOW2 = 3;
const int OUTFLOW3 = 4;
const int OUTFLOW4 = 5;
const int WALL     = 6;

LifeV::Real fZero( const LifeV::Real& /* t */,
                  const LifeV::Real& /* x */,
                  const LifeV::Real& /* y */,
                  const LifeV::Real& /* z */,
                  const LifeV::ID& /* i */ )
{
    return 0.0;
}

/*
   fluxes as fourier interpolation of physiological values, taken from
   Baek H, Jayaraman MV, Richardson PD, Karniadakis GE.
   Flow instability and wall shear stress variation in intracranial aneurysms.
   J R Soc Interface. 2010 Jun 6;7(47):967-88. Epub 2009 Dec 18.

*/
LifeV::Real aneurismFluxIn(const LifeV::Real&  t, const LifeV::Real& /*x*/, const LifeV::Real& /*y*/, const LifeV::Real& /*z*/, const LifeV::ID& /*i*/)
{
    // We change the flux for our geometry
    const LifeV::Real pi         = 3.141592653589793;
    const LifeV::Real area       = 0.0907122;
    const LifeV::Real areaFactor = area/((0.6/2)*(0.6/2)*pi);

    // Unit conversion from ml/min to cm^3/s
    const LifeV::Real unitFactor = 1./ 60.;

    // T is the period of the cardiac cycle
    const LifeV::Real T          = 1.0;

    // a0 is the average VFR (the value is taken from Karniadakis p970)
    const LifeV::Real a0         = 255;

    // Fourrier
    const LifeV::Int M(7);
    const LifeV::Real a[M] = {-0.152001,-0.111619, 0.043304,0.028871,0.002098,-0.027237,-0.000557};
    const LifeV::Real b[M] = { 0.129013,-0.031435,-0.086106,0.028263,0.010177, 0.012160,-0.026303};

    LifeV::Real flux(0);
    const LifeV::Real xi(2*pi*t/T);

    flux = a0;
    int k(1);
    for (; k<=M ; ++k)
        flux += a0*(a[k-1]*cos(k*xi) + b[k-1]*sin(k*xi));

    return - flux * areaFactor * unitFactor;
}

LifeV::Real aneurismFluxIn2(const LifeV::Real&  t, const LifeV::Real& x, const LifeV::Real& y, const LifeV::Real& z, const LifeV::ID& i)
{
    LifeV::Real n1(-0.000019768882940);
    LifeV::Real n2(-0.978289616345544);
    LifeV::Real n3( 0.207242433298975);
    LifeV::Real flux(aneurismFluxIn(t,x,y,z,i));
    LifeV::Real area(0.0907122); // Computed with the triangle on the INLET boundary

    switch(i) {
    case 0:
        return n1*flux/area;
    case 1:
        return n2*flux/area;
    case 2:
        return n3*flux/area;
    default:
        return 0.0;
    }
}


int main(int argc, char** argv)
{

    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

    boost::shared_ptr<Epetra_Comm>   comm;
#ifdef EPETRA_MPI
    comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
    comm.reset( new Epetra_SerialComm() );
#endif

    bool verbose(false);
    if(comm->MyPID() == 0){
        verbose = true;
        std::cout
            << " +-----------------------------------------------+" << std::endl
            << " |           Fluid benchmark for LifeV           |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl
            << " +-----------------------------------------------+" << std::endl
            << " |           Author: Gwenol Grandperrin          |" << std::endl
            << " |             Date: 2010-10-28                  |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl;

        std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
        std::cout << "Using MPI (" << nproc << " proc.)" << std::endl;
#else
        std::cout << "Using serial version" << std::endl;
#endif
    }

    LifeV::LifeChrono globalChrono;
    LifeV::LifeChrono initChrono;
    LifeV::LifeChrono iterChrono;
    globalChrono.start();
    initChrono.start();

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the data]" << std::endl;
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data", 2, "-f","--file");
    GetPot dataFile(dataFileName);

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the mesh]" << std::endl;
    LifeV::MeshData dataMesh;
    dataMesh.setup(dataFile, "fluid/space_discretization");
    if (verbose) std::cout << "Mesh file: " << dataMesh.meshDir() << dataMesh.meshFile() << std::endl;

    boost::shared_ptr< LifeV::RegionMesh3D<LifeV::LinearTetra> > fullMeshPtr(new LifeV::RegionMesh3D<LifeV::LinearTetra>);
    LifeV::readMesh(*fullMeshPtr, dataMesh);

    LifeV::MeshPartitioner< LifeV::RegionMesh3D<LifeV::LinearTetra> >   meshPart(fullMeshPtr, comm);

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order",   "P2");
    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");
    if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                           << "FE for the pressure: " << pOrder << std::endl;

    if (verbose) std::cout << "Building the velocity FE space... " << std::flush;
    fespace_ptr uFESpace (new fespace_type(meshPart, uOrder, 3, comm));
    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Building the pressure FE space... " << std::flush;
    fespace_ptr pFESpace(new fespace_type(meshPart,pOrder,1,comm));
    if (verbose) std::cout << "ok." << std::endl;

    LifeV::UInt totalVelDof   = uFESpace->map().map(LifeV::Unique)->NumGlobalElements();
    LifeV::UInt totalPressDof = pFESpace->map().map(LifeV::Unique)->NumGlobalElements();

    if (verbose) std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Boundary conditions]" << std::endl;

    LifeV::BCFunctionBase uZero(fZero);
    LifeV::BCFunctionBase uFlux(aneurismFluxIn2);

    LifeV::BCHandler bcH;
    bcH.addBC( "Inflow"  , INFLOW  , LifeV::Essential, LifeV::Full, uFlux, 3 );
    bcH.addBC( "Outflow1", OUTFLOW1, LifeV::Natural  , LifeV::Full, uZero, 3 );
    bcH.addBC( "Outflow2", OUTFLOW2, LifeV::Natural  , LifeV::Full, uZero, 3 );
    bcH.addBC( "Outflow3", OUTFLOW3, LifeV::Natural  , LifeV::Full, uZero, 3 );
    bcH.addBC( "Outflow4", OUTFLOW4, LifeV::Natural  , LifeV::Full, uZero, 3 );
    bcH.addBC( "Wall"    , WALL    , LifeV::Essential, LifeV::Full, uZero, 3 );

    // Getting the number of Lagrange Multiplyers (LM)
    // and setting the offsets
    std::vector<LifeV::bcName_Type> fluxVector = bcH.findAllBCWithType( LifeV::Flux );
    LifeV::UInt numLM = static_cast<LifeV::UInt>( fluxVector.size() );

    LifeV::UInt offset = uFESpace->map().map(LifeV::Unique)->NumGlobalElements()
                            + pFESpace->map().map(LifeV::Unique)->NumGlobalElements();

    for ( LifeV::UInt i = 0; i < numLM; ++i )
        bcH.setOffset( fluxVector[i], offset + i );

    // +-----------------------------------------------+
    // |             Creating the problem              |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Creating the problem]" << std::endl;
    boost::shared_ptr<LifeV::OseenData> oseenData(new LifeV::OseenData());
    oseenData->setup( dataFile );

    if (verbose) std::cout << "Time discretization order " << oseenData->dataTime()->orderBDF() << std::endl;

    LifeV::OseenSolver< LifeV::RegionMesh3D<LifeV::LinearTetra> > fluid (oseenData,
                                                                       *uFESpace,
                                                                       *pFESpace,
                                                                       comm,
                                                                       numLM);
    fluid.setUp(dataFile);
    // ==>Call setupPreconditioner from SolverAztecOO
    //    ==> Call the factory to build the preconditioner
    //        Call setSolver from the preconditioner
    //        Call setDataFromGetPot from the preconditioner

    fluid.buildSystem();

    LifeV::MapEpetra fullMap(fluid.getMap());

    MPI_Barrier(MPI_COMM_WORLD);

    // +-----------------------------------------------+
    // |       Initialization of the simulation        |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Initialization of the simulation]" << std::endl;
    LifeV::Real dt     = oseenData->dataTime()->timeStep();
    LifeV::Real t0     = oseenData->dataTime()->initialTime();
    LifeV::Real tFinal = oseenData->dataTime()->endTime ();

    // bdf object to store the previous solutions
    LifeV::TimeAdvanceBDFNavierStokes<vector_type> bdf;
    bdf.setup(oseenData->dataTime()->orderBDF());

    // initialization with exact solution: either interpolation or "L2-NS"-projection
    t0 -= dt * bdf.bdfVelocity().order();

    vector_type beta( fullMap );
    vector_type rhs ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    oseenData->dataTime()->setTime(t0);

    beta *= 0.;
    rhs  *= 0.;
    fluid.updateSystem(0.0,beta,rhs);
    fluid.setTolMaxIteration(1e-6,400); // Set the tolerance of the solver
    fluid.iterate(bcH);
    bdf.bdfVelocity().setInitialCondition(*fluid.solution());

    LifeV::Real time = t0 + dt;
    for (  ; time <=  oseenData->dataTime()->initialTime() + dt/2.; time += dt)
    {
        oseenData->dataTime()->setTime(time);

        fluid.updateSystem(0.0,beta,rhs);
        fluid.iterate(bcH);
        bdf.bdfVelocity().shiftRight( *fluid.solution() );
    }

    fluid.resetPreconditioner();

    boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh3D<LifeV::LinearTetra> > > exporter;

    vector_ptrtype velAndPressure;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
        exporter.reset( new LifeV::ExporterHDF5<LifeV::RegionMesh3D<LifeV::LinearTetra> > ( dataFile, "benchmark_HPCNavierStokes" ) );
        exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
        exporter->setMeshProcId( meshPart.meshPartition(), comm->MyPID() );
    }
    else
#endif
    {
        if (exporterType.compare("none") == 0)
        {
            exporter.reset( new LifeV::ExporterEmpty<LifeV::RegionMesh3D<LifeV::LinearTetra> > ( dataFile, meshPart.meshPartition(), "benchmark_HPCNavierStokes", comm->MyPID()) );
        } else {
            exporter.reset( new LifeV::ExporterEnsight<LifeV::RegionMesh3D<LifeV::LinearTetra> > ( dataFile, meshPart.meshPartition(), "benchmark_HPCNavierStokes", comm->MyPID()) );
        }
    }

    velAndPressure.reset( new vector_type(*fluid.solution(), exporter->mapType() ) );

    exporter->addVariable( LifeV::ExporterData::Vector, "velocity", velAndPressure,
                           LifeV::UInt(0), uFESpace->dof().numTotalDof() );

    exporter->addVariable( LifeV::ExporterData::Scalar, "pressure", velAndPressure,
                           LifeV::UInt(3*uFESpace->dof().numTotalDof()),
                           LifeV::UInt(pFESpace->dof().numTotalDof()) );
    exporter->postProcess( 0 );

    initChrono.stop();
    if (verbose) std::cout << "Initialization time:  " << initChrono.diff() << " s." << std::endl;

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    if (verbose) std::cout<< std::endl << "[Solving the problem]" << std::endl;

    // We now change the tolerance of the solver
    fluid.setTolMaxIteration(1e-10,400);

    int iter = 1;

    for ( ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        oseenData->dataTime()->setTime(time);

        if (verbose) std::cout << "[t = "<< oseenData->dataTime()->time() << " s.]" << std::endl;

        iterChrono.start();

        double alpha = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / oseenData->dataTime()->timeStep();

        beta = bdf.bdfVelocity().extrapolation(); // Extrapolation for the convective term

        bdf.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
        rhs  = fluid.matrixMass()*bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.getDisplayer().leaderPrint("alpha ", alpha);
        fluid.getDisplayer().leaderPrint("\n");
        fluid.getDisplayer().leaderPrint("norm beta ", beta.norm2());
        fluid.getDisplayer().leaderPrint("\n");
        fluid.getDisplayer().leaderPrint("norm rhs  ", rhs.norm2());
        fluid.getDisplayer().leaderPrint("\n");

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdfVelocity().shiftRight( *fluid.solution() );

        // Computation of the error
        vector_type vel  (uFESpace->map(), LifeV::Repeated);
        vector_type press(pFESpace->map(), LifeV::Repeated);
        vector_type velpressure ( *fluid.solution(), LifeV::Repeated );

        velpressure = *fluid.solution();
        vel.subset(velpressure);
        press.subset(velpressure, uFESpace->dim()*uFESpace->fieldDim());


        bool verbose = (comm->MyPID() == 0);


        // Exporting the solution
        *velAndPressure = *fluid.solution();
        exporter->postProcess( time );


        MPI_Barrier(MPI_COMM_WORLD);

        iterChrono.stop();
        if (verbose) std::cout << "Iteration time: " << iterChrono.diff() << " s." << std::endl << std::endl;
    }

    globalChrono.stop();
    if (verbose) std::cout << "Total simulation time:  " << globalChrono.diff() << " s." << std::endl;

#ifdef HAVE_HDF5
    exporter->closeFile();
#endif

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
