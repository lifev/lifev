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
 *  @file
 *  @brief File containing the Monolithic Test
 *
 *  @date 2009-04-09
 *  @author Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @contributor Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 * Monolithic problem. Features:
 * - fullMonolithic (CE):
 *  -# solution with exact Newton (semiImplicit = false, useShapeDerivatives = true, domainVelImplicit = false, convectiveTermDer = false)
 *  -# solution with quasi Newton (semiImplicit = false, useShapeDerivatives = false, domainVelImplicit = false, convectiveTermDer = false)
 *  -# preconditioner choice: see the classes Monolithic and fullMonolithic
 * - Monolithic (GCE):
 *  -# solution extrapolating the fluid domain (semiImplicit = true, useShapeDerivatives = false, domainVelImplicit = false, convectiveTermDer = false)
 *  -# preconditioner choice: see the classes Monolithic and fullMonolithic
 * - boundary conditions:
 *  -# Neumann
 *  -# Dirichlet
 *  -# Robin
 *  -# Fluxes (defective)
 *  -# absorbing \cite BadiaNobileVergara2008 :
 *   through the class flowConditions.
 * - optional: computation of wall shear stress (not properly tested in parallel)
 * - optional: computation of the largest singular values of the preconditioned matrix
 *
 * \b Features:
 * This test by default solves the FSI probem discretized in time using the GCE or CE methods, implemented respectively
 * in the files monolithicGE.hpp and monolithicGI.hpp . The geometry is that of a tube (benchmark test introduced in \cite Gerbeau2003).
 * In this test the boundary conditions assigned are of type:
 * - flux (defective b.c.) at the inlet
 * - absorbing (see \cite BadiaNobileVergara2008) at the outlet
 * - Robin b.c. on the solid external wall
 * - Dirichlet homogeneous at the solid rings on the inlet-outlet (clamped tube).
 *
 * The output is written at every timestep, in HDF5 (if available) or ensight format.
 * This test implements an inlet flux bundary condition for the first three time steps, then at the fourth time step
 * the inlet boundary condition is replaced by a Neumann one (this mechanism is useful to implement rudimental valves).
 * The outflow boundary condition is of absorbing type. At the outer wall for the structure a Robin condition is imposed.
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <cassert>
#include <cstdlib>

#include <boost/timer.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV includes
#include <life/lifefem/FESpace.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifefem/BCHandler.hpp>
#include <life/lifecore/LifeV.hpp>

#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>

#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/StructuralSolver.hpp>
#include <life/lifesolver/FSIMonolithicGI.hpp>

#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/ExporterHDF5.hpp>
#endif

#include "boundaryConditions.hpp"

class Problem
{
public:

    typedef boost::shared_ptr<LifeV::FSISolver> FSISolverPtr_Type;

    typedef LifeV::FSIOperator::data_Type                          data_Type;
    typedef LifeV::FSIOperator::dataPtr_Type                       dataPtr_Type;

    typedef LifeV::FSIOperator::vector_Type        vector_Type;
    typedef LifeV::FSIOperator::vectorPtr_Type     vectorPtr_Type;

    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;

    typedef LifeV::ExporterEnsight<LifeV::FSIOperator::mesh_Type>  ensightFilter_Type;
    typedef boost::shared_ptr<ensightFilter_Type>                 ensightFilterPtr_Type;
#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<LifeV::FSIOperator::mesh_Type>  hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                  hdf5FilterPtr_Type;
#endif
    typedef LifeV::FactorySingleton<LifeV::Factory<LifeV::FSIOperator,  std::string> > FSIFactory_Type;
  typedef LifeV::RegionMesh<LifeV::LinearTetra> mesh_Type;
  typedef LifeV::MapEpetra map_Type;
  typedef LifeV::FESpace<mesh_Type, map_Type> fespace_Type;
    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */

    Problem( GetPot const& data_file ):
        M_Tstart(0.),
        M_saveEvery(1)
    {
        using namespace LifeV;

        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhof", &FSIOperator::createVenantKirchhoffLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "exponential", &FSIOperator::createExponentialMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "neoHookean", &FSIOperator::createNeoHookeanMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "nonLinearVenantKirchhof", &FSIOperator::createVenantKirchhoffNonLinear );

        std::cout<<"register MonolithicGE : "<<FSIMonolithicGE::S_register<<std::endl;
        std::cout<<"register MonolithicGI : "<<FSIMonolithicGI::S_register<<std::endl;

        M_data = dataPtr_Type( new data_Type() );
        M_data->setup( data_file );

        M_fsi = FSISolverPtr_Type( new FSISolver( ) );
        MPI_Barrier( MPI_COMM_WORLD );

        M_fsi->setData( M_data );
        M_fsi->FSIOper()->setDataFile( data_file ); //TO BE REMOVED!

        // Setting FESpace and DOF

        std::string  fluidMeshPartitioned    =  data_file( "problem/fluidMeshPartitioned", "none" );
        std::string  solidMeshPartitioned    =  data_file( "problem/solidMeshPartitioned", "none" );
#ifdef HAVE_HDF5
        if ( fluidMeshPartitioned.compare( "none" ) )
        {
            FSIOperator::meshFilter_Type fluidMeshFilter( data_file, fluidMeshPartitioned );
            fluidMeshFilter.setComm( M_fsi->FSIOper()->worldComm() );
            FSIOperator::meshFilter_Type solidMeshFilter( data_file, solidMeshPartitioned );
            solidMeshFilter.setComm( M_fsi->FSIOper( )->worldComm( ) );
            M_fsi->FSIOper( )->partitionMeshes( fluidMeshFilter, solidMeshFilter );
            M_fsi->FSIOper( )->setupFEspace( );
            M_fsi->FSIOper( )->setupDOF( fluidMeshFilter );
            fluidMeshFilter.closeFile( );
            solidMeshFilter.closeFile( );
        }
        else
#endif
        {
            M_fsi->FSIOper( )->partitionMeshes( );
            M_fsi->FSIOper( )->setupFEspace( );
            M_fsi->FSIOper( )->setupDOF( );
        }

        MPI_Barrier( MPI_COMM_WORLD );

        M_fsi->setFluidBC( BCh_monolithicFlux( true ) );
        M_fsi->setSolidBC( BCh_monolithicSolid( *M_fsi->FSIOper( ) ) );

        M_fsi->setup();

        M_fsi->setFluidBC( BCh_monolithicFluid( *M_fsi->FSIOper( )) );
        M_fsi->setHarmonicExtensionBC( BCh_harmonicExtension( *M_fsi->FSIOper( ) ) );

        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->mergeBCHandlers();

        std::string const exporterType =  data_file( "exporter/type", "ensight" );
        std::string const fluidName    =  data_file( "exporter/fluid/filename", "fluid" );
        std::string const solidName    =  data_file( "exporter/solid/filename", "solid" );

#ifdef HAVE_HDF5
        if (exporterType.compare("hdf5") == 0)
        {
            M_exporterFluid.reset( new  hdf5Filter_Type( data_file, fluidName) );
            M_exporterSolid.reset( new  hdf5Filter_Type ( data_file,solidName));
        }
        else
#endif
        {
            if (exporterType.compare("none") == 0)
            {
                M_exporterFluid.reset( new ExporterEmpty<RegionMesh<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), fluidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
                M_exporterSolid.reset( new ExporterEmpty<RegionMesh<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), solidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
            }
            else
            {
                M_exporterFluid.reset( new  ensightFilter_Type( data_file, fluidName) );
                M_exporterSolid.reset( new  ensightFilter_Type ( data_file, solidName) );
            }
        }


        // load using ensight/hdf5
        M_saveEvery=data_file("exporter/saveEvery",1);

        // load using ensight/hdf5
        std::string restartType(data_file("importer/restartType","newSimulation"));
        std::cout << "The load state is: "<< restartType << std::endl;

        if (!restartType.compare("Stokes"))
            initializeStokes(data_file);
        else if (!restartType.compare("restartFSI"))
            restartFSI(data_file);
        else
        {
            M_fsi->initialize();

            M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluid->mapType() ));

            M_fluidDisp.reset     ( new vector_Type( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ));

            M_solidDisp.reset( new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        }


        M_exporterFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());
        M_exporterSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID());
        M_exporterFluid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-velocity",
                                      M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressure, UInt(0) );
        M_exporterFluid->addVariable( ExporterData<FSIOperator::mesh_Type>::ScalarField, "f-pressure",
                                      M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressure,
                                      UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()) );

        M_exporterFluid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-displacement",
                                      M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDisp, UInt(0) );




        M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-displacement",
                                      M_fsi->FSIOper()->dFESpacePtr(), M_solidDisp, UInt(0) );

	M_fsi->FSIOper()->fluid().setupPostProc(); //this has to be called if we want to initialize the postProcess

        FC0.initParameters( *M_fsi->FSIOper(), 3);

        M_data->dataFluid()->dataTime()->setInitialTime( M_Tstart );
        M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setInitialTime( M_Tstart );
        M_data->dataSolid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataALE()->setInitialTime( M_Tstart );
        M_data->dataALE()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
    }

  void errors(GetPot const& dataFile);

    /*!
      This routine runs the temporal loop
    */
    void
    run()
    {
        boost::timer _overall_timer;

        LifeV::UInt iter = 1;
        LifeV::UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->enableStressComputation(1);

// #ifdef HAVE_HDF5
//         if (M_exporterFluid->mapType() == LifeV::Unique)
//         {
//             M_exporterFluid->postProcess( M_Tstart-M_data->dataFluid()->dataTime()->getTimeStep() );//ugly way to avoid that hdf5 starts with a deformed mesh
//             M_exporterSolid->postProcess( M_Tstart-M_data->dataSolid()->dataTime()->getTimeStep() );//ugly way to avoid that hdf5 starts with a deformed mesh
//         }
// #endif

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(),M_data->dataSolid()->dataTime()->updateTime(), ++iter)
        {
            boost::timer _timer;

            int flag =2;
            FC0.renewParameters( *M_fsi, 3 );

            if(iter%M_saveEvery==0)
            {
                M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDisp);//    displacement(), M_offset);
                //M_fsi->FSIOper()->exportSolidVelocity(*M_solidVel);//    displacement(), M_offset);

                M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressure);
                M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );

                //*M_WS= *(dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->/*WS());//*/computeStress());

                *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();
                M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );
            }

            M_fsi->iterate();


	    M_fsi->FSIOper()->displayer().leaderPrintMax("[fsi_run] Iteration ", iter);
	    M_fsi->FSIOper()->displayer().leaderPrintMax(" was done in : ", _timer.elapsed());

//             std::cout << "solution norm " << iter << " : "
//                       << M_fsi->displacement().norm2() << "\n";

        }
    //     if (M_data->method().compare("monolithicGI"))
    //     {
    //         M_fsi->FSIOper()->iterateMesh(M_fsi->displacement());

    //         M_solidDisp->subset(M_fsi->displacement(), offset);

    //         *M_velAndPressure = M_fsi->displacement();
    //         M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );
    //         *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
    //         M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );
    //     }

    //     std::cout << "Total computation time = "
    //               << _overall_timer.elapsed() << "s" << "\n";

    }

private:

    void initializeStokes(  GetPot const& data_file);
    void restartFSI(  GetPot const& data_file);

    FSISolverPtr_Type M_fsi;
    dataPtr_Type   M_data;

    filterPtr_Type M_exporterSolid;
    filterPtr_Type M_exporterFluid;
    filterPtr_Type M_importerSolid;
    filterPtr_Type M_importerFluid;
    vectorPtr_Type M_velAndPressure;
    vectorPtr_Type M_fluidDisp;
    vectorPtr_Type M_solidDisp;

    std::vector<vectorPtr_Type> M_solidStencil;
    std::vector<vectorPtr_Type> M_fluidStencil;
    std::vector<vectorPtr_Type> M_ALEStencil;

    LifeV::FlowConditions FC0;

    LifeV::Real    M_Tstart;
    // vectorPtr_Type M_WS;
    LifeV::UInt           M_saveEvery;

    struct RESULT_CHANGED_EXCEPTION
    {
public:
        RESULT_CHANGED_EXCEPTION(LifeV::Real time)
        {
            std::cout<<"Some modifications led to changes in the l2 norm of the solution at time"<<time<<std::endl;
        }
    };
};

struct FSIChecker
{
    FSIChecker( GetPot const& _data_file ):
            data_file( _data_file )
    {}

    void operator()()
    {
        boost::shared_ptr<Problem> fsip;

        try
        {
            fsip = boost::shared_ptr<Problem>( new Problem( data_file ) );

	    if(data_file( "problem/computeErrors", false))
	      fsip->errors(data_file);
	    else
	      fsip->run();

        }
        catch ( std::exception const& _ex )
        {
            std::cout << "caught exception :  " << _ex.what() << "\n";
        }

    }

    GetPot                data_file;
    LifeV::Vector         disp;
};



namespace LifeV
{
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}
}


int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#else
    std::cout << "% using serial Version" << std::endl;
#endif

    GetPot command_line(argc,argv);
    const bool check = command_line.search(2, "-c", "--check");

    if (check)
    {
        GetPot data_fileGCE("dataGCE");
        FSIChecker _GCE_check( data_fileGCE );
        _GCE_check();

        GetPot data_fileCE("dataCE");
        FSIChecker _CE_check(data_fileCE);
        _CE_check();

#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return 0;
    }
    else
    {
        const std::string data_file_name = command_line.follow("data", 2, "-f","--file");
        GetPot data_file(data_file_name);
        FSIChecker _sp_check( data_file );
        _sp_check();
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;

}

void Problem::restartFSI(  GetPot const& data_file)
{

  using namespace LifeV;

  typedef FSIOperator::mesh_Type        mesh_Type;

  //Creating the pointer to the filter
  std::string const importerType =  data_file( "importer/type", "ensight");
  std::string const fluidName    =  data_file( "importer/fluid/filename", "fluid");
  std::string const solidName    =  data_file( "importer/solid/filename", "solid");

  std::string const loadInitSol      = data_file( "importer/initSol", "00000");
  std::string const loadInitSolFD    = data_file("importer/initSolFD","-1");
  std::string iterationString;

  M_Tstart  = data_file( "fluid/time_discretization/initialtime", 0.);

  std::cout << "The file for fluid is    : " << fluidName << std::endl;
  std::cout << "The file for solid is    : " << solidName << std::endl;
  std::cout << "The importerType is      : " << importerType << std::endl;
  std::cout << "The iteration is         : " << loadInitSol << std::endl;
  std::cout << "For the fluid disp is    : " << loadInitSolFD << std::endl;
  std::cout << "Starting time            : " << M_Tstart << std::endl;

#ifdef HAVE_HDF5
  if (importerType.compare("hdf5") == 0)
    {
      M_importerFluid.reset( new  hdf5Filter_Type( data_file, fluidName) );
      M_importerSolid.reset( new  hdf5Filter_Type( data_file, solidName) );
    }
  else
#endif
    {
      if (importerType.compare("none") == 0)
	{
	  M_importerFluid.reset( new ExporterEmpty<mesh_Type > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), "fluid", M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
	  M_importerSolid.reset( new ExporterEmpty<mesh_Type > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), "solid", M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
	}
      else
	{
	  M_importerFluid.reset( new  ensightFilter_Type( data_file, fluidName) );
	  M_importerSolid.reset( new  ensightFilter_Type ( data_file, solidName) );
	}
    }

  M_importerFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());
  M_importerSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID());

  //Each of the vectors of the stencils has the dimension of the big vector
  //Performing a cycle of the size of the timeAdvance classes for each problem
  //The size of TimeAdvanceClass depends on the order of the BDF (data file)

  //It should work just initializing the timeAdvance classes
  //Three stencils are needed (Fluid-Structure-Geometric)


  vectorPtr_Type fluidSol         (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique));
  vectorPtr_Type initFluid         (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
  vectorPtr_Type HarmonicSol      (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
  vectorPtr_Type structureSol      (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
  vectorPtr_Type temporarySol     (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));

  vectorPtr_Type vel           (new vector_Type(M_fsi->FSIOper()->uFESpace().map(), M_importerFluid->mapType()));
  vectorPtr_Type pressure      (new vector_Type(M_fsi->FSIOper()->pFESpace().map(), M_importerFluid->mapType()));
  vectorPtr_Type solidDisp    (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), M_importerSolid->mapType()));
  vectorPtr_Type fluidDisp(new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), M_importerFluid->mapType()));

  vectorPtr_Type firstFluidDisp(new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), M_importerFluid->mapType()));

  UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

    iterationString = loadInitSol;
    std::cout << "Fluid size TimeAdvance:" << M_fsi->FSIOper()->fluidTimeAdvance()->size() << std::endl;
    for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->fluidTimeAdvance()->size(); iterInit++ )
    {

      /*!
	definition of the vector to fill with the initialization.
      */
      temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
      fluidSol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique));
      initFluid.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));


      LifeV::ExporterData<mesh_Type> initSolFluidVel   (LifeV::ExporterData<mesh_Type>::VectorField, std::string("f-velocity."+iterationString), M_fsi->FSIOper()->uFESpacePtr(), vel, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
      LifeV::ExporterData<mesh_Type> initSolFluidPress (LifeV::ExporterData<mesh_Type>::ScalarField, std::string("f-pressure."+iterationString), M_fsi->FSIOper()->pFESpacePtr(), pressure, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

      /*!load of the solutions*/
      M_importerFluid->readVariable(initSolFluidVel);
      M_importerFluid->readVariable(initSolFluidPress);

      fluidSol.reset( new vector_Type(*pressure, Unique, Zero));
      initFluid->subset(*fluidSol, fluidSol->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());
      vector_Type tmpVec(*initFluid);
      //temporarySol->subset(*fluidSol, fluidSol->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());
      fluidSol.reset( new vector_Type(*vel, Unique, Zero));
      tmpVec=*fluidSol;
      *initFluid += tmpVec;

      // std::string firstFl="firstFluid";
      // initFluid->spy(firstFl);

      std::cout << "Norm of first Fluid sol: "<< initFluid->norm2() << std::endl;
      *temporarySol = *initFluid;

    //   if(!iterInit)
	// {
	//   *un += *temporarySol;
    //     }

      M_fluidStencil.push_back(temporarySol);
      //Updating string name
      int iterations = std::atoi(iterationString.c_str());
      iterations--;

      std::ostringstream iter;
      iter.fill( '0' );
      iter << std::setw(5) << ( iterations );
      iterationString=iter.str();

    }

  iterationString = loadInitSol;
  for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->solidTimeAdvance()->size(); iterInit++ )
    {
        /*!
          definition of the vector to fill with the initialization.
        */
      temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
      structureSol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));


      /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
      LifeV::ExporterData<mesh_Type> initSolSolidDisp  (LifeV::ExporterData<mesh_Type>::VectorField,"s-displacement."+iterationString, M_fsi->FSIOper()->dFESpacePtr(), solidDisp, UInt(0)/*offset*/, LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

        /*!load of the solutions*/
        M_importerSolid->readVariable(initSolSolidDisp);

        structureSol->subset(*solidDisp, solidDisp->map(), (UInt)0, offset);

	*temporarySol = *structureSol/(M_fsi->FSIOper()->solid().rescaleFactor());

	M_solidStencil.push_back(temporarySol);

	//Updating the string name for the next iteration
        UInt iterations = std::atoi(iterationString.c_str());
        iterations--;

        std::ostringstream iter;
        iter.fill( '0' );
        iter << std::setw(5) << ( iterations );
        iterationString=iter.str();

    }

  Int convectiveTerm=0; //loadInitSol;
  if(!M_data->dataFluid()->domainVelImplicit())
    convectiveTerm=1;
  HarmonicSol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));

  iterationString = loadInitSolFD;

  for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->ALETimeAdvance()->size()+convectiveTerm+1; iterInit++ )
  {
      temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV:: Unique, Zero));
      /*!
        definition of the vector to fill with the initialization.
      */

        LifeV::ExporterData<mesh_Type> initSolFluidDisp  (LifeV::ExporterData<mesh_Type>::VectorField, "f-displacement."+iterationString, M_fsi->FSIOper()->mmFESpacePtr(), fluidDisp, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

        /*!load of the solutions*/
        M_importerFluid->readVariable(initSolFluidDisp);

	std::cout << "Reloaded Harmonic sol norm: "<< fluidDisp->norm2() << std::endl;
	if(iterInit == 0)
	  {
	    HarmonicSol->subset(*fluidDisp, fluidDisp->map(), (UInt)0, dynamic_cast<LifeV::FSIMonolithicGI*>(M_fsi->FSIOper().get())->mapWithoutMesh().map(Unique)->NumGlobalElements());
	    if(!convectiveTerm)
	      {
		*firstFluidDisp = *fluidDisp;
	      }
	  }
	else
	  if(convectiveTerm && iterInit == 1)
	    *firstFluidDisp = *fluidDisp;
	  else
	    M_ALEStencil.push_back(/*temporarySol*/fluidDisp);
	//Updating the iteration String name
	int iterations = std::atoi(iterationString.c_str());
	iterations--;

	std::ostringstream iter;
	iter.fill( '0' );
	iter << std::setw(5) << ( iterations );
	iterationString=iter.str();
  }

    *M_fluidStencil[0]+=*M_solidStencil[0];
    *M_fluidStencil[0]+=*HarmonicSol;
    //this is going to be the global solutions returned by the method solution()

    M_fsi->initialize(M_fluidStencil, M_solidStencil, M_ALEStencil);

    if(!M_data->dataFluid()->domainVelImplicit())
      {
	//The following is needed because (and if) of the extrapolation of the fluid domain velocity is used, i.e. M_domainVelImplicit
	M_fsi->FSIOper()->ALETimeAdvance()->updateRHSFirstDerivative( M_data->dataSolid()->dataTime()->timeStep() );
	M_fsi->FSIOper()->ALETimeAdvance()->shiftRight(*firstFluidDisp);
      }

    M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(), M_importerFluid->mapType() ));
    M_velAndPressure->subset(*pressure, pressure->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());
    *M_velAndPressure += *vel;

    M_fluidDisp.reset     ( new vector_Type( *fluidDisp, M_importerFluid->mapType() ));

    M_solidDisp.reset     ( new vector_Type( *solidDisp, M_importerSolid->mapType() ));

    M_data->dataFluid()->dataTime()->updateTime(),M_data->dataSolid()->dataTime()->updateTime();
}


void Problem::initializeStokes( GetPot const& data_file)
{

  using namespace LifeV;

  typedef FSIOperator::mesh_Type        mesh_Type;


  //Creating the pointer to the filter
  filterPtr_Type importer;
  std::string const importerType =  data_file( "importer/type", "hdf5");
  std::string const filename    = data_file( "importer/fluid/filename", "fluid");

  std::string const loadInitSol      = data_file( "importer/initSol", "00000");
  std::string const loadInitSolFD    = data_file("importer/initSolFD","-1");
  std::string iterationString;

  std::cout << "The filename is    : " << filename << std::endl;
  std::cout << "The importerType is: " << importerType << std::endl;

#ifdef HAVE_HDF5
  if (importerType.compare("hdf5") == 0)
    {
      importer.reset( new  hdf5Filter_Type( data_file, filename) );
    }
  else
#endif
    {
      if (importerType.compare("none") == 0)
	{
	  importer.reset( new ExporterEmpty<RegionMesh<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), "fluid", M_fsi->FSIOper()->uFESpace().map().comm().MyPID()));
	}
      else
	{
	  importer.reset( new  ensightFilter_Type( data_file, filename) );
	}
    }

  importer->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());


  std::vector<vectorPtr_Type> fluidStencil;;

  vectorPtr_Type temporarySol     (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique));

  UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

    iterationString = loadInitSol;
    for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->fluidTimeAdvance()->size(); ++iterInit )
    {

        /*!
          definition of the vector to fill with the initialization.
        */
        temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), Unique));

        vectorPtr_Type vel           (new vector_Type(M_fsi->FSIOper()->uFESpace().map(), M_importerFluid->mapType()));
        vectorPtr_Type pressure      (new vector_Type(M_fsi->FSIOper()->pFESpace().map(), M_importerFluid->mapType()));

	*vel *= 0.0;
	*pressure *= 0.0;

        LifeV::ExporterData<mesh_Type> initSolFluidVel   (LifeV::ExporterData<mesh_Type>::VectorField, std::string("f-velocity."+iterationString), M_fsi->FSIOper()->uFESpacePtr(), vel, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
        LifeV::ExporterData<mesh_Type> initSolFluidPress (LifeV::ExporterData<mesh_Type>::ScalarField, std::string("f-pressure."+iterationString), M_fsi->FSIOper()->pFESpacePtr(), pressure, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

        /*!load of the solutions*/
        importer->readVariable(initSolFluidVel);
        importer->readVariable(initSolFluidPress);

        int iterations = std::atoi(iterationString.c_str());
        iterations--;

        std::ostringstream iter;
        iter.fill( '0' );
        iter << std::setw(5) << ( iterations );
        iterationString=iter.str();

        *temporarySol=*vel+*pressure;
        fluidStencil.push_back(temporarySol);
    }

    //    M_fsi->initialize(fluidStencil);

}

void
Problem::errors( GetPot const& dataFile)
{
  using namespace LifeV;
            vectorPtr_Type velAndPressure;
            vectorPtr_Type velAndPressure2;
            // only for export -->
            vectorPtr_Type exactPressPtr;
            vector_Type exactPress(M_fsi->FSIOper()->pFESpace().map(), LifeV::Repeated);
            vectorPtr_Type exactVelPtr;
            vector_Type exactVel(M_fsi->FSIOper()->uFESpace().map(), LifeV::Repeated);
            // <--

            std::string importerType =  dataFile( "importer/type", "ensight");
	    std::string name =  dataFile( "importer/name", "ethiersteinman");

            std::string importerType2 =  dataFile( "importer/type2", "ensight");
	    std::string name2 =  dataFile( "importer/name2", "ethiersteinman");

	    filterPtr_Type importer;
	    filterPtr_Type importer2;

	     boost::shared_ptr<fespace_Type> pFESpace(M_fsi->FSIOper()->pFESpacePtr());
	     boost::shared_ptr<fespace_Type> uFESpace(M_fsi->FSIOper()->uFESpacePtr());
	    UInt pressureOffset =UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());

#ifdef HAVE_HDF5
  if (importerType.compare("hdf5") == 0)
    {
      importer.reset( new  hdf5Filter_Type( dataFile, name) );
    }
  else
#endif
            {
                if (importerType.compare("none") == 0)
                {
                    importer.reset( new ExporterEmpty<mesh_Type > ( dataFile, M_fsi->FSIOper()->uFESpace().mesh(), name.c_str(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
                }
                else
                {
                    importer.reset( new ExporterEnsight<mesh_Type > ( dataFile, M_fsi->FSIOper()->uFESpace().mesh(), name.c_str(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
                }
            }

  velAndPressure.reset(  new vector_Type( M_fsi->FSIOper()->fluid().getMap(), importer->mapType()  ));

            importer->addVariable( ExporterData<mesh_Type>::VectorField, "f-velocity", uFESpace,
                                  velAndPressure, UInt(0));
            importer->addVariable( ExporterData<mesh_Type>::ScalarField, "f-pressure", pFESpace,
                                  velAndPressure, pressureOffset );


#ifdef HAVE_HDF5
  if (importerType2.compare("hdf5") == 0)
    {
      importer2.reset( new  hdf5Filter_Type( dataFile, name2) );
    }
  else
#endif
            {
                if (importerType2.compare("none") == 0)
                {
                    importer2.reset( new ExporterEmpty<mesh_Type > ( dataFile, M_fsi->FSIOper()->uFESpace().mesh(), name2.c_str(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
                }
                else
                {
		  importer2.reset( new ExporterEnsight<mesh_Type > ( dataFile, M_fsi->FSIOper()->uFESpace().mesh(), name2.c_str(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
                }
            }

  velAndPressure2.reset(  new vector_Type( M_fsi->FSIOper()->fluid().getMap(), importer2->mapType()  ));

            importer2->addVariable( ExporterData<mesh_Type>::VectorField, "f-velocity", uFESpace,
                                  velAndPressure2, UInt(0));
            importer2->addVariable( ExporterData<mesh_Type>::ScalarField, "f-pressure", pFESpace,
                                  velAndPressure2, pressureOffset );

  	    std::string tstart1 =  dataFile( "importer/tstart1", "00000");
  	    std::string tstart2 =  dataFile( "importer/tstart2", "00000");


  	    LifeV::ExporterData<mesh_Type> initSolVel(LifeV::ExporterData<mesh_Type>::VectorField,
  						      std::string("f-velocity." + tstart1),
  						      uFESpace,
  						      velAndPressure,
  						      UInt(0),
  						      LifeV::ExporterData<mesh_Type>::UnsteadyRegime,
  						      LifeV::ExporterData<mesh_Type>::Node );

        importer->readVariable(initSolVel  );


        LifeV::ExporterData<mesh_Type> initSolVel2(LifeV::ExporterData<mesh_Type>::VectorField,
  						   std::string("f-velocity." + tstart2),
  						   uFESpace,
  						   velAndPressure2,
  						   UInt(0),
  						   LifeV::ExporterData<mesh_Type>::UnsteadyRegime,
  						   LifeV::ExporterData<mesh_Type>::Node );

        importer2->readVariable(initSolVel2  );


  	    std::cout<< "norm 2 vel1 \n " <<velAndPressure->norm2()<<std::endl;

  	    vector_Type vel  (uFESpace->map(), Repeated);
  	    vector_Type velpressure ( *velAndPressure, Repeated );
  	    vel.subset(velpressure);

    vector_Type vel2  (uFESpace->map(), Repeated);
    vector_Type velpressure2 ( *velAndPressure2, Repeated );
    vel2.subset(velpressure2);

    std::cout<< "norm 2 vel2 \n " <<velAndPressure2->norm2()<<std::endl;

 Real error = M_fsi->FSIOper()->uFESpacePtr()->l2Norm (vel-vel2);
std::cout<<error<<std::endl;
}
