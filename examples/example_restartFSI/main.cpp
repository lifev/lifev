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
#undef HAVE_HDF5
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

#include "ud_functions.hpp"
#include "flowConditions.hpp"
#include "lumpedHeart.hpp"
#include "boundaryConditions.hpp"

class Problem
{
public:

  typedef boost::shared_ptr<LifeV::FSISolver> fsi_solver_ptr;

  typedef LifeV::FSIOperator::data_Type                          data_Type;
  typedef LifeV::FSIOperator::dataPtr_Type                       dataPtr_Type;

  typedef LifeV::FSIOperator::vector_Type        vector_Type;
  typedef LifeV::FSIOperator::vectorPtr_Type     vectorPtr_Type;

  typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh3D<LifeV::LinearTetra> > > filterPtr_Type;

  typedef LifeV::ExporterEnsight<LifeV::FSIOperator::mesh_Type>  ensightFilter_Type;
  typedef boost::shared_ptr<ensightFilter_Type>                 ensightFilterPtr_Type;
#ifdef HAVE_HDF5
  typedef LifeV::ExporterHDF5<LifeV::FSIOperator::mesh_Type>  hdf5Filter_Type;
  typedef boost::shared_ptr<hdf5Filter_Type>                  hdf5FilterPtr_Type;
#endif
  typedef LifeV::FactorySingleton<LifeV::Factory<LifeV::FSIOperator,  std::string> > FSIFactory_Type;
  /*!
    This routine sets up the problem:

    -# create the standard boundary conditions for the fluid and
    structure problems.

    -# initialize and setup the FSIsolver
  */

  Problem( GetPot const& data_file ):
    M_Tstart(0.)
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
    //M_data->dataSolid()->setTimeData( M_data->dataFluid()->dataTime() ); //Same TimeData for fluid & solid
    //M_data->showMe();

    M_fsi = fsi_solver_ptr( new FSISolver( ) );
    MPI_Barrier( MPI_COMM_WORLD );

    M_fsi->setData( M_data );
    M_fsi->FSIOper()->setDataFile( data_file ); //TO BE REMOVED!

    MPI_Barrier( MPI_COMM_WORLD );

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


    Debug( 10000 ) << "Setting up the FESpace and DOF \n";

    MPI_Barrier( MPI_COMM_WORLD );

#ifdef DEBUG
    Debug( 10000 ) << "Setting up the BC \n";
#endif
    M_fsi->setFluidBC( BCh_monolithicFlux( true ) );
    M_fsi->setSolidBC( BCh_monolithicSolid( *M_fsi->FSIOper( ) ) );

    M_fsi->setup(/*data_file*/);

    M_fsi->FSIOper()->fluid().setupPostProc();
    M_fsi->setFluidBC( BCh_monolithicFluid( *M_fsi->FSIOper( ), true ) );
    M_fsi->setHarmonicExtensionBC( BCh_harmonicExtension( *M_fsi->FSIOper( ) ) );

#ifdef DEBUG
    Debug( 10000 ) << "BC set\n";
#endif

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
	    M_exporterFluid.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), fluidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
	    M_exporterSolid.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), solidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
	  }
	else
	  {
	    M_exporterFluid.reset( new  ensightFilter_Type( data_file, fluidName) );
	    M_exporterSolid.reset( new  ensightFilter_Type ( data_file, solidName) );
	  }
      }
    //Creating the vector to export the solution
    M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluid->mapType() ));
    M_fluidDisp.reset     ( new vector_Type( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ));
    M_solidDisp.reset( new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
    M_solidVel.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
    M_solidAcc.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        
    //M_WS.reset           ( new vector_Type(  M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));

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
    M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-velocity",
				  M_fsi->FSIOper()->dFESpacePtr(), M_solidVel, UInt(0) );
    M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-acceleration",
				  M_fsi->FSIOper()->dFESpacePtr(), M_solidAcc, UInt(0) );

    //M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-ws",
    //M_fsi->FSIOper()->dFESpacePtr(), M_WS, UInt(0) );


    /*Initialization*/

    // load using ensight/hdf5
    std::string restartType(data_file("importer/restartType","newSimulation"));
    std::cout << "The load state is: "<< restartType << std::endl;

    if (!restartType.compare("restartFSI"))
      restartFSI(restartType,data_file);
    else
      M_fsi->initialize();

    dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->mergeBCHandlers();

    FC0.initParameters( *M_fsi->FSIOper(), 3);

    M_data->dataFluid()->dataTime()->setInitialTime( M_data->dataFluid()->dataTime()->initialTime()  ); //+ M_data->dataFluid()->dataTime()->timeStep()
    M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
    M_data->dataSolid()->dataTime()->setInitialTime( M_data->dataFluid()->dataTime()->initialTime() ); // + M_data->dataFluid()->dataTime()->timeStep()
    M_data->dataSolid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
    M_data->dataALE()->setInitialTime( M_data->dataFluid()->dataTime()->initialTime() ); //+ M_data->dataFluid()->dataTime()->timeStep()
    M_data->dataALE()->setTime( M_data->dataFluid()->dataTime()->initialTime() );

    M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDisp);//    displacement(), M_offset);
    M_fsi->FSIOper()->exportSolidVelocity(*M_solidVel);//    displacement(), M_offset);
    M_fsi->FSIOper()->exportSolidAcceleration(*M_solidAcc);//    displacement(), M_offset);
    M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressure);
    
    *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();

    M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->initialTime() );
    M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->initialTime() );
  }

  fsi_solver_ptr fsiSolver() { return M_fsi; }

  dataPtr_Type fsiData() { return M_data; }

  /*!
    This routine runs the temporal loop
  */
  void
  run()
  {
    using namespace LifeV;
    boost::timer _overall_timer;

    int iter = 0;
    LifeV::UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

    dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->enableStressComputation(1);

    /*
      #ifdef HAVE_HDF5
      if (M_exporterFluid->mapType() == LifeV::Unique)
      {
      M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->initialTime() );//ugly way to avoid that hdf5 starts with a deformed mesh
      M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->initialTime() );//ugly way to avoid that hdf5 starts with a deformed mesh
      }
      #endif
    */
    bool valveIsOpen=true;
    Real dt = M_data->dataFluid()->dataTime()->timeStep();
    Real T = M_data->dataFluid()->dataTime()->endTime();
    Real t0 = M_data->dataFluid()->dataTime()->initialTime();

    for ( Real time=t0+dt; time<=T; time+=dt)
      {

	boost::timer _timer;

	M_fsi->iterate();

	//*M_WS= *(dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->/*WS());//*/computeStress());

	M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDisp);//    displacement(), M_offset);
	M_fsi->FSIOper()->exportSolidVelocity(*M_solidVel);//    displacement(), M_offset);
	M_fsi->FSIOper()->exportSolidAcceleration(*M_solidAcc);//    displacement(), M_offset);
	M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressure);

	*M_fluidDisp      = M_fsi->FSIOper()->meshDisp();
	M_exporterSolid->postProcess( time );
	M_exporterFluid->postProcess( time );

	std::cout << "[fsi_run] Iteration " << ++iter << " was done in : "
		  << _timer.elapsed() << "\n";

	std::cout << "solution norm " << iter << " : "
		  << M_fsi->displacement().norm2() << "\n";

	// if ( M_data->dataFluid()->dataTime()->time() == 0.004 )
	//   {
	// std::string sol="solutionGlobal";
	// M_fsi->FSIOper()->solution().spy(sol);
	// std::string vAndP="velAndP";
	// M_velAndPressure->spy(vAndP);
	//   }
      }
    /*
      if (M_data->method().compare("monolithicGI"))
      {
      M_fsi->FSIOper()->iterateMesh(M_fsi->displacement());

      M_solidDisp->subset(M_fsi->displacement(), offset);
      //M_solidVel->subset(M_fsi->FSIOper()->solid().getVelocity(), offset);
      //             *M_solidDisp *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataFluid()->dataTime()->getTimeStep());
      //             *M_solidVel  *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataFluid()->dataTime()->getTimeStep());

      *M_velAndPressure = M_fsi->displacement();
      M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );
      *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
      M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );
      }
    */

    FC0.renewParameters( *M_fsi, 3 );
    std::cout << "Total computation time = "
	      << _overall_timer.elapsed() << "s" << "\n";

  }

private:

  void restartFSI(std::string& restartType,  GetPot const& data_file);

  fsi_solver_ptr M_fsi;
  dataPtr_Type   M_data;

  filterPtr_Type M_exporterSolid;
  filterPtr_Type M_exporterFluid;
  filterPtr_Type M_importerSolid;
  filterPtr_Type M_importerFluid;

  filterPtr_Type M_exporterSolidCheck;
  filterPtr_Type M_exporterFluidCheck;

  vectorPtr_Type M_velAndPressure;
  vectorPtr_Type M_fluidDisp;
  vectorPtr_Type M_solidDisp;
  vectorPtr_Type M_solidVel;
  vectorPtr_Type M_solidAcc;

  vectorPtr_Type M_velAndPressureCheck;
  vectorPtr_Type M_fluidDispCheck;
  vectorPtr_Type M_solidDispCheck;
  vectorPtr_Type M_solidVelCheck;
  vectorPtr_Type M_solidAccCheck;


  LifeV::FlowConditions FC0;
  
  LifeV::Real    M_Tstart;
  vectorPtr_Type M_WS;

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

	fsip->run();
      }
    catch ( std::exception const& _ex )
      {
	std::cout << "caught exception :  " << _ex.what() << "\n";
      }

    //@disp = fsip->fsiSolver()->FSIOper()->displacementOnInterface();
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
      //This is not a test in the testsuite!!
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


void Problem::restartFSI(std::string& restartType,  GetPot const& data_file)
{

  using namespace LifeV;

  typedef FSIOperator::mesh_Type        mesh_Type;

  //Creating the pointer to the filter
  std::string const importerType =  data_file( "importer/type", "ensight");
  std::string const fluidName     = data_file( "importer/fluid/filename", "fluid");
  std::string const solidName    =  data_file( "importer/solid/filename", "solid");

  std::string const loadInitSol      = data_file( "importer/initSol", "00000");
  std::string const loadInitSolFD    = data_file("importer/initSolFD","-1");
  std::string iterationString;

  M_Tstart=data_file( "fluid/time_discretization/initialtime", 0.);

  if(!loadInitSolFD.compare("-1"))
    std::cout << "You have to set the iteration of the previous fluid displacement"<< std::endl;

  std::cout << "The file for fluid is    : " << fluidName << std::endl;
  std::cout << "The file for solid is    : " << solidName << std::endl;
  std::cout << "The importerType is      : " << importerType << std::endl;
  std::cout << "The iteration is         : " << loadInitSol << std::endl;
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

  std::vector<vectorPtr_Type> solidStencil;
  std::vector<vectorPtr_Type> fluidStencil;
  std::vector<vectorPtr_Type> ALEStencil;

  vectorPtr_Type fluidSol         (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique));
  vectorPtr_Type initFluid         (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
  vectorPtr_Type HarmonicSol      (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
  vectorPtr_Type structureSol      (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
  vectorPtr_Type temporarySol     (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
  //  vectorPtr_Type un               (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap()));

  iterationString = loadInitSol;
  std::cout << "Fluid size TimeAdvance:" << M_fsi->FSIOper()->fluidTimeAdvance()->size() << std::endl;
  for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->fluidTimeAdvance()->size(); iterInit++ )
    {

      /*!
	definition of the vector to fill with the initialization.
      */
      temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
      initFluid.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));

      vectorPtr_Type vel           (new vector_Type(M_fsi->FSIOper()->uFESpace().map(), M_importerFluid->mapType()));
      vectorPtr_Type pressure      (new vector_Type(M_fsi->FSIOper()->pFESpace().map(), M_importerFluid->mapType()));

      *vel *= 0.0;
      *pressure *= 0.0;

      LifeV::ExporterData<mesh_Type> initSolFluidVel   (LifeV::ExporterData<mesh_Type>::VectorField, std::string("f-velocity."+iterationString), M_fsi->FSIOper()->uFESpacePtr(), vel, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
      LifeV::ExporterData<mesh_Type> initSolFluidPress (LifeV::ExporterData<mesh_Type>::ScalarField, std::string("f-pressure."+iterationString), M_fsi->FSIOper()->pFESpacePtr(), pressure, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

      /*!load of the solutions*/
      M_importerFluid->readVariable(initSolFluidVel);
      M_importerFluid->readVariable(initSolFluidPress);

      initFluid->subset(*vel, vel->map(), UInt(0), UInt(0));
      initFluid->subset(*pressure, pressure->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());

      *temporarySol = *initFluid;

      fluidStencil.push_back(temporarySol);

      //Updating string name
      int iterations = std::atoi(iterationString.c_str());
      iterations--;

      std::ostringstream iter;
      iter.fill( '0' );
      iter << std::setw(5) << ( iterations );
      iterationString=iter.str();

      std::cout << "Iterinit: " << iterInit << "Norm: "<< fluidStencil[iterInit]->norm2() << std::endl;

    }

  std::cout << "solid size TimeAdvance:" << M_fsi->FSIOper()->solidTimeAdvance()->size() << std::endl;
  UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();
  iterationString = loadInitSol;
  for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->solidTimeAdvance()->size(); iterInit++ )
    {
      /*!
	definition of the vector to fill with the initialization.
      */
      temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
      structureSol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
      vectorPtr_Type solidDisp    (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), M_importerSolid->mapType()));
      *solidDisp *= 0.0;

      /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
      LifeV::ExporterData<mesh_Type> initSolSolidDisp  (LifeV::ExporterData<mesh_Type>::VectorField,"s-displacement."+iterationString, M_fsi->FSIOper()->dFESpacePtr(), solidDisp, UInt(0)/*offset*/, LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

      /*!load of the solutions*/
      M_importerSolid->readVariable(initSolSolidDisp);

      structureSol->subset(*solidDisp, solidDisp->map(), (UInt)0, offset);
      *structureSol *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataSolid()->dataTime()->timeStep());

      *temporarySol = *structureSol;

      solidStencil.push_back(temporarySol);
      std::cout << "Norm2 structure"<< temporarySol->norm2() << std::endl;

      //Updating the string name for the next iteration
      UInt iterations = std::atoi(iterationString.c_str());
      iterations--;

      std::ostringstream iter;
      iter.fill( '0' );
      iter << std::setw(5) << ( iterations );
      iterationString=iter.str();

      std::cout << "Iterinit: " << iterInit << "Norm: "<< solidStencil[iterInit]->norm2() << std::endl;

    }

  iterationString = loadInitSol;
  std::cout << "ALE size TimeAdvance:" << M_fsi->FSIOper()->ALETimeAdvance()->size() << std::endl;
  for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->ALETimeAdvance()->size(); iterInit++ )
    {
      temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV:: Unique, Zero));
      /*!
	definition of the vector to fill with the initialization.
      */
      vectorPtr_Type fluidDisp(new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), M_importerFluid->mapType()));
      *fluidDisp *= 0.0;

      LifeV::ExporterData<mesh_Type> initSolFluidDisp  (LifeV::ExporterData<mesh_Type>::VectorField, "f-displacement."+iterationString, M_fsi->FSIOper()->mmFESpacePtr(), fluidDisp, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

      /*!load of the solutions*/
      M_importerFluid->readVariable(initSolFluidDisp);

      HarmonicSol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
      HarmonicSol->subset(*fluidDisp, fluidDisp->map(), (UInt)0, dynamic_cast<LifeV::FSIMonolithicGI*>(M_fsi->FSIOper().get())->mapWithoutMesh().map(Unique)->NumGlobalElements());
      *temporarySol = *HarmonicSol;

      ALEStencil.push_back(temporarySol);
      //Updating the iteration String name

      std::cout << "Norm2 fluid"<< temporarySol->norm2() << std::endl;
      int iterations = std::atoi(iterationString.c_str());
      iterations--;

      std::ostringstream iter;
      iter.fill( '0' );
      iter << std::setw(5) << ( iterations );
      iterationString=iter.str();
      
      std::cout << "Iterinit: " << iterInit << "Norm: "<< ALEStencil[iterInit]->norm2() << std::endl;
    }

  M_fsi->initialize(fluidStencil, solidStencil, ALEStencil);

  
  //Exporting the solution to see if it is correct
  std::string fluidTry = "fluidCheckRestarted";
  std::string solidTry = "solidCheckRestarted";
  std::string const exporterType =  data_file( "exporter/type", "ensight");

#ifdef HAVE_HDF5
  if (exporterType.compare("hdf5") == 0)
    {
      M_exporterFluidCheck.reset( new  hdf5Filter_Type( data_file, fluidTry) );
      M_exporterSolidCheck.reset( new  hdf5Filter_Type ( data_file,solidTry) );
    }
  else
#endif
    {
      if (exporterType.compare("none") == 0)
	{
	  M_exporterFluidCheck.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), fluidTry, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
	  M_exporterSolidCheck.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), solidTry, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
	}
      else
	{
	  M_exporterFluidCheck.reset( new  ensightFilter_Type( data_file, fluidTry) );
	  M_exporterSolidCheck.reset( new  ensightFilter_Type ( data_file, solidTry) );
	}
    }


  //Creating the vector to export the solution
  M_velAndPressureCheck.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluidCheck->mapType() ));
  M_fluidDispCheck.reset     ( new vector_Type( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluidCheck->mapType() ));
  M_solidDispCheck.reset( new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolidCheck->mapType() ));
  M_solidVelCheck.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolidCheck->mapType() ));
  M_solidAccCheck.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolidCheck->mapType() ));
        
  //M_WS.reset           ( new vector_Type(  M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));

  M_exporterFluidCheck->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());
  M_exporterSolidCheck->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID());
  M_exporterFluidCheck->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-velocity",
				     M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressureCheck, UInt(0) );
  M_exporterFluidCheck->addVariable( ExporterData<FSIOperator::mesh_Type>::ScalarField, "f-pressure",
				     M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressureCheck,
				     UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()) );

  M_exporterFluidCheck->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-displacement",
				     M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDispCheck, UInt(0) );


  M_exporterSolidCheck->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-displacement",
				     M_fsi->FSIOper()->dFESpacePtr(), M_solidDispCheck, UInt(0) );
  M_exporterSolidCheck->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-velocity",
				     M_fsi->FSIOper()->dFESpacePtr(), M_solidVelCheck, UInt(0) );
  M_exporterSolidCheck->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-acceleration",
				     M_fsi->FSIOper()->dFESpacePtr(), M_solidAccCheck, UInt(0) );


  //Saving the solution
  M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDispCheck);//    displacement(), M_offset);
  M_fsi->FSIOper()->exportSolidVelocity(*M_solidVelCheck);//    displacement(), M_offset);
  M_fsi->FSIOper()->exportSolidAcceleration(*M_solidAccCheck);//    displacement(), M_offset);
  M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressureCheck);

  *M_fluidDispCheck      = M_fsi->FSIOper()->meshDisp();
  M_exporterSolidCheck->postProcess( M_data->dataFluid()->dataTime()->initialTime() );
  M_exporterFluidCheck->postProcess( M_data->dataFluid()->dataTime()->initialTime() );

  
}


