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
#include "boundaryConditions.hpp"
#include "boundaryConditionsSteady.hpp"

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
      M_Tstart(0.),
      M_solution( )
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
        M_fsi->setFluidBC( BCh_steadyMonolithicFlux( true ) );
        M_fsi->setSolidBC( BCh_steadyMonolithicSolid( *M_fsi->FSIOper( ) ) );

        M_fsi->setup(/*data_file*/);
        M_fsi->setFluidBC( BCh_steadyMonolithicFluid( *M_fsi->FSIOper( ), true ) );
        M_fsi->setHarmonicExtensionBC( BCh_steadyHarmonicExtension( *M_fsi->FSIOper( ) ) );


#ifdef DEBUG
        Debug( 10000 ) << "BC set\n";
#endif

        std::string const exporterType =  data_file( "exporter/type", "ensight" );
        std::string fluidName    =  data_file( "exporter/fluid/filename", "fluid" );
        std::string solidName    =  data_file( "exporter/solid/filename", "solid" );

	fluidName += "Stokes";
	solidName += "Stokes";

#ifdef HAVE_HDF5
        if (exporterType.compare("hdf5") == 0)
        {
            M_exporterFluidStokes.reset( new  hdf5Filter_Type( data_file, fluidName) );
            M_exporterSolidStokes.reset( new  hdf5Filter_Type ( data_file,solidName));
        }
        else
#endif
        {
            if (exporterType.compare("none") == 0)
            {
                M_exporterFluidStokes.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), fluidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
                M_exporterSolidStokes.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), solidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID()) );
            }
            else
            {
                M_exporterFluidStokes.reset( new  ensightFilter_Type( data_file, fluidName) );
                M_exporterSolidStokes.reset( new  ensightFilter_Type ( data_file, solidName) );
            }
        }
        M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluidStokes->mapType() ));
        M_fluidDisp.reset     ( new vector_Type( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluidStokes->mapType() ));

        M_exporterFluidStokes->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());
        M_exporterSolidStokes->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID());
        M_exporterFluidStokes->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-velocity",
                                      M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressure, UInt(0) );
        M_exporterFluidStokes->addVariable( ExporterData<FSIOperator::mesh_Type>::ScalarField, "f-pressure",
                                      M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressure,
                                      UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()) );

        M_exporterFluidStokes->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-displacement",
                                      M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDisp, UInt(0) );



        M_solidDisp.reset( new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolidStokes->mapType() ));
        M_solidVel.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolidStokes->mapType() ));
        M_solidAcc.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolidStokes->mapType() ));
        //        M_solidVel.reset ( new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        //M_WS.reset           ( new vector_Type(  M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));

        M_exporterSolidStokes->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-displacement",
                                      M_fsi->FSIOper()->dFESpacePtr(), M_solidDisp, UInt(0) );
        M_exporterSolidStokes->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-velocity",
                                      M_fsi->FSIOper()->dFESpacePtr(), M_solidVel, UInt(0) );
        M_exporterSolidStokes->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-acceleration",
                                      M_fsi->FSIOper()->dFESpacePtr(), M_solidAcc, UInt(0) );

        //M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-ws",
	//M_fsi->FSIOper()->dFESpacePtr(), M_WS, UInt(0) );


        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->mergeBCHandlers();
	/*Initialization*/
	M_fsi->FSIOper()->displayer().leaderPrint( " The First Simulation is Steady! \n" );
	
	M_fsi->initialize();

        FC0.initParameters( *M_fsi->FSIOper(), 3);


        M_data->dataFluid()->dataTime()->setInitialTime( M_Tstart  ); //+ M_data->dataFluid()->dataTime()->timeStep()
        M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setInitialTime( M_Tstart); // + M_data->dataFluid()->dataTime()->timeStep() 
        M_data->dataSolid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataALE()->setInitialTime( M_Tstart ); //+ M_data->dataFluid()->dataTime()->timeStep() 
        M_data->dataALE()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
    }

  void
  runStokes()
    {
        using namespace LifeV;
        boost::timer _overall_timer;

	M_solution.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap() ,LifeV::Unique));

        int iter = 1;
        LifeV::UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->enableStressComputation(1);

#ifdef HAVE_HDF5
        if (M_exporterFluidStokes->mapType() == LifeV::Unique)
        {
            M_exporterFluidStokes->postProcess( M_Tstart );//ugly way to avoid that hdf5 starts with a deformed mesh
            M_exporterSolidStokes->postProcess( M_Tstart );//ugly way to avoid that hdf5 starts with a deformed mesh
        }
#endif

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(),M_data->dataSolid()->dataTime()->updateTime(), ++iter)
        {
             FC0.renewParameters( *M_fsi, 3 );

            boost::timer _timer;

            M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDisp);//    displacement(), M_offset);
	    M_fsi->FSIOper()->exportSolidVelocity(*M_solidVel);//    displacement(), M_offset);
	    M_fsi->FSIOper()->exportSolidAcceleration(*M_solidAcc);//    displacement(), M_offset);
	    //M_fsi->FSIOper()->exportSolidVelocity(*M_solidVel);//    displacement(), M_offset);

            M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressure);

            M_exporterSolidStokes->postProcess( M_data->dataFluid()->dataTime()->time() );

            M_fsi->iterate();

            //*M_WS= *(dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->/*WS());//*/computeStress());


            *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();

            M_exporterFluidStokes->postProcess( M_data->dataFluid()->dataTime()->time() );


            std::cout << "[fsi_run] Iteration " << iter << " was done in : "
                      << _timer.elapsed() << "\n";

        }

	FC0.renewParameters( *M_fsi, 3 );
        std::cout << "Total computation time = "
                  << _overall_timer.elapsed() << "s" << "\n";
	
	//resetting the pointer to the right type vector
	M_solution = M_fsi->FSIOper()->solutionPtr();
    }


  Problem( GetPot const& data_file, boost::shared_ptr<LifeV::VectorEpetra> soluzione ):
    M_Tstart(0.),
    M_solution( )
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

        M_fsi->setFluidBC( BCh_monolithicFluid( *M_fsi->FSIOper( ), true ) );
        M_fsi->setHarmonicExtensionBC( BCh_harmonicExtension( *M_fsi->FSIOper( ) ) );

#ifdef DEBUG
        Debug( 10000 ) << "BC set\n";
#endif

	std::cout << "Finito con le BC" << std::endl;
	int n;
	std::cin >> n;
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
	//Saving the solution
	std::cout << "saving the solution "<< std::endl;
	M_solution.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap() ,LifeV::Unique));
	M_solution = soluzione;
	std::cin >> n;
        M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluid->mapType() ));
        M_fluidDisp.reset     ( new vector_Type( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ));

        M_exporterFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());
        M_exporterSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID());
        M_exporterFluid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-velocity",
                                      M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressure, UInt(0) );
        M_exporterFluid->addVariable( ExporterData<FSIOperator::mesh_Type>::ScalarField, "f-pressure",
                                      M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressure,
                                      UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()) );

        M_exporterFluid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-displacement",
                                      M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDisp, UInt(0) );

	std::cout << "sono qui!!!"<< std::endl;
	std::cin >> n;

        M_solidDisp.reset( new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        M_solidVel.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        M_solidAcc.reset(  new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        //        M_solidVel.reset ( new vector_Type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        //M_WS.reset           ( new vector_Type(  M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));

        M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-displacement",
                                      M_fsi->FSIOper()->dFESpacePtr(), M_solidDisp, UInt(0) );
        M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-velocity",
                                      M_fsi->FSIOper()->dFESpacePtr(), M_solidVel, UInt(0) );
        M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-acceleration",
                                      M_fsi->FSIOper()->dFESpacePtr(), M_solidAcc, UInt(0) );

        //M_exporterSolid->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-ws",
	//M_fsi->FSIOper()->dFESpacePtr(), M_WS, UInt(0) );


        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->mergeBCHandlers();

	/*Initialization*/
	M_fsi->FSIOper()->displayer().leaderPrint( " The Simulation is Unsteady! \n" );
	//Getting the offset
	UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();
	UInt offsetLambda = offset + M_fsi->FSIOper()->dFESpace().dof().numTotalDof();

	//Initialization of the stencils
	std::vector<vectorPtr_Type> fluidStencil;
	std::vector<vectorPtr_Type> solidStencil;
	std::vector<vectorPtr_Type> ALEStencil;

	M_fsi->FSIOper()->displayer().leaderPrint( " \n " );

	//These vectors are created to put zeros in the right places in the global monolithic vector
	vector_Type vel      (M_fsi->FSIOper()->uFESpace().map(), M_exporterFluid->mapType());
	vector_Type pressure (M_fsi->FSIOper()->pFESpace().map(), M_exporterFluid->mapType());
	vector_Type solidDisp(M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType());
	vector_Type fluidDisp(M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType());
	vel *= 0.0;
	pressure *= 0.0;
	solidDisp *= 0.0;
	fluidDisp *= 0.0;
	//This vector is to put zero in the Lagrange multiplyer zone
	vector_Type zeros(*M_fsi->FSIOper()->solidInterfaceMap(), LifeV::Unique);
	zeros *= 0.0;
	//Create the proper vector for the fluid Stencil
	vectorPtr_Type copySol(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
	copySol = M_solution;
	std::cout << "The norm of the copied solution is: "<< copySol->norm2() << std::endl;
	std::cout << "sono qui1!!!"<< std::endl;
	std::cin >> n; 

	copySol->subset(solidDisp, solidDisp.map(),UInt(0),offset);
	std::cout << "sono qui2!!!"<< std::endl;
	std::cin >> n;

	copySol->subset(fluidDisp, fluidDisp.map(),UInt(0),dynamic_cast<LifeV::FSIMonolithicGI*>(M_fsi->FSIOper().get())->mapWithoutMesh().map(Unique)->NumGlobalElements());
	std::cout << "sono qui2!!!"<< std::endl;
	std::cin >> n;
	copySol->subset(zeros, zeros.map(),UInt(0),offsetLambda);
	std::cout << "sono qui2!!!"<< std::endl;
	std::cin >> n;

	
	std::string fluidPart="fluidPart";
	copySol->spy(fluidPart);

	for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->fluidTimeAdvance()->size(); ++iterInit )
	  {
	    fluidStencil.push_back(copySol);
	  }

	//This part is to create the right vector for the solid stencil
	copySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
	copySol = M_solution;	
	    
	copySol->subset(vel, vel.map(),UInt(0),UInt(0));
	copySol->subset(pressure, pressure.map(),UInt(0),UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()));
	copySol->subset(fluidDisp, fluidDisp.map(),UInt(0),dynamic_cast<LifeV::FSIMonolithicGI*>(M_fsi->FSIOper().get())->mapWithoutMesh().map(Unique)->NumGlobalElements());
	copySol->subset(zeros, zeros.map(),UInt(0),offsetLambda);

	 std::string solidPart="solidPart";
	 copySol->spy(solidPart);

	for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->solidTimeAdvance()->size(); ++iterInit )
	  {
	    solidStencil.push_back(copySol);
	  }
     
	//Create the proper vectors for ALE stencil
	copySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
	copySol = M_solution;	

	copySol->subset(vel, vel.map(),UInt(0),UInt(0));
	copySol->subset(pressure, pressure.map(),UInt(0),UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()));
	copySol->subset(solidDisp, solidDisp.map(),UInt(0),offset);
	copySol->subset(zeros, zeros.map(),UInt(0),offsetLambda);


	std::string ALEPart="ALEPart";
	copySol->spy(ALEPart);

	for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->ALETimeAdvance()->size(); ++iterInit )
	  {
	    ALEStencil.push_back(copySol);	    
	  }


	M_fsi->initialize(fluidStencil,solidStencil,ALEStencil);
	//Setting the solution
	M_fsi->FSIOper()->setSolution(*M_solution);

        FC0.initParameters( *M_fsi->FSIOper(), 3);

        //LH.initParameters( *M_fsi->FSIOper(), "dataHM");
        M_data->dataFluid()->dataTime()->setInitialTime( M_Tstart  ); //+ M_data->dataFluid()->dataTime()->timeStep()
        M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setInitialTime( M_Tstart); // + M_data->dataFluid()->dataTime()->timeStep() 
        M_data->dataSolid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataALE()->setInitialTime( M_Tstart ); //+ M_data->dataFluid()->dataTime()->timeStep() 
        M_data->dataALE()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
    }


    /*!
      This routine runs the temporal loop
    */
    void
    run()
    {
        boost::timer _overall_timer;

        int iter = 1;
        LifeV::UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->enableStressComputation(1);

#ifdef HAVE_HDF5
        if (M_exporterFluid->mapType() == LifeV::Unique)
        {
            M_exporterFluid->postProcess( M_Tstart );//ugly way to avoid that hdf5 starts with a deformed mesh
            M_exporterSolid->postProcess( M_Tstart );//ugly way to avoid that hdf5 starts with a deformed mesh
        }
#endif

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(),M_data->dataSolid()->dataTime()->updateTime(), ++iter)
        {
             FC0.renewParameters( *M_fsi, 3 );

            boost::timer _timer;

            M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDisp);//    displacement(), M_offset);
	    M_fsi->FSIOper()->exportSolidVelocity(*M_solidVel);//    displacement(), M_offset);
	    M_fsi->FSIOper()->exportSolidAcceleration(*M_solidAcc);//    displacement(), M_offset);
	    //M_fsi->FSIOper()->exportSolidVelocity(*M_solidVel);//    displacement(), M_offset);

            M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressure);

            M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );

            M_fsi->iterate();

            //*M_WS= *(dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->/*WS());//*/computeStress());


            *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();

            M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );

            std::cout << "solution norm " << iter << " : "
                      << M_fsi->displacement().norm2() << "\n";

	    
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
        std::cout << "Total computation time steadyFSI= "
                  << _overall_timer.elapsed() << "s" << "\n";

    }

    fsi_solver_ptr fsiSolver() { return M_fsi; }

    dataPtr_Type fsiData() { return M_data; }

    vectorPtr_Type solutionStokes() {return M_solution;}


private:

    fsi_solver_ptr M_fsi;
    dataPtr_Type   M_data;

    filterPtr_Type M_exporterSolid;
    filterPtr_Type M_exporterFluid;

    filterPtr_Type M_exporterSolidStokes;
    filterPtr_Type M_exporterFluidStokes;

    filterPtr_Type M_importerSolid;
    filterPtr_Type M_importerFluid;
    vectorPtr_Type M_velAndPressure;
    vectorPtr_Type M_fluidDisp;
    vectorPtr_Type M_solidDisp;
    vectorPtr_Type M_solidVel;
    vectorPtr_Type M_solidAcc;

    vectorPtr_Type M_solution;
    //    vectorPtr_Type M_solidVel;
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
  FSIChecker( GetPot const& _data_file, GetPot const& data_steady ):
    data_file( _data_file ),
    dataSteady( data_steady )
    {}

    void operator()()
    {
        boost::shared_ptr<Problem>             fsip;
	boost::shared_ptr<Problem>             fsipSteady;
	boost::shared_ptr<LifeV::VectorEpetra> soluzione;

        try
        {
	  fsipSteady = boost::shared_ptr<Problem>( new Problem( dataSteady ) );
	  fsipSteady->runStokes();
	 
	  std::cout << "Finito la parte sopra "<< std::endl;
	  int n;
	  std::cin >> n;
	  soluzione.reset(new LifeV::VectorEpetra(*(fsipSteady->fsiSolver()->FSIOper()->couplingVariableMap()), LifeV::Unique));
	  soluzione = fsipSteady->solutionStokes();

	  fsip = boost::shared_ptr<Problem>( new Problem( data_file,soluzione ) );
	  fsip->run();
        }
        catch ( std::exception const& _ex )
        {
            std::cout << "caught exception :  " << _ex.what() << "\n";
        }

        //@disp = fsip->fsiSolver()->FSIOper()->displacementOnInterface();
    }

    GetPot                data_file;
    GetPot                dataSteady;
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

  using namespace LifeV;

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#else
    std::cout << "% using serial Version" << std::endl;
#endif

    GetPot command_line(argc,argv);
    const bool check = command_line.search(2, "-c", "--check");

    if (check)
    {
      //This is not a test in the testsuite!!!

#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return 0;
    }
    else
    {
        const std::string data_file_name = command_line.follow("data", 2, "-f","--file");
	std::cout << "The first data file is: "<< data_file_name << std::endl;
	const std::string data_steady = "dataSteady";
	std::cout << "The second data file is: "<< data_steady << std::endl;
        GetPot data_file(data_file_name);
        GetPot data_Steady(data_steady);
        FSIChecker _sp_check( data_file, data_Steady);       
        _sp_check();
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return 0;

}

