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

        // load using ensight/hdf5
        std::string initializationType(data_file("importer/initType","newSimulation"));
	std::cout << "The load state is: "<< initializationType << std::endl;
	
	if (!initializationType.compare("stokes"))
	  {

	    vectorPtr_Type un (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
	    *un *= 0.0;

	    vectorPtr_Type initVec (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
	    *initVec *= 0.0;
	    //This method takes the stokesSolution, which will be updated and the FSISolver
	    initializeStokes(un, data_file, M_fsi, M_data);

	    //loading the stencil for the fluid (Pay Attention that this works with BDF 1)
	    std::vector<vectorPtr_Type> fluidStencil;
	    std::vector<vectorPtr_Type> solidStencil;
	    std::vector<vectorPtr_Type> ALEStencil;
	    for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->fluidTimeAdvance()->size(); iterInit++ )
	      {
		fluidStencil.push_back(un);
	      }
	    for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->solidTimeAdvance()->size(); iterInit++ )
	      {
		solidStencil.push_back(initVec);
	      }
	    for(UInt iterInit=0; iterInit<M_fsi->FSIOper()->ALETimeAdvance()->size(); iterInit++ )
	      {
		ALEStencil.push_back(initVec);
	      }

	    M_fsi->initialize(fluidStencil,solidStencil,ALEStencil );
	    //loading the initialSolution for monolithic
	    M_fsi->FSIOper()->setSolution(*un);
	    /*std::cout << "back in the main" << std::endl;
	    int n;
	    std::cin >> n;*/
	    std::cout << "norm of un" << un->norm2() << std::endl;
	    std::cout << "norm of gotten solution" << M_fsi->FSIOper()->solutionPtr()->norm2() << std::endl;
	  }
	else
	  M_fsi->initialize();

        FC0.initParameters( *M_fsi->FSIOper(), 3);

        //LH.initParameters( *M_fsi->FSIOper(), "dataHM");
        M_data->dataFluid()->dataTime()->setInitialTime( M_Tstart  ); //+ M_data->dataFluid()->dataTime()->timeStep()
        M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setInitialTime( M_Tstart); // + M_data->dataFluid()->dataTime()->timeStep() 
        M_data->dataSolid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataALE()->setInitialTime( M_Tstart ); //+ M_data->dataFluid()->dataTime()->timeStep() 
        M_data->dataALE()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
    }

    fsi_solver_ptr fsiSolver() { return M_fsi; }

    dataPtr_Type fsiData() { return M_data; }

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

        bool valveIsOpen=true;

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(),M_data->dataSolid()->dataTime()->updateTime(), ++iter)
        {
             FC0.renewParameters( *M_fsi, 3 );

            //                 FC0.renewParameters( *M_fsi, 6 );
            //                 FC1.renewParameters( *M_fsi, 3, 4 );
            //                 FC2.renewParameters( *M_fsi, 3, 5 );
            //                 FC3.renewParameters( *M_fsi, 3, 6 );
            //                 FC4.renewParameters( *M_fsi, 3, 7 );

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


            //M_fsi->FSIOper()->displayer().leaderPrint( "average inlet pressure  = ", M_fsi->FSIOper()->fluid().pressure(2, *M_velAndPressure));
            //M_fsi->FSIOper()->displayer().leaderPrint( "average outlet pressure = ", M_fsi->FSIOper()->fluid().pressure(3, *M_velAndPressure));
            //M_fsi->FSIOper()->displayer().leaderPrint( "inlet flux              = ", M_fsi->FSIOper()->fluid().flux(2, *M_velAndPressure));
            //M_fsi->FSIOper()->displayer().leaderPrint( "outlet flux             = ", M_fsi->FSIOper()->fluid().flux(3, *M_velAndPressure));

            std::cout << "[fsi_run] Iteration " << iter << " was done in : "
                      << _timer.elapsed() << "\n";

            std::cout << "solution norm " << iter << " : "
                      << M_fsi->displacement().norm2() << "\n";

            std::cout << "solution norm " << iter << " : "
                      << M_solidDisp->norm2() << "\n";


            std::cout << "Vel norm " << iter << " : "
                      << M_solidVel->norm2() << "\n";


            std::cout << "Acc norm " << iter << " : "
                      << M_solidAcc->norm2() << "\n";


            std::cout << "fluidDisp norm " << iter << " : "
                      << M_fluidDisp->norm2() << "\n";

            std::cout << "VelAndPressure norm " << iter << " : "
                      << M_velAndPressure->norm2() << "\n";

            std::cout << "VelAndPressure norm " << iter << " : "
                      << M_velAndPressure->size() << "\n";
	    

	    if ( M_data->dataFluid()->dataTime()->time() == 0.004 )
	      {
		std::string sol="solutionGlobal";
		M_fsi->FSIOper()->solution().spy(sol);
		std::string vAndP="velAndP";
		M_velAndPressure->spy(vAndP);
	      }
            ///////// CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            try
            {
                if (!M_data->method().compare("monolithicGI"))
                    checkCEResult(M_data->dataFluid()->dataTime()->time());
                else
                    checkGCEResult(M_data->dataFluid()->dataTime()->time());
            }
            catch (Problem::RESULT_CHANGED_EXCEPTION) {std::cout<<"res. changed"<<std::endl;}
            ///////// END OF CHECK
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

  void initializeStokes(vectorPtr_Type& un, GetPot const& data_file, fsi_solver_ptr  fsiSolver, dataPtr_Type fsiData);

    void checkCEResult(const LifeV::Real& time);
    void checkGCEResult(const LifeV::Real& time);

    fsi_solver_ptr M_fsi;
    dataPtr_Type   M_data;

    filterPtr_Type M_exporterSolid;
    filterPtr_Type M_exporterFluid;
    filterPtr_Type M_importerSolid;
    filterPtr_Type M_importerFluid;
    vectorPtr_Type M_velAndPressure;
    vectorPtr_Type M_fluidDisp;
    vectorPtr_Type M_solidDisp;
    vectorPtr_Type M_solidVel;
    vectorPtr_Type M_solidAcc;
    //    vectorPtr_Type M_solidVel;
    LifeV::FlowConditions FC0;
    //LifeV::LumpedHeart LH;
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



void Problem::initializeStokes(vectorPtr_Type& un,  GetPot const& data_file, fsi_solver_ptr fsiSolver, dataPtr_Type fsiData)
{

  //This is a really ugly way to implement this.
  using namespace LifeV;

  typedef RegionMesh3D<LinearTetra>                        Mesh;

  typedef OseenSolver< RegionMesh3D<LinearTetra> >         fluid_Type;
  typedef fluid_Type::vector_Type                          vector_Type;
  typedef boost::shared_ptr<vector_Type>                   vectorPtr_Type;
  typedef FESpace< Mesh, MapEpetra >                       feSpace_Type;
  typedef boost::shared_ptr<feSpace_Type>                  feSpacePtr_Type;
  typedef boost::shared_ptr<Epetra_Comm>                   commPtr_Type; 
  
  commPtr_Type comunicator=fsiSolver->FSIOper()->fluid().comm();

  fsiSolver->FSIOper()->displayer().leaderPrint( " The FSI simulation is going to be initialized by the Stokes solution of the fluid domain \n" );

  int numLM = 0; //this is the number of fluxes
  
  MeshData meshData;
  meshData.setup(data_file, "fluid/space_discretization");

  boost::shared_ptr<Mesh> fullMeshPtr(new Mesh);
  readMesh(*fullMeshPtr, meshData);

  MeshPartitioner< Mesh >   meshPart(fullMeshPtr, comunicator);

  std::string uOrder =  data_file( "fluid/space_discretization/vel_order", "P1");
  std::string pOrder =  data_file( "fluid/space_discretization/press_order", "P1");
  //Building the functional spaces
  feSpacePtr_Type uFESpacePtr( new feSpace_Type(meshPart,uOrder,3,comunicator) );
  feSpacePtr_Type pFESpacePtr( new feSpace_Type(meshPart,pOrder,1,comunicator) );

  fsiSolver->FSIOper()->displayer().leaderPrint( " Calling the fluid Solver constructor \n" );
  fluid_Type stokes(fsiData->dataFluid(),*uFESpacePtr,*pFESpacePtr,comunicator, numLM);

  //Boundary conditions for the fluid problem
  BCHandler bcH;
  BCFunctionBase uZero( fZero );
  BCFunctionBase prInlet( pressure );

  bcH.addBC( "Inlet",    2,    Natural,  Full,  prInlet, 3 );
  bcH.addBC( "Ringin",   20,   Essential,Full,  uZero, 3 );
  bcH.addBC( "Ringout",  30,   Essential,Full,  uZero, 3 );
  bcH.addBC( "Outlet",   3,    Natural,  Full,  uZero, 3 );
  bcH.addBC( "Wall",     1,    Essential,Full,  uZero, 3 );

  UInt totalVelDof   = uFESpacePtr->map().map(Unique)->NumGlobalElements();
  UInt totalPressDof = pFESpacePtr->map().map(Unique)->NumGlobalElements();

  bcH.setOffset("Inlet", totalVelDof + totalPressDof);

  MapEpetra fullMap(stokes.getMap());

  stokes.setUp(data_file);
  stokes.buildSystem();

  vector_Type beta( fullMap ); 
  vector_Type rhs ( fullMap );

  MPI_Barrier( MPI_COMM_WORLD );

  beta *= 0.;
  rhs  *= 0.;

  LifeChrono chrono;
  chrono.start();

  stokes.updateSystem(0, beta, rhs );
  stokes.iterate( bcH );

  MPI_Barrier(MPI_COMM_WORLD);

  chrono.stop();

  std::cout << "Norm of the Stokes Solution: "<< stokes.solution()->norm2() << std::endl;
  un->subset(*stokes.solution(), stokes.solution()->map(), UInt(0), UInt(0));

  std::string vector="un";
  un->spy(vector);

  std::string solStokes="stokesSolution";
  stokes.solution()->spy(solStokes);
  

  fsiSolver->FSIOper()->displayer().leaderPrint( "Total computational time of Stokes solution  =  ", chrono.diff());
  fsiSolver->FSIOper()->displayer().leaderPrint( "\n" );
  /*
  std::cout << "I am here" << std::endl;
  int n;
  std::cin >> n;
  */

  std::string exporterType =  data_file( "exporter/type", "ensight" );
  std::string fluidName    =  data_file( "exporter/fluid/filename", "fluid" );

  fluidName +="Initialization"; 

  filterPtr_Type exportInitialize;

#ifdef HAVE_HDF5
  if (exporterType.compare("hdf5") == 0)
    {
      exportInitialize.reset( new  hdf5Filter_Type( data_file, fluidName) );
    }
  else
#endif
    {
      if (exporterType.compare("none") == 0)
	{
	  exportInitialize.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, uFESpacePtr->mesh(), fluidName, uFESpacePtr->map().comm().MyPID()) );
	}
      else
	{
	  exportInitialize.reset( new  ensightFilter_Type( data_file, fluidName) );	 
	}
    }
    vectorPtr_Type velAndPressure ( new vector_Type(*stokes.solution(), exportInitialize->mapType() ) );
    velAndPressure.reset( new vector_Type( stokes.getMap(), exportInitialize->mapType() ));

    exportInitialize->setMeshProcId(uFESpacePtr->mesh(), uFESpacePtr->map().comm().MyPID());

    exportInitialize->addVariable( ExporterData<FSIOperator::mesh_Type>::VectorField, "Stokes-velocity",
				   uFESpacePtr, velAndPressure, UInt(0) );
    exportInitialize->addVariable( ExporterData<FSIOperator::mesh_Type>::ScalarField, "Stokes-pressure",
				   pFESpacePtr, velAndPressure,
				   UInt(3*uFESpacePtr->dof().numTotalDof()) );

    *velAndPressure = *stokes.solution();
    exportInitialize->postProcess(0);


}



void Problem::checkGCEResult(const LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().norm2();
    if (time==0.000 && (dispNorm-834634)     /dispNorm*(dispNorm-834634)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.001 && (dispNorm-1.15328e+06)     /dispNorm*(dispNorm-1.15328e+06)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.002 && (dispNorm-943681)/dispNorm*(dispNorm-943681)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.003 && (dispNorm-825984)     /dispNorm*(dispNorm-825984)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.004 && (dispNorm-815018)     /dispNorm*(dispNorm-815018)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}


void Problem::checkCEResult(const LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().norm2();
    if (time==0.000 && (dispNorm-772280)/dispNorm*(dispNorm-772280)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.001 && (dispNorm-1.12286e+06)/dispNorm*(dispNorm-1.12286e+06)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.002 && (dispNorm-943697)/dispNorm*(dispNorm-943697)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.003 && (dispNorm-836363)/dispNorm*(dispNorm-836363)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time);
    else if (time==0.004 && (dispNorm-819303)/dispNorm*(dispNorm-819303)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}
