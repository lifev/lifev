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
//#include <life/lifesolver/VenantKirchhoffSolverLinear.hpp>
#include <life/lifesolver/StructuralSolver.hpp>
#include <life/lifesolver/FSIMonolithicGI.hpp>

#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/ExporterHDF5.hpp>
#endif

#include "ud_functions.hpp"
#include "boundaryConditions.hpp"
#include "flowConditions.hpp"
#include "lumpedHeart.hpp"

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


	/*Initialization*/

        // load using ensight/hdf5
        std::string restartType(data_file("importer/restartType","newSimulation"));
	std::cout << "The load state is: "<< restartType << std::endl;

	if (!restartType.compare("Stokes"))
	  initializeStokes(restartType, data_file);
        else if (!restartType.compare("restartFSI"))
	  restartFSI(restartType,data_file);
	else
	  M_fsi->initialize();

        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->mergeBCHandlers();

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
	  //	  LifeV::Real flux=M_fsi->FSIOper()->fluid().flux(2, M_fsi->displacement());
	  //if ( false && valveIsOpen)
          //   {
          //       if ( iter == 3 /*flux < -100*/)
          //       {
          //           valveIsOpen=false;
          //           M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->FSIOper(), valveIsOpen));
                     //M_fsi->FSIOper()->BCh_fluid()->substituteBC( (const LifeV::bcFlag_Type) 2, bcf,  LifeV::Essential, LifeV::Full, (const LifeV::UInt) 3);
          //       }
          //   }
          //   // close the valve
          //   else
          //   {
          //       if (false && M_fsi->FSIOper()->fluid().pressure(2, M_fsi->displacement()) < LifeV::LumpedHeart::M_pressure )
	  //                 {
          //           valveIsOpen=true;
          //           M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->FSIOper(), valveIsOpen));
          //           //M_fsi->FSIOper()->BCh_fluid()->substituteBC( (const LifeV::bcFlag_Type) 2, bcf,  LifeV::Natural, LifeV::Full, 3);
          //       }
          //   }

          //   int flag =2;
             FC0.renewParameters( *M_fsi, 3 );
	     //   LH.renewParameters( *M_fsi->FSIOper(), flag, M_data->dataFluid()->dataTime()->time(), flux );

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
	    }*/

             FC0.renewParameters( *M_fsi, 3 );
        std::cout << "Total computation time = "
                  << _overall_timer.elapsed() << "s" << "\n";

    }

private:

    void initializeStokes(std::string& restartType,  GetPot const& data_file);
    void restartFSI(std::string& restartType,  GetPot const& data_file);

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

  
  /*!These are shared_ptr to the vectors in order to pass them to the initialization */ 
  vectorPtr_Type velocityAndPressure    (new vector_Type(M_fsi->FSIOper()->fluid().getMap()  ));
  vectorPtr_Type fluidSolution          (new vector_Type(M_fsi->FSIOper()->fluid().getMap()  ));
  vectorPtr_Type solidDisplacement      (new vector_Type(M_fsi->FSIOper()->solid().map() ));
  vectorPtr_Type solidVelocity          (new vector_Type(M_fsi->FSIOper()->solid().map() ));
  vectorPtr_Type solidAcceleration      (new vector_Type(M_fsi->FSIOper()->solid().map() ));
  vectorPtr_Type fluidDisplacement      (new vector_Type(M_fsi->FSIOper()->meshMotion().getMap() ));
  vectorPtr_Type fluidDisplacementOld   (new vector_Type(M_fsi->FSIOper()->meshMotion().getMap() ));

  vectorPtr_Type initialSolution        (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap()));
  vectorPtr_Type temporarySol           (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap()));

  UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

  /*! These are the vectors wich will contain the loaded data*/
  /*! Fluid */
  vectorPtr_Type vel           (new vector_Type(M_fsi->FSIOper()->uFESpace().map(), Unique ));
  vectorPtr_Type pressure      (new vector_Type(M_fsi->FSIOper()->pFESpace().map(), Unique ));
  /*! Solid */
  vectorPtr_Type solidDisp     (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), Unique ));
  vectorPtr_Type solidVel      (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), Unique ));
  vectorPtr_Type solidAcc      (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), Unique ));
  /*! Geometric */
  vectorPtr_Type fluidDisp     (new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), Unique));
  vectorPtr_Type fluidDispPrev (new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), M_importerFluid->mapType()));

  /*! Creating ExporterData object to be used to read later */
  LifeV::ExporterData<mesh_Type> initSolFluidVel   (LifeV::ExporterData<mesh_Type>::VectorField, std::string("f-velocity."+loadInitSol), M_fsi->FSIOper()->uFESpacePtr(), vel, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
  LifeV::ExporterData<mesh_Type> initSolFluidPress (LifeV::ExporterData<mesh_Type>::ScalarField, "f-pressure."+loadInitSol, M_fsi->FSIOper()->pFESpacePtr(), pressure, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
  LifeV::ExporterData<mesh_Type> initSolSolidDisp  (LifeV::ExporterData<mesh_Type>::VectorField,"s-displacement."+loadInitSol, M_fsi->FSIOper()->dFESpacePtr(), solidDisp, UInt(0)/*offset*/, LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
  LifeV::ExporterData<mesh_Type> initSolSolidVel   (LifeV::ExporterData<mesh_Type>::VectorField,"s-velocity."+loadInitSol, M_fsi->FSIOper()->dFESpacePtr(), solidVel, UInt(0)/*offset*/, LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
  LifeV::ExporterData<mesh_Type> initSolSolidAcc   (LifeV::ExporterData<mesh_Type>::VectorField,"s-acceleration."+loadInitSol, M_fsi->FSIOper()->dFESpacePtr(), solidAcc, UInt(0)/*offset*/, LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
  LifeV::ExporterData<mesh_Type> initSolFluidDisp  (LifeV::ExporterData<mesh_Type>::VectorField, "f-displacement."+loadInitSol, M_fsi->FSIOper()->mmFESpacePtr(), fluidDisp, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
  LifeV::ExporterData<mesh_Type> initSolFluidDispPrev(LifeV::ExporterData<mesh_Type>::VectorField, "f-displacement."+loadInitSolFD, M_fsi->FSIOper()->mmFESpacePtr(), fluidDispPrev, 0, LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
  
  /*!load of the solutions*/
  M_importerFluid->readVariable(initSolFluidVel);
  M_importerFluid->readVariable(initSolFluidPress);
  M_importerSolid->readVariable(initSolSolidDisp);
  M_importerSolid->readVariable(initSolSolidVel);
  M_importerSolid->readVariable(initSolSolidAcc);
  M_importerFluid->readVariable(initSolFluidDisp);
  M_importerFluid->readVariable(initSolFluidDispPrev);


  /*! Importing the loaded solution to Unique vectors, to avoid that the repeated part is summed:*/
  /*! initializing the mesh displacement (GCE case) */
  fluidDisplacementOld.reset(new vector_Type(*fluidDispPrev, Unique, Zero));

  velocityAndPressure.reset( new vector_Type(*pressure, Unique, Zero));
  fluidSolution->subset(*velocityAndPressure, velocityAndPressure->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());
  initialSolution->subset(*velocityAndPressure, velocityAndPressure->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());

 
  vector_Type tmpVec(*fluidSolution);
  velocityAndPressure.reset( new vector_Type(*vel, Unique, Zero));
  tmpVec=*velocityAndPressure;//possibly an Import here
  *fluidSolution += tmpVec;
  *initialSolution += tmpVec;

  std::vector<vectorPtr_Type> u0;
  u0.push_back(fluidSolution);

  //Structure (setting of the pointers + rescaling)
  solidDisplacement.reset( new vector_Type(*solidDisp, Unique, Zero));
  *solidDisplacement *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataSolid()->dataTime()->timeStep());

  temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), Unique, Zero)); 
  temporarySol->subset(*solidDisplacement, solidDisplacement->map(), (UInt)0, offset);
  *initialSolution += *temporarySol;

  solidVelocity.reset(new vector_Type(*solidVel, Unique, Zero));
  *solidVelocity *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataSolid()->dataTime()->timeStep());

  solidAcceleration.reset(new vector_Type(*solidAcc, Unique, Zero));
  *solidAcceleration *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataSolid()->dataTime()->timeStep());

  std::vector<vectorPtr_Type> ds0;
  ds0.push_back(solidDisplacement);
  ds0.push_back(solidVelocity);
  ds0.push_back(solidAcceleration);

  //Geometric Problem

  if (!M_data->method().compare("monolithicGI"))
    {
      fluidDisplacement.reset( new vector_Type(*fluidDisp, Unique, Zero));
      temporarySol.reset(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), Unique, Zero)); 
      temporarySol->subset(*fluidDisp, fluidDisp->map(), (UInt)0, dynamic_cast<LifeV::FSIMonolithicGI*>(M_fsi->FSIOper().get())->mapWithoutMesh().map(Unique)->NumGlobalElements());
      *initialSolution += *temporarySol;
    }
  else
    {
      *fluidDisplacement=*fluidDisplacementOld;
    }
  std::vector<vectorPtr_Type> df0;
  df0.push_back(fluidDisplacementOld);

  //Initialize the timeAdvances Procedures
  M_fsi->initialize(u0,ds0,df0);

  // Initialize the solution of the monolithic system
  M_fsi->FSIOper()->setSolution(*initialSolution);

  std::cout << std::endl;
  std::cout << "The norm inf of the solution for the fluid is              : " << M_fsi->displacement().normInf() << std::endl;
  std::cout << "The norm 2 of the solution for the fluid is                : " << M_fsi->displacement().norm2() << std::endl;
  std::cout << "The norm of the Fluid solution for the fluid is            : " << fluidSolution->norm2() << std::endl;
  std::cout << "The norm of the Solid displacement for the fluid is        : " << solidDisplacement->norm2() << std::endl;
  std::cout << "The norm of the Solid velocity for the fluid is            : " << solidVelocity->norm2() << std::endl;
  std::cout << "The norm of the Solid acceleration for the fluid is        : " << solidAcceleration->norm2() << std::endl;
  std::cout << "The norm of the HarmonicExtension solution for the fluid is: " << fluidDisplacement->norm2() << std::endl;

}

void Problem::initializeStokes(std::string& restartType,  GetPot const& data_file)
{

  using namespace LifeV;

  typedef FSIOperator::mesh_Type        mesh_Type;


  //Creating the pointer to the filter
  filterPtr_Type importer;
  std::string const importerType =  data_file( "importer/type", "hdf5");
  std::string const filename    = data_file( "importer/fluid/filename", "fluid");

  std::string const loadInitSol      = data_file( "importer/initSol", "00000");
  std::string const loadInitSolFD    = data_file("importer/initSolFD","-1");

  std::cout << "The filename is    : " << filename << std::endl;
  std::cout << "The importerType is: " << importerType << std::endl;

  if (importerType.compare("hdf5") == 0)
    {
      importer.reset( new  hdf5Filter_Type( data_file, filename) );
    }
  else if (importerType.compare("none") == 0)
    {
      importer.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), "fluid", M_fsi->FSIOper()->uFESpace().map().comm().MyPID()));
    }
  else
    {
      importer.reset( new  ensightFilter_Type( data_file, filename) );
    }
   
  importer->setMeshProcId( M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());

  
  vectorPtr_Type initSol (new vector_Type(M_fsi->FSIOper()->fluid().getMap() ) );
  vectorPtr_Type UniqueV (new vector_Type(M_fsi->FSIOper()->fluid().getMap() ) );
  vectorPtr_Type vel     (new vector_Type(M_fsi->FSIOper()->uFESpace().map(), Unique ));
  vectorPtr_Type pressure(new vector_Type(M_fsi->FSIOper()->pFESpace().map(), Unique ));

  //In this case, according to the interface of FSISolver, two nonempty zero vectors are needed for S and HE problems
  vectorPtr_Type solidDisp     (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), Unique));
  vectorPtr_Type solidVel      (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), Unique));
  vectorPtr_Type solidAcc      (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), Unique));
  vectorPtr_Type fluidDisp     (new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), Unique));


  LifeV::ExporterData<mesh_Type> initSolFluidVel   (LifeV::ExporterData<mesh_Type>::VectorField, std::string("velocity."+loadInitSol), M_fsi->FSIOper()->uFESpacePtr(), vel, UInt(0), LifeV::ExporterData<mesh_Type>::SteadyRegime );
  LifeV::ExporterData<mesh_Type> initSolFluidPress (LifeV::ExporterData<mesh_Type>::ScalarField, std::string("pressure."+loadInitSol), M_fsi->FSIOper()->pFESpacePtr(), pressure, UInt(0), LifeV::ExporterData<mesh_Type>::SteadyRegime );
  
  importer->readVariable(initSolFluidVel);
  importer->readVariable(initSolFluidPress);

  UniqueV.reset( new vector_Type(*pressure, Unique, Zero));
  initSol->subset(*UniqueV, UniqueV->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());

  vector_Type tmpVec(*initSol);
  UniqueV.reset( new vector_Type(*vel, Unique, Zero));
  tmpVec=*UniqueV;//possibly an Import here
  *initSol+=tmpVec;

  std::vector<vectorPtr_Type> u0;
  u0.push_back(initSol);

  //Structure
  *solidDisp *= 0.0;
  *solidVel *= 0.0;
  *solidAcc *= 0.0;
  std::vector<vectorPtr_Type> ds0;
  ds0.push_back(solidDisp);
  ds0.push_back(solidVel);
  ds0.push_back(solidAcc);


  //Fluid Problme
  *fluidDisp *= 0.0;
  std::vector<vectorPtr_Type> df0;
  df0.push_back(fluidDisp);

  M_fsi->initialize(u0,ds0,df0);
  std::cout << std::endl;
  std::cout << "The norm of the Fluid solution for the fluid is: " << initSol->norm2() << std::endl;  

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
