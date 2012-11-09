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
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @maintainer Paolo Crosetto <crosetto@iacspc70.epfl.ch>
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
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/fsi/solver/FSISolver.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/fsi/solver/FSIMonolithicGI.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

#include "ud_functions.hpp"
#include "boundaryConditions.hpp"
#include "flowConditions.hpp"

class Problem
{
public:

    typedef boost::shared_ptr<LifeV::FSISolver> fsi_solver_ptr;

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
    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */

    Problem( GetPot const& data_file ):
        M_Tstart(0.),
        M_saveEvery(1),
        M_returnValue(EXIT_FAILURE)
    {
        using namespace LifeV;

        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhoff", &FSIOperator::createVenantKirchhoffLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "exponential", &FSIOperator::createExponentialMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "neoHookean", &FSIOperator::createNeoHookeanMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "nonLinearVenantKirchhoff", &FSIOperator::createVenantKirchhoffNonLinear );

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

        M_fsi->setup();

        M_fsi->setFluidBC( BCh_monolithicFluid( *M_fsi->FSIOper( ), true ) );
        M_fsi->setHarmonicExtensionBC( BCh_harmonicExtension( *M_fsi->FSIOper( ) ) );

        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->mergeBCHandlers();

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
	std::string restartType(data_file("importer/restartFSI", "false" ));
        std::cout << "The load state is: "<< restartType << std::endl;

	if ( !restartType.compare("true") )
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

    //M_fsi->FSIOper()->fluid().setupPostProc(); //this has to be called if we want to initialize the postProcess

        FC0.initParameters( *M_fsi->FSIOper(), 3);

        M_data->dataFluid()->dataTime()->setInitialTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setInitialTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->timeDataALE()->setInitialTime( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->timeDataALE()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
    }

    /*!
      This routine runs the temporal loop
    */
    int
    run()
    {
        boost::timer _overall_timer;

        LifeV::UInt iter = 1;

        dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->enableStressComputation(0);

        vectorPtr_Type solution ( new vector_Type( (*M_fsi->FSIOper()->couplingVariableMap()) ) );

        M_fsi->FSIOper()->extrapolation( *solution );

        //Initialize the exporters
        M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDisp);//    displacement(), M_offset);
        M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressure);
        *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();
        M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );
        M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );

        for ( ; M_data->dataFluid()->dataTime()->canAdvance();  ++iter)
        {
            M_returnValue = EXIT_FAILURE;
            M_data->dataFluid()->dataTime()->updateTime();
            M_data->dataSolid()->dataTime()->updateTime();
            M_data->timeDataALE()->updateTime();


            boost::timer _timer;
            FC0.renewParameters( *M_fsi, 3 );

            M_fsi->FSIOper()->extrapolation( *solution );

            M_fsi->iterate( solution );

            // Saving the solution
            if( M_data->method().compare("monolithicGI") == 0 )
            {
                M_fsi->FSIOper()->updateSolution( *solution );
            }
            else
            {
                M_fsi->FSIOper()->updateSolution( *solution );
            }

            if(iter%M_saveEvery==0)
            {
                M_fsi->FSIOper()->exportSolidDisplacement(*M_solidDisp);
                M_fsi->FSIOper()->exportFluidVelocityAndPressure(*M_velAndPressure);
                *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();

                M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );
                M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );
            }

            std::cout << "solution norm at time: " <<  M_data->dataFluid()->dataTime()->time() << "(iter" << iter << ") : "
                      << M_fsi->displacement().norm2() << "\n";

            checkResult(M_data->dataFluid()->dataTime()->time());


        }

        std::cout << "Total computation time = "
                  << _overall_timer.elapsed() << "s" << "\n";

        return M_returnValue;

    }

private:

    void restartFSI(  GetPot const& data_file);
    void checkResult(const LifeV::Real& time);

    fsi_solver_ptr M_fsi;
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
    LifeV::UInt    M_saveEvery;
    LifeV::UInt    M_returnValue;
public:
  void resultCorrect(LifeV::Real time)
  {
    std::cout<<"Result correct at time: " << time << std::endl;
    M_returnValue = EXIT_SUCCESS;
  }

};

struct FSIChecker
{
    FSIChecker( GetPot const& _data_file ):
            data_file( _data_file )
    {}

    int operator()()
    {
        boost::shared_ptr<Problem> fsip;
        int returnVariable;
        returnVariable = EXIT_FAILURE;
        try
        {
            fsip = boost::shared_ptr<Problem>( new Problem( data_file ) );
            returnVariable = fsip->run();

        }
        catch ( std::exception const& _ex )
        {
            std::cout << "caught exception :  " << _ex.what() << "\n";
        }

        return returnVariable;
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
        const std::string data_file_name = command_line.follow("data", 2, "-f","--file");
        GetPot data_file(data_file_name);
        FSIChecker _sp_check( data_file );
        int returnValue = _sp_check();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return returnValue;

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
    std::string iterationStringCopy ;

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
    vectorPtr_Type initFluid        (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
    vectorPtr_Type HarmonicSol      (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
    vectorPtr_Type structureSol     (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));
    vectorPtr_Type temporarySol     (new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero));

    vectorPtr_Type vel           (new vector_Type(M_fsi->FSIOper()->uFESpace().map(), LifeV::Unique));
    vectorPtr_Type pressure      (new vector_Type(M_fsi->FSIOper()->pFESpace().map(), LifeV::Unique));
    vectorPtr_Type solidDisp     (new vector_Type(M_fsi->FSIOper()->dFESpace().map(), LifeV::Unique));
    vectorPtr_Type fluidDisp     (new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), LifeV::Unique));

    vectorPtr_Type firstFluidDisp(new vector_Type(M_fsi->FSIOper()->mmFESpace().map(), LifeV::Unique));

    //The hypothesis used for this method is that the three TimeAdvance classes have the same size
    iterationString = loadInitSol;

    UInt iterInit;

    for(iterInit = 0; iterInit<M_fsi->FSIOper()->fluidTimeAdvance()->size(); iterInit++ )
      {

	//First the three fields are read at the same time
	//Creating the exporter data for the fields
	//Fluid
	LifeV::ExporterData<mesh_Type> initSolFluidVel   (LifeV::ExporterData<mesh_Type>::VectorField, std::string("f-velocity."+iterationString), M_fsi->FSIOper()->uFESpacePtr(), vel, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
	LifeV::ExporterData<mesh_Type> initSolFluidPress (LifeV::ExporterData<mesh_Type>::ScalarField, std::string("f-pressure."+iterationString), M_fsi->FSIOper()->pFESpacePtr(), pressure, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

	//Structure
	LifeV::ExporterData<mesh_Type> initSolSolidDisp  (LifeV::ExporterData<mesh_Type>::VectorField,"s-displacement."+iterationString, M_fsi->FSIOper()->dFESpacePtr(), solidDisp, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );
	//ALE
	LifeV::ExporterData<mesh_Type> initSolFluidDisp  (LifeV::ExporterData<mesh_Type>::VectorField, "f-displacement."+iterationString, M_fsi->FSIOper()->mmFESpacePtr(), fluidDisp, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

	//Initialization of the vectors used to read
	*vel *= 0.0;
	*pressure *= 0.0;
	*solidDisp *= 0.0;
	*fluidDisp *= 0.0;

	//load of the solutions
	M_importerFluid->readVariable(initSolFluidVel);   //Fluid u
	M_importerFluid->readVariable(initSolFluidPress); //Fluid p
	M_importerSolid->readVariable(initSolSolidDisp);  //Solid d
	M_importerFluid->readVariable(initSolFluidDisp);  //Fluid df

	// std::cout << "Norm of the vel " << vel->norm2() << std::endl;
	// std::cout << "Norm of the pressure " << pressure->norm2() << std::endl;
	// std::cout << "Norm of the solid " << solidDisp->norm2() << std::endl;
	// std::cout << "Norm of the df " << fluidDisp->norm2() << std::endl;

	//We send the vectors to the FSIMonolithic class using the interface of FSIOper
	M_fsi->FSIOper()->setVectorInStencils(vel, pressure, solidDisp, fluidDisp, iterInit )
;
	//Updating string name
	int iterations = std::atoi(iterationString.c_str());
	iterations--;

	std::ostringstream iter;
	iter.fill( '0' );
	iter << std::setw(5) << ( iterations );
	iterationString=iter.str();

	iterationStringCopy = iterationString;

      }

    //Reading another vector for the solidTimeAdvance since its BDF has the same order
    //as the other ones but since the orderDerivative = 2, the size of the stencil is
    //orderBDF + 1

    *solidDisp *= 0.0;
    LifeV::ExporterData<mesh_Type> initSolSolidDisp  (LifeV::ExporterData<mesh_Type>::VectorField,"s-displacement."+iterationString, M_fsi->FSIOper()->dFESpacePtr(), solidDisp, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    M_importerSolid->readVariable(initSolSolidDisp);  //Solid d

    //Solid problem
    vectorPtr_Type vectorMonolithicSolidDisplacement(new vector_Type(*M_fsi->FSIOper()->couplingVariableMap(), Unique, Zero) );
    *vectorMonolithicSolidDisplacement *= 0.0;

    //This dynamic cast is ugly but it's the only wat to understand where the read vector can be put.
    UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();
    vectorMonolithicSolidDisplacement->subset( *solidDisp, solidDisp->map(), (UInt)0, offset);
    *vectorMonolithicSolidDisplacement *= 1.0 / (M_fsi->FSIOper()->solid().rescaleFactor());

    vector_Type* normalPointerToSolidVector( new vector_Type(*vectorMonolithicSolidDisplacement) );
    (M_fsi->FSIOper()->solidTimeAdvance()->stencil()).push_back( normalPointerToSolidVector );

    //Set the initialRHS for the TimeAdvance classes
    vector_Type zeroFluidSolid(*M_fsi->FSIOper()->couplingVariableMap(), LifeV::Unique, Zero);
    vector_Type zeroALE(M_fsi->FSIOper()->mmFESpace().map(), LifeV::Unique, Zero);

    zeroFluidSolid *= 0.0;
    zeroALE *= 0.0;

    M_fsi->FSIOper()->fluidTimeAdvance()->setInitialRHS(zeroFluidSolid);
    M_fsi->FSIOper()->solidTimeAdvance()->setInitialRHS(zeroFluidSolid);
    M_fsi->FSIOper()->ALETimeAdvance()->setInitialRHS(zeroALE);

    //This updates at the current value (the one when the simulation was stopped) the RHScontribution
    //of the first derivative which is use to compute the velocity in TimeAdvance::velocity().
    //Please note that, even if it is ugly, at this stage, the fluidTimeAdvance is leading the Time Discretization
    //and this is why there  is the dataFluid class to get the dt.
    M_fsi->FSIOper()->ALETimeAdvance()->updateRHSFirstDerivative( M_data->dataFluid()->dataTime()->timeStep() );

    //This are used to export the loaded solution to check it is correct.
    M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(), M_importerFluid->mapType() ));
    M_velAndPressure->subset(*pressure, pressure->map(), UInt(0), (UInt)3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof());
    *M_velAndPressure += *vel;

    M_fluidDisp.reset     ( new vector_Type( *fluidDisp, M_importerFluid->mapType() ));

    M_solidDisp.reset     ( new vector_Type( *solidDisp, M_importerSolid->mapType() ));
}


void Problem::checkResult(const LifeV::Real& time)
{

  //Extract the previous solution
  LifeV::Real dispNorm=M_fsi->displacement().norm2();
  if (time==0.006 &&      (dispNorm-100221)/dispNorm * (dispNorm-100221)/dispNorm < 1e-5) resultCorrect(time);
  else if (time==0.007 && (dispNorm-94940.1)/dispNorm * (dispNorm-94940.1)/dispNorm < 1e-5) resultCorrect(time);
  else if (time==0.008 && (dispNorm-91910.8)/dispNorm * (dispNorm-91910.8)/dispNorm < 1e-5) resultCorrect(time);
}
