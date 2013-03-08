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
 *  -# solution with exact Newton (semiImplicit = false, useShapeDerivatives = true, conservativeFormulation = false)
 *  -# solution with quasi Newton (semiImplicit = false, useShapeDerivatives = true, conservativeFormulation = false)
 *  -# preconditioner choice: see the classes Monolithic and fullMonolithic
 * - Monolithic (GCE):
 *  -# solution extrapolating the fluid domain (semiImplicit = false, useShapeDerivatives = true, conservativeFormulation = false)
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
 * The time discretization is carried out using BDF methods of order 2. At the moment, even is the Newmark method is available
 * for the temporal discretization of the single problems( e.g. in test_structuralsolver), it cannot be used in the FSI framework
 * since the class TimeAdvanceNewmark is not registered as one of the possible instances of the abstrac class TimeAdvance.
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
#include "lumpedHeart.hpp"


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

    Problem ( GetPot const& data_file ) :
        M_Tstart (0.),
        M_saveEvery (1),
        M_returnValue (EXIT_FAILURE)

    {
        using namespace LifeV;

        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "linearVenantKirchhoff", &FSIOperator::createVenantKirchhoffLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "exponential", &FSIOperator::createExponentialMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "neoHookean", &FSIOperator::createNeoHookeanMaterialNonLinear );
        FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "nonLinearVenantKirchhoff", &FSIOperator::createVenantKirchhoffNonLinear );

        std::cout << "register MonolithicGE : " << FSIMonolithicGE::S_register << std::endl;
        std::cout << "register MonolithicGI : " << FSIMonolithicGI::S_register << std::endl;

        M_data = dataPtr_Type ( new data_Type() );
        M_data->setup ( data_file );
        //M_data->dataSolid()->setTimeData( M_data->dataFluid()->dataTime() ); //Same TimeData for fluid & solid
        //M_data->showMe();

        M_fsi = fsi_solver_ptr ( new FSISolver( ) );
        MPI_Barrier ( MPI_COMM_WORLD );

        M_fsi->setData ( M_data );
        M_fsi->FSIOper()->setDataFile ( data_file ); //TO BE REMOVED!

        // Setting FESpace and DOF

        std::string  fluidMeshPartitioned    =  data_file ( "problem/fluidMeshPartitioned", "none" );
        std::string  solidMeshPartitioned    =  data_file ( "problem/solidMeshPartitioned", "none" );
#ifdef HAVE_HDF5
        if ( fluidMeshPartitioned.compare ( "none" ) )
        {
            FSIOperator::meshFilter_Type fluidMeshFilter ( data_file, fluidMeshPartitioned );
            fluidMeshFilter.setComm ( M_fsi->FSIOper()->worldComm() );
            FSIOperator::meshFilter_Type solidMeshFilter ( data_file, solidMeshPartitioned );
            solidMeshFilter.setComm ( M_fsi->FSIOper( )->worldComm( ) );
            M_fsi->FSIOper( )->partitionMeshes ( fluidMeshFilter, solidMeshFilter );
            M_fsi->FSIOper( )->setupFEspace( );
            M_fsi->FSIOper( )->setupDOF ( fluidMeshFilter );
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

        debugStream ( 10000 ) << "Setting up the FESpace and DOF \n";

        MPI_Barrier ( MPI_COMM_WORLD );

#ifdef DEBUG
        debugStream ( 10000 ) << "Setting up the BC \n";
#endif
        M_fsi->setFluidBC ( BCh_monolithicFlux ( true ) );
        M_fsi->setSolidBC ( BCh_monolithicSolid ( *M_fsi->FSIOper( ) ) );

        M_fsi->setup();

        M_fsi->setFluidBC ( BCh_monolithicFluid ( *M_fsi->FSIOper( ), true ) );
        M_fsi->setHarmonicExtensionBC ( BCh_harmonicExtension ( *M_fsi->FSIOper( ) ) );

        dynamic_cast<LifeV::FSIMonolithic*> (M_fsi->FSIOper().get() )->mergeBCHandlers();

#ifdef DEBUG
        debugStream ( 10000 ) << "BC set\n";
#endif

        std::string const exporterType =  data_file ( "exporter/type", "ensight" );
        std::string const fluidName    =  data_file ( "exporter/fluid/filename", "fluid" );
        std::string const solidName    =  data_file ( "exporter/solid/filename", "solid" );

#ifdef HAVE_HDF5
        if (exporterType.compare ("hdf5") == 0)
        {
            M_exporterFluid.reset ( new  hdf5Filter_Type ( data_file, fluidName) );
            M_exporterSolid.reset ( new  hdf5Filter_Type ( data_file, solidName) );
        }
        else
#endif
        {
            if (exporterType.compare ("none") == 0)
            {
                M_exporterFluid.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), fluidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID() ) );
                M_exporterSolid.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), solidName, M_fsi->FSIOper()->uFESpace().map().comm().MyPID() ) );
            }
            else
            {
                M_exporterFluid.reset ( new  ensightFilter_Type ( data_file, fluidName) );
                M_exporterSolid.reset ( new  ensightFilter_Type ( data_file, solidName) );
            }
        }


        // load using ensight/hdf5
        M_saveEvery = data_file ("exporter/saveEvery", 1);

        M_fsi->initializeMonolithicOperator();

        M_velAndPressure.reset ( new vector_Type ( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluid->mapType() ) );

        M_fluidDisp.reset     ( new vector_Type ( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ) );

        M_solidDisp.reset ( new vector_Type ( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ) );

        M_exporterFluid->setMeshProcId (M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID() );
        M_exporterSolid->setMeshProcId (M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID() );
        M_exporterFluid->addVariable ( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-velocity",
                                       M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressure, UInt (0) );
        M_exporterFluid->addVariable ( ExporterData<FSIOperator::mesh_Type>::ScalarField, "f-pressure",
                                       M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressure,
                                       UInt (3 * M_fsi->FSIOper()->uFESpace().dof().numTotalDof() ) );

        M_exporterFluid->addVariable ( ExporterData<FSIOperator::mesh_Type>::VectorField, "f-displacement",
                                       M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDisp, UInt (0) );




        M_exporterSolid->addVariable ( ExporterData<FSIOperator::mesh_Type>::VectorField, "s-displacement",
                                       M_fsi->FSIOper()->dFESpacePtr(), M_solidDisp, UInt (0) );

        //M_fsi->FSIOper()->fluid().setupPostProc(); //this has to be called if we want to initialize the postProcess

        FC0.initParameters ( *M_fsi->FSIOper(), 3);
        LH.initParameters ( *M_fsi->FSIOper(), "dataHM");

        M_data->dataFluid()->dataTime()->setInitialTime ( M_Tstart );
        M_data->dataFluid()->dataTime()->setTime ( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->dataSolid()->dataTime()->setInitialTime ( M_Tstart );
        M_data->dataSolid()->dataTime()->setTime ( M_data->dataFluid()->dataTime()->initialTime() );
        M_data->timeDataALE()->setInitialTime ( M_Tstart );
        M_data->timeDataALE()->setTime ( M_data->dataFluid()->dataTime()->initialTime() );
    }

    /*!
      This routine runs the temporal loop
     */
    int
    run()
    {
        boost::timer _overall_timer;

        LifeV::UInt iter = 1;
        //LifeV::UInt offset=dynamic_cast<LifeV::FSIMonolithic*>(M_fsi->FSIOper().get())->offset();

        dynamic_cast<LifeV::FSIMonolithic*> (M_fsi->FSIOper().get() )->enableStressComputation (1);

        bool valveIsOpen = true;

        vectorPtr_Type solution ( new vector_Type ( (*M_fsi->FSIOper()->couplingVariableMap() ) ) );

        M_fsi->FSIOper()->extrapolation ( *solution );

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(), M_data->dataSolid()->dataTime()->updateTime(), ++iter)
        {
            //Return value for the testsuite
            M_returnValue = EXIT_FAILURE;

            LifeV::Real flux = M_fsi->FSIOper()->fluid().flux (2, M_fsi->displacement() );
            if ( valveIsOpen)
            {
                if ( iter == 3 /*flux < -100*/)
                {
                    valveIsOpen = false;
                    M_fsi->setFluidBC (BCh_monolithicFluid (*M_fsi->FSIOper(), valveIsOpen) );
                    //M_fsi->FSIOper()->BCh_fluid()->substituteBC( (const LifeV::bcFlag_Type) 2, bcf,  LifeV::Essential, LifeV::Full, (const LifeV::UInt) 3);
                }
            }
            // close the valve
            else
            {
                if (false && M_fsi->FSIOper()->fluid().pressure (2, M_fsi->displacement() ) < LifeV::LumpedHeart::M_pressure )
                {
                    valveIsOpen = true;
                    M_fsi->setFluidBC (BCh_monolithicFluid (*M_fsi->FSIOper(), valveIsOpen) );
                    //M_fsi->FSIOper()->BCh_fluid()->substituteBC( (const LifeV::bcFlag_Type) 2, bcf,  LifeV::Natural, LifeV::Full, 3);
                }
            }

            int flag = 2;
            FC0.renewParameters ( *M_fsi, 3 );
            LH.renewParameters ( *M_fsi->FSIOper(), flag, M_data->dataFluid()->dataTime()->time(), flux );

            boost::timer _timer;

            if (iter % M_saveEvery == 0)
            {
                M_fsi->FSIOper()->exportSolidDisplacement (*M_solidDisp);

                M_fsi->FSIOper()->exportFluidVelocityAndPressure (*M_velAndPressure);
                M_exporterSolid->postProcess ( M_data->dataFluid()->dataTime()->time() );

                *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();
                M_exporterFluid->postProcess ( M_data->dataFluid()->dataTime()->time() );
            }

            // This is just the previous solution. Should use the extrapolation from time advance
            M_fsi->FSIOper()->extrapolation ( *solution );

            M_fsi->iterate ( solution );

            // shift_right of the solution of all the time advance classes in the FSIOperator
            M_fsi->FSIOper()->updateSolution ( *solution );

            M_fsi->FSIOper()->displayer().leaderPrintMax ("[fsi_run] Iteration ", iter);
            M_fsi->FSIOper()->displayer().leaderPrintMax (" was done in : ", _timer.elapsed() );

            std::cout << "solution norm " << iter << " : "
                      << M_fsi->displacement().norm2() << "\n";

            //     ///////// CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            if (!M_data->method().compare ("monolithicGI") )
            {
                checkCEResult (M_data->dataFluid()->dataTime()->time() );
            }
            else
            {
                checkGCEResult (M_data->dataFluid()->dataTime()->time() );
            }

        }

        std::cout << "Total computation time = "
                  << _overall_timer.elapsed() << "s" << "\n";

        return M_returnValue;

    }

private:

    void checkCEResult (const LifeV::Real& time);
    void checkGCEResult (const LifeV::Real& time);

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
    LifeV::LumpedHeart LH;

    LifeV::Real    M_Tstart;
    // vectorPtr_Type M_WS;
    LifeV::UInt           M_saveEvery;
    LifeV::UInt M_returnValue;

public:

    void resultCorrect (LifeV::Real time)
    {
        std::cout << "Result correct at time: " << time << std::endl;
        M_returnValue = EXIT_SUCCESS;
    }
};


struct FSIChecker
{
    FSIChecker ( GetPot const& _data_file ) :
        data_file ( _data_file )
    {}

    int operator() ()
    {
        boost::shared_ptr<Problem> fsip;

        try
        {
            fsip = boost::shared_ptr<Problem> ( new Problem ( data_file ) );
            return fsip->run();
        }
        catch ( std::exception const& _ex )
        {
            std::cout << "caught exception :  " << _ex.what() << "\n";
        }
        return EXIT_FAILURE;

    }

    GetPot                data_file;
    LifeV::Vector         disp;
};


namespace LifeV
{

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}
}


int main (int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#else
    std::cout << "% using serial Version" << std::endl;
#endif

    GetPot command_line (argc, argv);

    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot data_file (data_file_name);
    FSIChecker _sp_check ( data_file );
    int returnValue = _sp_check();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return returnValue;

}

void Problem::checkCEResult (const LifeV::Real& time)
{
    LifeV::Real dispNorm = M_fsi->displacement().norm2();
    if (time == 0.000 && (dispNorm - 106344) / dispNorm * (dispNorm - 106344) / dispNorm < 1e-3)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.001 && (dispNorm - 147532) / dispNorm * (dispNorm - 147532) / dispNorm < 1e-3)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.002 && (dispNorm - 108216) / dispNorm * (dispNorm - 108216) / dispNorm < 1e-3)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.003 && (dispNorm - 105437) / dispNorm * (dispNorm - 105437) / dispNorm < 1e-3)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.004 && (dispNorm - 104585) / dispNorm * (dispNorm - 104585) / dispNorm < 1e-3)
    {
        Problem::resultCorrect (time);
    }

}


void Problem::checkGCEResult (const LifeV::Real& time)
{
    LifeV::Real dispNorm = M_fsi->displacement().norm2();
    if (time == 0.000 && (dispNorm - 110316) / dispNorm * (dispNorm - 110316) / dispNorm < 1e-5)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.001 && (dispNorm - 99468.8) / dispNorm * (dispNorm - 99468.8) / dispNorm < 1e-5)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.002 && (dispNorm - 91003.1) / dispNorm * (dispNorm - 91003.1) / dispNorm < 1e-5)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.003 && (dispNorm - 90179.9) / dispNorm * (dispNorm - 90179.9) / dispNorm < 1e-5)
    {
        Problem::resultCorrect (time);
    }
    else if (time == 0.004 && (dispNorm - 88319.3) / dispNorm * (dispNorm - 88319.3) / dispNorm < 1e-5)
    {
        Problem::resultCorrect (time);
    }
}
