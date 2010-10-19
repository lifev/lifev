/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/**
   \file main.cpp
   \author Paolo Crosetto <paolo.crosetto@epfl.ch>
   \date 2009-04-09
*/

/**
 *\include fluidstructure.dox
 * @file
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
 *  -# mixte
 *  -# fluxes (defective)
 *  -# absorbing \ref BNV08 :
 *   through the class flowConditions.
 * - optional: computation of wall shear stress (not properly tested in parallel)
 * - optional: computation of the largest singular values of the preconditioned matrix
 *
 * \b Features:
 * This test by default solves the FSI probem discretized in time using the GCE or CE methods, implemented respectively
 in the files monolithicGE.hpp and monolithicGI.hpp . The geometry is that of a tube (benchmark test introduced in \ref GV03).
 In this test the boundary conditions assigned are of type:
 - flux (defective b.c.) at the inlet
 - absorbing (see \ref BNV08) at the outlet
 - Robin b.c. on the solid external wall
 - Dirichlet homogeneous at the solid rings on the inlet-outlet (clamped tube).

 The output is written at every timestep, in both ensight and HDF5 (if available) formats.
*/

#ifdef TWODIM
#error test_monolithic cannot be compiled in 2D
#endif

#include <cassert>
#include <cstdlib>

#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/life.hpp>

#include <boost/timer.hpp>

#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/NonLinearVenantKirchhofSolver.hpp>
#include <life/lifesolver/LinearVenantKirchhofSolver.hpp>
#include <lifemc/lifesolver/MonolithicGI.hpp>

#include <life/lifesolver/DataFSI.hpp>

#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/noexport.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/hdf5exporter.hpp>
#endif

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
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
    typedef LifeV::FSIOperator::data_PtrType                       data_PtrType;

    typedef LifeV::FSIOperator::vector_type        vector_type;
    typedef LifeV::FSIOperator::vector_ptrtype     vector_ptrtype;

    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh3D<LifeV::LinearTetra> > > filter_ptrtype;

    typedef LifeV::Ensight<LifeV::FSIOperator::mesh_type>  ensightfilter_type;
    typedef boost::shared_ptr<ensightfilter_type>                 ensightfilter_ptrtype;
#ifdef HAVE_HDF5
    typedef LifeV::Hdf5exporter<LifeV::FSIOperator::mesh_type>  hdf5filter_type;
    typedef boost::shared_ptr<hdf5filter_type>                  hdf5filter_ptrtype;
#endif
    typedef LifeV::singleton<LifeV::factory<LifeV::FSIOperator,  std::string> > FSIFactory;
    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */

    Problem( GetPot const& data_file, std::string _oper = "" )
    {
        using namespace LifeV;

        VenantKirchhofSolver< FSIOperator::mesh_type, SolverTrilinos >::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhof", &FSIOperator::createLinearStructure );

        VenantKirchhofSolver< FSIOperator::mesh_type, SolverTrilinos >::StructureSolverFactory::instance().registerProduct( "nonLinearVenantKirchhof", &FSIOperator::createNonLinearStructure );

        M_data = data_PtrType( new data_Type() );
        M_data->setup( data_file );
        M_data->dataSolid()->setDataTime( M_data->dataFluid()->dataTime() ); //Same dataTime for fluid & solid
        //M_data->showMe();

#ifdef DEBUG
        Debug( 10000 ) << "creating FSISolver with operator :  " << method << "\n";
#endif
        M_fsi = fsi_solver_ptr( new FSISolver( ) );
        MPI_Barrier( MPI_COMM_WORLD );

#ifdef DEBUG
        Debug( 10000 ) << "Setting up data from GetPot \n";
#endif
        M_fsi->setData( M_data );
        M_fsi->FSIOper()->setDataFile( data_file ); //TO BE REMOVED!
        MPI_Barrier( MPI_COMM_WORLD );

        // Setting FESpace and DOF

        std::string  fluidMeshPartitioned    =  data_file( "problem/fluidMeshPartitioned", "none" );
        std::string  solidMeshPartitioned    =  data_file( "problem/solidMeshPartitioned", "none" );
        if( fluidMeshPartitioned.compare( "none" ) )
        {
            FSIOperator::mesh_filtertype fluidMeshFilter( data_file, fluidMeshPartitioned );
            fluidMeshFilter.setComm( M_fsi->FSIOper()->worldComm() );
            FSIOperator::mesh_filtertype solidMeshFilter( data_file, solidMeshPartitioned );
            solidMeshFilter.setComm( M_fsi->FSIOper( )->worldComm( ) );
            M_fsi->FSIOper( )->partitionMeshes( fluidMeshFilter, solidMeshFilter );
            M_fsi->FSIOper( )->setupFEspace( );
            M_fsi->FSIOper( )->setupDOF( fluidMeshFilter );
            fluidMeshFilter.CloseFile( );
            solidMeshFilter.CloseFile( );
        }
        else
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
	    M_fsi->setFluidBC( BCh_monolithicFlux( ) );
	    M_fsi->setSolidBC( BCh_monolithicRobin( *M_fsi->FSIOper( ) ) );

	    M_fsi->setup(/*data_file*/);

	    M_fsi->setFluidBC( BCh_monolithicFluid( *M_fsi->FSIOper( ) ) );
        M_fsi->setHarmonicExtensionBC( BCh_harmonicExtension( *M_fsi->FSIOper( ) ) );
        M_fsi->setSolidBC( BCh_monolithicSolid( *M_fsi->FSIOper( ) ) );
#ifdef DEBUG
        Debug( 10000 ) << "BC set\n";
#endif

        std::string const exporterType =  data_file( "exporter/type", "ensight" );
        std::string const fluidName    =  data_file( "exporter/fluid/filename", "fluid" );
        std::string const solidName    =  data_file( "exporter/solid/filename", "solid" );

#ifdef HAVE_HDF5
        if (exporterType.compare("hdf5") == 0)
        {
            M_exporterFluid.reset( new  hdf5filter_type( data_file, fluidName) );
            M_exporterSolid.reset( new  hdf5filter_type ( data_file,solidName));
        }
        else
#endif
        {
            if (exporterType.compare("none") == 0)
            {
                M_exporterFluid.reset( new NoExport<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), fluidName, M_fsi->FSIOper()->uFESpace().map().Comm().MyPID()) );
                M_exporterSolid.reset( new NoExport<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), solidName, M_fsi->FSIOper()->uFESpace().map().Comm().MyPID()) );
            } else {
                M_exporterFluid.reset( new  ensightfilter_type( data_file, fluidName) );
                M_exporterSolid.reset( new  ensightfilter_type ( data_file, solidName) );
            }
        }
        M_velAndPressure.reset( new vector_type( M_fsi->FSIOper()->fluid().getMap(), M_exporterFluid->mapType() ));
        M_fluidDisp.reset     ( new vector_type( M_fsi->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ));

        M_exporterFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());
        M_exporterSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().Comm().MyPID());
        M_exporterFluid->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                      UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );
        M_exporterFluid->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
                                      UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()),
                                      UInt(M_fsi->FSIOper()->pFESpace().dof().numTotalDof()) );

        M_exporterFluid->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
                                      UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );



        M_solidDisp.reset( new vector_type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        M_solidVel.reset ( new vector_type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));

        M_exporterSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                                      UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
        M_exporterSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
                                      UInt(0),
                                      M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );


		// load using ensight/hdf5
		std::string loadInitSol(data_file("importer/initSol","-1"));

		if(loadInitSol.compare("-1"))
        {
            initialize(loadInitSol, data_file);
        }
        else
        {
            M_fsi->initialize();
        }
        M_fsi->FSIOper()->mergeBCHandlers();

        FC0.initParameters( *M_fsi->FSIOper(), 3);
        LH.initParameters( *M_fsi->FSIOper(), "dataHM");
//         M_data->dataFluid()->dataTime()->setInitialTime( M_data->dataFluid()->dataTime()->getInitialTime() + M_data->dataFluid()->dataTime()->getTimeStep() );
//         M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->getInitialTime() + M_data->dataFluid()->dataTime()->getTimeStep());
    }

    fsi_solver_ptr fsiSolver() { return M_fsi; }

    data_PtrType fsiData() { return M_data; }

    /*!
      This routine runs the temporal loop
    */
    void
    run()
    {
        boost::timer _overall_timer;
        M_Tstart=M_fsi->FSIOper()->dataFluid()->dataTime()->getInitialTime();
        int _i = 1;
        LifeV::UInt offset=dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->getOffset();

#ifdef HAVE_HDF5
        if (M_exporterFluid->mapType() == LifeV::Unique)
        {
            M_exporterFluid->postProcess( M_Tstart );//ugly way to avoid that hdf5 starts with a deformed mesh
            M_exporterSolid->postProcess( M_Tstart );//ugly way to avoid that hdf5 starts with a deformed mesh
        }
#endif

        bool valveIsOpen=true;

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(), ++_i)
        {
            LifeV::Real flux=M_fsi->FSIOper()->fluid().flux(2, M_fsi->displacement());
            //std::cout<<"flux : "<<flux<<std::endl;
            if ( valveIsOpen)
            {
                if (false && flux > 0.5)
                {
                    valveIsOpen=false;
                    M_fsi->setFluxBC(LifeV::BCh_monolithicFlux(valveIsOpen));
                    M_fsi->setup(/*data_file*/);
                    M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->FSIOper(), valveIsOpen));
                }
            }
            // close the valve
            else{
                if (M_fsi->FSIOper()->fluid().pressure(2, M_fsi->displacement()) < LifeV::LumpedHeart::M_pressure )
                {
                    valveIsOpen=true;
                    M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->FSIOper(), valveIsOpen));
                    M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->FSIOper(), valveIsOpen));
                }
            }

            int flag =2;
            FC0.renewParameters( *M_fsi, 3 );
            LH.renewParameters( *M_fsi->FSIOper(), flag, M_data->dataFluid()->dataTime()->getTime(), flux );

            //                 FC0.renewParameters( *M_fsi, 6 );
            //                 FC1.renewParameters( *M_fsi, 3, 4 );
            //                 FC2.renewParameters( *M_fsi, 3, 5 );
            //                 FC3.renewParameters( *M_fsi, 3, 6 );
            //                 FC4.renewParameters( *M_fsi, 3, 7 );

            boost::timer _timer;

            M_fsi->FSIOper()->getSolidDisp(*M_solidDisp);//    displacement(), M_offset);
            M_fsi->FSIOper()->getSolidVel(*M_solidVel);//    displacement(), M_offset);
//             *M_solidDisp *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataFluid()->dataTime()->getTimeStep());
//             *M_solidVel  *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataFluid()->dataTime()->getTimeStep());

            M_fsi->FSIOper()->getFluidVelAndPres(*M_velAndPressure);

            M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->getTime() );

            M_fsi->iterate();

            *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();

            M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->getTime() );


            M_fsi->FSIOper()->displayer().leaderPrint( "average inlet pressure  = ", M_fsi->FSIOper()->fluid().pressure(2, *M_velAndPressure));
            M_fsi->FSIOper()->displayer().leaderPrint( "average outlet pressure = ", M_fsi->FSIOper()->fluid().pressure(3, *M_velAndPressure));
            M_fsi->FSIOper()->displayer().leaderPrint( "inlet flux              = ", M_fsi->FSIOper()->fluid().flux(2, *M_velAndPressure));
            M_fsi->FSIOper()->displayer().leaderPrint( "outlet flux             = ", M_fsi->FSIOper()->fluid().flux(3, *M_velAndPressure));

            std::cout << "[fsi_run] Iteration " << _i << " was done in : "
                      << _timer.elapsed() << "\n";

            std::cout << "solution norm " << _i << " : "
                      << M_fsi->displacement().Norm2() << "\n";

            ///////// CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            try{
                if(!M_data->method().compare("monolithicGI"))
                 checkCEResult(M_data->dataFluid()->dataTime()->getTime());
             else
                 checkGCEResult(M_data->dataFluid()->dataTime()->getTime());
            }catch(Problem::RESULT_CHANGED_EXCEPTION){std::cout<<"res. changed"<<std::endl;}
            ///////// END OF CHECK
        }
        if(M_data->method().compare("monolithicGI"))
        {
            M_fsi->FSIOper()->iterateMesh(M_fsi->displacement());

            M_solidDisp->subset(M_fsi->displacement(), offset);
            M_solidVel->subset(M_fsi->FSIOper()->solid().vel(), offset);
//             *M_solidDisp *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataFluid()->dataTime()->getTimeStep());
//             *M_solidVel  *= 1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataFluid()->dataTime()->getTimeStep());


            *M_velAndPressure = M_fsi->displacement();
            M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->getTime() );
            *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
            M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->getTime() );
        }

        std::cout << "Total computation time = "
                  << _overall_timer.elapsed() << "s" << "\n";

    }

private:

    void initialize(std::string& loadInitSol,  GetPot const& data_file);

    void checkCEResult(const LifeV::Real& time);
    void checkGCEResult(const LifeV::Real& time);

    fsi_solver_ptr M_fsi;
    data_PtrType   M_data;

    filter_ptrtype M_exporterSolid;
    filter_ptrtype M_exporterFluid;
    filter_ptrtype M_importerSolid;
    filter_ptrtype M_importerFluid;
    vector_ptrtype M_velAndPressure;
    vector_ptrtype M_fluidDisp;
    vector_ptrtype M_solidDisp;
    vector_ptrtype M_solidVel;
    LifeV::FlowConditions FC0;
    LifeV::LumpedHeart LH;
    LifeV::Real    M_Tstart;

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
        data_file( _data_file ),
        oper     ( _data_file( "problem/method", "exactJacobian" ) ),
        prec     ( ( LifeV::Preconditioner )_data_file( "problem/precond", LifeV::NEUMANN_NEUMANN ) )
    {}

    FSIChecker( GetPot const& _data_file,
                std::string _oper,
                LifeV::Preconditioner _prec = LifeV::NO_PRECONDITIONER ):
        data_file( _data_file ),
        oper     ( _oper ),
        prec     ( ( LifeV::Preconditioner )_data_file( "problem/precond", LifeV::NEUMANN_NEUMANN ) )
    {}

    void operator()()
    {
        boost::shared_ptr<Problem> fsip;

        try
        {
            fsip = boost::shared_ptr<Problem>( new Problem( data_file, oper ) );

            fsip->fsiData()->setPreconditioner( prec );

            fsip->run();
        }
        catch ( std::exception const& _ex )
        {
            std::cout << "caught exception :  " << _ex.what() << "\n";
        }

        //@disp = fsip->fsiSolver()->FSIOper()->displacementOnInterface();
    }

    GetPot                data_file;
    std::string           oper;
    LifeV::Preconditioner prec;
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

void Problem::initialize(std::string& loadInitSol,  GetPot const& data_file)
{

    using namespace LifeV;
    std::string const importerType =  data_file( "importer/type", "ensight");
    std::string const fluidName    =  data_file( "importer/fluid/filename", "fluid");
    std::string const solidName    =  data_file( "importer/solid/filename", "solid");


#ifdef HAVE_HDF5
    if (importerType.compare("hdf5") == 0)
    {
        M_importerFluid.reset( new  hdf5filter_type( data_file, fluidName) );
        M_importerSolid.reset( new  hdf5filter_type ( data_file,solidName));// M_fsi->FSIOper()->solidMesh().mesh(), "solid", M_fsi->FSIOper()->dFESpace().map().Comm().MyPID()) );

    }
    else
#endif
    {
        if (importerType.compare("none") == 0)
        {
            M_importerFluid.reset( new NoExport<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->uFESpace().mesh(), "fluid", M_fsi->FSIOper()->uFESpace().map().Comm().MyPID()) );
            M_importerSolid.reset( new NoExport<RegionMesh3D<LinearTetra> > ( data_file, M_fsi->FSIOper()->dFESpace().mesh(), "solid", M_fsi->FSIOper()->uFESpace().map().Comm().MyPID()) );
        } else {
            M_importerFluid.reset( new  ensightfilter_type( data_file, fluidName) );
            M_importerSolid.reset( new  ensightfilter_type ( data_file, solidName) );
        }
    }

    M_importerFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());
    M_importerSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().Comm().MyPID());

    M_importerFluid->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                  UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );

    M_importerFluid->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
                                  UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()),
                                  UInt(M_fsi->FSIOper()->pFESpace().dof().numTotalDof()) );

    M_importerFluid->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
                                  UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );



    M_importerSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                                  UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
    M_importerSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
                                  UInt(0),
                                  M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );


    using namespace LifeV;
    typedef EpetraVector vector_type;

    std::string loadInitSolPrev(data_file("problem/initSolPrev","-1"));


    boost::shared_ptr<LifeV::EpetraVector> initSol(new LifeV::EpetraVector(*M_fsi->FSIOper()->getCouplingVariableMap()));
    boost::shared_ptr<LifeV::EpetraVector> initSolSVel(new LifeV::EpetraVector(*M_fsi->FSIOper()->getCouplingVariableMap()));
    boost::shared_ptr<LifeV::EpetraVector> UniqueV(new LifeV::EpetraVector(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique));
    boost::shared_ptr<LifeV::EpetraVector> UniqueVFD;
    boost::shared_ptr<LifeV::EpetraVector> UniqueVFDOld;


    UInt offset=dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->getOffset();

    Real dt= M_fsi->FSIOper()->dataFluid()->dataTime()->getTimeStep();//data_file("problem/Tstart"   ,0.);
    M_fsi->FSIOper()->displayer().leaderPrint( "Starting time = " ,M_Tstart);

    M_importerFluid->import(M_Tstart-dt, dt);
    M_importerSolid->import(M_Tstart-dt, dt);

    UniqueVFDOld.reset(new vector_type(*M_fluidDisp, Unique, Zero));
    dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->initializeMesh(UniqueVFDOld);

    M_importerFluid->import(M_Tstart);
    M_importerSolid->import(M_Tstart);


    UniqueV.reset( new vector_type(*M_velAndPressure, Unique, Zero));
    *initSol=*UniqueV;
    M_fsi->FSIOper()->fluid().initialize(*initSol);



    UniqueV.reset(new vector_type(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique, Zero));
    UniqueV->subset(*M_solidDisp, M_solidDisp->getMap(), (UInt)0, offset);
    *UniqueV*=1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataFluid()->dataTime()->getTimeStep());
    M_fsi->FSIOper()->solid().initialize(UniqueV);
    *initSol+=*UniqueV;

    if(!M_data->method().compare("monolithicGI"))
    {
        UniqueVFD.reset(new vector_type(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique, Zero));
        UniqueVFD->subset(*M_fluidDisp, M_fluidDisp->getMap(), (UInt)0, dynamic_cast<LifeV::MonolithicGI*>(M_fsi->FSIOper().get())->mapWithoutMesh().getMap(Unique)->NumGlobalElements());
        *initSol+=*UniqueVFD;
    }

    initSolSVel.reset(new vector_type(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique, Zero));
    initSolSVel->subset(*M_solidVel,M_solidVel->getMap(), (UInt)0, offset);
    *initSolSVel*=1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_data->dataSolid()->dataTime()->getTimeStep());
    M_fsi->FSIOper()->solid().initializeVel(*initSolSVel);
    M_fsi->initialize(initSol);
}

void Problem::checkGCEResult(const LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().Norm2();
    if(time==0.001 && (dispNorm-684898)     /dispNorm*(dispNorm-684898)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.002 && (dispNorm-854345)     /dispNorm*(dispNorm-850537)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.003 && (dispNorm-1.11118e+06)/dispNorm*(dispNorm-1.10523e+06)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.004 && (dispNorm-802296)     /dispNorm*(dispNorm-807697)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.005 && (dispNorm-869612)     /dispNorm*(dispNorm-869367)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.006 && (dispNorm-799188)     /dispNorm*(dispNorm-794390)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.007 && (dispNorm-795947)     /dispNorm*(dispNorm-794135)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.008 && (dispNorm-756083)     /dispNorm*(dispNorm-752333)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.009 && (dispNorm-765216)     /dispNorm*(dispNorm-762949)     /dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}


void Problem::checkCEResult(const LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().Norm2();
    if(time==0.001 && (dispNorm-615015)/dispNorm*(dispNorm-615015)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.002 && (dispNorm-1.00181e+06)/dispNorm*(dispNorm-1.00181e+06)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.003 && (dispNorm-1.01128e+06)/dispNorm*(dispNorm-1.01128e+06)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time);else
    if(time==0.004 && (dispNorm-644936)/dispNorm*(dispNorm-644936)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.005 && (dispNorm-652025)/dispNorm*(dispNorm-652025)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.006 && (dispNorm-555216)/dispNorm*(dispNorm-555216)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.007 && (dispNorm-538934)/dispNorm*(dispNorm-538934)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.008 && (dispNorm-520004)/dispNorm*(dispNorm-520004)/dispNorm>1e-3) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}
