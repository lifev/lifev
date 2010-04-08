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
 in the files Monolithic.hpp and fullMonolithic.hpp . The geometry is that of a tube (benchmark test introduced in \ref GV03).
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

#include <life/lifesolver/dataNavierStokes.hpp>

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


// LifeV::FlowConditions FC0;
// LifeV::FlowConditions FC1;
// LifeV::FlowConditions FC2;
// LifeV::FlowConditions FC3;
// LifeV::FlowConditions FC4;

LifeV::FSIOperator* createFM(){ return new LifeV::fullMonolithic(); }
LifeV::FSIOperator* createM(){ return new LifeV::Monolithic(); }

class Problem
{
public:
    typedef boost::shared_ptr<LifeV::FSISolver> fsi_solver_ptr;

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
        if(!_oper.compare("monolithic"))
            FSIFactory::instance().registerProduct( "monolithic", &createM );
        else if(!_oper.compare("fullMonolithic"))
            {
                FSIFactory::instance().registerProduct( "fullMonolithic", &createFM );
            }

#ifdef DEBUG
			Debug( 10000 ) << "creating FSISolver with operator :  " << method << "\n";
#endif
			M_fsi = fsi_solver_ptr( new FSISolver( _oper ) );

			MPI_Barrier( MPI_COMM_WORLD );

#ifdef DEBUG
			Debug( 10000 ) << "Setting up data from GetPot \n";
#endif
			M_fsi->setDataFromGetPot( data_file );

			MPI_Barrier( MPI_COMM_WORLD );

            MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
            Debug( 10000 ) << "Setting up the BC \n";
#endif
	    M_fsi->setFluxBC(BCh_monolithicFlux());
	    M_fsi->setup(/*data_file*/);
	    M_fsi->setRobinBC(BCh_monolithicRobin(*M_fsi->FSIOper()));
	    M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->FSIOper()));
        M_fsi->setHarmonicExtensionBC (BCh_harmonicExtension(*M_fsi->FSIOper()));
        M_fsi->setSolidBC(BCh_monolithicSolid(*M_fsi->FSIOper()));
#ifdef DEBUG
        Debug( 10000 ) << "BC set\n";
#endif

        std::string const exporterType =  data_file( "exporter/type", "ensight");
        std::string const fluidName    =  data_file( "exporter/fluid/filename", "fluid");
        std::string const solidName    =  data_file( "exporter/solid/filename", "solid");

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

        M_WSS.reset      ( new vector_type( M_fsi->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
        M_exporterSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                                      UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
         M_exporterSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
                                       UInt(0),
                                       M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );

        M_exporterSolid->addVariable( ExporterData::Vector, "s-wss", M_WSS,
                                      UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );


        M_Tstart = data_file("fluid/time_discretization/initialtime"   ,0.);

		//to load using ensight/hdf5
		std::string loadInitSol(data_file("problem/initSol","-1"));

		if(loadInitSol.compare("-1"))
            {
                initialize(loadInitSol, data_file);
            }

        FC0.initParameters( *M_fsi->FSIOper(), 3);
    }

    fsi_solver_ptr fsiSolver() { return M_fsi; }


    /*!
      This routine runs the temporal loop
    */
    void
    run( LifeV::Real dt, LifeV::Real T)
    {
        boost::timer _overall_timer;
        int _i = 1;
        LifeV::Real time=M_Tstart + dt;
        LifeV::UInt offset=3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()+M_fsi->FSIOper()->pFESpace().dof().numTotalDof()+M_fsi->FSIOper()->BCh_flux()->size();

	    dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->enableWssComputation(1);
#ifdef HAVE_HDF5
        M_exporterFluid->postProcess( 0 );//ugly way to avoid that hdf5 starts with a deformed mesh
#endif

        for ( ; time <= T; time += dt, ++_i)
            {
                FC0.renewParameters( *M_fsi, 3 );
                //                 FC0.renewParameters( *M_fsi, 6 );
                //                 FC1.renewParameters( *M_fsi, 3, 4 );
                //                 FC2.renewParameters( *M_fsi, 3, 5 );
                //                 FC3.renewParameters( *M_fsi, 3, 6 );
                //                 FC4.renewParameters( *M_fsi, 3, 7 );

                boost::timer _timer;



                M_solidDisp->subset(M_fsi->displacement(), offset);
                *M_solidDisp *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();
                M_solidVel->subset(M_fsi->FSIOper()->solid().vel(), offset);
                *M_solidVel *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();

                *M_velAndPressure = M_fsi->displacement();
                *M_WSS= *(dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->/*WS());//*/computeWS());

                M_exporterSolid->postProcess( time );

                M_fsi->iterate( time );

                dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->computeMaxSingularValue();

                *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();

                M_exporterFluid->postProcess( time );

                //                    }

                std::cout << "[fsi_run] Iteration " << _i << " was done in : "
                          << _timer.elapsed() << "\n";

            std::cout << "solution norm " << _i << " : "
                      << M_fsi->displacement().Norm2() << "\n";

                ///////// CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            //try{
            if(dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->isFullMonolithic())
                checkCEResult(time);
            else
                checkGCEResult(time);
            //}catch(Problem::RESULT_CHANGED_EXCEPTION){std::cout<<"res. changed"<<std::endl;}
                ///////// END OF CHECK
            }
        if(!M_fsi->FSIOper()->dataFluid().useShapeDerivatives())
            {
                M_fsi->FSIOper()->iterateMesh(M_fsi->displacement());

                M_solidDisp->subset(M_fsi->displacement(), offset);
                *M_solidDisp *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();
                M_solidVel->subset(M_fsi->FSIOper()->solid().vel(), offset);
                *M_solidVel *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();

                *M_velAndPressure = M_fsi->displacement();
                M_exporterSolid->postProcess( time );
                *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
                M_exporterFluid->postProcess( time );
            }

        std::cout << "Total computation time = "
                  << _overall_timer.elapsed() << "s" << "\n";

    }

private:

    void initialize(std::string& loadInitSol,  GetPot const& data_file);

    void checkCEResult(LifeV::Real& time);
    void checkGCEResult(LifeV::Real& time);

    fsi_solver_ptr M_fsi;
    double         M_Tstart;

    MPI_Comm*      M_comm;

    filter_ptrtype M_exporterSolid;
    filter_ptrtype M_exporterFluid;
    filter_ptrtype M_importerSolid;
    filter_ptrtype M_importerFluid;
    vector_ptrtype M_velAndPressure;
    vector_ptrtype M_fluidDisp;
    vector_ptrtype M_solidDisp;
    vector_ptrtype M_solidVel;
    vector_ptrtype M_WSS;
    LifeV::FlowConditions FC0;
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

                fsip->fsiSolver()->FSIOper()->setPreconditioner( prec );

                fsip->run( fsip->fsiSolver()->timeStep(), fsip->fsiSolver()->timeEnd() );
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

        M_importerSolid->addVariable( ExporterData::Vector, "s-wss", M_WSS,
                                      UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );



    using namespace LifeV;
    typedef EpetraVector vector_type;

    std::string loadInitSolPrev(data_file("problem/initSolPrev","-1"));


    boost::shared_ptr<LifeV::EpetraVector> initSol(new LifeV::EpetraVector(*M_fsi->FSIOper()->getCouplingVariableMap()));
    boost::shared_ptr<LifeV::EpetraVector> initSolSVel(new LifeV::EpetraVector(*M_fsi->FSIOper()->getCouplingVariableMap()));
    boost::shared_ptr<LifeV::EpetraVector> UniqueV(new LifeV::EpetraVector(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique));
    boost::shared_ptr<LifeV::EpetraVector> UniqueVFD;
    boost::shared_ptr<LifeV::EpetraVector> UniqueVFDOld;


         UInt offset=3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()+M_fsi->FSIOper()->pFESpace().dof().numTotalDof() +M_fsi->FSIOper()->BCh_flux()->size();
    //                        }
    //upcast necessary

    Real init= data_file("problem/Tstart"   ,0.);
    Real dt= M_fsi->FSIOper()->dataFluid().dataTime()->getTimeStep();//data_file("problem/Tstart"   ,0.);
    M_fsi->FSIOper()->displayer().leaderPrint( "Starting time = " ,init);

    M_importerFluid->import(init-dt, dt);
    M_importerSolid->import(init-dt, dt);

    UniqueVFDOld.reset(new vector_type(*M_fluidDisp, Unique, Zero));
    dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->initializeMesh(UniqueVFDOld);

    M_importerFluid->import(init);
    M_importerSolid->import(init);


    UniqueV.reset( new vector_type(*M_velAndPressure, Unique, Zero));
    *initSol=*UniqueV;
    M_fsi->FSIOper()->fluid().initialize(*initSol);



    UniqueV.reset(new vector_type(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique, Zero));
    UniqueV->subset(*M_solidDisp, M_solidDisp->getMap(), (UInt)0, offset);
    *UniqueV*=1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_fsi->timeStep());
    M_fsi->FSIOper()->solid().initialize(UniqueV);
    *initSol+=*UniqueV;

    if(dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->isFullMonolithic())
    {
        UniqueVFD.reset(new vector_type(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique, Zero));
        UniqueVFD->subset(*M_fluidDisp, M_fluidDisp->getMap(), (UInt)0, dynamic_cast<LifeV::fullMonolithic*>(M_fsi->FSIOper().get())->mapWithoutMesh().getMap(Unique)->NumGlobalElements());
        *initSol+=*UniqueVFD;
    }


    initSolSVel.reset(new vector_type(*M_fsi->FSIOper()->getCouplingVariableMap(), Unique, Zero));
    initSolSVel->subset(*M_solidVel,M_solidVel->getMap(), (UInt)0, offset);
    *initSolSVel*=1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_fsi->timeStep());
    M_fsi->FSIOper()->solid().initializeVel(*initSolSVel);
    M_fsi->initialize(initSol);
}

void Problem::checkGCEResult(LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().Norm2();
    if(time==0.001 && (dispNorm-684898)/dispNorm*(dispNorm-504035)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.002 && (dispNorm-850537)/dispNorm*(dispNorm-914879)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.003 && (dispNorm-1.10523e+06)/dispNorm*(dispNorm-1.26589e+06)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.004 && (dispNorm-807697)/dispNorm*(dispNorm-777523)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.005 && (dispNorm-869367)/dispNorm*(dispNorm-938905)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.006 && (dispNorm-794390)/dispNorm*(dispNorm-753279)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.007 && (dispNorm-794135)/dispNorm*(dispNorm-835605)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.008 && (dispNorm-752333)/dispNorm*(dispNorm-717966)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.009 && (dispNorm-762949)/dispNorm*(dispNorm-779217)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}


void Problem::checkCEResult(LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().Norm2();
    if(time==0.001 && (dispNorm-615015)/dispNorm*(dispNorm-472128)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.002 && (dispNorm-787299)/dispNorm*(dispNorm-913593)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.003 && (dispNorm-773835)/dispNorm*(dispNorm-1.28524e+06)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);else
    if(time==0.004 && (dispNorm-654622)/dispNorm*(dispNorm-788936)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.005 && (dispNorm-543915)/dispNorm*(dispNorm-946694)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.006 && (dispNorm-517692)/dispNorm*(dispNorm-755566)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.007 && (dispNorm-497380)/dispNorm*(dispNorm-825823)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.008 && (dispNorm-486165)/dispNorm*(dispNorm-695446)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.009 && (dispNorm-478644)/dispNorm*(dispNorm-765094)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}
