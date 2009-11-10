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
 *  -# solution with exact Newton (useShapeDerivatives = true)
 *  -# solution with quasi Newton (useShapeDerivatives = false)
 *  -# preconditioner choice: see the classes Monolithic and fullMonolithic
 * - Monolithic (GCE):
 *  -# solution extrapolating the fluid domain
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


// namespace flowCond
// {
//  LifeV::FlowConditions FC0;
// }
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

    typedef LifeV::Ensight<LifeV::FSIOperator::mesh_type>  filter_type;
    typedef boost::shared_ptr<filter_type>                 filter_ptrtype;
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



        M_ensightMesh.reset( new  filter_type( data_file, "fluid") );

        M_ensightMesh->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());
        M_velAndPressure.reset( new vector_type( M_fsi->FSIOper()->fluid().getMap(), Repeated ));
        M_fluidDisp.reset     ( new vector_type( M_fsi->FSIOper()->meshMotion().getMap(), Repeated ));
        M_WSS.reset           ( new vector_type(  M_fsi->FSIOper()->uFESpace().map(), Repeated ));

        M_ensightMesh->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                    UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );

        M_ensightMesh->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
                                    UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()),
                                    UInt(M_fsi->FSIOper()->pFESpace().dof().numTotalDof()) );

        M_ensightMesh->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
                                    UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );

        M_ensightMesh->addVariable( ExporterData::Vector, "f-wss", M_WSS,
                                    UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );



        UInt offset=3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()+M_fsi->FSIOper()->pFESpace().dof().numTotalDof()+M_fsi->FSIOper()->BCh_flux()->size();

        M_ensightSolid.reset( new  filter_type ( data_file, "solid") );

        M_ensightSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().Comm().MyPID());

        M_solidDisp.reset( new vector_type( M_fsi->FSIOper()->solid().getMap(), Repeated ));
        M_solidVel.reset ( new vector_type( M_fsi->FSIOper()->solid().getMap(), Repeated ));
        M_ensightSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                                     UInt(offset), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
        M_ensightSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
                                     UInt(offset), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );



#ifdef HAVE_HDF5
        M_hdf5Mesh.reset( new  hdf5filter_type( data_file, "fluid") );

        M_hdf5Mesh->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());

        M_hdf5Mesh->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                 UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );

        M_hdf5Mesh->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
                                 UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()),
                                 UInt(M_fsi->FSIOper()->pFESpace().dof().numTotalDof()) );

        M_hdf5Mesh->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
                                 UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );

        M_hdf5Mesh->addVariable( ExporterData::Vector, "f-wss", M_WSS,
                                 UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );

        M_hdf5Solid.reset( new  hdf5filter_type ( data_file, "solid") );
        M_hdf5Solid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().Comm().MyPID());
        M_hdf5Solid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                                  UInt(offset), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
        M_hdf5Solid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
                                  UInt(offset), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
#endif


        M_Tstart = 0.;

		//to load using ensight
		std::string loadInitSol(data_file("problem/initSol","-1"));

		if(loadInitSol.compare("-1"))
            {
                initialize(loadInitSol, data_file);
            }
        else
            {
            }


        FC0.initParameters( *M_fsi->FSIOper(), 3);
    }

    fsi_solver_ptr fsiSolver() { return M_fsi; }


    /*!
      This routine runs the temporal loop
    */
    void
    run( double dt, double T)
    {
        boost::timer _overall_timer;
        int _i = 1;
        double time=M_Tstart + dt;

        dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->enableWssComputation(1);

        for (time=M_Tstart + dt; time <= T; time += dt, ++_i)
            {
                FC0.renewParameters( *M_fsi, 3 );
                //                 FC0.renewParameters( *M_fsi, 6 );
                //                 FC1.renewParameters( *M_fsi, 3, 4 );
                //                 FC2.renewParameters( *M_fsi, 3, 5 );
                //                 FC3.renewParameters( *M_fsi, 3, 6 );
                //                 FC4.renewParameters( *M_fsi, 3, 7 );

                boost::timer _timer;



                *M_solidDisp = M_fsi->displacement();
                *M_solidDisp *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();
                *M_solidVel = M_fsi->FSIOper()->solid().vel();
                *M_solidVel *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();
                *M_velAndPressure = M_fsi->displacement();
                M_ensightSolid->postProcess( time );
#ifdef HAVE_HDF5
                M_hdf5Solid->postProcess( time );
#endif

                M_fsi->iterate( time );


                *M_WSS= *(dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->/*WS());//*/computeWS());

                *M_fluidDisp      = M_fsi->FSIOper()->meshDisp();

                M_ensightMesh->postProcess( time );
#ifdef HAVE_HDF5
                M_hdf5Mesh->postProcess( time );
#endif

                //                    }

                std::cout << "[fsi_run] Iteration " << _i << " was done in : "
                          << _timer.elapsed() << "\n";

            std::cout << "solution norm " << _i << " : "
                      << M_fsi->displacement().Norm2() << "\n";

                ///////// CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            if(dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->isFullMonolithic())
                checkCEResult(time);
            else
                checkGCEResult(time);
                ///////// END OF CHECK
            }
        if(!M_fsi->FSIOper()->dataFluid().useShapeDerivatives())
            {
                M_fsi->FSIOper()->iterateMesh(M_fsi->displacement());
                *M_solidDisp = M_fsi->displacement();
                *M_solidDisp *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();
                *M_solidVel = M_fsi->FSIOper()->solid().vel();
                *M_solidVel *= M_fsi->timeStep()*M_fsi->FSIOper()->solid().rescaleFactor();
                *M_velAndPressure = M_fsi->displacement();
                M_ensightSolid->postProcess( time );
                *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
                M_ensightMesh->postProcess( time );
#ifdef HAVE_HDF5
                M_hdf5Solid->postProcess( time );
                M_hdf5Mesh->postProcess( time );
#endif
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

    filter_ptrtype M_ensightSolid;
    filter_ptrtype M_ensightMesh;
#ifdef HAVE_HDF5
    hdf5filter_ptrtype M_hdf5Solid;
    hdf5filter_ptrtype M_hdf5Mesh;
#endif
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
                std::cout << "calling problem constructor ... " << std::flush;
                fsip = boost::shared_ptr<Problem>( new Problem( data_file, oper ) );
                std::cout << "problem set" << std::endl;


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
    std::cout << "% using MPI" << std::endl;
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
    typedef EpetraVector vector_type;

    std::string loadInitSolPrev(data_file("problem/initSolPrev","-1"));


    boost::shared_ptr<LifeV::EpetraVector> initSol(new LifeV::EpetraVector(*M_fsi->FSIOper()->couplingVariableMap()));
    boost::shared_ptr<LifeV::EpetraVector> initSolSVel(new LifeV::EpetraVector(*M_fsi->FSIOper()->couplingVariableMap()));

    boost::shared_ptr<LifeV::EpetraVector> initSolV(new LifeV::EpetraVector(*M_fsi->FSIOper()->couplingVariableMap(), Repeated));
    boost::shared_ptr<LifeV::EpetraVector> initSolP(new LifeV::EpetraVector(*M_fsi->FSIOper()->couplingVariableMap(), Repeated));
    boost::shared_ptr<LifeV::EpetraVector> initSolS(new LifeV::EpetraVector(*M_fsi->FSIOper()->couplingVariableMap(), Repeated));
    boost::shared_ptr<LifeV::EpetraVector> initSolFD;
    boost::shared_ptr<LifeV::EpetraVector> initSolidV(new LifeV::EpetraVector(*M_fsi->FSIOper()->couplingVariableMap(), Repeated));
    boost::shared_ptr<LifeV::EpetraVector> initSolFDPrev;
    boost::shared_ptr<LifeV::EpetraVector> UniqueV(new LifeV::EpetraVector(*M_fsi->FSIOper()->couplingVariableMap(), Unique));
    boost::shared_ptr<LifeV::EpetraVector> UniqueVFD;
    boost::shared_ptr<LifeV::EpetraVector> UniqueVFDOld;

    boost::shared_ptr<LifeV::ExporterData> initSolFluidDisp;
    boost::shared_ptr<LifeV::ExporterData> initSolFluidDispPrev;

    UInt offset=3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()+M_fsi->FSIOper()->pFESpace().dof().numTotalDof() +M_fsi->FSIOper()->BCh_flux()->size();
    LifeV::ExporterData initSolFluidVel(LifeV::ExporterData::Vector, std::string("f-velocity."+loadInitSol), initSolV, UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof(), UInt(0) );
    LifeV::ExporterData initSolFluidPress(LifeV::ExporterData::Scalar, "f-pressure."+loadInitSol, initSolP,3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof(), M_fsi->FSIOper()->pFESpace().dof().numTotalDof(), UInt(0) );
    LifeV::ExporterData initSolSolidDisp(LifeV::ExporterData::Vector,"s-displacement."+loadInitSol, initSolS, offset, M_fsi->FSIOper()->dFESpace().dof().numTotalDof() , UInt(0) );
    LifeV::ExporterData initSolSolidVel(LifeV::ExporterData::Vector,"s-velocity."+loadInitSol, initSolidV, offset, M_fsi->FSIOper()->dFESpace().dof().numTotalDof() , UInt(0) );

    initSolFD.reset(new LifeV::EpetraVector(M_fsi->FSIOper()->mmFESpace().map(), Repeated));
    initSolFDPrev.reset(new LifeV::EpetraVector(M_fsi->FSIOper()->mmFESpace().map(), Repeated));
    initSolFluidDisp.reset(new LifeV::ExporterData(LifeV::ExporterData::Vector, "f-displacement."+loadInitSol, initSolFD, 0, M_fsi->FSIOper()->mmFESpace().dof().numTotalDof(), UInt(0)));
    initSolFluidDispPrev.reset(new LifeV::ExporterData(LifeV::ExporterData::Vector, "f-displacement."+loadInitSolPrev, initSolFDPrev, 0, M_fsi->FSIOper()->mmFESpace().dof().numTotalDof(), UInt(0)));
    //                        }
    //upcast necessary

    M_Tstart= M_fsi->FSIOper()->dataFluid().getInitialTime();//data_file("problem/Tstart"   ,0.);
    Real dt= M_fsi->FSIOper()->dataFluid().getTimeStep();//data_file("problem/Tstart"   ,0.);
    M_fsi->FSIOper()->displayer().leaderPrint( "Starting time = " ,M_Tstart);

    //M_ensightMesh->import(M_Tstart);
    //M_ensightSolid->import(M_Tstart);
    M_ensightMesh->M_rd_ascii(initSolFluidVel);
    M_ensightMesh->M_rd_ascii(initSolFluidPress);
    M_ensightSolid->M_rd_ascii(initSolSolidDisp);
    M_ensightSolid->M_rd_ascii(initSolSolidVel);
    M_ensightMesh->M_rd_ascii(*initSolFluidDisp);
    M_ensightMesh->M_rd_ascii(*initSolFluidDispPrev);

    //                     M_hdf5Mesh->M_rd_ascii(initSolFluidVel);
    //                     M_hdf5Mesh->M_rd_ascii(initSolFluidPress);
    //                     M_hdf5Solid->M_rd_ascii(initSolSolidDisp);
    //                     M_hdf5Solid->M_rd_ascii(initSolSolidVel);
    //                     M_hdf5Mesh->M_rd_ascii(*initSolFluidDisp);
    //                     M_hdf5Mesh->M_rd_ascii(*initSolFluidDispPrev);


    UniqueV.reset( new vector_type(*initSolV, Unique, Zero));
    *initSol=*UniqueV;
    UniqueV.reset( new vector_type(*initSolP, Unique, Zero));
    *initSol+=*UniqueV;
    M_fsi->FSIOper()->fluid().initialize(*initSol);

    UniqueV.reset(new vector_type(*initSolS, Unique, Zero));
    *UniqueV*=1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_fsi->timeStep());
    M_fsi->FSIOper()->solid().initialize(UniqueV);
    *initSol+=*UniqueV;
    UniqueVFD.reset(new vector_type(*initSolFD, Unique, Zero));

    if(dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->isFullMonolithic())
        initSol->add(*UniqueVFD, offset + 3*M_fsi->FSIOper()->dFESpace().dof().numTotalDof() + (dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get()))->dimInterface());

    UniqueVFDOld.reset(new vector_type(*initSolFDPrev, Unique, Zero));

    dynamic_cast<LifeV::Monolithic*>(M_fsi->FSIOper().get())->initializeMesh(UniqueVFDOld);

    initSolSVel.reset(new vector_type(*initSolidV, Unique, Zero));
    *initSolSVel*=1/(M_fsi->FSIOper()->solid().rescaleFactor()*M_fsi->timeStep());
    M_fsi->FSIOper()->solid().initializeVel(*initSolSVel);
    M_fsi->initialize(initSol);
}

void Problem::checkGCEResult(LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().Norm2();
    if(time==0.001 && (dispNorm-504035)/dispNorm*(dispNorm-504035)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.002 && (dispNorm-914879)/dispNorm*(dispNorm-914879)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.003 && (dispNorm-1.26589e+06)/dispNorm*(dispNorm-1.26589e+06)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.004 && (dispNorm-777523)/dispNorm*(dispNorm-777523)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.005 && (dispNorm-938905)/dispNorm*(dispNorm-938905)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.006 && (dispNorm-753279)/dispNorm*(dispNorm-753279)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.007 && (dispNorm-835605)/dispNorm*(dispNorm-835605)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.008 && (dispNorm-717966)/dispNorm*(dispNorm-717966)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.009 && (dispNorm-779217)/dispNorm*(dispNorm-779217)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}


void Problem::checkCEResult(LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().Norm2();
    if(time==0.001 && (dispNorm-472128)/dispNorm*(dispNorm-472128)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.002 && (dispNorm-913593)/dispNorm*(dispNorm-913593)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.003 && (dispNorm-1.28524e+06)/dispNorm*(dispNorm-1.28524e+06)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);else
    if(time==0.004 && (dispNorm-788936)/dispNorm*(dispNorm-788936)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.005 && (dispNorm-946694)/dispNorm*(dispNorm-946694)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.006 && (dispNorm-755566)/dispNorm*(dispNorm-755566)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.007 && (dispNorm-825823)/dispNorm*(dispNorm-825823)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.008 && (dispNorm-695446)/dispNorm*(dispNorm-695446)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time); else
    if(time==0.009 && (dispNorm-765094)/dispNorm*(dispNorm-765094)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}
