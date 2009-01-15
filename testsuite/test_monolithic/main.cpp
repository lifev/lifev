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
#include <cassert>

#include <life/lifecore/life.hpp>

#include <boost/timer.hpp>

#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSIOperator.hpp>
//#include <life/lifesolver/Monolithic.hpp>
#include <life/lifesolver/fixedPointBase.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefilters/ensight.hpp>

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





//LifeV::FSIOperator* createM(){ return new LifeV::Monolithic(); }
class Problem
{
public:
    typedef boost::shared_ptr<LifeV::FSISolver> fsi_solver_ptr;

    typedef LifeV::FSIOperator::vector_type        vector_type;
    typedef LifeV::FSIOperator::vector_ptrtype     vector_ptrtype;

    typedef LifeV::Ensight<LifeV::FSIOperator::mesh_type>  filter_type;
    typedef boost::shared_ptr<filter_type>                 filter_ptrtype;
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
            //            bool reg = FSIFactory::instance().registerProduct( "monolithic", &createM );

            Debug( 10000 ) << "creating FSISolver with operator :  " << _oper << "\n";
//            DataNavierStokes< RegionMesh3D_ALE<LinearTetra> > M_dataNS(data_file);

            M_fsi = fsi_solver_ptr(  new FSISolver( data_file, _oper ) );
            Debug( 10000 ) << _oper << " set \n";

            MPI_Barrier(MPI_COMM_WORLD);

//            M_fsi->setSourceTerms( fZero, fZero );

            Debug( 10000 ) << "Setting up the BC \n";
            M_fsi->setFluidBC(BCh_monolithicFluid(*M_fsi->operFSI()));
            M_fsi->setHarmonicExtensionBC (BCh_harmonicExtension(*M_fsi->operFSI()));
            M_fsi->setSolidBC(BCh_monolithicSolid(*M_fsi->operFSI()));
            //            M_fsi->setLinFluidBC(BCh_fluidLin(*M_fsi->operFSI()));
            //            M_fsi->setLinSolidBC(BCh_solidLin(*M_fsi->operFSI()));
            Debug( 10000 ) << "BC set\n";


            int restart = data_file("problem/restart",0);
            M_Tstart = 0.;

            if (restart)
            {
                std::string velFName  = data_file("fluid/miscellaneous/velname"  ,"vel");
                std::string pressName = data_file("fluid/miscellaneous/pressname","press");
                std::string velwName  = data_file("fluid/miscellaneous/velwname", "velw");
                std::string depName   = data_file("solid/miscellaneous/depname"  ,"dep");
                std::string velSName  = data_file("solid/miscellaneous/velname"  ,"velw");
                M_Tstart= data_file("problem/Tstart"   ,0.);
                std::cout << "Starting time = " << M_Tstart << std::endl;
                M_fsi->initialize(velFName, pressName, velwName, depName, velSName, M_Tstart);
            }
            else
            {
//                 M_fsi->initialize( u0, p0, d0, w0 );
            }

            if (M_fsi->isFluid())
                {
                    M_ensightFluid.reset( new  filter_type( data_file, "fixedPtFluid") );
                    M_ensightMesh.reset( new  filter_type( data_file, "fixedPtMesh") );

                    //assert( M_fsi->FSIOper()->uFESpace().get() );
                    //assert( M_fsi->FSIOper()->uFESpace().mesh().get() );

                    M_ensightFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());
                    M_ensightMesh->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());

                    M_velAndPressure.reset( new vector_type( M_fsi->FSIOper()->fluid().getMap(), Repeated ));
                    M_fluidDisp.reset     ( new vector_type( M_fsi->FSIOper()->meshMotion().getMap(), Repeated ));
                    M_ensightMesh->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                                 UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );

                    M_ensightMesh->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
                                                 UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()),
                                                 UInt(M_fsi->FSIOper()->pFESpace().dof().numTotalDof()) );

                    M_ensightMesh->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
                                                 UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );


                }
            if (M_fsi->isSolid())
                {

                    UInt offset=3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()+M_fsi->FSIOper()->pFESpace().dof().numTotalDof();

                    M_ensightSolid.reset( new  filter_type ( data_file, "fixedPtSolid") );

                    //assert( M_fsi->FSIOper()->uFESpace().get() );
                    //assert( M_fsi->FSIOper()->uFESpace().mesh().get() );

                    M_ensightSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().Comm().MyPID());

                    M_solidDisp.reset( new vector_type( M_fsi->FSIOper()->solid().getMap(), Repeated ));
                    M_solidVel.reset( new vector_type( M_fsi->FSIOper()->solid().getMap(), Repeated ));
                    M_ensightSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                                                 offset, M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
                    M_ensightSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
                                                                     UInt(offset), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );

                }
//            std::cout << "in problem" << std::endl;
//            M_fsi->FSIOper()->fluid().postProcess();
        }

    fsi_solver_ptr fsiSolver() { return M_fsi; }


    /*!
      This routine runs the temporal loop
    */
    void
    run( double dt, double T)
        {
//            std::ofstream ofile( "fluxes.res" );

            boost::timer _overall_timer;

            if (M_Tstart != 0.) M_Tstart -= dt;

            int _i = 1;
            double time=M_Tstart + dt;

            for (time=M_Tstart + dt; time <= T; time += dt, ++_i)
            {
                boost::timer _timer;

                *M_velAndPressure = M_fsi->displacement();


                *M_solidDisp = M_fsi->displacement();
                *M_solidDisp *= M_fsi->timeStep();
                *M_solidVel = M_fsi->FSIOper()->solid().vel();
                *M_solidVel *= M_fsi->timeStep();
                *M_velAndPressure = M_fsi->displacement();
                M_ensightSolid->postProcess( time );
                M_ensightFluid->postProcess( time );


                M_fsi->iterate( time );

//                 ofile << time << " ";
//                 ofile << M_fsi->operFSI()->fluid().flux(2) << " ";
//                 ofile << M_fsi->operFSI()->fluid().flux(3) << " ";
//                 ofile << std::endl;
                        *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();

                        M_ensightMesh->postProcess( time );

                        //                    }

                std::cout << "[fsi_run] Iteration " << _i << " was done in : "
                          << _timer.elapsed() << "\n";
            }

            *M_solidDisp = M_fsi->displacement();
            *M_solidDisp *= M_fsi->timeStep();
            *M_solidVel = M_fsi->FSIOper()->solid().vel();
            *M_solidVel *= M_fsi->timeStep();
            *M_velAndPressure = M_fsi->displacement();
            M_ensightSolid->postProcess( time );
            M_ensightFluid->postProcess( time );

            M_fsi->FSIOper()->iterateMesh(M_fsi->displacement());
            *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
            M_ensightMesh->postProcess( time );

            std::cout << "Total computation time = "
                      << _overall_timer.elapsed() << "s" << "\n";

//            ofile.close();
        }

private:

    fsi_solver_ptr M_fsi;
    double         M_Tstart;
//    LifeV::DataNavierStokes< LifeV::RegionMesh3D<LifeV::LinearTetra> > M_dataNS;

    MPI_Comm*      M_comm;

    filter_ptrtype M_ensightFluid;
    filter_ptrtype M_ensightSolid;
    filter_ptrtype M_ensightMesh;

    vector_ptrtype M_velAndPressure;
    vector_ptrtype M_fluidDisp;
    vector_ptrtype M_solidDisp;
    vector_ptrtype M_solidVel;

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

//                fsip->fsiSolver()->FSIOperator()->setDataFromGetPot( data_file );
//                 std::cout << "in operator 1" << std::endl;
//                 fsip->fsiSolver()->operFSI()->fluid().postProcess();
//                 fsip->fsiSolver()->operFSI()->solid().postProcess();

                fsip->fsiSolver()->operFSI()->setPreconditioner( prec );

//                std::cout << "in operator 2" << std::endl;
//                fsip->fsiSolver()->operFSI()->fluid().postProcess();
                fsip->run( fsip->fsiSolver()->timeStep(), fsip->fsiSolver()->timeEnd() );
            }
            catch ( std::exception const& _ex )
            {
                std::cout << "caught exception :  " << _ex.what() << "\n";
            }

            //@disp = fsip->fsiSolver()->operFSI()->displacementOnInterface();
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
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);


    const bool check = command_line.search(2, "-c", "--check");

    if (check)
    {
        LifeV::Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        FSIChecker _ej_check( data_file, "exactJacobian" );
        _ej_check();

        LifeV::Debug( 10000 ) << "_ej_disp size : "  << _ej_check.disp.size() << "\n";
        LifeV::Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        LifeV::Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        FSIChecker _sp_check( data_file, "steklovPoincare" );
        _sp_check();


        LifeV::Debug( 10000 ) << "_fp_disp size : "  << _sp_check.disp.size() << "\n";
        LifeV::Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        double norm1 = LifeV::norm_2( _ej_check.disp - _sp_check.disp );

        std::cout << "norm_2(EJ displacement)          = " << LifeV::norm_2( _ej_check.disp ) << " \n"
                  << "norm_2(SP displacement)          = " << LifeV::norm_2( _sp_check.disp ) << " \n"
                  << "norm_2(displacement error EJ/SP) = " << norm1 << "\n";

#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        if (norm1 < 1e-04) return 0;
        else return -1;
    }
    else
    {
        FSIChecker _sp_check( data_file );
        _sp_check();
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;

}

