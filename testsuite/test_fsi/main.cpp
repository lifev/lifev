/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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

#ifdef TWODIM
#error test_fsi cannot be compiled in 2D
#endif

#include <boost/timer.hpp>

#include <cassert>
#include <iomanip>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/chrono.hpp>

#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSIOperator.hpp>
#include "life/lifesolver/exactJacobianBase.hpp"
#include "life/lifesolver/fixedPointBase.hpp"
#include <life/lifesolver/dataNavierStokes.hpp>

#include <life/lifefilters/ensight.hpp>

#include "Epetra_config.h"
#ifdef HAVE_MPI
	#include "Epetra_MpiComm.h"
#else
	#include "Epetra_SerialComm.h"
#endif

#include "ud_functions.hpp"
#include "boundaryConditions.hpp"

namespace LifeV
{
namespace
{
	static bool regIF = PRECFactory::instance().registerProduct( "Ifpack", &createIfpack );
	static bool regML = PRECFactory::instance().registerProduct( "ML", &createML );
	static bool regFP = FSIFactory::instance().registerProduct( "fixedPoint", &createFP );
	static bool regEJ = FSIFactory::instance().registerProduct( "exactJacobian", &createEJ );
}
}


#define PI 3.141592653589793
class bc_adaptor
{
public:

    bc_adaptor( LifeV::FSIOperator &_oper ):
            M_oper   ( _oper )
    {
    	LifeV::Real area0	= 0.78540;
    	//LifeV::Real area0	= M_oper.fluid().area(3);
        LifeV::Real area	= area0;

        LifeV::Real beta	= M_oper.solid().thickness()*M_oper.solid().young() /
							(1 - M_oper.solid().poisson()*M_oper.solid().poisson()) * PI/area0;

        LifeV::Real qn		= M_oper.fluid().flux(3);

    	M_outflow			= std::pow(std::sqrt(M_oper.solid().rho())/(2*std::sqrt(2))*qn/area + std::sqrt(beta*std::sqrt(area0)), 2)
							- beta*std::sqrt(area0);

        std::cout << "--------------- Absorbing boundary condition ---------------" << std::endl;
        std::cout << "  Outflow BC : density   = " << M_oper.solid().rho() << std::endl;
        std::cout << "  Outflow BC : thickness = " << M_oper.solid().thickness() << std::endl;
        std::cout << "  Outflow BC : young     = " << M_oper.solid().young() << std::endl;
        std::cout << "  Outflow BC : poisson   = " << M_oper.solid().poisson() << std::endl;
        std::cout << "  Outflow BC : area0     = " << area0 << std::endl;
        std::cout << "  Outflow BC : area      = " << M_oper.fluid().area(3) << std::endl;
        std::cout << "  Outflow BC : radius    = " << std::sqrt(area0/PI) << std::endl;
        std::cout << "  Outflow BC : beta      = " << beta << std::endl;
        std::cout << "  Outflow BC : Flow rate = " << qn << std::endl;
        std::cout << "  Outflow BC : outflow   = " << M_outflow << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
    }

    LifeV::Real operator()(LifeV::Real /*__time*/,
                           LifeV::Real /*__x*/,
                           LifeV::Real /*__y*/,
                           LifeV::Real /*__z*/,
                           LifeV::ID   __i)
    {
        switch(__i)
        {
            case 1:
                return 0.;
                break;
            case 2:
            	return 0.;
                break;
            case 3:
                //return 0.;
                return -M_outflow;
                break;
            default:
                ERROR_MSG("This entrie is not allowed: disp_adatptor");
                break;
        }

        return 0.;
    }
private:

    LifeV::FSIOperator&		M_oper;
    LifeV::Real				M_outflow;
};

class Problem
{
public:

	typedef boost::shared_ptr<LifeV::FSISolver>				fsi_solver_ptr;

    typedef LifeV::FSIOperator::vector_type					vector_type;
    typedef LifeV::FSIOperator::vector_ptrtype				vector_ptrtype;

    typedef LifeV::Ensight<LifeV::FSIOperator::mesh_type>	filter_type;
    typedef boost::shared_ptr<filter_type>					filter_ptrtype;

    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */
    Problem( GetPot const& data_file, std::string _oper = "" )
    {
    	using namespace LifeV;

    	Debug( 10000 ) << "creating FSISolver with operator :  " << _oper << "\n";
		//DataNavierStokes< RegionMesh3D_ALE<LinearTetra> > M_dataNS(data_file);

    	M_fsi = fsi_solver_ptr(  new FSISolver( data_file, _oper ) );
    	//M_fsi = fsi_solver_ptr(  new FSISolver( _oper ) );

    	Debug( 10000 ) << _oper << " set \n";
		//M_fsi->setup( data_file );

    	MPI_Barrier(MPI_COMM_WORLD);

		//M_fsi->setSourceTerms( fZero, fZero );
    	Debug( 10000 ) << "Setting up the BC \n";

    	M_fsi->setFluidBC(BCh_fluid(*M_fsi->FSIOper()));
    	M_fsi->setHarmonicExtensionBC (BCh_harmonicExtension(*M_fsi->FSIOper()));
    	M_fsi->setSolidBC(BCh_solid(*M_fsi->FSIOper()));
    	M_fsi->setLinFluidBC(BCh_fluidLin(*M_fsi->FSIOper()));
    	M_fsi->setLinSolidBC(BCh_solidLin(*M_fsi->FSIOper()));

    	M_absorbingBC = data_file("fluid/absorbing_bc", false);
    	std::cout << "   absorbing BC = " << M_absorbingBC << std::endl;

    	Debug( 10000 ) << "BC set\n";

    	//Debug( 10000 ) << "Setting up the FSISolver \n";

    	bool restart = data_file("problem/restart",false);
    	M_Tstart = 0.;

    	//M_fsi->resetFSISolvers();

    	if (M_fsi->isFluid())
    	{
    		M_ensightFluid.reset( new  filter_type( data_file, "fixedPtFluid") );

    		M_ensightFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());

    		M_velAndPressure.reset( new vector_type( M_fsi->FSIOper()->fluid().getMap(), Repeated ));
    		M_fluidDisp.reset     ( new vector_type( M_fsi->FSIOper()->meshMotion().getMap(), Repeated ));
    		M_ensightFluid->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
										UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );

    		M_ensightFluid->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
										UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()),
										UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()+M_fsi->FSIOper()->pFESpace().dof().numTotalDof()) );

    		M_ensightFluid->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
										UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );
    	}

    	if (M_fsi->isSolid())
    	{
    		M_ensightSolid.reset( new  filter_type ( data_file, "fixedPtSolid") );

    		//assert( M_fsi->FSIOper()->uFESpace().get() );
    		//assert( M_fsi->FSIOper()->uFESpace().mesh().get() );

    		M_ensightSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().Comm().MyPID());

    		M_solidDisp.reset( new vector_type( M_fsi->FSIOper()->solid().getMap(), Repeated ));
    		M_solidVel.reset( new vector_type( M_fsi->FSIOper()->solid().getMap(), Repeated ));
    		M_ensightSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
										UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
    		M_ensightSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
										UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
    	}

    	if (restart)
    	{
    		std::string velFName  = data_file("fluid/miscellaneous/velname"  ,"vel");
    		std::string pressName = data_file("fluid/miscellaneous/pressname","press");
    		std::string velwName  = data_file("fluid/miscellaneous/velwname", "velw");
    		std::string depName   = data_file("solid/miscellaneous/depname"  ,"dep");
    		std::string velSName  = data_file("solid/miscellaneous/velname"  ,"velw");
    		M_Tstart= data_file("problem/Tstart"   ,0.);
    		std::cout << "Starting time = " << M_Tstart << std::endl;

    		//M_fsi->initialize(velFName, pressName, velwName, depName, velSName, M_Tstart);

    		if (M_fsi->isFluid())
    		{
    			M_ensightFluid->import(M_Tstart, M_fsi->timeStep());
    			M_fsi->FSIOper()->initializeFluid( *M_velAndPressure, *M_fluidDisp );
    		}
    		if (M_fsi->isSolid())
    		{
    			M_ensightSolid->import(M_Tstart, M_fsi->timeStep());
				M_fsi->FSIOper()->initializeSolid( *M_solidDisp, *M_solidVel );
    		}

    	}
    	//std::cout << "in problem" << std::endl;
    	//M_fsi->FSIOper()->fluid().postProcess();
	}

    fsi_solver_ptr fsiSolver() { return M_fsi; }

    /*!
      This routine runs the temporal loop
    */
    void
    run( double dt, double T)
	{
    	std::ofstream ofile;

    	bool const isFluidLeader( M_fsi->FSIOper()->isFluid() && M_fsi->FSIOper()->isLeader() );
    	if (isFluidLeader) ofile.open("fluxes.res");

    	boost::timer _overall_timer;

    	double flux;
    	int _i = 1;

    	//double time=M_Tstart;

		//if ( M_fsi->isFluid() )
    	//{
			//if (isFluidLeader) ofile << time << " ";
			//flux = M_fsi->FSIOper()->fluid().flux(2);
			//if (isFluidLeader) ofile << flux << " ";
			//flux = M_fsi->FSIOper()->fluid().flux(3);
			//if (isFluidLeader) ofile << flux << " ";
			//flux = M_fsi->FSIOper()->fluid().pressure(2);
			//if (isFluidLeader) ofile << flux << " ";
			//flux = M_fsi->FSIOper()->fluid().pressure(3);
			//if (isFluidLeader) ofile << flux << " " << std::endl;

			//*M_velAndPressure = M_fsi->FSIOper()->fluid().solution();
			//*M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();

			//M_ensightFluid->postProcess( time );
		//}

    	for (double time = M_Tstart + dt; time <= T; time += dt, ++_i)
    	{
    		boost::timer _timer;

    		if (M_absorbingBC && M_fsi->isFluid())
			{
				LifeV::BCFunctionBase outFlow;
				outFlow.setFunction(bc_adaptor(*M_fsi->FSIOper()));
				M_fsi->FSIOper()->BCh_fluid()->modifyBC(3, outFlow);
			}

    		M_fsi->iterate( time );

    		if ( M_fsi->isFluid() )
            {
				if (isFluidLeader) ofile << time << " ";
    			flux = M_fsi->FSIOper()->fluid().flux(2);
				if (isFluidLeader) ofile << flux << " ";
				flux = M_fsi->FSIOper()->fluid().flux(3);
				if (isFluidLeader) ofile << flux << " ";
				flux = M_fsi->FSIOper()->fluid().pressure(2);
				if (isFluidLeader) ofile << flux << " ";
				flux = M_fsi->FSIOper()->fluid().pressure(3);
				if (isFluidLeader) ofile << flux << " " << std::endl;

				*M_velAndPressure = M_fsi->FSIOper()->fluid().solution();
				*M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
				M_ensightFluid->postProcess( time );
            }

    		if ( M_fsi->isSolid() )
    		{
    			*M_solidDisp = M_fsi->FSIOper()->solid().disp();
    			*M_solidVel = M_fsi->FSIOper()->solid().vel();
    			M_ensightSolid->postProcess( time );
    		}

				std::cout << "[fsi_run] Iteration " << _i << " was done in : " << _timer.elapsed() << "\n";
    	}

		std::cout << "Total computation time = " << _overall_timer.elapsed() << "s" << "\n";
		ofile.close();
	}

private:

	fsi_solver_ptr M_fsi;
	double         M_Tstart;
    //LifeV::DataNavierStokes< LifeV::RegionMesh3D<LifeV::LinearTetra> > M_dataNS;

	MPI_Comm*      M_comm;

	bool           M_absorbingBC;

	filter_ptrtype M_ensightFluid;
	filter_ptrtype M_ensightSolid;

	vector_ptrtype M_velAndPressure;
	vector_ptrtype M_fluidDisp;
	vector_ptrtype M_solidDisp;
	vector_ptrtype M_solidVel;
};

struct FSIChecker
{
    FSIChecker( GetPot const& _data_file ):
        data_file( _data_file ),
        Oper     ( _data_file( "problem/method", "exactJacobian" ) ),
        prec     ( ( LifeV::Preconditioner )_data_file( "problem/precond", LifeV::NEUMANN_NEUMANN ) )
        {}

    FSIChecker( GetPot const& _data_file,
                std::string _oper):
        data_file( _data_file ),
        Oper     ( _oper ),
        prec     ( ( LifeV::Preconditioner )_data_file( "problem/precond", LifeV::NEUMANN_NEUMANN ) )
        {}

    void operator()()
	{
    	boost::shared_ptr<Problem> fsip;

    	try
    	{
    		std::cout << "calling problem constructor ... " << std::flush;
    		fsip = boost::shared_ptr<Problem>( new Problem( data_file, Oper ) );
    		std::cout << "problem set" << std::endl;

    		//fsip->fsiSolver()->FSIOperator()->setDataFromGetPot( data_file );
    		//std::cout << "in operator 1" << std::endl;
    		//fsip->fsiSolver()->FSIOper()->fluid().postProcess();
			//fsip->fsiSolver()->FSIOper()->solid().postProcess();

    		fsip->fsiSolver()->FSIOper()->setPreconditioner( prec );

			//std::cout << "in operator 2" << std::endl;
			//fsip->fsiSolver()->FSIOper()->fluid().postProcess();
    		fsip->run( fsip->fsiSolver()->timeStep(), fsip->fsiSolver()->timeEnd() );
    	}
    	catch ( std::exception const& _ex )
    	{
    		std::cout << "caught exception :  " << _ex.what() << "\n";
    	}

    	//@disp = fsip->fsiSolver()->FSIOper()->displacementOnInterface();
	}

	GetPot                data_file;
	std::string           Oper;
	LifeV::Preconditioner prec;
	LifeV::Vector         disp;
};

int main(int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    std::cout << "% using MPI" << std::endl;
#else
    std::cout << "% using serial Version" << std::endl;
#endif

    LifeV::Chrono chrono;
    chrono.start();

    GetPot command_line(argc,argv);
    string data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);

    const bool check = command_line.search(2, "-c", "--check");

    if (check)
    {
        LifeV::Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        FSIChecker _ej_check( data_file, "exactJacobian" );
        _ej_check();

        LifeV::Debug( 10000 ) << "_ej_disp size : "  << static_cast<double> (_ej_check.disp.size()) << "\n";
        LifeV::Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        LifeV::Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        FSIChecker _sp_check( data_file, "steklovPoincare" );
        _sp_check();


        LifeV::Debug( 10000 ) << "_fp_disp size : "  << static_cast<double> (_sp_check.disp.size()) << "\n";
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

    std::cout << "Total sum up " << chrono.diff_cumul() << " s." << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
