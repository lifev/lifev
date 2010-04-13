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

//#include <life/lifealg/IfpackPreconditioner.hpp>
//#include <life/lifealg/MLPreconditioner.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/chrono.hpp>

#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSIOperator.hpp>
//#include "life/lifesolver/exactJacobianBase.hpp"
//#include "life/lifesolver/fixedPointBase.hpp"
#include <life/lifesolver/dataNavierStokes.hpp>

#ifdef HAVE_HDF5
	#include <life/lifefilters/hdf5exporter.hpp>
#else
	#include <life/lifefilters/ensight.hpp>	
#endif

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
	#include <Epetra_MpiComm.h>
#else
	#include <Epetra_SerialComm.h>
#endif

#include "ud_functions.hpp"
#include "boundaryConditions.hpp"

/*
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
*/
using namespace LifeV;




#define PI 3.141592653589793
class bc_adaptor
{
public:

    bc_adaptor( FSIOperator& Operator ):
            M_oper   ( Operator )
    {
    	Real area0	= 0.7854;
    	//Real area0	= M_oper.fluid().area(3);
        Real area	= area0;

        Real beta	= M_oper.solid().thickness()*M_oper.solid().young() /
							(1 - M_oper.solid().poisson()*M_oper.solid().poisson()) * PI/area0;

        Real qn		= M_oper.fluid().flux(3);

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

    Real operator()( Real /*t*/, Real /*x*/, Real /*y*/, Real /*z*/, ID id)
    {
        switch( id )
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

    FSIOperator& M_oper;
    Real         M_outflow;
};





class Problem
{
public:

	typedef boost::shared_ptr<FSISolver>                    fsi_solver_ptr;

    typedef FSIOperator::vector_type                        vector_type;
    typedef FSIOperator::vector_ptrtype                     vector_ptrtype;

#ifdef HAVE_HDF5
    typedef Hdf5exporter<FSIOperator::mesh_type>            filter_type;
#else
    typedef Ensight<FSIOperator::mesh_type>                 filter_type;
#endif

    typedef boost::shared_ptr<filter_type>                  filter_ptrtype;

    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */
    Problem( const GetPot& dataFile, std::string method = "" )
    {
    	Debug( 10000 ) << "creating FSISolver with operator :  " << method << "\n";
    	M_fsi = fsi_solver_ptr( new FSISolver( method ) );

    	MPI_Barrier( MPI_COMM_WORLD );

		Debug( 10000 ) << "Setting up data from GetPot \n";
    	M_fsi->setDataFromGetPot( dataFile );

    	MPI_Barrier( MPI_COMM_WORLD );

    	Debug( 10000 ) << "Setting up the BC \n";

		M_fsi->setFluidBC(BCh_fluid(*M_fsi->FSIOper()));
		M_fsi->setHarmonicExtensionBC( BCh_harmonicExtension(*M_fsi->FSIOper()));
		M_fsi->setSolidBC(BCh_solid(*M_fsi->FSIOper()));
		M_fsi->setLinFluidBC(BCh_fluidLin(*M_fsi->FSIOper()));
		M_fsi->setLinSolidBC(BCh_solidLin(*M_fsi->FSIOper()));

		M_absorbingBC = dataFile("fluid/absorbing_bc", false);
		std::cout << "   absorbing BC = " << M_absorbingBC << std::endl;

    	MPI_Barrier( MPI_COMM_WORLD );

		Debug( 10000 ) << "Setting up the problem \n";
    	M_fsi->setup( );

    	//M_fsi->resetFSISolvers();

    	MPI_Barrier( MPI_COMM_WORLD );

		Debug( 10000 ) << "Setting up Ensight \n";
    	if ( M_fsi->isFluid() )
    	{
    		M_ensightFluid.reset( new  filter_type( dataFile, "fixedPtFluid") );

    		M_ensightFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().Comm().MyPID());

    		M_velAndPressure.reset( new vector_type( M_fsi->FSIOper()->fluid().getMap(),      M_ensightFluid->mapType() ));
    		M_fluidDisp.reset     ( new vector_type( M_fsi->FSIOper()->meshMotion().getMap(), M_ensightFluid->mapType() ));

    		M_ensightFluid->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                         UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );

    		M_ensightFluid->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
										UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()),
										UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof()+M_fsi->FSIOper()->pFESpace().dof().numTotalDof()) );

    		M_ensightFluid->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
										UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );
    	}
    	if ( M_fsi->isSolid() )
    	{
    		M_ensightSolid.reset( new  filter_type ( dataFile, "fixedPtSolid") );

    		M_ensightSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().Comm().MyPID());

    		M_solidDisp.reset( new vector_type( M_fsi->FSIOper()->solid().getMap(), M_ensightSolid->mapType() ));
    		M_solidVel.reset ( new vector_type( M_fsi->FSIOper()->solid().getMap(), M_ensightSolid->mapType() ));
    		M_ensightSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
										UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
    		M_ensightSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
										UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
    	}

    	bool restart = dataFile("problem/restart",false);
    	M_Tstart = 0.;
    	if ( restart )
    	{
    		std::string velFName  = dataFile("fluid/miscellaneous/velname"  ,"vel");
    		std::string pressName = dataFile("fluid/miscellaneous/pressname","press");
    		std::string velwName  = dataFile("fluid/miscellaneous/velwname", "velw");
    		std::string depName   = dataFile("solid/miscellaneous/depname"  ,"dep");
    		std::string velSName  = dataFile("solid/miscellaneous/velname"  ,"velw");
    		M_Tstart = dataFile("problem/Tstart"   ,0.);
    		std::cout << "Starting time = " << M_Tstart << std::endl;

    		//M_fsi->initialize(velFName, pressName, velwName, depName, velSName, M_Tstart);

    		if ( M_fsi->isFluid() )
    		{
    			M_ensightFluid->import(M_Tstart, M_fsi->timeStep());
    			M_fsi->FSIOper()->initializeFluid( *M_velAndPressure, *M_fluidDisp );
    		}
    		if ( M_fsi->isSolid() )
    		{
    			M_ensightSolid->import(M_Tstart, M_fsi->timeStep());
				M_fsi->FSIOper()->initializeSolid( M_solidDisp, M_solidVel );
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
    run( Real dt, Real T)
	{
    	std::ofstream ofile;

    	bool const isFluidLeader( M_fsi->FSIOper()->isFluid() && M_fsi->FSIOper()->isLeader() );
    	if ( isFluidLeader )
    		ofile.open( "fluxes.res" );

    	boost::timer _overall_timer;

    	Real flux;
    	int _i = 1;

    	for ( Real time = M_Tstart + dt; time <= T; time += dt, ++_i )
    	{
    		boost::timer _timer;

    		if ( M_absorbingBC && M_fsi->isFluid() )
			{
				BCFunctionBase outFlow;
				outFlow.setFunction(bc_adaptor(*M_fsi->FSIOper()));
				M_fsi->FSIOper()->BCh_fluid()->modifyBC(3, outFlow);
                std::cout << "   F- pressure = " << outFlow(0., 0., 0., 0., 3) << std::endl;
			}

    		M_fsi->iterate( time );

    		if ( M_fsi->isFluid() )
            {
				if ( isFluidLeader )
					ofile << time << " ";

				flux = M_fsi->FSIOper()->fluid().flux(2);
				if ( isFluidLeader )
					ofile << flux << " ";

				flux = M_fsi->FSIOper()->fluid().flux(3);
				if ( isFluidLeader )
					ofile << flux << " ";

				flux = M_fsi->FSIOper()->fluid().pressure(2);
				if ( isFluidLeader )
					ofile << flux << " ";

				flux = M_fsi->FSIOper()->fluid().pressure(3);
				if ( isFluidLeader )
					ofile << flux << " " << std::endl;

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

            std::cout << "solution norm " << _i << " : "
                      << M_fsi->displacement().Norm2() << "\n";

            // CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            checkResult( time );
    	}
		std::cout << "Total computation time = " << _overall_timer.elapsed() << "s" << "\n";
		ofile.close();
	}

private:

    void checkResult( Real& time );
	fsi_solver_ptr M_fsi;
	Real           M_Tstart;

	bool           M_absorbingBC;

	filter_ptrtype M_ensightFluid;
	filter_ptrtype M_ensightSolid;

	vector_ptrtype M_velAndPressure;
	vector_ptrtype M_fluidDisp;
	vector_ptrtype M_solidDisp;
	vector_ptrtype M_solidVel;

    struct RESULT_CHANGED_EXCEPTION
    {
    public:
        RESULT_CHANGED_EXCEPTION(LifeV::Real time)
        {
            std::cout << "Some modifications led to changes in the l2 norm of the solution at time " << time << std::endl;
        }
    };
};





struct FSIChecker
{
    FSIChecker( const GetPot&      dataFile ):
        M_dataFile( dataFile ),
        M_method  ( dataFile( "problem/method", "exactJacobian" ) )
        {}

    FSIChecker( const GetPot&      dataFile,
                const std::string& method   ):
        M_dataFile ( dataFile ),
        M_method   ( method )
        {}

    void operator()()
	{
    	boost::shared_ptr<Problem> FSIproblem;

    	try
    	{
    		FSIproblem = boost::shared_ptr<Problem> ( new Problem( M_dataFile, M_method ) );
    		FSIproblem->run( FSIproblem->fsiSolver()->timeStep(), FSIproblem->fsiSolver()->timeEnd() );
    	}
    	catch ( const std::exception& _ex )
    	{
    		std::cout << "caught exception :  " << _ex.what() << "\n";
    	}

    	//disp = fsip->fsiSolver()->FSIOper()->displacementOnInterface();
	}

	GetPot                M_dataFile;
	std::string           M_method;
	//Vector                 disp;
};


void Problem::checkResult(LifeV::Real& time)
{
    LifeV::Real dispNorm=M_fsi->displacement().Norm2();
    if(time==0.0001 && (dispNorm-0.0621691)/dispNorm*(dispNorm-4.43565e-5)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);  else
    if(time==0.0002 && (dispNorm-0.10668)/dispNorm*(dispNorm-0.000129848)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);else
    if(time==0.0003 && (dispNorm-0.113252)/dispNorm*(dispNorm-0.000251846)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);else
    if(time==0.0004 && (dispNorm-0.107976)/dispNorm*(dispNorm-0.000402367)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);else
    if(time==0.0005 && (dispNorm-0.0995921)/dispNorm*(dispNorm-0.000579832)/dispNorm>1e-5) throw Problem::RESULT_CHANGED_EXCEPTION(time);
}


int main( int argc, char** argv )
{
#ifdef HAVE_MPI
	std::cout << "MPI Initialization" << std::endl;
	MPI_Init( &argc, &argv );
#else
    std::cout << "Serial Version" << std::endl;
#endif

    Chrono chrono;
    chrono.start();

    GetPot command_line( argc, argv );
    std::string dataFileName = command_line.follow( "data", 2, "-f", "--file" );
    GetPot dataFile( dataFileName );

    int returnValue = EXIT_FAILURE;

/*
    const bool check = command_line.search(2, "-c", "--check");
    if (check)
    {
        Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        FSIChecker _ej_check( dataFile, "exactJacobian" );

        _ej_check();

        Debug( 10000 ) << "_ej_disp size : "  << static_cast<Real> (_ej_check.disp.size()) << "\n";
        Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        FSIChecker _sp_check( dataFile, "steklovPoincare" );

        _sp_check();

        Debug( 10000 ) << "_fp_disp size : "  << static_cast<Real> (_sp_check.disp.size()) << "\n";
        Debug( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        Real norm1 = norm_2( _ej_check.disp - _sp_check.disp );

        std::cout << "norm_2(EJ displacement)          = " << norm_2( _ej_check.disp ) << " \n"
                  << "norm_2(SP displacement)          = " << norm_2( _sp_check.disp ) << " \n"
                  << "norm_2(displacement error EJ/SP) = " << norm1 << "\n";

        if (norm1 < 1e-04)
        	returnValue = EXIT_SUCCESS;
        else
        	returnValue = EXIT_FAILURE;
    }
    else
    {

        FSIChecker _sp_check( dataFile );
        _sp_check();

        returnValue = EXIT_SUCCESS;
    }
*/

    FSIChecker FSIProblem( dataFile );
    FSIProblem();

    returnValue = EXIT_SUCCESS;

    std::cout << "Total sum up " << chrono.diff_cumul() << " s." << std::endl;

#ifdef HAVE_MPI
	std::cout << "MPI Finalization" << std::endl;
	MPI_Finalize();
#endif

    return returnValue;
}
