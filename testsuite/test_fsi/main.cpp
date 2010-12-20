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
    @file
    @brief

    @author
    @date 00-00-0000
*/


#ifdef TWODIM
#error test_fsi cannot be compiled in 2D
#endif

#include <boost/timer.hpp>

#include <cassert>
#include <iomanip>
#include <cmath>

//#include <life/lifealg/IfpackPreconditioner.hpp>
//#include <life/lifealg/MLPreconditioner.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/chrono.hpp>

#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSI.hpp>
//#include "life/lifesolver/exactJacobianBase.hpp"
//#include "life/lifesolver/fixedPointBase.hpp"
#include <life/lifesolver/DataFSI.hpp>
#include <life/lifesolver/VenantKirchhoffSolverLinear.hpp>

#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifefilters/ensight.hpp>

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

namespace LifeV
{
namespace
{
LifeV::VenantKirchhoffSolver< LifeV::FSI::mesh_Type, LifeV::SolverTrilinos >*    createLinearStructure() { return new VenantKirchhoffSolverLinear< LifeV::FSI::mesh_Type, LifeV::SolverTrilinos >(); }

//NOTE: the nonlinear structure solver is still in development in the FSI framework
//LifeV::VenantKirchhofSolver< LifeV::FSI::mesh_Type, LifeV::SolverTrilinos >*    createNonLinearStructure(){ return new NonLinearVenantKirchhofSolver< LifeV::FSI::mesh_Type, LifeV::SolverTrilinos >(); }
}
}

using namespace LifeV;

int returnValue = EXIT_SUCCESS; // For the final check



#define PI 3.141592653589793
class bc_adaptor
{
public:

    bc_adaptor( FSI& Operator ):
            M_oper   ( Operator )
    {
        Real area0	= 0.7854;
        //Real area0	= M_oper.fluid().area(3);
        Real area	= area0;

        Real beta	= M_oper.solid().getThickness()*M_oper.solid().getYoung() /
                    (1 - M_oper.solid().getPoisson()*M_oper.solid().getPoisson()) * PI/area0;

        Real qn		= M_oper.fluid().flux(3);

        M_outflow			= std::pow(std::sqrt(M_oper.solid().getRho())/(2*std::sqrt(2.))*qn/area + std::sqrt(beta*std::sqrt(area0)), 2)
                      - beta*std::sqrt(area0);

        std::cout << "--------------- Absorbing boundary condition ---------------" << std::endl;
        std::cout << "  Outflow BC : density   = " << M_oper.solid().getRho() << std::endl;
        std::cout << "  Outflow BC : thickness = " << M_oper.solid().getThickness() << std::endl;
        std::cout << "  Outflow BC : young     = " << M_oper.solid().getYoung() << std::endl;
        std::cout << "  Outflow BC : poisson   = " << M_oper.solid().getPoisson() << std::endl;
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
        switch ( id )
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

    FSI& M_oper;
    Real         M_outflow;
};





class Problem
{
public:

    typedef boost::shared_ptr<FSISolver>                    fsi_solver_ptr;
    typedef FSI::data_Type                          data_Type;
    typedef FSI::dataPtr_Type                       dataPtr_Type;

    typedef FSI::vector_Type                        vector_Type;
    typedef FSI::vectorPtr_Type                     vectorPtr_Type;

    typedef Exporter<FSI::mesh_Type>                filter_type;
    typedef boost::shared_ptr<filter_type>                  filter_ptrtype;

    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */
    Problem( const std::string& dataFileName, std::string method = "" )
    {

        VenantKirchhoffSolver< FSI::mesh_Type, SolverTrilinos >::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhof", &createLinearStructure );
        //        VenantKirchhofSolver< FSIOperator::mesh_Type, SolverTrilinos >::StructureSolverFactory::instance().registerProduct( "nonLinearVenantKirchhof", &createNonLinearStructure );

        Debug( 10000 ) << "Setting up data from GetPot \n";
        GetPot dataFile( dataFileName );
        M_data = dataPtr_Type( new data_Type() );
        M_data->setup( dataFile );
        M_data->dataSolid()->setDataTime( M_data->dataFluid()->dataTime() ); //Same dataTime for fluid & solid
        //M_data->showMe();
        MPI_Barrier( MPI_COMM_WORLD );

        Debug( 10000 ) << "creating FSISolver with operator :  " << method << "\n";
        M_fsi = fsi_solver_ptr( new FSISolver( ) );
        M_fsi->setData( M_data );
        M_fsi->FSIOper()->setDataFile( dataFile ); //TO BE REMOVED!
        MPI_Barrier( MPI_COMM_WORLD );

        // Setting FESpace and DOF
        Debug( 10000 ) << "Setting up the FESpace and DOF \n";
        M_fsi->FSIOper( )->partitionMeshes( );
        M_fsi->FSIOper()->setupFEspace();
        M_fsi->FSIOper()->setupDOF();
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

        std::string const exporterType =  dataFile( "exporter/type", "hdf5");
        std::string const exporterName =  dataFile( "exporter/name", "fixedPt");

        Debug( 10000 ) << "Setting up Ensight \n";
        if ( M_fsi->isFluid() )
        {
#ifdef HAVE_HDF5
            if (exporterType.compare("hdf5") == 0)
                M_exporterFluid.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, exporterName+"Fluid" ) );
            else
#endif
                M_exporterFluid.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, exporterName+"Fluid" ) );


            M_exporterFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());

            M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(),      M_exporterFluid->mapType() ));
            M_fluidDisp.reset     ( new vector_Type( M_fsi->FSIOper()->meshMotion().getMap(), M_exporterFluid->mapType() ));

            M_exporterFluid->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                          UInt(0), M_fsi->FSIOper()->uFESpace().dof().numTotalDof() );

            M_exporterFluid->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
                                          UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof() ),
                                          UInt(  M_fsi->FSIOper()->pFESpace().dof().numTotalDof() ) );

            M_exporterFluid->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
                                          UInt(0), M_fsi->FSIOper()->mmFESpace().dof().numTotalDof() );
        }
        if ( M_fsi->isSolid() )
        {
#ifdef HAVE_HDF5
            if (exporterType.compare("hdf5") == 0)
                M_exporterSolid.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, exporterName+"Solid" ) );
            else
#endif
                M_exporterSolid.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, exporterName+"Solid" ) );

            M_exporterSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID());

            M_solidDisp.reset( new vector_Type( M_fsi->FSIOper()->solid().getMap(), M_exporterSolid->mapType() ));
            M_solidVel.reset ( new vector_Type( M_fsi->FSIOper()->solid().getMap(), M_exporterSolid->mapType() ));
            M_exporterSolid->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                                          UInt(0), M_fsi->FSIOper()->dFESpace().dof().numTotalDof() );
            M_exporterSolid->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
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
                M_exporterFluid->import(M_Tstart, M_data->dataFluid()->dataTime()->timeStep());
                M_fsi->FSIOper()->initializeFluid( *M_velAndPressure, *M_fluidDisp );
            }
            if ( M_fsi->isSolid() )
            {
               M_exporterSolid->import(M_Tstart, M_data->dataSolid()->getDataTime()->timeStep());
               M_fsi->FSIOper()->initializeSolid( M_solidDisp, M_solidVel );
            }
        }
        else
        {
            M_fsi->initialize();
        }
        M_data->dataFluid()->dataTime()->setInitialTime( M_Tstart + M_data->dataFluid()->dataTime()->timeStep() );
        M_data->dataFluid()->dataTime()->setTime( M_data->dataFluid()->dataTime()->initialTime() );
        //std::cout << "in problem" << std::endl;
        //M_fsi->FSIOper()->fluid().postProcess();
    }

    fsi_solver_ptr fsiSolver() { return M_fsi; }

    dataPtr_Type fsiData() { return M_data; }

    /*!
      This routine runs the temporal loop
    */
    void run()
    {
        std::ofstream ofile;

        bool const isFluidLeader( M_fsi->FSIOper()->isFluid() && M_fsi->FSIOper()->isLeader() );
        if ( isFluidLeader )
            ofile.open( "fluxes.res" );

        boost::timer _overall_timer;

        Real flux;
        int _i = 1;

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(), ++_i )
        {
            boost::timer _timer;

            if ( M_absorbingBC && M_fsi->isFluid() )
            {
                BCFunctionBase outFlow;
                outFlow.setFunction(bc_adaptor(*M_fsi->FSIOper()));
                M_fsi->FSIOper()->BCh_fluid()->modifyBC(3, outFlow);
                //std::cout << "  F-  Pressure = " << outFlow(0., 0., 0., 0., 3) << std::endl;
            }

            M_fsi->iterate();

            if ( M_fsi->isFluid() )
            {
                if ( isFluidLeader )
                    ofile << M_data->dataFluid()->dataTime()->time() << " ";

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

                *M_velAndPressure = *M_fsi->FSIOper()->fluid().solution();
                *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
                M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );
            }

            if ( M_fsi->isSolid() )
	      {
                *M_solidDisp = M_fsi->FSIOper()->solid().getDisplacement();
                *M_solidVel = M_fsi->FSIOper()->solid().getVelocity();
                M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );
	      }

            std::cout << "[fsi_run] Iteration " << _i << " was done in : " << _timer.elapsed() << "\n";

            std::cout << "solution norm " << _i << " : "
                      << M_fsi->displacement().norm2() << "\n";

            // CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            checkResult( M_data->dataFluid()->dataTime()->time() );
        }
        std::cout << "Total computation time = " << _overall_timer.elapsed() << "s" << "\n";
        ofile.close();
    }

private:

    void checkResult(const LifeV::Real& time)
    {
        assert(M_data->dataFluid()->dataTime()->timeStep()==0.001);
        double dispNorm(M_fsi->displacement().norm2());

        const LifeV::Real relTol(5e-3);

        if ( sameAs(time,0.001) && sameAs(dispNorm, 0.0621691, relTol) ) return;
        if ( sameAs(time,0.002) && sameAs(dispNorm, 0.10668,   relTol) )  return;
        if ( sameAs(time,0.003) && sameAs(dispNorm, 0.113252,  relTol) )  return;
        if ( sameAs(time,0.004) && sameAs(dispNorm, 0.107976,  relTol) )  return;
        if ( sameAs(time,0.005) && sameAs(dispNorm, 0.0995918, relTol) )  return;
        if ( sameAs(time,0.006) && sameAs(dispNorm, 0.0751478, relTol) ) return;

        throw Problem::RESULT_CHANGED_EXCEPTION(time);

    }

    bool sameAs(const LifeV::Real& a, const LifeV::Real& b, const LifeV::Real& relTol = 1e-6)
    {
        const LifeV::Real maxAbs (std::max(std::abs(a),std::abs(b)));
        if (maxAbs < relTol*relTol) return true;

        return std::abs(a-b) < relTol*maxAbs;
    }

    fsi_solver_ptr M_fsi;
    dataPtr_Type   M_data;
    Real           M_Tstart;

    bool           M_absorbingBC;

    filter_ptrtype M_exporterFluid;
    filter_ptrtype M_exporterSolid;

    vectorPtr_Type M_velAndPressure;
    vectorPtr_Type M_fluidDisp;
    vectorPtr_Type M_solidDisp;
    vectorPtr_Type M_solidVel;

    struct RESULT_CHANGED_EXCEPTION
    {
public:
        RESULT_CHANGED_EXCEPTION(LifeV::Real time)
        {
            std::cout << "Some modifications led to changes in the l2 norm of the solution at time " << time << std::endl;
            returnValue = EXIT_FAILURE;
        }
    };
};





struct FSIChecker
{
    FSIChecker( const std::string& dataFileName ):
            M_dataFileName ( dataFileName ),
            M_method       ()
    {
        GetPot dataFile( dataFileName );
        M_method = dataFile( "problem/method", "exactJacobian" );
    }

    FSIChecker( const std::string& dataFileName,
                const std::string& method   ):
            M_dataFileName ( dataFileName ),
            M_method       ( method )
    {}

    void operator()()
    {
        boost::shared_ptr<Problem> FSIproblem;

        try
        {
            FSIproblem = boost::shared_ptr<Problem> ( new Problem( M_dataFileName, M_method ) );
            FSIproblem->run();
        }
        catch ( const std::exception& _ex )
        {
            std::cout << "caught exception :  " << _ex.what() << "\n";
        }
    }

    std::string           M_dataFileName;
    std::string           M_method;
};

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

    FSIChecker FSIProblem( dataFileName );
    FSIProblem();

    std::cout << "Total sum up " << chrono.diffCumul() << " s." << std::endl;

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    return returnValue;
}
