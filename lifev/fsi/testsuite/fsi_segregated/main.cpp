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

    * The time discretization is carried out using BDF methods of order 2. At the moment, even is the Newmark method is available
    * for the temporal discretization of the single problems( e.g. in test_structuralsolver), it cannot be used in the FSI framework
    * since the class TimeAdvanceNewmark is not registered as one of the possible instances of the abstrac class TimeAdvance.
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

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/fsi/solver/FSISolver.hpp>
#include <lifev/fsi/solver/FSIOperator.hpp>
#include <lifev/fsi/solver/FSIExactJacobian.hpp>
#include <lifev/fsi/solver/FSIFixedPoint.hpp>
#include <lifev/fsi/solver/FSIData.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>

#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include "ud_functions.hpp"
#include "boundaryConditions.hpp"


namespace LifeV
{
namespace
{
static bool regIF = PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack );
static bool regML = PRECFactory::instance().registerProduct ( "ML", &createML );
static bool regFP = FSIOperator::FSIFactory_Type::instance().registerProduct ( "fixedPoint", &createFP );
static bool regEJ = FSIOperator::FSIFactory_Type::instance().registerProduct ( "exactJacobian", &createEJ );
}
}


namespace LifeV
{
namespace
{

//LifeV::VenantKirchhoffSolver< LifeV::FSIOperator::mesh_Type, LifeV::SolverAztecOO >*    createLinearStructure() { return new VenantKirchhoffSolverLinear< LifeV::FSIOperator::mesh_Type, LifeV::SolverAztecOO >(); }



//LifeV::StructuralMaterial< LifeV::FSIOperator::mesh_Type >*    createVenantKirchhoffLinear() { return new VenantKirchhoffMaterialLinear< LifeV::FSIOperator::mesh_Type >(); }




//LifeV::StructuralMaterial< LifeV::FSIOperator::mesh_Type >*    createVenantKirchhoffLinear() { return new VenantKirchhoffMaterialLinear< LifeV::FSIOperator::mesh_Type >(); }


//LifeV::VenantKirchhoffSolver< LifeV::FSIOperator::mesh_Type, LifeV::SolverAztecOO >*    createLinearStructure() {
//return new VenantKirchhoffSolverLinear< LifeV::FSIOperator::mesh_Type, LifeV::SolverAztecOO >();
//}


//NOTE: the nonlinear structure solver is still in development in the FSI framework
//LifeV::VenantKirchhofSolver< LifeV::FSI::mesh_Type, LifeV::SolverAztecOO >*    createNonLinearStructure(){ return new NonLinearVenantKirchhofSolver< LifeV::FSI::mesh_Type, LifeV::SolverAztecOO >(); }
}
}

using namespace LifeV;

int returnValue = EXIT_SUCCESS; // For the final check



#define PI 3.141592653589793
class bc_adaptor
{
public:

    bc_adaptor ( FSIOperator& Operator ) :
        M_oper   ( Operator )
    {
        Real area0    = 0.7854;
        //Real area0    = M_oper.fluid().area(3);
        Real area    = area0;
        UInt flag       = M_oper.dataSolid()->vectorFlags() [0];

        Real beta    = M_oper.solid().thickness() * M_oper.solid().young ( flag ) /
                       (1 - M_oper.solid().poisson (  flag ) * M_oper.solid().poisson ( flag ) ) * PI / area0;

        Real qn        = M_oper.fluid().flux (3);

        M_outflow    = std::pow (std::sqrt (M_oper.solid().rho() ) / (2 * std::sqrt (2.) ) * qn / area + std::sqrt (beta * std::sqrt (area0) ), 2) - beta * std::sqrt (area0);

        std::cout << "--------------- Absorbing boundary condition ---------------" << std::endl;
        std::cout << "  Outflow BC : density   = " << M_oper.solid().rho() << std::endl;
        std::cout << "  Outflow BC : thickness = " << M_oper.solid().thickness() << std::endl;
        std::cout << "  Outflow BC : young     = " << M_oper.solid().young ( flag ) << std::endl;
        std::cout << "  Outflow BC : poisson   = " << M_oper.solid().poisson ( flag ) << std::endl;
        std::cout << "  Outflow BC : area0     = " << area0 << std::endl;
        std::cout << "  Outflow BC : area      = " << M_oper.fluid().area (3) << std::endl;
        std::cout << "  Outflow BC : radius    = " << std::sqrt (area0 / PI) << std::endl;
        std::cout << "  Outflow BC : beta      = " << beta << std::endl;
        std::cout << "  Outflow BC : Flow rate = " << qn << std::endl;
        std::cout << "  Outflow BC : outflow   = " << M_outflow << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
    }

    Real operator() ( Real /*t*/, Real /*x*/, Real /*y*/, Real /*z*/, ID id)
    {
        switch ( id )
        {
            case 0:
                return 0.;
                break;

            case 1:
                return 0.;
                break;

            case 2:
                //return 0.;
                return -M_outflow;
                break;

            default:
                ERROR_MSG ("This entrie is not allowed: disp_adatptor");
                break;
        }

        return 0.;
    }

private:

    FSIOperator& M_oper;
    Real         M_outflow;
};

class bc_adaptorFace
{
public:

    bc_adaptorFace ( FSIOperator& Operator ) :
        M_oper   ( Operator )
    {
        Real area0    = 0.147439938;
        //Real area    = area0;
        UInt flag       = 1;

        //Real beta    = M_oper.solid().thickness()*M_oper.solid().young() / (1 - M_oper.solid().poisson()*M_oper.solid().poisson()) * PI/area0;
        //Alexandra's Abc
        Real exp  = 5 / 4;
        Real beta = ( std::sqrt (PI) * M_oper.solid().thickness() * M_oper.solid().young ( flag ) ) / (1 - M_oper.solid().poisson ( flag ) * M_oper.solid().poisson ( flag ) );
        Real R    = ( std::sqrt (M_oper.solid().rho( ) * beta ) ) / ( std::sqrt (2.0) * std::pow (area0, exp) );

        Real qn        = M_oper.fluid().flux (3);

        M_outflowFace    = R * qn;
        //M_outflowFace    = std::pow(std::sqrt(M_oper.solid().rho())/(2*std::sqrt(2.))*qn/area + std::sqrt(beta*std::sqrt(area0)), 2)- beta*std::sqrt(area0);

        std::cout << "--------------- Absorbing boundary condition for Face-------" << std::endl;
        std::cout << "  Outflow BC : density   = " << M_oper.solid().rho() << std::endl;
        std::cout << "  Outflow BC : thickness = " << M_oper.solid().thickness() << std::endl;
        std::cout << "  Outflow BC : young     = " << M_oper.solid().young ( flag ) << std::endl;
        std::cout << "  Outflow BC : poisson   = " << M_oper.solid().poisson ( flag ) << std::endl;
        std::cout << "  Outflow BC : area0     = " << area0 << std::endl;
        std::cout << "  Outflow BC : area      = " << M_oper.fluid().area (3) << std::endl;
        std::cout << "  Outflow BC : radius    = " << std::sqrt (area0 / PI) << std::endl;
        std::cout << "  Outflow BC : beta      = " << beta << std::endl;
        std::cout << "  Outflow BC : Flow rate = " << qn << std::endl;
        std::cout << "  Outflow BC : outflow   = " << M_outflowFace << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;

    }

    Real operator() ( Real /*t*/, Real /*x*/, Real /*y*/, Real /*z*/, ID id)
    {
        switch ( id )
        {
            case 0:
                return 0.155851 * M_outflowFace;
                break;

            case 1:
                return -0.987781 * M_outflowFace;
                break;

            case 2:
                //return 0.;
                return 9.02676e-06 * M_outflowFace;
                break;

            default:
                ERROR_MSG ("This entrie is not allowed: disp_adatptor");
                break;
        }

        return 0.;
    }

private:

    FSIOperator& M_oper;
    Real         M_outflowFace;
};


class bc_adaptorBrain
{
public:

    bc_adaptorBrain ( FSIOperator& Operator ) :
        M_oper   ( Operator )
    {
        Real area0    = 0.191155176;
        //Real area    = area0;
        UInt flag       = 1 ;

        //Real beta    = M_oper.solid().thickness()*M_oper.solid().young() / (1 - M_oper.solid().poisson()*M_oper.solid().poisson()) * PI/area0;

        //Alexandra's Abc
        Real exp  = 5 / 4;
        Real beta = ( std::sqrt (PI) * M_oper.solid().thickness() * M_oper.solid().young ( flag ) ) / (1 - M_oper.solid().poisson ( flag ) * M_oper.solid().poisson ( flag ) );
        Real R    = ( std::sqrt (M_oper.solid().rho( ) * beta ) ) / ( std::sqrt (2.0) * std::pow (area0, exp) );

        Real qn        = M_oper.fluid().flux (4);

        M_outflowBrain = R * qn;
        //M_outflowBrain  = std::pow(std::sqrt(M_oper.solid().rho())/(2*std::sqrt(2.))*qn/area + std::sqrt(beta*std::sqrt(area0)), 2) - beta*std::sqrt(area0);


        std::cout << "--------------- Absorbing boundary condition for Brain-------" << std::endl;
        std::cout << "  Outflow BC : density   = " << M_oper.solid().rho() << std::endl;
        std::cout << "  Outflow BC : thickness = " << M_oper.solid().thickness() << std::endl;
        std::cout << "  Outflow BC : young     = " << M_oper.solid().young ( flag ) << std::endl;
        std::cout << "  Outflow BC : poisson   = " << M_oper.solid().poisson ( flag ) << std::endl;
        std::cout << "  Outflow BC : area0     = " << area0 << std::endl;
        std::cout << "  Outflow BC : area      = " << M_oper.fluid().area (4) << std::endl;
        std::cout << "  Outflow BC : radius    = " << std::sqrt (area0 / PI) << std::endl;
        std::cout << "  Outflow BC : beta      = " << beta << std::endl;
        std::cout << "  Outflow BC : Flow rate = " << qn << std::endl;
        std::cout << "  Outflow BC : outflow   = " << M_outflowBrain << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
    }

    Real operator() ( Real /*t*/, Real /*x*/, Real /*y*/, Real /*z*/, ID id)
    {
        switch ( id )
        {
            case 0:
                return -0.0979856 * M_outflowBrain;
                break;

            case 1:
                return -0.995188 * M_outflowBrain;
                break;

            case 2:
                return -3.13375e-05 * M_outflowBrain;
                break;

            default:
                ERROR_MSG ("This entrie is not allowed: disp_adatptor");
                break;
        }

        return 0.;
    }

private:

    FSIOperator& M_oper;
    Real         M_outflowBrain;
};




class Problem
{
public:

    typedef boost::shared_ptr<FSISolver>                    fsi_solver_ptr;
    typedef FSIOperator::data_Type                          data_Type;
    typedef FSIOperator::dataPtr_Type                       dataPtr_Type;

    typedef FSIOperator::vector_Type                        vector_Type;
    typedef FSIOperator::vectorPtr_Type                     vectorPtr_Type;
    typedef FSIOperator::mesh_Type                          mesh_Type;

    typedef Exporter<FSIOperator::mesh_Type>                filter_type;
    typedef boost::shared_ptr<filter_type>                  filter_ptrtype;

    typedef FESpace<FSIOperator::mesh_Type, MapEpetra>      feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>                 feSpacePtr_Type;

    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */
    Problem ( const std::string& dataFileName, std::string method = "" )
    {

        //VenantKirchhoffSolver< FSIOperator::mesh_Type, SolverAztecOO >::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhoff", &createLinearStructure );


        //StructuralSolver< FSIOperator::mesh_Type, SolverAztecOO >::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhoff", &createVenantKirchhoffLinear );

        StructuralOperator< FSIOperator::mesh_Type>();

        //StructuralSolver< FSIOperator::mesh_Type, SolverAztecOO >::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhoff", &createVenantKirchhoffLinear );

        //        VenantKirchhofSolver< FSIOperator::mesh_Type, SolverAztecOO >::StructureSolverFactory::instance().registerProduct( "nonLinearVenantKirchhoff", &createNonLinearStructure );

        debugStream ( 10000 ) << "Setting up data from GetPot \n";
        GetPot dataFile ( dataFileName );
        M_data = dataPtr_Type ( new data_Type() );
        M_data->setup ( dataFile );
        M_data->dataSolid()->setTimeData ( M_data->dataFluid()->dataTime() ); //Same TimeData for fluid & solid
        M_data->showMe();
        //    M_data->dataSolid()->showMe();
        MPI_Barrier ( MPI_COMM_WORLD );

        debugStream ( 10000 ) << "creating FSISolver with operator :  " << method << "\n";
        M_fsi = fsi_solver_ptr ( new FSISolver( ) );
        M_fsi->setData ( M_data );
        M_fsi->FSIOper()->setDataFile ( dataFile ); //TO BE REMOVED!
        MPI_Barrier ( MPI_COMM_WORLD );

        // Setting FESpace and DOF
        debugStream ( 10000 ) << "Setting up the FESpace and DOF \n";
        M_fsi->FSIOper( )->partitionMeshes( );
        M_fsi->FSIOper()->setupFEspace();
        M_fsi->FSIOper()->setupDOF();
        MPI_Barrier ( MPI_COMM_WORLD );

        debugStream ( 10000 ) << "Setting up the BC \n";
        M_fsi->setFluidBC (BCh_fluid (*M_fsi->FSIOper() ) );
        M_fsi->setHarmonicExtensionBC ( BCh_harmonicExtension (*M_fsi->FSIOper() ) );
        M_fsi->setSolidBC (BCh_solid (*M_fsi->FSIOper() ) );
        M_fsi->setLinFluidBC (BCh_fluidLin (*M_fsi->FSIOper() ) );
        M_fsi->setLinSolidBC (BCh_solidLin (*M_fsi->FSIOper() ) );

        M_absorbingBC = dataFile ("fluid/absorbing_bc", false);
        std::cout << "   absorbing BC = " << M_absorbingBC << std::endl;

        MPI_Barrier ( MPI_COMM_WORLD );

        debugStream ( 10000 ) << "Setting up the problem \n";
        M_fsi->setup( );

        //M_fsi->resetFSISolvers();

        MPI_Barrier ( MPI_COMM_WORLD );

        std::string const exporterType =  dataFile ( "exporter/type", "hdf5");
        std::string const exporterName =  dataFile ( "exporter/name", "fixedPt");

        debugStream ( 10000 ) << "Setting up ExporterEnsight \n";
        if ( M_fsi->isFluid() )
        {
#ifdef HAVE_HDF5
            if (exporterType.compare ("hdf5") == 0)
            {
                M_exporterFluid.reset ( new ExporterHDF5<mesh_Type > ( dataFile, exporterName + "Fluid" ) );
            }
            else
#endif
                M_exporterFluid.reset ( new ExporterEnsight<mesh_Type > ( dataFile, exporterName + "Fluid" ) );


            M_exporterFluid->setMeshProcId (M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID() );

            M_velAndPressure.reset ( new vector_Type ( M_fsi->FSIOper()->fluid().getMap(),      M_exporterFluid->mapType() ) );
            M_fluidDisp.reset     ( new vector_Type ( M_fsi->FSIOper()->meshMotion().getMap(), M_exporterFluid->mapType() ) );

            M_exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f-velocity",
                                           M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressure, UInt (0) );

            M_exporterFluid->addVariable ( ExporterData<mesh_Type>::ScalarField, "f-pressure",
                                           M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressure,
                                           UInt (3 * M_fsi->FSIOper()->uFESpace().dof().numTotalDof() ) );

            M_exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f-displacement",
                                           M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDisp, UInt (0) );

        }
        if ( M_fsi->isSolid() )
        {
#ifdef HAVE_HDF5
            if (exporterType.compare ("hdf5") == 0)
            {
                M_exporterSolid.reset ( new ExporterHDF5<mesh_Type > ( dataFile, exporterName + "Solid" ) );
            }
            else
#endif
                M_exporterSolid.reset ( new ExporterEnsight<mesh_Type > ( dataFile, exporterName + "Solid" ) );

            M_exporterSolid->setMeshProcId (M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID() );
            M_solidDisp.reset ( new vector_Type ( M_fsi->FSIOper()->solid().map(), M_exporterSolid->mapType() ) );
            M_solidVel.reset ( new vector_Type ( M_fsi->FSIOper()->solid().map(), M_exporterSolid->mapType() ) );

            M_exporterSolid->addVariable ( ExporterData<mesh_Type>::VectorField, "s-displacement",
                                           M_fsi->FSIOper()->dFESpacePtr(), M_solidDisp, UInt (0) );
            M_exporterSolid->addVariable ( ExporterData<mesh_Type>::VectorField, "s-velocity",
                                           M_fsi->FSIOper()->dFESpacePtr(), M_solidVel, UInt (0) );

        }

        bool restart = dataFile ("problem/restart", false);
        M_Tstart = 0.;
        if ( restart )
        {
            std::string velFName  = dataFile ("fluid/miscellaneous/velname"  , "vel");
            std::string pressName = dataFile ("fluid/miscellaneous/pressname", "press");
            std::string velwName  = dataFile ("fluid/miscellaneous/velwname", "velw");
            std::string depName   = dataFile ("solid/miscellaneous/depname"  , "dep");
            std::string velSName  = dataFile ("solid/miscellaneous/velname"  , "velw");
            M_Tstart = dataFile ("problem/initialtime"   , 0.);
            std::cout << "Starting time = " << M_Tstart << std::endl;

            //M_fsi->initialize(velFName, pressName, velwName, depName, velSName, M_Tstart);

            if ( M_fsi->isFluid() )
            {
                M_exporterFluid->import (M_Tstart, M_data->dataFluid()->dataTime()->timeStep() );
                M_fsi->FSIOper()->initializeFluid ( *M_velAndPressure, *M_fluidDisp );
            }
            if ( M_fsi->isSolid() )
            {
                M_exporterSolid->import (M_Tstart, M_data->dataSolid()->dataTime()->timeStep() );
                M_fsi->FSIOper()->initializeSolid ( M_solidDisp, M_solidVel );
            }
        }
        else
        {
            M_fsi->initialize();
        }
        M_data->dataFluid()->dataTime()->setInitialTime ( M_Tstart ); //+ M_data->dataFluid()->dataTime()->timeStep()
        M_data->dataFluid()->dataTime()->setTime ( M_Tstart  );
        //std::cout << "in problem" << std::endl;
        //M_fsi->FSIOper()->fluid().postProcess();
    }

    fsi_solver_ptr fsiSolver()
    {
        return M_fsi;
    }

    dataPtr_Type fsiData()
    {
        return M_data;
    }

    /*!
      This routine runs the temporal loop
    */
    void run()
    {
        std::ofstream ofile;

        bool const isFluidLeader ( M_fsi->FSIOper()->isFluid() && M_fsi->FSIOper()->isLeader() );
        if ( isFluidLeader )
        {
            ofile.open ( "fluxes.res" );
        }

        boost::timer _overall_timer;

        Real flux;
        int _i = 1;

        for ( ; M_data->dataFluid()->dataTime()->canAdvance(); M_data->dataFluid()->dataTime()->updateTime(), ++_i )
        {
            boost::timer _timer;

            if ( M_absorbingBC && M_fsi->isFluid() )
            {

                BCFunctionBase outFlow;
                outFlow.setFunction (bc_adaptor (*M_fsi->FSIOper() ) );
                M_fsi->FSIOper()->BCh_fluid()->modifyBC (3, outFlow);

                /*
                  BCFunctionBase outFlowFace;
                      BCFunctionBase outFlowBrain;

                      outFlowFace.setFunction(bc_adaptorFace(*M_fsi->FSIOper()));
                  M_fsi->FSIOper()->BCh_fluid()->modifyBC(3, outFlowFace);

                  outFlowBrain.setFunction(bc_adaptorBrain(*M_fsi->FSIOper()));
                  M_fsi->FSIOper()->BCh_fluid()->modifyBC(4, outFlowBrain);
                */
                //std::cout << "  F-  Pressure = " << outFlow(0., 0., 0., 0., 3) << std::endl;
            }

            M_fsi->iterate();

            if ( M_fsi->isFluid() )
            {
                if ( isFluidLeader )
                {
                    ofile << M_data->dataFluid()->dataTime()->time() << " ";
                }

                flux = M_fsi->FSIOper()->fluid().flux (2);
                if ( isFluidLeader )
                {
                    ofile << flux << " ";
                }

                flux = M_fsi->FSIOper()->fluid().flux (3);
                if ( isFluidLeader )
                {
                    ofile << flux << " ";
                }

                flux = M_fsi->FSIOper()->fluid().pressure (2);
                if ( isFluidLeader )
                {
                    ofile << flux << " ";
                }

                flux = M_fsi->FSIOper()->fluid().pressure (3);
                if ( isFluidLeader )
                {
                    ofile << flux << " " << std::endl;
                }

                *M_velAndPressure = *M_fsi->FSIOper()->fluid().solution();
                *M_fluidDisp      = M_fsi->FSIOper()->meshMotion().disp();
                M_exporterFluid->postProcess ( M_data->dataFluid()->dataTime()->time() );
            }

            if ( M_fsi->isSolid() )
            {
                *M_solidDisp = M_fsi->FSIOper()->solid().displacement();
                // *M_solidVel = M_fsi->FSIOper()->solid().velocity();
                *M_solidVel = M_fsi->FSIOper()->solidTimeAdvance()->firstDerivative();

                M_exporterSolid->postProcess ( M_data->dataFluid()->dataTime()->time() );
            }

            std::cout << "[fsi_run] Iteration " << _i << " was done in : " << _timer.elapsed() << "\n";

            // CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            checkResult ( M_data->dataFluid()->dataTime()->time() );
        }
        std::cout << "Total computation time = " << _overall_timer.elapsed() << "s" << "\n";
        ofile.close();
    }

private:

    void checkResult (const LifeV::Real& time)
    {
        assert (M_data->dataFluid()->dataTime()->timeStep() == 0.001);
        double dispNorm (M_fsi->displacement().norm2() );

        std::cout << "Displ Norm: " << dispNorm << std::endl;

        const LifeV::Real relTol (5e-3);

        if ( sameAs (time, 0) && sameAs (dispNorm, 0.0474091, relTol) )
        {
            return;
        }
        if ( sameAs (time, 0.001) && sameAs (dispNorm, 0.0564487,   relTol) )
        {
            return;
        }
        if ( sameAs (time, 0.002) && sameAs (dispNorm, 0.0618011,  relTol) )
        {
            return;
        }

        throw RESULT_CHANGED_EXCEPTION (time);

    }

    bool sameAs (const LifeV::Real& a, const LifeV::Real& b, const LifeV::Real& relTol = 1e-6)
    {
        const LifeV::Real maxAbs (std::max (std::fabs (a), std::fabs (b) ) );
        if (maxAbs < relTol * relTol)
        {
            return true;
        }

        return std::fabs (a - b) < relTol * maxAbs;
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
        RESULT_CHANGED_EXCEPTION (LifeV::Real time)
        {
            std::cout << "Some modifications led to changes in the l2 norm of the solution at time " << time << std::endl;
            returnValue = EXIT_FAILURE;
        }
    };
};





struct FSIChecker
{
    FSIChecker ( const std::string& dataFileName ) :
        M_dataFileName ( dataFileName ),
        M_method       ()
    {
        GetPot dataFile ( dataFileName );
        M_method = dataFile ( "problem/method", "exactJacobian" );
    }

    FSIChecker ( const std::string& dataFileName,
                 const std::string& method   ) :
        M_dataFileName ( dataFileName ),
        M_method       ( method )
    {}

    void operator() ()
    {
        boost::shared_ptr<Problem> FSIproblem;

        try
        {
            FSIproblem = boost::shared_ptr<Problem> ( new Problem ( M_dataFileName, M_method ) );
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

int main ( int argc, char** argv )
{
#ifdef HAVE_MPI
    std::cout << "MPI Initialization" << std::endl;
    MPI_Init ( &argc, &argv );
#else
    std::cout << "Serial Version" << std::endl;
#endif

    LifeChrono chrono;
    chrono.start();

    GetPot command_line ( argc, argv );
    std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );

    /*
        const bool check = command_line.search(2, "-c", "--check");
        if (check)
        {
            debugStream( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            FSIChecker _ej_check( dataFile, "exactJacobian" );

            _ej_check();

            debugStream( 10000 ) << "_ej_disp size : "  << static_cast<Real> (_ej_check.disp.size()) << "\n";
            debugStream( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            debugStream( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            FSIChecker _sp_check( dataFile, "steklovPoincare" );

            _sp_check();

            debugStream( 10000 ) << "_fp_disp size : "  << static_cast<Real> (_sp_check.disp.size()) << "\n";
            debugStream( 10000 ) << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

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

    FSIChecker FSIProblem ( dataFileName );
    FSIProblem();

    std::cout << "Total sum up " << chrono.diffCumul() << " s." << std::endl;

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    return returnValue;
}
