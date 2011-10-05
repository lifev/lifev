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

#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/LifeChrono.hpp>

#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifesolver/FSIData.hpp>
#include <life/lifesolver/StructuralSolver.hpp>
#include <life/lifesolver/StructuralMaterial.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialNonLinear.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialLinear.hpp>

#include <life/lifesolver/StructuralMaterial.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialNonLinear.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialLinear.hpp>

#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

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
class analyticalSolution
{
public:

    analyticalSolution( FSIOperator& Operator ):
            M_oper   ( Operator )
    {
     flag = 1 ;
     rhoFluid =  M_oper.fluid().density();
     rhoSolid =  M_oper.solid().rho( );
     theta  =  1./5. * (1-cos(50*PI*  M_oper.data().dataFluid()->dataTime()->timeStep()));
     dtheta =  10. *PI *sin(50*PI* M_oper.data().dataFluid()->dataTime()->timeStep());
     ddtheta = 500.* PI*PI * cos(50*PI* M_oper.data().dataFluid()->dataTime()->timeStep());
     dtheta2 = dtheta*dtheta;
     c1 = 0,   c2 =0,   c3  = 0.0,
     dc1 = 0,  dc2 = 0,  dc3 = 0.0,
     ddc1 = 0, ddc2 = 0, ddc3 = 0.0;
    }

  typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_Type;
  typedef FSIOperator::vector_Type                        vector_Type;
  typedef FSIOperator::vectorPtr_Type                     vectorPtr_Type;

Real fluidAccelerate(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
  switch(i) {
    case 0:
      return  dtheta2*( c1-x ) + ddtheta*(c2-y) + ddc1 ;
        break;
    case 1:
      return dtheta2*(c2-y)  + ddtheta*(x-c1)+ ddc2 ;
      break;
    case 2:
      return ddc3;
      break;
    }
    return 0;
}

Real
fluidRHS(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
 {
   switch(i) {
   case 0:
     return rhoFluid*( dtheta2*( c1-x ) + ddtheta*(c2-y) + ddc1 );
     break;
   case 1:
     return rhoFluid*(dtheta2*(c2-y)  + ddtheta*(x-c1)+ ddc2 );
     break;
   case 2:
     return rhoFluid*ddc3;
     break;
   default:
     ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
     break;
   }
}

//Real
//pressure(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
//{
//  return 2*(M_oper.solid().mu(flag)+M_oper.solid().lambda(flag) )*(1-cos(theta));
//}


Real
solidRHS(const Real& /*t*/, const Real& X, const Real& Y, const Real& /*Z*/, const ID& i)
{

  switch(i) {
  case 0:
    return - rhoSolid*(ddtheta * ( X*sin(theta) + Y*cos(theta) )+ dtheta2*( X*cos(theta) - Y*sin(theta))-ddc1);
    break;
  case 1:
    return rhoSolid*( ddtheta* (X*cos(theta)-Y*sin(theta) ) - dtheta2 * ( X*sin(theta)+Y*cos(theta))+ddc2);
    break;
  case 2:
    return rhoSolid*ddc3;
    break;
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }
}

Real
fluidVelocity(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
  switch(i) {
  case 0:
    return ( c2- y )* dtheta + dc1;
    break;
  case 1:
    return (x - c1) * dtheta + dc2;
    break;
    case 2:
      return dc3;
      break;
  }
  return 0;
}

Real
solidAccelerate(const Real& /*t*/, const Real& X, const Real& Y, const Real& /*Z*/, const ID& i)
{
  switch(i) {
  case 0:
    return - ddtheta*( X*sin(theta)+Y*cos(theta) ) + ddc1
      - dtheta2 *( X*cos(theta) - Y*sin(theta) )	;
    break;
  case 1:
    return ddtheta * ( X*cos(theta) - Y*sin(theta) ) + ddc2
      - dtheta2 *( X*sin(theta) + Y*cos(theta) )	;
      break;
  case 2:
    return dc3;
    break;
  default:
        ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
        break;
  }
}

// Initial velocity
Real
solidDisplacement(const Real& /*t*/, const Real&  X, const Real& Y, const Real& /*Z*/, const ID& i)
 {
  switch(i) {
   case 0:
     return X * ( cos(theta) - 1 ) - Y * sin( theta ) + c1;
     break;
   case 1:
     return  X * sin(theta) + Y * ( cos( theta ) - 1 ) + c2;
     break;
   case 2:
     return c3;
     break;
   default:
     ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
     break;
   }
}

Real
solidVelocity(const Real& /*t*/, const Real& X, const Real& Y, const Real& /*Z*/, const ID& i)
{
  switch(i) {
  case 0:
    return - dtheta*( X*sin(theta) + Y*cos(theta) ) + dc1;
        break;
  case 1:
    return dtheta *( X*cos(theta) - Y*sin(theta) ) + dc2;
    break;
  case 2:
        return dc3;
        break;
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }
}

fct_Type FluidAccelerate()
    {
        fct_Type f;
        f = boost::bind(&analyticalSolution::fluidAccelerate, this, _1, _2, _3, _4, _5);
        return f;
    }

fct_Type FluidRHS()
    {
        fct_Type f;
        f = boost::bind(&analyticalSolution::fluidRHS, this, _1, _2, _3, _4, _5);
        return f;
    }

const fct_Type& FluidVelocity()
    {
        fct_Type f;
        f = boost::bind(&analyticalSolution::fluidVelocity, this, _1, _2, _3, _4, _5);
        return f;
    }

fct_Type SolidAccelerate()
    {
        fct_Type f;
        f = boost::bind(&analyticalSolution::solidAccelerate, this, _1, _2, _3, _4, _5);
        return f;
    }

fct_Type SolidDisplacement()
    {
        fct_Type f;
        f = boost::bind(&analyticalSolution::solidDisplacement, this, _1, _2, _3, _4, _5);
        return f;
    }

fct_Type SolidVelocity()
    {
        fct_Type f;
        f = boost::bind(&analyticalSolution::solidVelocity, this, _1, _2, _3, _4, _5);
        return f;
    }

fct_Type SolidRHS()
    {
        fct_Type f;
        f = boost::bind(&analyticalSolution::solidRHS, this, _1, _2, _3, _4, _5);
        return f;
    }

void initializeFluid( std::vector <vectorPtr_Type>& fluidStart)
  {
      vectorPtr_Type velocity;
      vectorPtr_Type accelerate ( new vector_Type(M_oper.fluid().getMap()) );

      if( M_oper.fluidTimeAdvanceMethod() == "Newmark")
        {
            velocity.reset ( new vector_Type(M_oper.fluid().getMap()) );
	  // M_oper.uFESpace().interpolate(*this->fluidVelocity,   velocity,  tStart);
          //M_oper.uFESpace().interpolate(this->fluidAccelerate(), accelerate,  tStart);
          fluidStart.push_back(velocity);
          fluidStart.push_back(accelerate);
        }

      if(M_oper.fluidTimeAdvanceMethod() == "BDF")
      {
          for ( UInt previousPass = 0; previousPass < M_oper.data().dataFluid()->dataTime()->orderBDF() ; previousPass++)
          {
              Real previousTimeStep = -previousPass*  M_oper.data().dataFluid()->dataTime()->timeStep() ;
              velocity.reset ( new vector_Type(M_oper.fluid().getMap()) );
              // M_oper.uFESpace().interpolate(this->FluidVelocity(), velocity, previousTimeStep);
              fluidStart.push_back(velocity);
          }
      }
  }

void initializeALE( std::vector <vectorPtr_Type>& ALEStart )
  {
      vectorPtr_Type displacement ( new vector_Type(M_oper.meshMotion().getMap()) );
      vectorPtr_Type velocity( new vector_Type(M_oper.meshMotion().getMap()) );
      //vectorPtr_Type accelerate      ( new vector_Type(M_oper.meshMotion().getMap()) );

      if( M_oper.ALETimeAdvanceMethod() == "Newmark")
      {
          displacement.reset( new vector_Type(M_oper.meshMotion().getMap()) );
          // M_oper.meshMotion().mFESpace().interpolate(this->SolidDisplacement(),     displacement,  tStart);
          // M_oper.meshMotion().mFESpace().interpolate(this->SolidVelocity(), velocity,  tStart);
          //M_oper.meshMotion().mFESpace().interpolate(this->SolidVelocity(), accelerate,  tStart);
          ALEStart.push_back(displacement);
          ALEStart.push_back(velocity);
          //ALEStart.push_back(accelerate);
      }

      if(M_oper.ALETimeAdvanceMethod() == "BDF")
	{
	  for ( UInt previousPass = 0; previousPass < M_oper.data().dataFluid()->dataTime()->orderBDF() ; previousPass++)
	    {
          displacement.reset( new vector_Type(M_oper.meshMotion().getMap()) );
	      Real previousTimeStep = -previousPass*  M_oper.data().dataFluid()->dataTime()->timeStep() ;
          //M_oper.meshMotion().mFESpace().interpolate(SolidDisplacement(), displacement, previousTimeStep);
          ALEStart.push_back( displacement) ;
	    }
        }
}

void initializeSolid( std::vector <vectorPtr_Type>& solidStart )
  {
      vectorPtr_Type displacement ;
      vectorPtr_Type velocity     ( new vector_Type(M_oper.solid().map()) );
      vectorPtr_Type accelerate   ( new vector_Type(M_oper.solid().map()) );

      if( M_oper.solidTimeAdvanceMethod() == "Newmark")
        {
            displacement.reset ( new vector_Type(M_oper.solid().map()) );

            //M_oper.meshMotion().mFESpace().interpolate(SolidDisplacement(),     displacement,  tStart);
            //M_oper.meshMotion().mFESpace().interpolate(SolidVelocity(), velocity,  tStart);
            //M_oper.meshMotion().mFESpace().interpolate(SolidVelocity(), accelerate,  tStart);
            solidStart.push_back(displacement);
            solidStart.push_back(velocity);
            solidStart.push_back(accelerate);
        }

      if(M_oper.solidTimeAdvanceMethod() == "BDF")
	{
	   for ( UInt previousPass ; previousPass < M_oper.data().dataSolid()->dataTime()->orderBDF() ; previousPass++)
	    {
            displacement.reset ( new vector_Type(M_oper.solid().map()) );
            Real previousTimeStep = -previousPass*  M_oper.data().dataSolid()->dataTime()->timeStep() ;
            //  M_oper.meshMotion().mFESpace().interpolate(this->SolidDisplacement(), displacement,previousTimeStep);
            solidStart.push_back( displacement) ;
	    }
	}

}


private:

    FSIOperator& M_oper;
    Real  flag;
    Real  rhoFluid, rhoSolid;
    Real  theta,  dtheta,  ddtheta,  dtheta2;
    Real  c1,   c2,   c3,
         dc1,  dc2,  dc3,
        ddc1, ddc2, ddc3;

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
    Problem( const std::string& dataFileName, std::string method = "" )
    {

      //VenantKirchhoffSolver< FSIOperator::mesh_Type, SolverAztecOO >::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhof", &createLinearStructure );


      //StructuralSolver< FSIOperator::mesh_Type, SolverAztecOO >::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhoff", &createVenantKirchhoffLinear );

      StructuralSolver< FSIOperator::mesh_Type, SolverAztecOO >();

      //StructuralSolver< FSIOperator::mesh_Type, SolverAztecOO >::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhoff", &createVenantKirchhoffLinear );

        //        VenantKirchhofSolver< FSIOperator::mesh_Type, SolverAztecOO >::StructureSolverFactory::instance().registerProduct( "nonLinearVenantKirchhof", &createNonLinearStructure );

        Debug( 10000 ) << "Setting up data from GetPot \n";
        GetPot dataFile( dataFileName );
        M_data = dataPtr_Type( new data_Type() );
        M_data->setup( dataFile );
        M_data->dataSolid()->setTimeData( M_data->dataFluid()->dataTime() ); //Same TimeData for fluid & solid
        M_data->showMe();
	M_data->dataSolid()->showMe();
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

        MPI_Barrier( MPI_COMM_WORLD );

        Debug( 10000 ) << "Setting up the problem \n";
        M_fsi->setup( );

        //M_fsi->resetFSISolvers();

        MPI_Barrier( MPI_COMM_WORLD );

        std::string const exporterType =  dataFile( "exporter/type", "hdf5");
        std::string const exporterName =  dataFile( "exporter/name", "fixedPt");

        Debug( 10000 ) << "Setting up ExporterEnsight \n";
        if ( M_fsi->isFluid() )
        {
#ifdef HAVE_HDF5
            if (exporterType.compare("hdf5") == 0)
                M_exporterFluid.reset( new ExporterHDF5<mesh_Type > ( dataFile, exporterName+"Fluid" ) );
            else
#endif
                M_exporterFluid.reset( new ExporterEnsight<mesh_Type > ( dataFile, exporterName+"Fluid" ) );


            M_exporterFluid->setMeshProcId(M_fsi->FSIOper()->uFESpace().mesh(), M_fsi->FSIOper()->uFESpace().map().comm().MyPID());

            M_velAndPressure.reset( new vector_Type( M_fsi->FSIOper()->fluid().getMap(),      M_exporterFluid->mapType() ));
            M_fluidDisp.reset     ( new vector_Type( M_fsi->FSIOper()->meshMotion().getMap(), M_exporterFluid->mapType() ));

            M_exporterFluid->addVariable( ExporterData<mesh_Type>::VectorField, "f-velocity",
                                          M_fsi->FSIOper()->uFESpacePtr(), M_velAndPressure, UInt(0) );

            M_exporterFluid->addVariable( ExporterData<mesh_Type>::ScalarField, "f-pressure",
                                          M_fsi->FSIOper()->pFESpacePtr(), M_velAndPressure,
                                          UInt(3*M_fsi->FSIOper()->uFESpace().dof().numTotalDof() ) );

            M_exporterFluid->addVariable( ExporterData<mesh_Type>::VectorField, "f-displacement",
                                          M_fsi->FSIOper()->mmFESpacePtr(), M_fluidDisp, UInt(0) );

        }
        if ( M_fsi->isSolid() )
        {
#ifdef HAVE_HDF5
            if (exporterType.compare("hdf5") == 0)
                M_exporterSolid.reset( new ExporterHDF5<mesh_Type > ( dataFile, exporterName+"Solid" ) );
            else
#endif
                M_exporterSolid.reset( new ExporterEnsight<mesh_Type > ( dataFile, exporterName+"Solid" ) );

            M_exporterSolid->setMeshProcId(M_fsi->FSIOper()->dFESpace().mesh(), M_fsi->FSIOper()->dFESpace().map().comm().MyPID());
	    M_solidDisp.reset( new vector_Type( M_fsi->FSIOper()->solid().map(), M_exporterSolid->mapType() ));
            M_solidVel.reset ( new vector_Type( M_fsi->FSIOper()->solid().map(), M_exporterSolid->mapType() ));

	    M_exporterSolid->addVariable( ExporterData<mesh_Type>::VectorField, "s-displacement",
                                          M_fsi->FSIOper()->dFESpacePtr(), M_solidDisp, UInt(0) );
            M_exporterSolid->addVariable( ExporterData<mesh_Type>::VectorField, "s-velocity",
                                          M_fsi->FSIOper()->dFESpacePtr(), M_solidVel, UInt(0) );

        }

        bool restart = dataFile("problem/restart",false);
        M_Tstart = 0.;

	analyticalSolution M_analyticalSolution(*M_fsi->FSIOper());

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
               M_exporterSolid->import(M_Tstart, M_data->dataSolid()->dataTime()->timeStep());
               M_fsi->FSIOper()->initializeSolid( M_solidDisp, M_solidVel );
            }
        }
        else
        {
            std::vector<vectorPtr_Type> fluidStart(0), ALEStart(0), solidStart(0);

//             M_analyticalSolution.initializeFluid(fluidStart);
//             M_analyticalSolution.initializeALE(ALEStart);
//             M_analyticalSolution.initializeSolid(solidStart);

            M_fsi->initialize(/*fluidStart, ALEStart, solidStart*/);
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

            if ( M_fsi->isFluid() )
            {
              BCFunctionBase meshDisplacement;
	      // meshDisplacement.setFunction(M_analyticalSolution.SolidDisplacement());

	      M_fsi->FSIOper()->BCh_harmonicExtension()->modifyBC(3,  meshDisplacement);
	      M_fsi->FSIOper()->BCh_harmonicExtension()->modifyBC(2,  meshDisplacement);
	      M_fsi->FSIOper()->BCh_harmonicExtension()->modifyBC(20, meshDisplacement);
	      M_fsi->FSIOper()->BCh_harmonicExtension()->modifyBC(30, meshDisplacement);

              BCFunctionBase fluidVelocity;
	      //fluidVelocity.setFunction(M_analyticalSolution.FluidVelocity());

	      M_fsi->FSIOper()->BCh_fluid()->modifyBC(3,  fluidVelocity);
	      M_fsi->FSIOper()->BCh_fluid()->modifyBC(2,  fluidVelocity);
	      M_fsi->FSIOper()->BCh_fluid()->modifyBC(20, fluidVelocity);
	      M_fsi->FSIOper()->BCh_fluid()->modifyBC(30, fluidVelocity);
            }

            if ( M_fsi->isFluid() )
            {
              BCFunctionBase solidDisplacement;
 	      //solidDisplacement.setFunction(M_analyticalSolution.SolidDisplacement());

	      M_fsi->FSIOper()->BCh_solid()->modifyBC(3,  solidDisplacement);
	      M_fsi->FSIOper()->BCh_solid()->modifyBC(2,  solidDisplacement);
	      M_fsi->FSIOper()->BCh_solid()->modifyBC(20, solidDisplacement);
	      M_fsi->FSIOper()->BCh_solid()->modifyBC(30, solidDisplacement);
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
                *M_solidDisp = M_fsi->FSIOper()->solid().displacement();
		// *M_solidVel = M_fsi->FSIOper()->solid().velocity();
		*M_solidVel = M_fsi->FSIOper()->solidTimeAdvance()->velocity();

                M_exporterSolid->postProcess( M_data->dataFluid()->dataTime()->time() );
            }

            std::cout << "[fsi_run] Iteration " << _i << " was done in : " << _timer.elapsed() << "\n";

            std::cout << "solution norm " << _i << " : "
                      << M_fsi->displacement().norm2() << "\n";

            // CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
            //checkResult( M_data->dataFluid()->dataTime()->time() );
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
        const LifeV::Real maxAbs (std::max(std::fabs(a),std::fabs(b)));
        if (maxAbs < relTol*relTol) return true;

        return std::fabs(a-b) < relTol*maxAbs;
    }

    fsi_solver_ptr M_fsi;

    dataPtr_Type   M_data;
    Real           M_Tstart;

    bool               M_absorbingBC;

  //analyticalSolution M_analyticalSolution (*M_fsi->FSIOper()) ;

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

    LifeChrono chrono;
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
