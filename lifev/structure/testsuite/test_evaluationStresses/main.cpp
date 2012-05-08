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
/**
   \file main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-16
 */
#undef HAVE_HDF5
#ifdef TWODIM
#error test_structure cannot be compiled in 2D
#endif

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

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/structure/solver/VenantKirchhoffElasticData.hpp>

#include <lifev/structure/solver/StructuralMaterial.hpp>
#include <lifev/structure/solver/StructuralSolver.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/ExponentialMaterialNonLinear.hpp>

#include <lifev/structure/solver/WallTensionEstimator.hpp>
#include <lifev/structure/solver/WallTensionEstimatorData.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <iostream>


using namespace LifeV;

int returnValue = EXIT_SUCCESS;

std::set<UInt> parseList( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    size_t commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );
        setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
        stringList = stringList.substr( commaPos+1 );
    }
    setList.insert( atoi( stringList.c_str() ) );
    return setList;
}


class Structure
{
public:

  typedef LifeV::RegionMesh<LinearTetra>                                       mesh_Type;

  // Filters
  typedef boost::shared_ptr< LifeV::Exporter<mesh_Type  > >                filterPtr_Type;
    
  typedef LifeV::ExporterEmpty<mesh_Type >    emptyFilter_Type;
  typedef boost::shared_ptr<emptyFilter_Type>                              emptyFilterPtr_Type;                         
  typedef LifeV::ExporterEnsight<mesh_Type >  ensightFilter_Type;
  typedef boost::shared_ptr<ensightFilter_Type>                           ensightFilterPtr_Type;

#ifdef HAVE_HDF5
  typedef LifeV::ExporterHDF5<mesh_Type >     hdf5Filter_Type;
  typedef boost::shared_ptr<hdf5Filter_Type>                              hdf5FilterPtr_Type;
#endif



/** @name Constructors, destructor
 */
//@{
    Structure( int                                   argc,
               char**                                argv,
               boost::shared_ptr<Epetra_Comm>        structComm );

    ~Structure()
    {}
//@}

//@{
    void run()
    {
        run3d();
    }
    void CheckResult(const Real& dispNorm, const Real& time);
    void resultChanged(Real time);
//@}

protected:

private:

    /**
     * run the 2D driven cylinder simulation
     */
    void run2d();

    /**
     * run the 3D driven cylinder simulation
     */
    void run3d();

private:
    struct Private;
    boost::shared_ptr<Private> parameters;
    filterPtr_Type M_importer;
    filterPtr_Type M_exporter;
};



struct Structure::Private
{
    Private() :
            rho(1), 
	    poisson(1), 
	    young(1), 
	    bulk(1), 
	    alpha(1), 
	    gamma(1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

};



Structure::Structure( int                                   argc,
                      char**                                argv,
                      boost::shared_ptr<Epetra_Comm>        structComm):
        	      parameters( new Private() )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->rho     = dataFile( "solid/physics/density", 1. );
    parameters->young   = dataFile( "solid/physics/young",   1. );
    parameters->poisson = dataFile( "solid/physics/poisson", 1. );
    parameters->bulk    = dataFile( "solid/physics/bulk",    1. );
    parameters->alpha   = dataFile( "solid/physics/alpha",   1. );
    parameters->gamma   = dataFile( "solid/physics/gamma",   1. );

    std::cout << "density = " << parameters->rho     << std::endl
              << "young   = " << parameters->young   << std::endl
              << "poisson = " << parameters->poisson << std::endl
              << "bulk    = " << parameters->bulk    << std::endl
              << "alpha   = " << parameters->alpha   << std::endl
              << "gamma   = " << parameters->gamma   << std::endl;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID()) std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
}



void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}



void
Structure::run3d()
{
    // General typedefs
    typedef WallTensionEvaluator< mesh_Type >::solutionVect_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                                vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >                 solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                          solidFESpacePtr_Type;    


    bool verbose = (parameters->comm->MyPID() == 0);

    //! dataElasticStructure for parameters
    GetPot dataFile( parameters->data_file_name.c_str() );

    boost::shared_ptr<VenantKirchhoffElasticData> dataStructure(new VenantKirchhoffElasticData( ));
    dataStructure->setup(dataFile);

    //! Parameters for the analysis
    boost::shared_ptr<WallTensionEstimatorData> tensionData(new WallTensionEstimatorData( ));
    tensionData->setup(dataFile);

    
    //! Read and partition mesh
    MeshData             meshData;
    meshData.setup(dataFile, "solid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr(new mesh_Type);
    readMesh(*fullMeshPtr, meshData);

    MeshPartitioner< mesh_Type > meshPart( fullMeshPtr, parameters->comm );

    //! Functional spaces - needed for the computations of the gradients
    std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace( new solidFESpace_Type(meshPart,dOrder,3,parameters->comm) );
    if (verbose) std::cout << std::endl;

    //! 1. Constructor of the structuralSolver
    WallTensionEvaluator< mesh_Type > solid;

    //! 2. Setup of the structuralSolver
    solid.setup(dataStructure,
    		tensionData,
                dFESpace,
                parameters->comm);

    //! 3. Creation of the importers to read the displacement field
    std::string const importerType =  data_file( "importer/type", "hdf5");
    std::string const filename    = data_file( "importer/fluid/filename", "solid");
    std::string iterationString; //useful to iterate over the strings

    if (verbose) 
      {
	std::cout << "The filename is    : " << filename << std::endl;
	std::cout << "The importerType is: " << importerType << std::endl;
      }

#ifdef HAVE_HDF5
    if (importerType.compare("hdf5") == 0)
      {
	M_importer.reset( new  hdf5Filter_Type( data_file, filename) );
      }
    else
#endif
      {
	if (importerType.compare("none") == 0)
	  {
	    M_importer.reset( new ExporterEmpty<mesh_Type > ( data_file, solid.dFESpace().mesh(), "solid", solid.dFESpace().map().comm().MyPID()) );
	  }
	else
	  {
	    M_importer.reset( new  ensightFilter_Type ( data_file, fileName) );
	  }
      }

    // The vector where the solution will be stored
    vectorPtr_Type solidDisp (new vector_Type(solid.dFESpace().map(),M_importer->mapType() ));

    //! 4. Creation of the expoters to save the tensions




    //! =================================================================================
    //! Analysis - Istant or Interval
    //! =================================================================================

    //! 5. For each interval, the analysis is performed

    //Cycle of the intervals
    for ( UInt i(0); i< tensionData->iterStart().size(); i++ )
      {
	std::string initial = tensionData->iterStart(i);
	std::string last    = tensionData->iterEnd(i);

	//! Todo: A check on the existence of the initial and last strings would be appreciated
	
	iterationString = initial;

	while( iterationString.compare(last)  )
	  {

	    //Read variable at iterationString

	    /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
	    LifeV::ExporterData<mesh_Type> solidDisp  (LifeV::ExporterData<mesh_Type>::VectorField,"s-displacement."+iterationString, solid.dFESpacePtr(), solidDisp, UInt(0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );



	    //Increment the iterationString
	    
	  }
	
      }


    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose ) std::cout << "ok." << std::endl;


    //! 6. Post-processing setting
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
      exporter.reset( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    }
    else
#endif
    {
        if (exporterType.compare("none") == 0)
	{
	    exporter.reset( new ExporterEmpty<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID()) );
	}

        else
        {
	    exporter.reset( new ExporterEnsight<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID()) );
	}
    }

    exporter->setPostDir( "./" );
    exporter->setMeshProcId( meshPart.meshPartition(), parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type(solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidVel  ( new vector_Type(solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidAcc  ( new vector_Type(solid.displacement(),  exporter->mapType() ) );

    exporter->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt(0) );
    exporter->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt(0) );
    exporter->addVariable( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt(0) );

    exporter->postProcess( 0 );


    
    //! 7. For each of the iterations (which correspond to a time) the analysis is performed.
    vectorPtr_Type rhs (new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type disp(new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type vel (new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type acc (new vector_Type(solid.displacement(), Unique));

    //! 8. Loop on the iterations

    exporter->postProcess( time ); //a time is needed!


	Real normVect;
	normVect =  solid.displacement().norm2();
	std::cout << "The norm 2 of the displacement field is: "<< normVect << std::endl;



	//!--------------------------------------------------------------------------------------------------

        MPI_Barrier(MPI_COMM_WORLD);
    }
}



void Structure::CheckResultLE(const Real& dispNorm,const Real& time)
{
    if ( time == 0.1  && std::fabs(dispNorm-0.276527)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.2  && std::fabs(dispNorm-0.276536)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.3  && std::fabs(dispNorm-0.276529)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.4  && std::fabs(dispNorm-0.276531)>1e-5 )
        this->resultChanged(time);
}

void Structure::CheckResultSVK(const Real& dispNorm,const Real& time)
{
    if ( time == 0.1  && std::fabs(dispNorm-0.263348)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.2  && std::fabs(dispNorm-0.263350)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.3  && std::fabs(dispNorm-0.263350)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.4  && std::fabs(dispNorm-0.263351)>1e-5 )
        this->resultChanged(time);
}
void Structure::CheckResultEXP(const Real& dispNorm,const Real& time)
{
    if ( time == 0.1  && std::fabs(dispNorm-0.284844)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.2  && std::fabs(dispNorm-0.284853)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.3  && std::fabs(dispNorm-0.284846)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.4  && std::fabs(dispNorm-0.284848)>1e-5 )
        this->resultChanged(time);
}

void Structure::CheckResultNH(const Real& dispNorm,const Real& time)
{
    if ( time == 0.1  && std::fabs(dispNorm-0.286120)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.2  && std::fabs(dispNorm-0.286129)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.3  && std::fabs(dispNorm-0.286122)>1e-5 )
        this->resultChanged(time);
    if ( time == 0.4  && std::fabs(dispNorm-0.286123)>1e-5 )
        this->resultChanged(time);
}



void Structure::resultChanged(Real time)
{
  std::cout << "Some modifications led to changes in the l2 norm of the solution at time " << time << std::endl;
  returnValue = EXIT_FAILURE;
}



int
main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> Comm(new Epetra_MpiComm( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
        cout << "% using MPI" << endl;
#else
    boost::shared_ptr<Epetra_SerialComm> Comm( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif

    Structure structure( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
