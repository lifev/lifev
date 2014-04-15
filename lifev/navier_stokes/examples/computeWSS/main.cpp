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
   \author Paolo Tricerri <paolo.tricerri@epfl.ch>
   \date 2012-05-01
 */
#undef HAVE_HDF5
#ifdef TWODIM
#error example_computeWSS cannot be compiled in 2D
#endif

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>

//Vectors and related things
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>

// Mesh infos
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

// FESpaces and ETFESpaces
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

//Oseen data to make use of the fields that are available there.
//No real need for this.
#include <lifev/navier_stokes/solver/OseenData.hpp>

// Evaluation operations
#include <lifev/eta/expression/Evaluate.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <iostream>

using namespace LifeV;

int returnValue = EXIT_SUCCESS;

std::set<UInt> parseList ( const std::string& list )
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
        commaPos = stringList.find ( "," );
        setList.insert ( atoi ( stringList.substr ( 0, commaPos ).c_str() ) );
        stringList = stringList.substr ( commaPos + 1 );
    }
    setList.insert ( atoi ( stringList.c_str() ) );
    return setList;
}

class WSS
{
public:

    typedef LifeV::RegionMesh<LinearTetra>                              mesh_Type;

    // Filters
    typedef typename LifeV::Exporter<mesh_Type  >                       filter_Type;
    typedef boost::shared_ptr< LifeV::Exporter<mesh_Type  > >           filterPtr_Type;

    typedef LifeV::ExporterEmpty<mesh_Type >                            emptyFilter_Type;
    typedef boost::shared_ptr<emptyFilter_Type>                         emptyFilterPtr_Type;
    typedef LifeV::ExporterEnsight<mesh_Type >                          ensightFilter_Type;
    typedef boost::shared_ptr<ensightFilter_Type>                       ensightFilterPtr_Type;

#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<mesh_Type >                             hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                          hdf5FilterPtr_Type;
#endif



    /** @name Constructors, destructor
     */
    //@{
    WSS ( int                                   argc,
                char**                                argv,
                boost::shared_ptr<Epetra_Comm>        structComm );

    ~WSS()
    {}
    //@}

    //@{
    void run()
    {
        run3d();
    }
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
};



struct WSS::Private
{
    Private()
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

};



WSS::WSS ( int                                   argc,
           char**                                argv,
           boost::shared_ptr<Epetra_Comm>        structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID() )
    {
        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}



void
WSS::run2d()
{
    std::cout << "2D WSS example is not available yet\n";
}



void
WSS::run3d()
{
    typedef VectorEpetra                                                 vector_Type;
    typedef boost::shared_ptr<VectorEpetra>                              vectorPtr_Type;

    typedef FESpace< mesh_Type, MapEpetra >                              wssFESpace_Type;
    typedef boost::shared_ptr<wssFESpace_Type>                           wssFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >        vectorialETFESpace_Type;
    typedef boost::shared_ptr<vectorialETFESpace_Type>                   vectorialETFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >        scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;

    bool verbose = (parameters->comm->MyPID() == 0);

    //! dataElasticStructure for parameters
    GetPot dataFile ( parameters->data_file_name.c_str() );
    boost::shared_ptr<OseenData> dataClass (new OseenData( ) );
    dataClass->setup (dataFile);

    dataClass->showMe();

    //! Read and partition mesh
    MeshData             meshData;
    meshData.setup (dataFile, "fluid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr (new RegionMesh<LinearTetra> ( ( parameters->comm ) ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, parameters->comm );

    wssFESpacePtr_Type velFESpace ( new wssFESpace_Type (meshPart, dataClass->uOrder() , 3, parameters->comm) );

    // This FESpace is to have the nice quadrature rule
    wssFESpacePtr_Type quadFESpace ( new wssFESpace_Type (meshPart, dataClass->uOrder() , 2, parameters->comm) );

    vectorialETFESpacePtr_Type uETFESpace ( new vectorialETFESpace_Type (meshPart, & (velFESpace->refFE() ),
                                                                         & (velFESpace->fe().geoMap() ), parameters->comm) );

    //! 3. Creation of the importers to read the solution (u,p)
    std::string const filename    = dataFile ( "importer/filename", "noNameFile");
    std::string const importerType = dataFile ( "importer/type", "hdf5");

    // How many solution do we have to read?
    // If you want to read a specific solution, choose instant
    // If you want to read a set of solutions, choose interval.
    std::string readType = dataFile ( "importer/analysis", "instant");
    UInt numberOfSol(0);
    UInt start(0);
    UInt end(0);

    if( !readType.compare("instant") )
      {
          numberOfSol = dataFile.vector_variable_size ( "importer/iteration"  );
          ASSERT( numberOfSol, "You did not set any solution to read!! ");
      }
    else
      {
          start = dataFile ( "importer/iterationStart" , 0 );
          end = dataFile ( "importer/iterationEnd" , 0 );
          numberOfSol = end - start;
      }
    ASSERT( numberOfSol, "You did not set any solution to read!! ");

    if (verbose)
    {
        std::cout << "The filename is    : " << filename << std::endl;
        std::cout << "The importerType is: " << importerType << std::endl;
    }

#ifdef HAVE_HDF5
    if (importerType.compare ("hdf5") == 0)
    {
        M_importer.reset ( new hdf5Filter_Type (dataFile, filename) );
    }
    else
#endif
    {
        if (importerType.compare ("none") == 0)
        {
            M_importer.reset ( new emptyFilter_Type ( dataFile, velFESpace->mesh(), "WSS", velFESpace->map().comm().MyPID() ) );
        }
        else
        {
            M_importer.reset ( new ensightFilter_Type ( dataFile, filename ) );
        }
    }
    M_importer->setMeshProcId (velFESpace->mesh(), velFESpace->map().comm().MyPID() );

    //! 6. Post-processing setting
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "hdf5");
    std::string const nameExporter =  dataFile ( "exporter/name", "jacobian");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new hdf5Filter_Type ( dataFile, nameExporter) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new emptyFilter_Type ( dataFile, meshPart.meshPartition(), nameExporter, parameters->comm->MyPID() ) ) ;
        }

        else
        {
            exporter.reset ( new ensightFilter_Type ( dataFile, meshPart.meshPartition(), nameExporter, parameters->comm->MyPID() ) ) ;
        }
    }
    exporter->setMeshProcId (velFESpace->mesh(), velFESpace->map().comm().MyPID() );

    // Scalar vector to have scalar quantities
    vectorPtr_Type patchAreaVector;

    vectorPtr_Type velocity;
    vectorPtr_Type velocityRead;

    vectorPtr_Type WSSvector;

    patchAreaVector.reset ( new vector_Type ( uETFESpace->map() ) );
    velocity.reset( new vector_Type( velFESpace->map() ) );
    velocityRead.reset( new vector_Type( velFESpace->map() ) );

    WSSvector.reset ( new vector_Type ( uETFESpace->map() ) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "WSS",
                            velFESpace, WSSvector, UInt (0) );

    // Debug purposes
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocityRead",
                              velFESpace, velocityRead, UInt (0) );

    exporter->postProcess ( 0.0 );

    //! =================================================================================
    //! Analysis - Istant or Interval
    //! =================================================================================

    MPI_Barrier (MPI_COMM_WORLD);

    QuadratureRule fakeQuadratureRule( quadRuleTria3pt );

    GeoVector firstPoint(2); Real firstWeight(0.0);
    firstPoint[0] = 0.0; firstPoint[1] = 0.0;; firstWeight = ( 0.5 ) / 3;
    QuadraturePoint first( firstPoint, firstWeight );

    GeoVector secondPoint(2); Real secondWeight(0.0);
    secondPoint[0] = 1.0; secondPoint[1] = 0.0; secondWeight = ( 0.5 ) / 3;
    QuadraturePoint second( secondPoint, secondWeight );

    GeoVector thirdPoint(2); Real thirdWeight(0.0);
    thirdPoint[0] = 0.0; thirdPoint[1] = 1.0; thirdWeight = ( 0.5 ) / 3;
    QuadraturePoint third( thirdPoint, thirdWeight );

    std::vector<GeoVector> coordinates( 3 );
    coordinates[0] = firstPoint;
    coordinates[1] = secondPoint;
    coordinates[2] = thirdPoint;
    std::vector<Real> weights ( 3 , 1.0 / 6.0 );

    fakeQuadratureRule.setPoints(coordinates, weights);
    QuadratureBoundary adaptedBDQuadRule( buildTetraBDQR( fakeQuadratureRule ) );

    using namespace ExpressionAssembly;

    // Trying to compute the Jacobian using ET
    MatrixSmall<3,3> identity;
    identity (0, 0) = 1.0;
    identity (0, 1) = 0.0;
    identity (0, 2) = 0.0;
    identity (1, 0) = 0.0;
    identity (1, 1) = 1.0;
    identity (1, 2) = 0.0;
    identity (2, 0) = 0.0;
    identity (2, 1) = 0.0;
    identity (2, 2) = 1.0;

    ExpressionVectorFromNonConstantScalar<ExpressionMeasBDCurrentFE, 3  > vMeas( meas_BDk );
    evaluateNode( boundary( uETFESpace->mesh(), 1 ),
                  adaptedBDQuadRule,
                  uETFESpace,
                  dot( vMeas , phi_i )
                  ) >> patchAreaVector;
    patchAreaVector->globalAssemble();

    std::string const nameField =  dataFile ( "importer/nameField", "f-velocity");
    UInt i,k;

    for( i=0,k =0; i < numberOfSol; i++, k++ )
    {

        // Reading the solution
        // resetting the element of the vector
        *velocity *= 0.0;
	*WSSvector *= 0.0;
        UInt current(0);
        if( !readType.compare("interval") )
        {
            current = i + start;
        }
        else
        {
            current = dataFile ( "importer/iteration" , 100000, i );
        }
        // Transform current ingot a string
        std::string iterationString;
        std::ostringstream number;
        number.fill ( '0' );
        number << std::setw (5) << ( current );
        iterationString = number.str();

        std::cout << "Current reading: " << iterationString << std::endl;

        /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
        LifeV::ExporterData<mesh_Type> solutionVel  (LifeV::ExporterData<mesh_Type>::VectorField,nameField + "." + iterationString,
                                                     velFESpace, velocity, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );


        //Read the variable
        M_importer->readVariable (solutionVel);
        *velocityRead = *velocity;

        Real dynamicViscosity = dataClass->viscosity( );

	std::cout << "Viscosity: " << dataClass->viscosity( ) << std::endl;
        // Defining expressions
        evaluateNode( boundary( uETFESpace->mesh(), 1 ),
                      adaptedBDQuadRule,
                      uETFESpace,
		      meas_BDk * dot ( value(dynamicViscosity) * ( grad(uETFESpace, *velocity) + transpose(grad(uETFESpace, *velocity)) ) * Nface - 
				       dot( value(dynamicViscosity) * ( grad(uETFESpace, *velocity) + transpose(grad(uETFESpace, *velocity)) ) * Nface , Nface ) * Nface, phi_i )
                      ) >> WSSvector;
        WSSvector->globalAssemble();
        *WSSvector = *WSSvector / *patchAreaVector;

        // // Defining expressions
        // integrate( boundary( uETFESpace->mesh(), 1 ),
	// 	   adaptedBDQuadRule,
	// 	   uETFESpace,
	// 	   dot ( Nface, phi_i )
	// 	   ) >> WSSvector;
        // WSSvector->globalAssemble();


        exporter->postProcess( dataClass->dataTime()->initialTime() + k * dataClass->dataTime()->timeStep() );
    }


    if (verbose )
    {
        std::cout << "Analysis Completed!" << std::endl;
    }

    //Closing files
    exporter->closeFile( );
    M_importer->closeFile( );

    MPI_Barrier (MPI_COMM_WORLD);
    //!---------------------------------------------.-----------------------------------------------------
}

int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }
#else
    boost::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif

    WSS wss ( argc, argv, Comm );
    wss.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
