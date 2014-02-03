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


Attention: At the moment the restart works only if the solution is saved at
each time step and with the BDF method!!

\date 2005-04-16
*/
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
//Include fils which were in the structure.cpp file
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <iostream>
#include <sstream>


using namespace LifeV;

int returnValue = EXIT_FAILURE;

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


class Structure
{
public:

    typedef RegionMesh<LifeV::LinearTetra >                       mesh_Type;
    typedef VectorEpetra                                          vector_Type;
    typedef boost::shared_ptr<vector_Type>                        vectorPtr_Type;

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    //Exporters Typedefs
    typedef LifeV::Exporter<mesh_Type >                           filter_Type;
    typedef boost::shared_ptr<filter_Type >                       filterPtr_Type;

    typedef LifeV::ExporterEmpty<mesh_Type >                      emptyExporter_Type;
    typedef boost::shared_ptr<emptyExporter_Type>                 emptyExporterPtr_Type;

    typedef LifeV::ExporterEnsight<mesh_Type >                    ensightFilter_Type;
    typedef boost::shared_ptr<ensightFilter_Type>                 ensightFilterPtr_Type;

#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<mesh_Type >                       hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                    hdf5FilterPtr_Type;
#endif

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >         solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                  solidFESpacePtr_Type;

    /** @name Constructors, destructor
     */
    //@{
    Structure ( int                                   argc,
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
    filterPtr_Type importerSolid;

};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

    static Real uexact (const Real& /*t*/, const Real&  X, const Real& Y, const Real& /*Z*/, const ID& i)
    {
        //Setting up the data of the problem
        Real E = 8e+6;
        Real poi = 0.45;
        Real mu = E / ( 2 * ( 1 + poi ) );
        Real lambda = ( E * poi ) / ( ( 1 + poi ) * ( 1 - 2 * poi ) );
        Real Rin = 0.5;
        Real Rout = 0.6;
        Real Pin = 1000.0;
        Real Pout = 0;


        // Defining the new variables
        Real radius = std::sqrt ( X * X + Y * Y  );
        Real theta = std::atan ( Y / X );

        switch (i)
        {
            case 0:
                //u_x = cos(theta) * u_r(r)
                return std::cos (theta) * ( ( ( radius / (2.0 * (mu + lambda) ) ) * ( ( Rin * Rin * Pin - Rout * Rout * Pout ) / ( Rout * Rout - Rin * Rin ) ) )
                                            + ( ( (Rin * Rin * Rout * Rout ) / ( 2 * mu * radius)  ) * ( ( Pin - Pout ) / ( Rout * Rout - Rin * Rin ) ) ) );

                break;
            case 1:
                //u_x = sin(theta) * u_r(r)
                return std::sin (theta) * ( ( ( radius / (2.0 * (mu + lambda) ) ) * ( ( Rin * Rin * Pin - Rout * Rout * Pout ) / ( Rout * Rout - Rin * Rin ) ) )
                                            + ( ( (Rin * Rin * Rout * Rout ) / ( 2 * mu * radius)  ) * ( ( Pin - Pout ) / ( Rout * Rout - Rin * Rin ) ) ) );

                break;
            case 2:
                return 0.0;
                break;
            default:
                ERROR_MSG ("This entry is not allowed!!");
                return 0;
                break;
        }

        return  0.;
    }

};



Structure::Structure ( int                                   argc,
                       char**                                argv,
                       boost::shared_ptr<Epetra_Comm>        structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    parameters->data_file_name = data_file_name;

    // parameters->rho     = dataFile( "solid/physics/density", 1. );
    // parameters->young   = dataFile( "solid/physics/young",   1. );
    // parameters->poisson = dataFile( "solid/physics/poisson", 1. );
    // parameters->bulk    = dataFile( "solid/physics/bulk",    1. );
    // parameters->alpha   = dataFile( "solid/physics/alpha",   1. );
    // parameters->gamma   = dataFile( "solid/physics/gamma",   1. );

    // std::cout << "density = " << parameters->rho     << std::endl
    //           << "young   = " << parameters->young   << std::endl
    //           << "poisson = " << parameters->poisson << std::endl
    //           << "bulk    = " << parameters->bulk    << std::endl
    //           << "alpha   = " << parameters->alpha   << std::endl
    //           << "gamma   = " << parameters->gamma   << std::endl;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID() )
    {
        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}



void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}



void
Structure::run3d()
{

    bool verbose = (parameters->comm->MyPID() == 0);

    //! Number of boundary conditions for the velocity and mesh motion
    boost::shared_ptr<BCHandler> BCh ( new BCHandler() );

    //! dataElasticStructure
    GetPot dataFile ( parameters->data_file_name.c_str() );

    //Loading a partitoned mesh or reading a new one
    const std::string partitioningMesh = dataFile ( "partitioningOffline/loadMesh", "no");

    //Creation of pointers
    boost::shared_ptr<MeshPartitioner<mesh_Type> > meshPart;
    boost::shared_ptr<mesh_Type> pointerToMesh;

#ifdef LIFEV_HAS_HDF5
    if ( ! (partitioningMesh.compare ("no") ) )
#endif
    {
        boost::shared_ptr<mesh_Type > fullMeshPtr (new mesh_Type ( ( parameters->comm ) ) );
        //Creating a new mesh from scratch
        MeshData             meshData;
        meshData.setup (dataFile, "solid/space_discretization");
        readMesh (*fullMeshPtr, meshData);

        meshPart.reset ( new MeshPartitioner<mesh_Type> (fullMeshPtr, parameters->comm ) );

        pointerToMesh = meshPart->meshPartition();
    }
#ifdef LIFEV_HAS_HDF5
    else
    {
        //Creating a mesh object from a partitioned mesh
        const std::string partsFileName (dataFile ("partitioningOffline/hdf5_file_name", "NO_DEFAULT_VALUE.h5") );

        boost::shared_ptr<Epetra_MpiComm> mpiComm =
            boost::dynamic_pointer_cast<Epetra_MpiComm>(parameters->comm);
        PartitionIO<mesh_Type> partitionIO (partsFileName, mpiComm);

        partitionIO.read (pointerToMesh);

    }
#endif

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");

    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 3, parameters->comm) );

    // setting precise quadrature rule for fine meshes
    const QuadratureRule fineQuadRule = quadRuleTetra15pt;
    QuadratureRule fineBdQuadRule = quadRuleTria4pt;

    dFESpace->setQuadRule ( fineQuadRule );
    dFESpace->setBdQuadRule ( fineBdQuadRule );
    dFESpace->qr().showMe();

    if (verbose)
    {
        std::cout << std::endl;
    }

    if (verbose)
    {
        std::cout << "Setting up the reader and the iterations!! " <<  std::endl;
    }

    //Reading fileNames - setting data for reading
    std::string const importerType =  dataFile ( "importer/type", "ensight");
    std::string const fileName     =  dataFile ( "importer/filename", "structure");
    std::string const initialLoaded     =  dataFile ( "importer/initialSol", "NO_DEFAULT_VALUE");
    LifeV::Real initialTime        =  dataFile ( "importer/initialTime", 0.0);

    //Creating the importer
#ifdef HAVE_HDF5
    if ( !importerType.compare ("hdf5") )
    {
        importerSolid.reset ( new  hdf5Filter_Type ( dataFile, fileName) );
    }
    else
#endif
    {
        if ( !importerType.compare ("none") )
        {
            importerSolid.reset ( new emptyExporter_Type ( dataFile, dFESpace->mesh(), "solid", dFESpace->map().comm().MyPID() ) );
        }
        else
        {
            importerSolid.reset ( new  ensightFilter_Type ( dataFile, fileName) );
        }
    }

    importerSolid->setMeshProcId (dFESpace->mesh(), dFESpace->map().comm().MyPID() );

    //Creation of Exporter to check the loaded solution (working only for HDF5)
    // std::string expVerFile = "verificationDisplExporter";
    // LifeV::ExporterHDF5<RegionMesh<LinearTetra> > exporter( dataFile, pointerToMesh, expVerFile, parameters->comm->MyPID());
    // vectorPtr_Type vectVer ( new vector_Type(solid.displacement(),  LifeV::Unique ) );

    // exporter.addVariable( ExporterData<mesh_Type >::VectorField, "displVer", dFESpace, vectVer, UInt(0) );

    // exporter.postProcess(0.0);

    //Reading the displacement field
    vectorPtr_Type solidDisp (new vector_Type (dFESpace->map(), importerSolid->mapType() ) );
    *solidDisp *= 0.0;

    std::string iterationString;

    //Loading the stencil
    iterationString = initialLoaded;

    //Reading
    LifeV::ExporterData<mesh_Type> solidDataReader (LifeV::ExporterData<mesh_Type>::VectorField, std::string ("displacement." + iterationString), dFESpace, solidDisp, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    importerSolid->readVariable (solidDataReader);

    //Exporting the just loaded solution (debug purposes)
    // Real currentLoading(iterInit + 1.0);
    // *vectVer = *solidDisp;
    // exporter.postProcess( currentLoading );

    //It is not useful anymore.
    importerSolid->closeFile();

    //Compute the error at the moment, the L2 error
    Real L2_Error, L2_RelError;

    vector_Type solution (*solidDisp, Repeated);


    //Creation of Exporter to check the loaded solution (working only for HDF5)
    // std::string expVerFile = "exactSolution";
    // LifeV::ExporterHDF5<RegionMesh<LinearTetra> > exporter( dataFile, pointerToMesh, expVerFile, parameters->comm->MyPID());
    // vectorPtr_Type vectVer ( new vector_Type(*solidDisp,  LifeV::Unique ) );

    // exporter.addVariable( ExporterData<mesh_Type >::VectorField, "exactDispl", dFESpace, vectVer, UInt(0) );

    // exporter.postProcess(0.0);

    // //interpolating the exact function on the FESpace
    // dFESpace->interpolate( static_cast<solidFESpace_Type::function_Type>( Private::uexact ), *vectVer , 0);

    // exporter.postProcess(1.0);

    L2_Error = dFESpace->l2Error (Private::uexact, solution, initialTime , &L2_RelError);

    std::ofstream out_norm;
    //save the norm
    if (verbose)
    {
        out_norm.open ("norm.txt", std::ios::app);
        out_norm << L2_Error       << "   "
                 << L2_RelError << "\n";
        out_norm.close();
    }


    MPI_Barrier (MPI_COMM_WORLD);

    std::cout << "Relative error L2: " << L2_RelError << std::endl;
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

    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
