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

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
// #include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
// #include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
// #include <lifev/structure/solver/ExponentialMaterialNonLinear.hpp>

#include <lifev/structure/solver/WallTensionEstimator.hpp>
#include <lifev/structure/solver/WallTensionEstimatorData.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <iostream>


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

    typedef LifeV::RegionMesh<LinearTetra>                              mesh_Type;

    // Filters
    typedef typename LifeV::Exporter<mesh_Type  >                                filter_Type;
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
    void CheckResultDisplacement (const Real tensNorm);
    void CheckResultEigenvalues (const Real tensNorm);
    void CheckResultTensions (const Real tensNorm);
    void resultChanged ( void );
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
        rho (1),
        poisson (1),
        young (1),
        bulk (1),
        alpha (1),
        gamma (1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

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

    parameters->rho     = dataFile ( "solid/physics/density", 1. );
    parameters->young   = dataFile ( "solid/physics/young",   1. );
    parameters->poisson = dataFile ( "solid/physics/poisson", 1. );
    parameters->bulk    = dataFile ( "solid/physics/bulk",    1. );
    parameters->alpha   = dataFile ( "solid/physics/alpha",   1. );
    parameters->gamma   = dataFile ( "solid/physics/gamma",   1. );

    std::cout << "density = " << parameters->rho     << std::endl
              << "young   = " << parameters->young   << std::endl
              << "poisson = " << parameters->poisson << std::endl
              << "bulk    = " << parameters->bulk    << std::endl
              << "alpha   = " << parameters->alpha   << std::endl
              << "gamma   = " << parameters->gamma   << std::endl;

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
    // General typedefs
    typedef WallTensionEstimator< mesh_Type >::solutionVect_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                        vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >                       solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                  solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;


    bool verbose = (parameters->comm->MyPID() == 0);

    //! dataElasticStructure for parameters
    GetPot dataFile ( parameters->data_file_name.c_str() );

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    //! Parameters for the analysis
    boost::shared_ptr<WallTensionEstimatorData> tensionData (new WallTensionEstimatorData( ) );
    tensionData->setup (dataFile);

    tensionData->showMe();
    //! Read and partition mesh
    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr (new RegionMesh<LinearTetra> ( ( parameters->comm ) ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, parameters->comm );

    //! Functional spaces - needed for the computations of the gradients
    std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace( new solidFESpace_Type(meshPart,dOrder,3,parameters->comm) );
    solidETFESpacePtr_Type dETFESpace( new solidETFESpace_Type(meshPart,&(dFESpace->refFE()),&(dFESpace->fe().geoMap()), parameters->comm) );


    if (verbose) std::cout << std::endl;


    //! Setting the marker for the volumes
    UInt marker = dataFile ( "solid/physics/material_flag", 1);

    //! 1. Constructor of the class to compute the tensions
    boost::shared_ptr<WallTensionEstimator< mesh_Type > >  solid ( new WallTensionEstimator< mesh_Type >() );

    //! 2. Its setup
    solid->setup(dataStructure,
                 tensionData,
                 dFESpace,
                 dETFESpace,
                 parameters->comm,
                 marker);

    //! 3. Creation of the importers to read the displacement field
    std::string const filename    = tensionData->nameFile();
    std::string const importerType =  tensionData->typeFile();

    std::string iterationString; //useful to iterate over the strings

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
            M_importer.reset ( new emptyFilter_Type ( dataFile, solid->dFESpace().mesh(), "solid", solid->dFESpace().map().comm().MyPID() ) );
        }
        else
        {
            M_importer.reset ( new ensightFilter_Type ( dataFile, filename ) );
        }
    }
    M_importer->setMeshProcId (solid->dFESpace().mesh(), solid->dFESpace().map().comm().MyPID() );

    // The vector where the solution will be stored
    vectorPtr_Type solidDisp (new vector_Type (solid->dFESpace().map(), M_importer->mapType() ) );


    //! 6. Post-processing setting
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "hdf5");
    std::string const nameExporter =  dataFile ( "exporter/name", "tensions");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        M_exporter.reset ( new hdf5Filter_Type ( dataFile, nameExporter ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            M_exporter.reset ( new emptyFilter_Type ( dataFile, meshPart.meshPartition(), nameExporter, parameters->comm->MyPID() ) ) ;
        }

        else
        {
            M_exporter.reset ( new ensightFilter_Type ( dataFile, meshPart.meshPartition(), nameExporter, parameters->comm->MyPID() ) ) ;
        }
    }

    M_exporter->setMeshProcId (solid->dFESpace().mesh(), solid->dFESpace().map().comm().MyPID() );

    vectorPtr_Type solidTensions ( new vector_Type (solid->principalStresses(),  M_exporter->mapType() ) );

    M_exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "vonMises", dFESpace, solidTensions, UInt (0) );
    M_exporter->postProcess ( 0.0 );


    //Post processing for the displacement gradient
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterX;
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterY;
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterZ;

    //Setting pointers
    exporterX.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "gradX" ) );
    exporterY.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "gradY" ) );
    exporterZ.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "gradZ" ) );

    exporterX->setMeshProcId (solid->dFESpace().mesh(), solid->dFESpace().map().comm().MyPID() );
    exporterY->setMeshProcId (solid->dFESpace().mesh(), solid->dFESpace().map().comm().MyPID() );
    exporterZ->setMeshProcId (solid->dFESpace().mesh(), solid->dFESpace().map().comm().MyPID() );

    //Defining the vectors
    vectorPtr_Type gradX ( new vector_Type (solid->gradientX(),  M_exporter->mapType() ) );
    vectorPtr_Type gradY ( new vector_Type (solid->gradientY(),  M_exporter->mapType() ) );
    vectorPtr_Type gradZ ( new vector_Type (solid->gradientZ(),  M_exporter->mapType() ) );

    //Adding variable
    exporterX->addVariable ( ExporterData<mesh_Type >::VectorField, "gradX", solid->dFESpacePtr(),
                             gradX, UInt (0) );
    exporterY->addVariable ( ExporterData<mesh_Type >::VectorField, "gradY", solid->dFESpacePtr(),
                             gradY, UInt (0) );
    exporterZ->addVariable ( ExporterData<mesh_Type >::VectorField, "gradZ", solid->dFESpacePtr(),
                             gradZ, UInt (0) );

    exporterX->postProcess ( 0.0 );
    exporterY->postProcess ( 0.0 );
    exporterZ->postProcess ( 0.0 );


    //! =================================================================================
    //! Analysis - Istant or Interval
    //! =================================================================================

    MPI_Barrier (MPI_COMM_WORLD);

    //! 5. For each interval, the analysis is performed
    LifeV::Real dt =  dataFile ( "solid/time_discretization/timestep", 0.0);
    std::string const nameField =  dataFile ( "solid/analysis/nameField", "NO_DEFAULT_VALUE");

    if ( !tensionData->analysisType().compare ("istant") )
    {
        //Get the iteration number
        iterationString = tensionData->iterStart (0);
        LifeV::Real startTime = tensionData->initialTime (0);

        /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
        LifeV::ExporterData<mesh_Type> solutionDispl  (LifeV::ExporterData<mesh_Type>::VectorField, nameField + "." + iterationString, solid->dFESpacePtr(), solidDisp, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

        //Read the variable
        M_importer->readVariable (solutionDispl);
        M_importer->closeFile();


        //Create and exporter to check importing
        std::string expVerFile = "verificationDisplExporter";
        LifeV::ExporterHDF5<RegionMesh<LinearTetra> > exporter ( dataFile, meshPart.meshPartition(), expVerFile, parameters->comm->MyPID() );
        vectorPtr_Type vectVer ( new vector_Type (solid->displacement(),  exporter.mapType() ) );

        exporter.addVariable ( ExporterData<mesh_Type >::VectorField, "displVer", solid->dFESpacePtr(),
                               vectVer, UInt (0) );

        exporter.postProcess (0.0);
        *vectVer = *solidDisp;
        exporter.postProcess (startTime);


        //Set the current solution as the displacement vector to use
        solid->setDisplacement (*solidDisp);

        std::cout << "The norm of the set displacement, at time " << startTime << ", is: " << solid->displacement().norm2() << std::endl;

        //Perform the analysis
        solid->analyzeTensions();

        //Extracting the gradient
        *gradX = solid->gradientX();
        *gradY = solid->gradientY();
        *gradZ = solid->gradientZ();

        exporterX->postProcess ( startTime );
        exporterY->postProcess ( startTime );
        exporterZ->postProcess ( startTime );

        //Extracting the tensions
        std::cout << std::endl;
        std::cout << "Norm of the tension vector: " << solid->principalStresses().norm2() << std::endl;

        *solidTensions = solid->principalStresses();
        M_exporter->postProcess ( startTime );

        if (verbose )
        {
            std::cout << "Analysis Completed!" << std::endl;
        }

        ///////// CHECKING THE RESULTS OF THE TEST AT EVERY TIMESTEP
        if ( !tensionData->recoveryVariable().compare ("displacement")  )
        {
            CheckResultDisplacement ( solid->principalStresses().norm2()  );
        }
        else if ( !tensionData->recoveryVariable().compare ("eigenvalues") )
        {
            CheckResultEigenvalues ( solid->principalStresses().norm2() );
        }
        else
        {
            CheckResultTensions ( solid->principalStresses().norm2()  );
        }

        ///////// END OF CHECK

        //Closing files
        M_exporter->closeFile();
        exporterX->closeFile();
        exporterY->closeFile();
        exporterZ->closeFile();
        exporter.closeFile();


    }
    else
    {
        std::cout << "we are still working idiot! " << std::endl;
    }

    if (verbose )
    {
        std::cout << "finished" << std::endl;
    }

    MPI_Barrier (MPI_COMM_WORLD);
    //!---------------------------------------------.-----------------------------------------------------
}

void Structure::CheckResultDisplacement (const Real tensNorm)
{
    if ( ( (std::fabs (tensNorm - 4.67086e6) / 4.67086e6) <= 1e-5 ) )
    {
        this->resultChanged( );
    }
}

void Structure::CheckResultEigenvalues (const Real tensNorm)
{
    if ( ( (std::fabs (tensNorm - 4.67086e6) / 4.67086e6) <= 1e-5 ) )
    {
        this->resultChanged( );
    }
}

void Structure::CheckResultTensions (const Real tensNorm)
{
    if ( ( ( std::fabs (tensNorm - 4.67086e6) / 4.67086e6) <= 1e-5 ) )
    {
        this->resultChanged( );
    }
}

void Structure::resultChanged ( void )
{
    std::cout << "Correct Result after the Analysis" << std::endl;
    returnValue = EXIT_SUCCESS;
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
