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

   This test is the case of traction of a cube. It does not use the symmetry BCs
   This test uses the FESpace which is the standard in LifeV and the ETFESpace
   The FESpace is used for the BCs of Neumann type since in the ET branch there
   is not the integration over the boundary faces.

   \author Paolo Tricerri <paolo.tricerri@epfl.ch>
   \date 1861-03-17
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
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/anisotropic/StructuralAnisotropicConstitutiveLaw.hpp>
#include <lifev/structure/solver/anisotropic/AnisotropicMultimechanismMaterialNonLinear.hpp>
#include <lifev/structure/solver/anisotropic/HolzapfelMaterialNonLinear.hpp>


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

//Includes for the Expression Template
#include <lifev/eta/fem/ETFESpace.hpp>

#include <iostream>
#include <fstream>

#include "ud_functions.hpp"


using namespace LifeV;

int returnValue = EXIT_FAILURE;
enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}


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

    // Public typedefs
    typedef RegionMesh<LinearTetra>                                     mesh_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                              vectorPtr_Type;
    typedef boost::shared_ptr< TimeAdvance< vector_Type > >             timeAdvance_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    // typedefs for fibers
    // Boost function for fiber direction
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fiberFunction_Type;
    typedef boost::shared_ptr<fiberFunction_Type> fiberFunctionPtr_Type;

    typedef std::vector<fiberFunctionPtr_Type>                          vectorFiberFunction_Type;
    typedef boost::shared_ptr<vectorFiberFunction_Type>                 vectorFiberFunctionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 listOfFiberDirections_Type;

    typedef LifeV::ExporterEmpty<mesh_Type >                      emptyExporter_Type;
    typedef boost::shared_ptr<emptyExporter_Type>                 emptyExporterPtr_Type;

    typedef LifeV::ExporterEnsight<mesh_Type >                    ensightFilter_Type;
    typedef boost::shared_ptr<ensightFilter_Type>                 ensightFilterPtr_Type;

#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<mesh_Type >                       hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                    hdf5FilterPtr_Type;
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
};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
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

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    //dataStructure->showMe();

    //Loading a partitoned mesh or reading a new one
    const std::string partitioningMesh = dataFile ( "partitioningOffline/loadMesh", "no");

    //Creation of pointers
    boost::shared_ptr<MeshPartitioner<mesh_Type> > meshPart;
    boost::shared_ptr<mesh_Type> pointerToMesh;

    if ( ! (partitioningMesh.compare ("no") ) )
    {
        boost::shared_ptr<mesh_Type > fullMeshPtr (new mesh_Type ( ( parameters->comm ) ) );
        //Creating a new mesh from scratch
        MeshData             meshData;
        meshData.setup (dataFile, "solid/space_discretization");
        readMesh (*fullMeshPtr, meshData);

        meshPart.reset ( new MeshPartitioner<mesh_Type> (fullMeshPtr, parameters->comm ) );

        pointerToMesh = meshPart->meshPartition();
    }
    else
    {
        //Creating a mesh object from a partitioned mesh
        const std::string partsFileName (dataFile ("partitioningOffline/hdf5_file_name", "NO_DEFAULT_VALUE.h5") );

        boost::shared_ptr<Epetra_MpiComm> mpiComm =
            boost::dynamic_pointer_cast<Epetra_MpiComm>(parameters->comm);
        PartitionIO<mesh_Type> partitionIO (partsFileName, mpiComm);

        partitionIO.read (pointerToMesh);
    }

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");

    //Mainly used for BCs assembling (Neumann type)
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 3, parameters->comm) );
    solidFESpacePtr_Type dScalarFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 1, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (pointerToMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    // // setting precise quadrature rule for fine meshes
    // const QuadratureRule fineQuadRule = quadRuleTetra15pt;
    // QuadratureRule fineBdQuadRule = quadRuleTria4pt;

    // dFESpace->setQuadRule ( fineQuadRule );
    // dFESpace->setBdQuadRule ( fineBdQuadRule );


    // Setting the fibers
    vectorFiberFunctionPtr_Type pointerToVectorOfFamilies( new vectorFiberFunction_Type( ) );
    (*pointerToVectorOfFamilies).resize( dataStructure->numberFibersFamilies( ) );

    listOfFiberDirections_Type fiberDirections;
    fiberDirections.resize( dataStructure->numberFibersFamilies( ) );

    listOfFiberDirections_Type selectionExport;
    listOfFiberDirections_Type activationDispl;

    selectionExport.resize( dataStructure->numberFibersFamilies( ) );
    activationDispl.resize( dataStructure->numberFibersFamilies( ) );

    if( verbose )
        std::cout << "Size of the number of families: " << (*pointerToVectorOfFamilies).size() << std::endl;

    fibersDirectionList setOfFiberFunctions;
    setOfFiberFunctions.setupFiberDefinitions( dataStructure->numberFibersFamilies( ) );

    if (verbose)
    {
        std::cout << std::endl;
    }

    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "BDF");

    timeAdvance_Type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );

    UInt OrderDev = 2;

    timeAdvance->setup (dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);

    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
    //timeAdvance->showMe();

    //! #################################################################################
    //! BOUNDARY CONDITIONS
    //! #################################################################################
    vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
    compx[0] = 0;
    compy[0] = 1, compz[0] = 2;
    compxy[0] = 0;
    compxy[1] = 1;
    compxz[0] = 0;
    compxz[1] = 2;
    compyz[0] = 1;
    compyz[1] = 2;

    BCFunctionBase zero (bcZero);
    BCFunctionBase nonZero;
    BCFunctionBase pressure;

    nonZero.setFunction (bcNonZero);
    pressure.setFunction (traction);

    //! =================================================================================
    //! BC for StructuredCube4_test_structuralsolver.mesh
    //! =================================================================================
    BCh->addBC ("EdgesIn",      20,  Natural,   Component, pressure, compy);
    BCh->addBC ("EdgesIn",      40,  Essential, Component, zero,    compy);

    //! Symmetry BC
    BCh->addBC ("EdgesIn",      50,   EssentialVertices, Component, zero, compxy);
    BCh->addBC ("EdgesIn",      30,   EssentialVertices, Component, zero, compyz);
    BCh->addBC ("EdgesIn",      80,   EssentialVertices, Component, zero, compxz);
    BCh->addBC ("EdgesIn",      100,  EssentialVertices,  Full, zero, 3);

    BCh->addBC ("EdgesIn",      7, Essential, Component , zero, compx);
    BCh->addBC ("EdgesIn",      3, Essential, Component , zero, compz);
    //! =================================================================================

    // Case of a tube
    //Condition for Inflation
    // BCh->addBC ("EdgesIn",      200, Natural,   Full, pressure, 3);
    // BCh->addBC ("EdgesIn",      40,  Natural,   Full, zero, 3);
    // BCh->addBC ("EdgesIn",      70,  Essential, Full, zero, 3);
    // BCh->addBC ("EdgesIn",      60,  Essential, Full, zero, 3);


    //! 1. Constructor of the structuralSolver
    StructuralOperator< RegionMesh<LinearTetra> > solid;

    //! 2. Setup of the structuralSolver
    solid.setup (dataStructure,
                 dFESpace,
                 dETFESpace,
                 BCh,
                 parameters->comm);

    //! 3. Setting data from getPot
    solid.setDataFromGetPot (dataFile);

    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
    {
        // Setting up the name of the function to define the family
        std::string family="Family";
        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = family + familyNumber;
        (*pointerToVectorOfFamilies)[ k-1 ].reset( new fiberFunction_Type() );
        (*pointerToVectorOfFamilies)[ k-1 ] = setOfFiberFunctions.fiberDefinition( creationString );

        fiberDirections[ k-1 ].reset( new vector_Type(solid.displacement(), Unique) );

	selectionExport[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
	activationDispl[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
    }

    if( !dataStructure->constitutiveLaw().compare("anisotropic") )
    {
        //! 3.b Setting the fibers in the abstract class of Anisotropic materials
        solid.material()->anisotropicLaw()->setupFiberDirections( pointerToVectorOfFamilies );
    }

    //! 4. Building system using TimeAdvance class
    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
    solid.buildSystem (timeAdvanceCoefficient);


    //dataStructure->showMe();
    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type vel (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type acc (new vector_Type (solid.displacement(), Unique) );

    if (verbose)
    {
        std::cout << "S- initialization ... ";
    }

    // Initialization
    //Initialization of TimeAdvance
    std::string const restart =  dataFile ( "importer/restart", "none");
    std::vector<vectorPtr_Type> solutionStencil;
    solutionStencil.resize ( timeAdvance->size() );

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > importerSolid;

    if ( restart.compare ( "none" ) )
    {
        //Reading fileNames - setting data for reading
        std::string const importerType =  dataFile ( "importer/type", "ensight");
        std::string const fileName     =  dataFile ( "importer/filename", "structure");
        std::string const initialLoaded     =  dataFile ( "importer/initialSol", "NO_DEFAULT_VALUE");
        //LifeV::Real initialTime        =  dataFile ( "importer/initialTime", 0.0);

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

        //Reading the displacement field and the timesteps need by the TimeAdvance class
        vectorPtr_Type solidDisp (new vector_Type (dFESpace->map(), importerSolid->mapType() ) );

        std::string iterationString;

        // if( verbose )
        //     std::cout << "size TimeAdvance:" << timeAdvance->size() << std::endl;

        //Loading the stencil
        iterationString = initialLoaded;
        for (UInt iterInit = 0; iterInit < timeAdvance->size(); iterInit++ )
        {
            // if( verbose )
            //     std::cout << "new iterationString" << iterationString << std::endl;

            solidDisp.reset ( new vector_Type (solid.displacement(),  Unique ) );
            *solidDisp *= 0.0;

            LifeV::ExporterData<mesh_Type> solidDataReader (LifeV::ExporterData<mesh_Type>::VectorField, std::string ("displacement." + iterationString), dFESpace, solidDisp, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

            importerSolid->readVariable (solidDataReader);

            // if( verbose )
            //     std::cout << "Norm of the " << iterInit + 1 << "-th solution : " << solidDisp->norm2() << std::endl;

            //Exporting the just loaded solution (debug purposes)
            // Real currentLoading(iterInit + 1.0);
            // *vectVer = *solidDisp;
            // exporter.postProcess( currentLoading );

            solutionStencil[ iterInit ] = solidDisp;

            //initializing the displacement field in the StructuralSolver class with the first solution
            if ( !iterInit )
            {
                solid.initialize ( solidDisp );
            }

            //Updating string name
            int iterations = std::atoi (iterationString.c_str() );
            iterations--;

            std::ostringstream iter;
            iter.fill ( '0' );
            iter << std::setw (5) << ( iterations );
            iterationString = iter.str();

        }

        importerSolid->closeFile();

        //Putting the vector in the TimeAdvance Stencil
        timeAdvance->setInitialCondition (solutionStencil);

    }
    else //Initialize with zero vectors
    {

        // if( verbose )
        //     std::cout << "Starting from scratch" << std::endl;

        vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );

        for ( UInt previousPass = 0; previousPass < timeAdvance->size() ; previousPass++)
        {
            solutionStencil[ previousPass ] = disp;
        }

        timeAdvance->setInitialCondition (solutionStencil);

        //It is automatically done actually. It's clearer if specified.
        solid.initialize ( disp );
    }


    timeAdvance->setTimeStep ( dataStructure->dataTime()->timeStep() );

    timeAdvance->updateRHSContribution ( dataStructure->dataTime()->timeStep() );

    MPI_Barrier (MPI_COMM_WORLD);

    // if (verbose )
    // {
    //     std::cout << "ok." << std::endl;
    // }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
    std::string const exportFileName = dataFile ( "exporter/nameFile", "structure");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, exportFileName ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( dataFile, pointerToMesh, exportFileName, parameters->comm->MyPID() ) );
        }

        else
        {
            exporter.reset ( new ExporterEnsight<RegionMesh<LinearTetra> > ( dataFile, pointerToMesh, exportFileName, parameters->comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( pointerToMesh, parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidVel  ( new vector_Type (solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidAcc  ( new vector_Type (solid.displacement(),  exporter->mapType() ) );
    //    vectorPtr_Type rhsVector ( new vector_Type (solid.displacement(),  exporter->mapType() ) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt (0) );
    //    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "forcingTerm", dFESpace, rhsVector,  UInt (0) );

    // Adding the fibers vectors
    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
    {
        // Setting up the name of the function to define the family
        std::string family="Family-";
        std::string familySel="FamilySel-";
        std::string familyAct="FamilyAct-";
        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = family + familyNumber;
        std::string creationStringSel = familySel + familyNumber;
        std::string creationStringAct = familyAct + familyNumber;
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationString, dFESpace, fiberDirections[ k-1 ], UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringSel, dFESpace, selectionExport[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringAct, dFESpace, activationDispl[ k-1 ], UInt (0) );

        if( !dataStructure->constitutiveLaw().compare("anisotropic") )
        {
            // Extracting the fibers vectors
            *(fiberDirections[ k-1 ]) = solid.material()->anisotropicLaw()->ithFiberVector( k );

            *(selectionExport[ k-1 ]) = *(solid.material()->anisotropicLaw()->selectionCriterion( k ));
            *(activationDispl[ k-1 ]) = *(solid.material()->anisotropicLaw()->activationDisplacement( k ));
        }
    }


    exporter->postProcess ( 0 );
    cout.precision(16);

    //! 5. Initial data
    Real initialTime = dataStructure->dataTime()->initialTime();
    Real dt = dataStructure->dataTime()->timeStep();
    Real T  = dataStructure->dataTime()->endTime();

    //! 6. Setting the pillow saving for restart
    UInt tol = dataStructure->dataTimeAdvance()->orderBDF() + 1;
    UInt saveEvery = dataFile( "exporter/saveEvery", 1);
    UInt r(0);
    UInt d(0);
    UInt iter(0);

    //! =============================================================================
    //! Temporal loop
    //! =============================================================================
    //    for (Real time = dt; time <= T; time += dt)
    for (Real time = initialTime + dt; time <= T; time += dt)
    {
        dataStructure->dataTime()->setTime( time );

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
        }

        //! 6. Updating right-hand side
        *rhs *= 0;
        timeAdvance->updateRHSContribution ( dt );
        *rhs += *solid.massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;
        solid.setRightHandSide ( *rhs );

        //! 7. Iterate --> Calling Newton
        solid.iterate ( BCh );

        timeAdvance->shiftRight ( solid.displacement() );

        *solidDisp = solid.displacement();
        *solidVel  = timeAdvance->firstDerivative();
        *solidAcc  = timeAdvance->secondDerivative();

        // This vector is to export the forcing term
        //*rhsVector = solid.rhsCopy();

        if( !dataStructure->constitutiveLaw().compare("anisotropic") )
        {
            for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
            {
                *(selectionExport[ k-1 ]) = *(solid.material()->anisotropicLaw()->selectionCriterion( k ));
                *(activationDispl[ k-1 ]) = *(solid.material()->anisotropicLaw()->activationDisplacement( k ));
            }
        }


        iter = iter + 1;

        r = iter % saveEvery;
        d = iter - r;

        if ( (iter - d) <= tol || ( (std::floor(d/saveEvery) + 1)*saveEvery - iter ) <= tol )
        {
            exporter->postProcess( time );
        }
	}

        //!--------------------------------------------------------------------------------------------------

        MPI_Barrier (MPI_COMM_WORLD);
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
