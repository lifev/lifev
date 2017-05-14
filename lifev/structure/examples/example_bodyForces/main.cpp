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

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/isotropic/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/anisotropic/StructuralAnisotropicConstitutiveLaw.hpp>
#include <lifev/structure/solver/anisotropic/HolzapfelMaterialNonLinear.hpp>


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

//Includes for the Expression Template
#include <lifev/eta/fem/ETFESpace.hpp>

#include <iostream>

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
    typedef StructuralOperator< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef std::shared_ptr<vector_Type>                              vectorPtr_Type;
    typedef std::shared_ptr< TimeAdvance< vector_Type > >             timeAdvance_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef std::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef std::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    // typedefs for fibers
    // std function for fiber direction
    typedef std::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > analyticalFunction_Type;
    typedef std::shared_ptr<analyticalFunction_Type> analyticalFunctionPtr_Type;

    typedef std::function<VectorSmall<3> ( Real const&, Real const&, Real const&, Real const& ) > bodyFunction_Type;
    typedef std::shared_ptr<bodyFunction_Type> bodyFunctionPtr_Type;


    typedef std::vector<analyticalFunctionPtr_Type>                          vectorFiberFunction_Type;
    typedef std::shared_ptr<vectorFiberFunction_Type>                 vectorFiberFunctionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 listOfFiberDirections_Type;

    /** @name Constructors, destructor
     */
    //@{
    Structure ( int                                   argc,
                char**                                argv,
                std::shared_ptr<Epetra_Comm>        structComm );

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
    std::shared_ptr<Private> parameters;
};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    std::shared_ptr<Epetra_Comm>     comm;
};



Structure::Structure ( int                                   argc,
                       char**                                argv,
                       std::shared_ptr<Epetra_Comm>        structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
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
    std::shared_ptr<BCHandler> BCh ( new BCHandler() );

    //! dataElasticStructure
    GetPot dataFile ( parameters->data_file_name.c_str() );

    std::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    dataStructure->showMe();

    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    std::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra> (  parameters->comm  ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< RegionMesh<LinearTetra> > meshPart ( fullMeshPtr, parameters->comm );
    std::shared_ptr<RegionMesh<LinearTetra> > localMeshPtr (new RegionMesh<LinearTetra> (  parameters->comm  ) );
    localMeshPtr = meshPart.meshPartition();

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");

    //Mainly used for BCs assembling (Neumann type)
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localMeshPtr, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    vectorFiberFunctionPtr_Type pointerToVectorOfFamilies;
    listOfFiberDirections_Type fiberDirections;
    fibersDirectionList setOfFiberFunctions;

    if( !dataStructure->constitutiveLaw().compare("anisotropic") )
    {
        // Setting the fibers
        pointerToVectorOfFamilies.reset( new vectorFiberFunction_Type( ) );
        (*pointerToVectorOfFamilies).resize( dataStructure->numberFibersFamilies( ) );

        fiberDirections.resize( dataStructure->numberFibersFamilies( ) );

        // if ( verbose )
        // {
        //     std::cout << "Size of the number of families: " << (*pointerToVectorOfFamilies).size() << std::endl;
        // }

        setOfFiberFunctions.setupFiberDefinitions( dataStructure->numberFibersFamilies( ) );
    }

    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "BDF");

    timeAdvance_Type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );

    UInt OrderDev = 2;

    timeAdvance->setup (dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);

    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
    //timeAdvance->showMe();

    //! #################################################################################
    //! BOUNDARY CONDITIONS + SOURCE TERMS
    //! #################################################################################
    std::vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
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
    BCFunctionBase displacementImposed;

    nonZero.setFunction (bcNonZero);
    displacementImposed.setFunction( analyticalDisplacement );

    //! =================================================================================
    //! BC for StructuredCube4_test_structuralsolver.mesh
    //! =================================================================================
    BCh->addBC ("EdgesIn",      20,  Natural, Component, nonZero, compy);
    BCh->addBC ("EdgesIn",      40,  Essential, Component, zero, compy);

    //! Symmetry BC
    BCh->addBC ("EdgesIn",      50,   EssentialVertices, Component, zero, compxy);
    BCh->addBC ("EdgesIn",      30,   EssentialVertices, Component, zero, compyz);
    BCh->addBC ("EdgesIn",      80,   EssentialVertices, Component, zero, compxz);
    BCh->addBC ("EdgesIn",      100,  EssentialVertices,  Full, zero, 3);

    BCh->addBC ("EdgesIn",      7, Essential, Component, zero, compx);
    BCh->addBC ("EdgesIn",      3, Essential, Component, zero, compz);

    //! =================================================================================

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

    // Comment this part in order not to have body force in the RHS of the equation
    bool bodyForce = dataFile ( "solid/physics/bodyForce" , false );
    if( bodyForce )
    {
        solid.setHavingSourceTerm(  bodyForce );
        bodyFunctionPtr_Type bodyTerm( new bodyFunction_Type( f ) );

        solid.setSourceTerm( bodyTerm );
    }

    if( !dataStructure->constitutiveLaw().compare("anisotropic") )
    {
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
            (*pointerToVectorOfFamilies)[ k-1 ].reset( new analyticalFunction_Type() );
            (*pointerToVectorOfFamilies)[ k-1 ] = setOfFiberFunctions.fiberDefinition( creationString );

            fiberDirections[ k-1 ].reset( new vector_Type(solid.displacement(), Unique) );
        }

        //! 3.b Setting the fibers in the abstract class of Anisotropic materials
        solid.material()->anisotropicLaw()->setupFiberDirections( pointerToVectorOfFamilies );
    }

    //! 4. Building system using TimeAdvance class
    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
    solid.buildSystem (timeAdvanceCoefficient);

    vectorPtr_Type analyticDispl (new vector_Type ( dFESpace->map() ) );
    // Interpolating the solution on the mesh
    dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type> ( analyticalDisplacement ), *analyticDispl, 0.0 );


    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    //! 5. Initial data
    Real dt = dataStructure->dataTime()->timeStep();

    vectorPtr_Type rhs (new vector_Type ( dFESpace->map() ) );
    vectorPtr_Type disp (new vector_Type ( dFESpace->map() ) );
    vectorPtr_Type vel (new vector_Type ( dFESpace->map() ) );
    vectorPtr_Type acc (new vector_Type ( dFESpace->map() ) );

    std::vector<vectorPtr_Type> uv0;

    vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );

    if (timeAdvanceMethod == "BDF")
    {
        Real tZero = dataStructure->dataTime()->initialTime();

        for ( UInt previousPass = 0; previousPass < timeAdvance->size() ; previousPass++)
        {
            Real previousTimeStep = tZero - previousPass * dt;
            std::cout << "BDF " << previousTimeStep << "\n";
            uv0.push_back (disp);
        }
    }

    timeAdvance->setInitialCondition (uv0);

    timeAdvance->setTimeStep ( dt );

    timeAdvance->updateRHSContribution ( dt );

    solid.initialize ( disp );

    MPI_Barrier (MPI_COMM_WORLD);

    std::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    std::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterCheck;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
    std::string const exporterNameFile =  dataFile ( "exporter/nameFile", "structure");
    std::string const exporterCheckName =  dataFile ( "exporter/nameCheck", "verifyVectors");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, exporterNameFile ) );
        exporterCheck.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, exporterCheckName ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }

        else
        {
            exporter.reset ( new ExporterEnsight<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( meshPart.meshPartition(), parameters->comm->MyPID() );
    exporterCheck->setPostDir ( "./" );
    exporterCheck->setMeshProcId ( meshPart.meshPartition(), parameters->comm->MyPID() );

    vectorPtr_Type imposedSolution ( new vector_Type ( dFESpace->map() ) );
    // vectorPtr_Type bodyForceVector ( new vector_Type ( dFESpace->map() ) );

    vectorPtr_Type solidDisp ( new vector_Type ( dFESpace->map() ) );
    vectorPtr_Type solidVel  ( new vector_Type ( dFESpace->map() ) );
    vectorPtr_Type solidAcc  ( new vector_Type ( dFESpace->map() ) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt (0) );

    //exporterCheck->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "forcingTerm", dFESpace, bodyForceVector,  UInt (0) );
    exporterCheck->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "imposedSolution", dFESpace, imposedSolution,  UInt (0) );

    *imposedSolution = *analyticDispl;

    if( !dataStructure->constitutiveLaw().compare("anisotropic") )
    {
        // Adding the fibers vectors
        // Setting the vector of fibers functions
        for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
        {
            // Setting up the name of the function to define the family
            std::string family="Family-";
            // adding the number of the family
            std::string familyNumber;
            std::ostringstream number;
            number << ( k );
            familyNumber = number.str();

            // Name of the function to create
            std::string creationString = family + familyNumber;
            exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationString, dFESpace, fiberDirections[ k-1 ], UInt (0) );

            // Extracting the fibers vectors
            *(fiberDirections[ k-1 ]) = solid.material()->anisotropicLaw()->ithFiberVector( k );
        }
    }

    exporter->postProcess ( 0 );
    exporterCheck->postProcess ( 0 );
    std::cout.precision(16);

    // Saving the analyitical displacement
    exporterCheck->postProcess ( 1.0 );

    Real normVect;

    // Computing the stiffness vector to set as RHS
    // solid.material()->computeStiffness( *analyticDispl, 0, 1.0, dataStructure,
    //                                     solid.mapMarkersVolumes(), solid.mapMarkersIndexes(),
    //                                     solid.displayerPtr() );

    // vectorPtr_Type stiffnessVector;
    // stiffnessVector.reset( new vector_Type( dFESpace->map() ) );

    // *stiffnessVector = *( solid.material()->stiffVector() );
    //! =============================================================================
    //! Temporal loop
    //! =============================================================================
    //    for (Real time = dt; time <= T; time += dt)
    for (dataStructure->dataTime()->setTime ( dt ) ; dataStructure->dataTime()->canAdvance( ); dataStructure->dataTime()->updateTime( ) )
    {

        //! 6. Updating right-hand side
        *rhs *= 0;
        timeAdvance->updateRHSContribution ( dt );
        *rhs += *solid.massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;

        // Summing the term coming from the imposed solution
        //*rhs += *stiffnessVector;

        if( !solid.havingSourceTerm() )
        {
            solid.setRightHandSide ( *rhs );
        }
        else
        {
            solid.updateRightHandSideWithBodyForce( dataStructure->dataTime()->time(), *rhs );
        }

        //debug
        //*bodyForceVector = solid.bodyForce();

        //! 7. Iterate --> Calling Newton
        solid.iterate ( BCh );

        timeAdvance->shiftRight ( solid.displacement() );

        *solidDisp = solid.displacement();
        *solidVel  = timeAdvance->firstDerivative();
        *solidAcc  = timeAdvance->secondDerivative();

        exporter->postProcess ( dataStructure->dataTime()->time() );
        //!--------------------------------------------------------------------------------------------------

        MPI_Barrier (MPI_COMM_WORLD);
    }
    Real errorAbs, error;

    // Computing the L2 error of the solution
    errorAbs = dFESpace->l2Error( static_cast<solidFESpace_Type::function_Type> ( analyticalDisplacement ),
                                  solid.displacement(), 0.0, &error);

    std::cout << "L2 relative error: " << error << std::endl;
    std::cout << "L2 abs error: " << errorAbs << std::endl;
}

int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "% using MPI" << std::endl;
    }
#else
    std::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
    std::cout << "% using serial Version" << std::endl;
#endif

    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
