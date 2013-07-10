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
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

//Include fils which were in the structure.cpp file
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinearPenalized.hpp>

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

    typedef RegionMesh<LifeV::LinearTetra >                       mesh_Type;
    typedef StructuralOperator<mesh_Type >::vector_Type           vector_Type;
    typedef boost::shared_ptr<vector_Type>                        vectorPtr_Type;
    typedef boost::shared_ptr< TimeAdvance< vector_Type > >       timeAdvance_type;

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    //Exporters Typedefs
    typedef typename LifeV::Exporter<mesh_Type >                  filter_Type;
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


    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;
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
    filterPtr_Type exporterSolid;

    //filterPtr_Type checkExporter;


};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

    static Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
    {
        return  0.;
    }

    static Real bcNonZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
    {
        return  20000;
    }

    static Real bcPressure (const Real& t, const Real&  x, const Real& y, const Real& /*Z*/, const ID& i)
    {
        Real radius = 0.5;
        Real pressure = 20000;

        if ( t >= 0.0 )
        {
            switch (i)
            {
                case 0:
                    return  pressure *  ( y / radius );
                    break;
                case 1:
                    return  pressure *  ( x / radius );
                    break;
                case 2:
                    return 0.0;
                    break;
            }
        }

        return 0;

    }

    static Real smoothPressure(const Real& t, const Real&  x, const Real& y, const Real& /*Z*/, const ID& i)
    {
        Real radius = 0.1;
        Real pressure( 14663 );

        //Real highestPressure(146630);
        // Real totalTime = 3.0;
        // Real halfTime = totalTime / 2.0;

        // Real a = ( highestPressure / 2 ) * ( 1/ ((totalTime/2)*(totalTime/2)) );

        // if ( t <= halfTime )
        //     pressure = a * t*t;

        // if ( t > halfTime )
        //     pressure = - a * (t - totalTime)*(t - totalTime) + highestPressure;

        switch (i)
        {
        case 0:
            return  pressure *  ( x / radius ) ;
            break;
        case 1:
            return  pressure *  ( y / radius ) ;
            break;
        case 2:
            return 0.0;
            break;

        }
        return 0;

    }

    static Real pressureUsingNormal (const Real& t, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
    {

        return -5000;
        // if( t < 15.0 )
        //   return  -(300000/(2*3.1415962*0.5*40))*(1/15)*t;
        // else
        //   return  -300000/(2*3.1415962*0.5*40);
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

    parameters->rho     = dataFile ( "solid/physics/density", 1. );
    parameters->young   = dataFile ( "solid/physics/young",   1. );
    parameters->poisson = dataFile ( "solid/physics/poisson", 1. );
    parameters->bulk    = dataFile ( "solid/physics/bulk",    1. );
    parameters->alpha   = dataFile ( "solid/physics/alpha",   1. );
    parameters->gamma   = dataFile ( "solid/physics/gamma",   1. );

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

    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (pointerToMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );


    //Setting a more complex quadrature rule
    const QuadratureRule fineQuadRule = quadRuleTetra15pt;
    QuadratureRule fineBdQuadRule = quadRuleTria7pt;

    dFESpace->setQuadRule ( fineQuadRule );
    dFESpace->setBdQuadRule ( fineBdQuadRule );
    dFESpace->qr().showMe();


    if (verbose)
    {
        std::cout << std::endl;
    }

    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "Newmark");

    timeAdvance_type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );

    UInt OrderDev = 2;

    //! initialization of parameters of time Advance method:
    if (timeAdvanceMethod == "Newmark")
    {
        timeAdvance->setup ( dataStructure->dataTimeAdvance()->coefficientsNewmark() , OrderDev);
    }

    if (timeAdvanceMethod == "BDF")
    {
        timeAdvance->setup (dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);
    }

    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
    timeAdvance->showMe();
    dataStructure->showMe();

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

    BCFunctionBase zero (Private::bcZero);
    BCFunctionBase nonZero (Private::bcNonZero);
    BCFunctionBase pressure (Private::smoothPressure);
    // BCFunctionBase pressureNormal(Private::pressureUsingNormal);


    //! =================================================================================
    //! BC for StructuredCube4_test_structuralsolver.mesh
    //! =================================================================================
    //Condition for Inflation
    BCh->addBC ("EdgesIn",      200, Natural,   Full, pressure, 3);
    BCh->addBC ("EdgesIn",      210,  Natural,   Full, zero, 3);
    BCh->addBC ("EdgesIn",      20,  EssentialVertices, Full, zero, 3);
    BCh->addBC ("EdgesIn",      30,  EssentialVertices, Full, zero, 3);
    BCh->addBC ("EdgesIn",      2,  Essential, Full, zero, 3);
    BCh->addBC ("EdgesIn",      3,  Essential, Full, zero, 3);


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

    //! 4. Building system using TimeAdvance class
    Real timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
    solid.buildSystem (timeAdvanceCoefficient);

    if (verbose)
    {
        std::cout << "S- initialization ... ";
    }

    //Initialization of TimeAdvance
    std::string const restart =  dataFile ( "importer/restart", "none");
    std::vector<vectorPtr_Type> solutionStencil;
    solutionStencil.resize ( timeAdvance->size() );

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

        std::cout << "size TimeAdvance:" << timeAdvance->size() << std::endl;

        //Loading the stencil
        iterationString = initialLoaded;
        for (UInt iterInit = 0; iterInit < timeAdvance->size(); iterInit++ )
        {
            std::cout << "new iterationString" << iterationString << std::endl;

            solidDisp.reset ( new vector_Type (solid.displacement(),  Unique ) );
            *solidDisp *= 0.0;

            LifeV::ExporterData<mesh_Type> solidDataReader (LifeV::ExporterData<mesh_Type>::VectorField, std::string ("displacement." + iterationString), dFESpace, solidDisp, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

            importerSolid->readVariable (solidDataReader);

            std::cout << "Norm of the " << iterInit + 1 << "-th solution : " << solidDisp->norm2() << std::endl;

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

        std::cout << "Starting from scratch" << std::endl;
        vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );

        for ( UInt previousPass = 0; previousPass < timeAdvance->size() ; previousPass++)
        {
            solutionStencil[ previousPass ] = disp;
        }

        timeAdvance->setInitialCondition (solutionStencil);

        //It is automatically done actually. It's clearer if specified.
        solid.initialize ( disp );
    }



    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );

    timeAdvance->updateRHSContribution (dataStructure->dataTime()->timeStep() );

    // timeAdvance->spyStateVector();
    // timeAdvance->spyRHS();

    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose )
    {
        std::cout << "ok." << std::endl;
    }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterSolid;

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterCheck;


    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
    std::string const exportFileName = dataFile ( "exporter/nameFile", "structure");
    std::string const exportCheckName = "checkExporter";

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporterSolid.reset ( new hdf5Filter_Type ( dataFile, exportFileName ) );
        exporterCheck.reset ( new hdf5Filter_Type ( dataFile, exportCheckName ) ); //The check is done in hdf5 for sake of brevity.
    }

    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporterSolid.reset ( new emptyExporter_Type ( dataFile, pointerToMesh, exportFileName, parameters->comm->MyPID() ) );
        }
        else
        {
            exporterSolid.reset ( new ensightFilter_Type ( dataFile, pointerToMesh, exportFileName, parameters->comm->MyPID() ) );
        }
    }

    exporterSolid->setPostDir ( "./" );
    exporterSolid->setMeshProcId ( pointerToMesh, parameters->comm->MyPID() );
    exporterCheck->setMeshProcId ( pointerToMesh, parameters->comm->MyPID() );

    vectorPtr_Type rhsVector ( new vector_Type (solid.rhsCopy(),  Unique ) );

    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(),  Unique ) );
    vectorPtr_Type solidVel  ( new vector_Type (solid.displacement(),  Unique ) );
    vectorPtr_Type solidAcc  ( new vector_Type (solid.displacement(),  Unique ) );

    exporterSolid->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporterSolid->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt (0) );
    exporterSolid->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt (0) );

    exporterCheck->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "rhs", dFESpace, rhsVector,  UInt (0) );

    * solidDisp = solid.displacement();
    *solidVel = timeAdvance->firstDerivative();
    *solidAcc = timeAdvance->secondDerivative();

    //!--------------------------------------------------------------------------------------------
    //!The update of the RHS is done by the TimeAdvance class
    //solid.updateSystem();
    //! =================================================================================

    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );

    //! 5. Initial data
    Real initialTime = dataStructure->dataTime()->initialTime();
    Real dt = dataStructure->dataTime()->timeStep();
    Real T  = dataStructure->dataTime()->endTime();

    exporterSolid->postProcess ( dataStructure->dataTime()->initialTime() );
    exporterCheck->postProcess ( 0 );

    //! 6. Setting the pillow saving for restart
    UInt tol = dataStructure->dataTimeAdvance()->orderBDF() + 1;
    UInt saveEvery = dataFile( "exporter/saveEvery", 1);
    UInt r(0);
    UInt d(0);
    UInt iter(0);

    //! =============================================================================
    //! Temporal loop
    //! =============================================================================
    for (Real time = initialTime + dt; time <= T; time += dt)
    {
        std::cout << "time: " << time << std::endl;
        dataStructure->dataTime()->setTime (time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
        }

        //! 6. Updating right-hand side
        *rhs *= 0;
        timeAdvance->updateRHSContribution ( dt );
        *rhs += *solid.massMatrix() * timeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;

        std::cout << "Norm of the rhsNoBC: " << (*rhs).norm2() << std::endl;
        solid.setRightHandSide ( *rhs );

        //! 7. Iterate --> Calling Newton
        solid.iterate ( BCh );

        timeAdvance->shiftRight ( solid.displacement() );

        *solidDisp = solid.displacement();
        *solidVel  = timeAdvance->firstDerivative();
        *solidAcc  = timeAdvance->secondDerivative();

        std::cout << "vel Inf: " << solidVel->normInf() << std::endl;
        std::cout << "vel 2: "   << solidVel->norm2() << std::endl;

        std::cout << "acc Inf: " << solidAcc->normInf() << std::endl;
        std::cout << "acc 2: "   << solidAcc->norm2() << std::endl;

        //timeAdvance->spyStateVector();

        // Real normVect;
        // normVect =  solid.displacement().norm2();
        // std::cout << "The norm 2 of the displacement field is: "<< normVect << std::endl;
        // std::cout << "The norm 2 of the velocity field is: "<< timeAdvance->velocity().norm2() << std::endl;
        // std::cout << "The norm 2 of the acceleration field is: "<< timeAdvance->acceleration().norm2() << std::endl;

        *rhsVector = solid.rhsCopy();

        iter = iter + 1;

        r = iter % saveEvery;
        d = iter - r;

        if ( (iter - d) <= tol || ( (std::floor(d/saveEvery) + 1)*saveEvery - iter ) <= tol )
        {
            exporterSolid->postProcess( time );
            exporterCheck->postProcess( time );
        }

        //!--------------------------------------------------------------------------------------------------

        MPI_Barrier (MPI_COMM_WORLD);
    }

    exporterSolid->closeFile();
    exporterCheck->closeFile();
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
