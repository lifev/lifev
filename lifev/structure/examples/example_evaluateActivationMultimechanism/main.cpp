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
#include <lifev/structure/solver/isotropic/SecondOrderExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/anisotropic/StructuralAnisotropicConstitutiveLaw.hpp>
#include <lifev/structure/solver/anisotropic/AnisotropicMultimechanismMaterialNonLinear.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

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

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                     scalarETFESpacePtr_Type;

    // typedefs for fibers
    // Boost function for fiber direction
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fiberFunction_Type;
    typedef boost::shared_ptr<fiberFunction_Type> fiberFunctionPtr_Type;

    typedef std::vector<fiberFunctionPtr_Type>                          vectorFiberFunction_Type;
    typedef boost::shared_ptr<vectorFiberFunction_Type>                 vectorFiberFunctionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 listOfFiberDirections_Type;

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

    GetPot dataFile ( parameters->data_file_name.c_str() );

    boost::shared_ptr<BCHandler> BCh ( new BCHandler() );

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
    //Mainly used for BCs assembling (Neumann type)
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 3, parameters->comm) );
    solidFESpacePtr_Type dScalarFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 1, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (pointerToMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );
    scalarETFESpacePtr_Type dScalarETFESpace ( new scalarETFESpace_Type (pointerToMesh, & (dScalarFESpace->refFE() ),
                                                                         & (dScalarFESpace->fe().geoMap() ), parameters->comm) );

    vectorPtr_Type patchAreaVector( new vector_Type( dETFESpace->map() ) );
    vectorPtr_Type patchAreaVectorScalar( new vector_Type( dScalarETFESpace->map() ) );
    QuadratureRule fakeQuadratureRule;

    Real refElemArea (0); //area of reference element
    //compute the area of reference element
    for (UInt iq = 0; iq < dFESpace->qr().nbQuadPt(); iq++)
    {
        refElemArea += dFESpace->qr().weight (iq);
    }

    Real wQuad (refElemArea / dFESpace->refFE().nbDof() );

    //Setting the quadrature Points = DOFs of the element and weight = 1
    std::vector<GeoVector> coords = dFESpace->refFE().refCoor();
    std::vector<Real> weights (dFESpace->fe().nbFEDof(), wQuad);
    fakeQuadratureRule.setDimensionShape ( shapeDimension (dFESpace->refFE().shape() ), dFESpace->refFE().shape() );
    fakeQuadratureRule.setPoints (coords, weights);

    //fakeQuadratureRule.showMe();
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

    // Computing areas
    evaluateNode( elements ( dScalarFESpace->mesh() ),
                  fakeQuadratureRule,
                  dScalarETFESpace,
                  meas_K * phi_i
                  ) >> patchAreaVectorScalar;
    patchAreaVectorScalar->globalAssemble();

    ExpressionVectorFromNonConstantScalar<ExpressionMeas, 3  > vMeas( meas_K );
    evaluateNode( elements ( dETFESpace->mesh() ),
                  fakeQuadratureRule,
                  dETFESpace,
                  dot( vMeas , phi_i )
                  ) >> patchAreaVector;
    patchAreaVector->globalAssemble();

    // Setting the fibers
    // vectorFiberFunctionPtr_Type pointerToVectorOfFamilies( new vectorFiberFunction_Type( ) );
    // (*pointerToVectorOfFamilies).resize( dataStructure->numberFibersFamilies( ) );

    vectorPtr_Type jacobian;
    vectorPtr_Type jacobianIsochoric;

    listOfFiberDirections_Type fiberDirections; // store vector of fibers
    listOfFiberDirections_Type deformedFiberDirections; // store vector of fibers
    listOfFiberDirections_Type normalizedDefFibers; // store vector of fibers
    listOfFiberDirections_Type selectionExport; // vector of the nodal satisfaction of the activation condition
    listOfFiberDirections_Type activationDispl; // active displacement for each fiber
    listOfFiberDirections_Type stretchFibers;   // stretch of the fiber at a certain displacement
    //listOfFiberDirections_Type activationFibers;// corresponding activation function

    fiberDirections.resize( dataStructure->numberFibersFamilies( ) );
    deformedFiberDirections.resize( dataStructure->numberFibersFamilies( ) );
    normalizedDefFibers.resize( dataStructure->numberFibersFamilies( ) );
    selectionExport.resize( dataStructure->numberFibersFamilies( ) );
    activationDispl.resize( dataStructure->numberFibersFamilies( ) );
    stretchFibers.resize( dataStructure->numberFibersFamilies( ) );
    //activationFibers.resize( dataStructure->numberFibersFamilies( ) );

    std::cout << "Size of the number of families: " << dataStructure->numberFibersFamilies( ) << std::endl;

    // fibersDirectionList setOfFiberFunctions;
    // setOfFiberFunctions.setupFiberDefinitions( dataStructure->numberFibersFamilies( ) );

    //! 1. Constructor of the structuralSolver
    StructuralOperator< RegionMesh<LinearTetra> > solid;

    //! 2. Setup of the structuralSolver
    solid.setup (dataStructure,
                 dFESpace,
                 dETFESpace,
                 BCh,
                 parameters->comm);

    solid.material()->anisotropicLaw()->setNumberOfFamilies( dataStructure->numberFibersFamilies( ) );

    // Setting the vector of fibers functions
    for( UInt k(1); k <= dataStructure->numberFibersFamilies( ); k++ )
    {
        // // Setting up the name of the function to define the family
        // std::string family="Family";
        // // adding the number of the family
        // std::string familyNumber;
        // std::ostringstream number;
        // number << ( k );
        // familyNumber = number.str();

        // // Name of the function to create
        // std::string creationString = family + familyNumber;
        // (*pointerToVectorOfFamilies)[ k-1 ].reset( new fiberFunction_Type() );
        // (*pointerToVectorOfFamilies)[ k-1 ] = setOfFiberFunctions.fiberDefinition( creationString );

        fiberDirections[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
        deformedFiberDirections[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
        normalizedDefFibers[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
        selectionExport[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
        activationDispl[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
        stretchFibers[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
        //activationFibers[ k-1 ].reset( new vector_Type( dFESpace->map() ) );
    }
    jacobian.reset( new vector_Type( dScalarFESpace->map() ) );
    jacobianIsochoric.reset( new vector_Type( dScalarFESpace->map() ) );

    //! 3.b Setting the fibers in the abstract class of Anisotropic materials
    //! This example is made only for anisotropic laws that means that if an isotropic
    //! law will be used, a seg fault will be obtained at this point.
    //! The method setupFiberDirections sets the fiber interpolating the functions.
    //! In this example the fibers are read from the paraview file.

    //! solid.material()->anisotropicLaw()->setupFiberDirections( pointerToVectorOfFamilies );

    // Set up the reading apparatus: importer, start and end iterations etc.
    //! 3. Creation of the importers to read the displacement field
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > importerSolid;

    std::string const filename    = dataFile ( "importer/filename", "structure");
    std::string const nameField = dataFile ( "importer/nameField", "displacement");
    std::string const importerType = dataFile ( "importer/type", "hdf5");

    ASSERT( !importerType.compare("hdf5"), "The example only works for HDF5 files!!" );

    // How many solution do we have to read?
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
        numberOfSol = end - start + 1;
    }
    ASSERT( numberOfSol, "You did not set any solution to read!! ");

    if (verbose)
    {
        std::cout << std::endl;
        std::cout << "The filename is    : " << filename << std::endl;
        std::cout << "The importerType is: " << importerType << std::endl;
    }

#ifdef HAVE_HDF5
    if (importerType.compare ("hdf5") == 0)
    {
        importerSolid.reset ( new hdf5Filter_Type (dataFile, filename) );
    }
#endif

    importerSolid->setMeshProcId (dFESpace->mesh(), dFESpace->map().comm().MyPID() );

    MPI_Barrier (MPI_COMM_WORLD);

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
    std::string const exportFileName = dataFile ( "exporter/nameFile", "activationFunction");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, exportFileName ) );
    }
#endif

    exporter->setMeshProcId ( pointerToMesh, parameters->comm->MyPID() );

    MPI_Barrier (MPI_COMM_WORLD);

    // Setting exporter quantities
    vectorPtr_Type solidDisp ( new vector_Type ( dFESpace->map() ) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "jacobian", dScalarFESpace, jacobian, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "jacobianIsochoric", dScalarFESpace, jacobianIsochoric, UInt (0) );

    // Adding the fibers vectors
    // Setting the vector of fibers functions
    for( UInt k(1); k <= dataStructure->numberFibersFamilies( ); k++ )
    {
        // Setting up the name of the function to define the family
        std::string family="Family-";
        std::string familyDef="FamilyDef-";
        std::string familyNorm="FamilyNorm-";
        std::string familySel="FamilySel-";
        std::string familyAct="FamilyAct-";
        std::string familyStretch="FamilyStretch-";
        //std::string familyAtan="FamilyAtan-";
        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = family + familyNumber;
        std::string creationStringDef = familyDef + familyNumber;
        std::string creationStringNorm = familyNorm + familyNumber;
        std::string creationStringSel = familySel + familyNumber;
        std::string creationStringAct = familyAct + familyNumber;
        std::string creationStringStretch = familyStretch + familyNumber;
        //std::string creationStringAtan = familyAtan + familyNumber;

        // Setting the vectors with the maps
        (fiberDirections[ k-1 ]).reset( new vector_Type( dFESpace->map() ) );
        (deformedFiberDirections[ k-1 ]).reset( new vector_Type( dFESpace->map() ) );
        (normalizedDefFibers[ k-1 ]).reset( new vector_Type( dFESpace->map() ) );
        (selectionExport[ k-1 ]).reset( new vector_Type( dFESpace->map() ) );
        (activationDispl[ k-1 ]).reset( new vector_Type( dFESpace->map() ) );
        (stretchFibers[ k-1 ]).reset( new vector_Type( dScalarFESpace->map() ) );
        //(activationFibers[ k-1 ]).reset( new vector_Type( dScalarFESpace->map() ) );

        // Adding them to the exporter
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationString, dFESpace, fiberDirections[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringDef, dFESpace, deformedFiberDirections[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringNorm, dFESpace, normalizedDefFibers[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringSel, dFESpace, selectionExport[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringAct, dFESpace, activationDispl[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, creationStringStretch, dScalarFESpace, stretchFibers[ k-1 ], UInt (0) );
        //exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, creationStringAtan, dScalarFESpace, activationFibers[ k-1 ], UInt (0) );
    }

    exporter->postProcess ( 0.0 );
    cout.precision(16);

    //! 5. Initial data
    Real initialTime = dataStructure->dataTime()->initialTime();
    Real dt = dataStructure->dataTime()->timeStep();
    Real T  = dataStructure->dataTime()->endTime();
    Real time;

    vectorPtr_Type readDispl( new vector_Type( dFESpace->map() ) );

    // Vector for reading fibers. Sorry it is hard coded!
    // TODO: make this more general for n fibers.
    vectorPtr_Type stFamily( new vector_Type( dFESpace->map() ) );
    vectorPtr_Type ndFamily( new vector_Type( dFESpace->map() ) );


    // Loop on the solutions
    for ( UInt i(start); i < numberOfSol; i++ )
    {
        time = i * dt + initialTime;
        *readDispl *= 0.0;

        UInt current(0);
        if( !readType.compare("interval") )
            current = i ;
        else
            current = dataFile ( "importer/iteration" , 100000, i );

        // Transform current ingot a string
        std::string iterationString;
        std::ostringstream number;
        number.fill ( '0' );
        number << std::setw (5) << ( current );
        iterationString = number.str();

        std::cout << "Current reading: " << iterationString << std::endl;

        /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
        LifeV::ExporterData<mesh_Type> solutionDispl  (LifeV::ExporterData<mesh_Type>::VectorField,nameField + "." + iterationString,
                                                       dFESpace, readDispl, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

        /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
        LifeV::ExporterData<mesh_Type> familyFirst  (LifeV::ExporterData<mesh_Type>::VectorField,"Family-1." + iterationString,
                                                     dFESpace, stFamily, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

        /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
        LifeV::ExporterData<mesh_Type> familySecond  (LifeV::ExporterData<mesh_Type>::VectorField,"Family-2." + iterationString,
                                                      dFESpace,ndFamily, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

        //Read the variables ( displacement and fibers )
        importerSolid->readVariable (solutionDispl);
        importerSolid->readVariable (familyFirst);
        importerSolid->readVariable (familySecond);

        // Setting fibers in the material class
        solid.material()->anisotropicLaw()->setIthFiberVector( *stFamily, 1 );
        solid.material()->anisotropicLaw()->setIthFiberVector( *ndFamily, 2 );

        solid.material()->anisotropicLaw()->computeReferenceConfigurations( *readDispl, dataStructure, solid.displayerPtr() );

        *solidDisp = *readDispl;

        // Defining general expressions
        ExpressionDefinitions::deformationGradient_Type  F =
            ExpressionDefinitions::deformationGradient ( dETFESpace,  *readDispl, 0, identity );

	ExpressionDefinitions::determinantTensorF_Type J = 
	  ExpressionDefinitions::determinantF( F );

	ExpressionDefinitions::isochoricChangeOfVariable_Type isoJ =
	  ExpressionDefinitions::isochoricDeterminant( J  );

	*jacobian *= 0.0;
	evaluateNode( elements ( dScalarETFESpace->mesh() ),
		      fakeQuadratureRule,
		      dScalarETFESpace,
		      meas_K * J  * phi_i
		      ) >> jacobian;
	jacobian->globalAssemble();
	*(jacobian) = *(jacobian) / *patchAreaVectorScalar;

	*jacobianIsochoric *= 0.0;
	evaluateNode( elements ( dScalarETFESpace->mesh() ),
		      fakeQuadratureRule,
		      dScalarETFESpace,
		      meas_K * isoJ  * phi_i
		      ) >> jacobianIsochoric;
	jacobianIsochoric->globalAssemble();
	*(jacobianIsochoric) = *(jacobianIsochoric) / *patchAreaVectorScalar;
	
        // // Definition of C = F^T F
        ExpressionDefinitions::rightCauchyGreenTensor_Type C =
             ExpressionDefinitions::tensorC( transpose(F), F );

        for( UInt k(1); k <= dataStructure->numberFibersFamilies( ) ; k++ )
        {
            // Fibers definitions
            ExpressionDefinitions::interpolatedValue_Type fiberIth =
                ExpressionDefinitions::interpolateFiber( dETFESpace,
                                                         solid.material()->anisotropicLaw()->ithFiberVector( k ) );

            ExpressionMultimechanism::activatedFiber_Type defFiberIth =
                ExpressionMultimechanism::activateFiberDirection( F, fiberIth);

            ExpressionMultimechanism::normActivatedFiber_Type normDef =
                ExpressionMultimechanism::normActivatedFiber( defFiberIth );

            ExpressionMultimechanism::normalizedFiber_Type normalizedFiber =
                ExpressionMultimechanism::normalizedFiberDirection( defFiberIth, normDef );

            // f /otimes f
            ExpressionDefinitions::outerProduct_Type Mith = ExpressionDefinitions::fiberTensor( fiberIth );

            // Definition of the fourth invariant : I_4^i = C:Mith
            ExpressionDefinitions::stretch_Type IVith = ExpressionDefinitions::fiberStretch( C, Mith );

            *stretchFibers[ k-1 ] *= 0.0;
            evaluateNode( elements ( dScalarETFESpace->mesh() ),
                          fakeQuadratureRule,
                          dScalarETFESpace,
                          meas_K * IVith  * phi_i
                          ) >> stretchFibers[ k-1 ];
            stretchFibers[ k-1 ]->globalAssemble();
            *(stretchFibers[ k-1 ]) = *(stretchFibers[ k-1 ]) / *patchAreaVectorScalar;

            *deformedFiberDirections[ k-1 ] *= 0.0;
            evaluateNode( elements ( dETFESpace->mesh() ),
                          fakeQuadratureRule,
                          dETFESpace,
                          meas_K * dot( defFiberIth ,phi_i )
                          ) >> deformedFiberDirections[ k-1 ];
            deformedFiberDirections[ k-1 ]->globalAssemble();
            *(deformedFiberDirections[ k-1 ]) = *(deformedFiberDirections[ k-1 ]) / *patchAreaVector;

            *normalizedDefFibers[ k-1 ] *= 0.0;
            evaluateNode( elements ( dETFESpace->mesh() ),
                          fakeQuadratureRule,
                          dETFESpace,
                          meas_K * dot ( normalizedFiber , phi_i )
                          ) >> normalizedDefFibers[ k-1 ];
            normalizedDefFibers[ k-1 ]->globalAssemble();
            *(normalizedDefFibers[ k-1 ]) = *(normalizedDefFibers[ k-1 ]) / *patchAreaVector;

            // Savint to Export
            *(fiberDirections[ k-1 ]) = solid.material()->anisotropicLaw()->ithFiberVector( k );
            *(selectionExport[ k-1 ]) = *(solid.material()->anisotropicLaw()->selectionCriterion( k ));
            *(activationDispl[ k-1 ]) = *(solid.material()->anisotropicLaw()->activationDisplacement( k ));
        }
        exporter->postProcess( time );

	}

    //!--------------------------------------------------------------------------------------------------

    MPI_Barrier (MPI_COMM_WORLD);

    // Closing files
    importerSolid->closeFile();
    exporter->closeFile();
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
