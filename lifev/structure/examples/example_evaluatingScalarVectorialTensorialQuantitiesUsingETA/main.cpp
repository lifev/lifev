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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/fem/ExpressionDefinitions.hpp>

#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>
// #include <lifev/structure/solver/anisotropic/AnisotropicMultimechanismMaterialNonLinear.hpp>

// Evaluation operations
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/eta/expression/Evaluate.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <iostream>
#include "ud_functions.hpp"

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

class Structure
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
    filterPtr_Type M_importer;
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
    typedef boost::shared_ptr<VectorEpetra>                              vectorPtr_Type;

    // typedefs for fibers
    // Boost function for fiber direction
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fiberFunction_Type;
    typedef boost::shared_ptr<fiberFunction_Type> fiberFunctionPtr_Type;

    typedef std::vector<fiberFunctionPtr_Type>                          vectorFiberFunction_Type;
    typedef boost::shared_ptr<vectorFiberFunction_Type>                 vectorFiberFunctionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 listOfFiberDirections_Type;

    // General typedefs
    typedef typename StructuralOperator<mesh_Type >::vector_Type vector_Type;
    typedef boost::shared_ptr<vector_Type>                        vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >                       solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                  solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                     scalarETFESpacePtr_Type;

    bool verbose = (parameters->comm->MyPID() == 0);

    //! dataElasticStructure for parameters
    GetPot dataFile ( parameters->data_file_name.c_str() );
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    dataStructure->showMe();

    //! Read and partition mesh
    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

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


    //! Functional spaces - needed for the computations of the gradients
    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (pointerToMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    //! Scalar ETFEspace to evaluate scalar quantities
    solidFESpacePtr_Type dScalarFESpace ( new solidFESpace_Type (pointerToMesh, dOrder, 1, parameters->comm) );
    scalarETFESpacePtr_Type dScalarETFESpace ( new scalarETFESpace_Type (pointerToMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    // Setting the fibers
    vectorFiberFunctionPtr_Type pointerToVectorOfFamilies( new vectorFiberFunction_Type( ) );
    (*pointerToVectorOfFamilies).resize( dataStructure->numberFibersFamilies( ) );

    // Interpolate fibers
    std::vector<vectorPtr_Type> vectorInterpolated(0);
    vectorPtr_Type              referenceDirectionVector;

    // Interpolating fiber fields
    vectorInterpolated.resize( (*pointerToVectorOfFamilies).size() );

    fibersDirectionList setOfFiberFunctions;
    setOfFiberFunctions.setupFiberDefinitions( dataStructure->numberFibersFamilies( ) );

    if( !dataStructure->constitutiveLaw().compare("anisotropic") )
    {

        if( verbose )
            std::cout << "Number of families: " << (*pointerToVectorOfFamilies).size() << std::endl;

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
        }


        for( UInt k(0); k < (*pointerToVectorOfFamilies).size(); k++ )
        {
            vectorInterpolated[ k ].reset( new vector_Type( dFESpace->map() ) );
            dFESpace->interpolate ( *( ( *(pointerToVectorOfFamilies) )[ k ] ) ,
                                    * ( ( vectorInterpolated )[ k ] ),
                                    0.0 );
        }

        referenceDirectionVector.reset( new vector_Type( dFESpace->map() ) );
        dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type> ( referenceDirection ),
                                    * ( referenceDirectionVector ),
                                    0.0 );
    }


    //! 3. Creation of the importers to read the displacement field
    std::string const filename    = dataFile ( "importer/filename", "structure");
    std::string const importerType = dataFile ( "importer/type", "hdf5");

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
            M_importer.reset ( new emptyFilter_Type ( dataFile, dFESpace->mesh(), "solid", dFESpace->map().comm().MyPID() ) );
        }
        else
        {
            M_importer.reset ( new ensightFilter_Type ( dataFile, filename ) );
        }
    }
    M_importer->setMeshProcId (dFESpace->mesh(), dFESpace->map().comm().MyPID() );

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
            exporter.reset ( new emptyFilter_Type ( dataFile, pointerToMesh, nameExporter, parameters->comm->MyPID() ) ) ;
        }

        else
        {
            exporter.reset ( new ensightFilter_Type ( dataFile, pointerToMesh, nameExporter, parameters->comm->MyPID() ) ) ;
        }
    }
    exporter->setMeshProcId (dFESpace->mesh(), dFESpace->map().comm().MyPID() );

    // Scalar vector to have scalar quantities
    vectorPtr_Type patchAreaVector;
    vectorPtr_Type patchAreaVectorScalar;
    vectorPtr_Type jacobianF;

    vectorPtr_Type disp;
    vectorPtr_Type dispRead;

    vectorPtr_Type family1;
    vectorPtr_Type family2;
    vectorPtr_Type family1Read;
    vectorPtr_Type family2Read;
    vectorPtr_Type family1Computed;
    vectorPtr_Type family2Computed;

    std::vector< vectorPtr_Type > activationFunction(0);
    std::vector< vectorPtr_Type > stretch(0);
    std::vector< vectorPtr_Type > deformedVector(0);
    std::vector< vectorPtr_Type > characteristicAngle(0);
    std::vector< vectorPtr_Type > referenceAngle(0);

    // vectorPtr_Type deformationGradient_1;
    // vectorPtr_Type deformationGradient_2;
    // vectorPtr_Type deformationGradient_3;

    // // To check
    // vectorPtr_Type reconstructed_1;
    // vectorPtr_Type reconstructed_2;
    // vectorPtr_Type reconstructed_3;

    if( !dataStructure->constitutiveLaw().compare("anisotropic") )
    {
        activationFunction.resize( dataStructure->numberFibersFamilies( ) );
        stretch.resize( dataStructure->numberFibersFamilies( ) );
        deformedVector.resize( dataStructure->numberFibersFamilies( ) );
        characteristicAngle.resize( dataStructure->numberFibersFamilies( ) );
        referenceAngle.resize( dataStructure->numberFibersFamilies( ) );

        family1.reset( new vector_Type( dFESpace->map() ) );
        family2.reset( new vector_Type( dFESpace->map() ) );
        family1Computed.reset( new vector_Type( dFESpace->map() ) );
        family2Computed.reset( new vector_Type( dFESpace->map() ) );

        family1Read.reset( new vector_Type( dFESpace->map() ) );
        family2Read.reset( new vector_Type( dFESpace->map() ) );
    }

    patchAreaVector.reset ( new vector_Type ( dETFESpace->map() ) );
    patchAreaVectorScalar.reset ( new vector_Type ( dScalarETFESpace->map() ) );
    jacobianF.reset ( new vector_Type ( dScalarETFESpace->map() ) );

    disp.reset( new vector_Type( dFESpace->map() ) );
    dispRead.reset( new vector_Type( dFESpace->map() ) );

    // deformationGradient_1.reset( new vector_Type( dFESpace->map() ) );
    // deformationGradient_2.reset( new vector_Type( dFESpace->map() ) );
    // deformationGradient_3.reset( new vector_Type( dFESpace->map() ) );

    // reconstructed_1.reset( new vector_Type( dFESpace->map() ) );
    // reconstructed_2.reset( new vector_Type( dFESpace->map() ) );
    // reconstructed_3.reset( new vector_Type( dFESpace->map() ) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement",
                              dFESpace, dispRead, UInt (0) );

    // exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "deformationGradient_1",
    //                           dFESpace, deformationGradient_1, UInt (0) );

    // exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "deformationGradient_2",
    //                           dFESpace, deformationGradient_2, UInt (0) );

    // exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "deformationGradient_3",
    //                           dFESpace, deformationGradient_3, UInt (0) );

    // exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "reconstructed_1",
    //                           dFESpace, reconstructed_1, UInt (0) );

    // exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "reconstructed_2",
    //                           dFESpace, reconstructed_2, UInt (0) );

    // exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "reconstructed_3",
    //                           dFESpace, reconstructed_3, UInt (0) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "jacobianF",
				dScalarFESpace, jacobianF, UInt (0) );

    if( !dataStructure->constitutiveLaw().compare("anisotropic") )
    {
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "family1",
                                dFESpace, family1Read, UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "family2",
                                dFESpace, family2Read, UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "family1Computed",
                                dFESpace, family1Computed, UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "family2Computed",
                                dFESpace, family2Computed, UInt (0) );

        for( UInt i(0); i < dataStructure->numberFibersFamilies( ); i++ )
        {
            std::string familyNumber;
            std::ostringstream number;
            number << ( i + 1 );
            familyNumber = number.str();

            activationFunction[ i ].reset( new vector_Type( dScalarFESpace->map() ) );
            stretch[ i ].reset( new vector_Type( dScalarFESpace->map() ) );
            deformedVector[ i ].reset( new vector_Type( dFESpace->map() ) );
            characteristicAngle[ i ].reset( new vector_Type( dScalarFESpace->map() ) );
            referenceAngle[ i ].reset( new vector_Type( dScalarFESpace->map() ) );

            std::string stringNameA("activation");
            std::string stringNameS("stretch");
            std::string deformedName("deformed");
            std::string angleName("angle");
            std::string refAngleName("refAngle");

            stringNameA  += "-"; stringNameA  += familyNumber;
            stringNameS  += "-"; stringNameS  += familyNumber;
            deformedName += "-"; deformedName += familyNumber;
            angleName += "-"; angleName += familyNumber;
            refAngleName += "-"; refAngleName += familyNumber;

            exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, stringNameA,
                                    dScalarFESpace, activationFunction[ i ], UInt (0) );

            exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, stringNameS,
                                    dScalarFESpace, stretch[ i ], UInt (0) );

            exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, deformedName,
                                dFESpace, deformedVector[ i ], UInt (0) );

            exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, angleName,
                                    dScalarFESpace, characteristicAngle[ i ], UInt (0) );

            exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, refAngleName,
                                    dScalarFESpace, referenceAngle[ i ], UInt (0) );

        }

        *family1Computed = *(vectorInterpolated[ 0 ] );
        *family2Computed = *(vectorInterpolated[ 1 ] );

    }

    exporter->postProcess ( 0.0 );

    //! =================================================================================
    //! Analysis - Istant or Interval
    //! =================================================================================

    MPI_Barrier (MPI_COMM_WORLD);

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
    evaluateNode( elements ( dScalarETFESpace->mesh() ),
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

    std::string const nameField =  dataFile ( "importer/nameField", "displacement");
    UInt i,k;

    for( i=0,k =0; i < numberOfSol; i++, k++ )
    {

        // *reconstructed_1 *=0.0; *reconstructed_2 *=0.0; *reconstructed_3 *=0.0;
        // *deformationGradient_1 *=0.0; *deformationGradient_2 *=0.0; *deformationGradient_3 *=0.0;
        *jacobianF *= 0.0;

        // Reading the solution
        // resetting the element of the vector
        *disp *= 0.0;
        *dispRead *= 0.0;

        if( !dataStructure->constitutiveLaw().compare("anisotropic") )
        {
            *family1 *= 0.0;
            *family2 *= 0.0;

            *family1Read *= 0.0;
            *family2Read *= 0.0;
        }

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
        std::string fam1 = "Family-1";
        std::string fam2 = "Family-2";

        std::cout << "Current reading: " << iterationString << std::endl;

        /*!Definition of the ExporterData, used to load the solution inside the previously defined vectors*/
        LifeV::ExporterData<mesh_Type> solutionDispl  (LifeV::ExporterData<mesh_Type>::VectorField,nameField + "." + iterationString,
                                                       dFESpace, disp, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );


        //Read the variable
        M_importer->readVariable (solutionDispl);
        *dispRead = *disp;

        // Computing references
        // *reconstructed_1 = dFESpace->gradientRecovery (*disp, 0);
        // *reconstructed_2 = dFESpace->gradientRecovery (*disp, 1);
        // *reconstructed_3 = dFESpace->gradientRecovery (*disp, 2);

        // Defining expressions
        ExpressionDefinitions::deformationGradient_Type  F =
            ExpressionDefinitions::deformationGradient ( dETFESpace,  *disp, 0, identity );

        // Definition of C = F^T F
        ExpressionDefinitions::rightCauchyGreenTensor_Type C =
            ExpressionDefinitions::tensorC( transpose(F), F );

        // Definition of J
        ExpressionDefinitions::determinantTensorF_Type J =
            ExpressionDefinitions::determinantF( F );

        // Definition of J^-(2/3) = det( C ) using the isochoric/volumetric decomposition
        ExpressionDefinitions::powerExpression_Type  Jel =
            ExpressionDefinitions::powerExpression( J , (-2.0/3.0) );

        // ExpressionVectorFromNonConstantMatrix< ExpressionDefinitions::deformationGradient_Type,3 ,3 > Fi1( F, 0 );
        // ExpressionVectorFromNonConstantMatrix< ExpressionDefinitions::deformationGradient_Type,3 ,3 > Fi2( F, 1 );
        // ExpressionVectorFromNonConstantMatrix< ExpressionDefinitions::deformationGradient_Type,3 ,3 > Fi3( F, 2 );

        // // Deformation gradient
        // evaluateNode( elements ( dETFESpace->mesh() ),
        //               fakeQuadratureRule,
        //               dETFESpace,
        //               meas_K *  dot ( Fi1, phi_i)
        //               ) >> deformationGradient_1;
        // deformationGradient_1->globalAssemble();
        // *deformationGradient_1 = *deformationGradient_1 / *patchAreaVector;

        // evaluateNode( elements ( dETFESpace->mesh() ),
        //               fakeQuadratureRule,
        //               dETFESpace,
        //               meas_K *  dot ( Fi2, phi_i)
        //               ) >> deformationGradient_2;
        // deformationGradient_2->globalAssemble();
        // *deformationGradient_2 = *deformationGradient_2 / *patchAreaVector;

        // evaluateNode( elements ( dETFESpace->mesh() ),
        //               fakeQuadratureRule,
        //               dETFESpace,
        //               meas_K *  dot ( Fi3, phi_i)
        //               ) >> deformationGradient_3;
        // deformationGradient_3->globalAssemble();
        // *deformationGradient_3 = *deformationGradient_3 / *patchAreaVector;


        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dScalarETFESpace,
                      meas_K *  J * phi_i
                      ) >> jacobianF;
        jacobianF->globalAssemble();
        *jacobianF = *jacobianF / *patchAreaVectorScalar;

        if( !dataStructure->constitutiveLaw().compare("anisotropic") )
        {
            LifeV::ExporterData<mesh_Type> family1Data  (LifeV::ExporterData<mesh_Type>::VectorField,fam1 + "." + iterationString,
                                                         dFESpace, family1, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

            LifeV::ExporterData<mesh_Type> family2Data  (LifeV::ExporterData<mesh_Type>::VectorField,fam2 + "." + iterationString,
                                                         dFESpace, family2, UInt (0), LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

            M_importer->readVariable (family2Data);
            M_importer->readVariable (family1Data);

            *family1Read = *family1;
            *family2Read = *family2;

            // Evaluating quantities that are related to fibers
            for( UInt j(0); j < dataStructure->numberFibersFamilies( ); j++ )
            {
                // Fibers definitions
                ExpressionDefinitions::interpolatedValue_Type fiberIth =
                    ExpressionDefinitions::interpolateFiber( dETFESpace, *(vectorInterpolated[ j ] ) );

                // f /otimes f
                ExpressionDefinitions::outerProduct_Type Mith = ExpressionDefinitions::fiberTensor( fiberIth );

                // Definition of the fourth invariant : I_4^i = C:Mith
                ExpressionDefinitions::stretch_Type IVith = ExpressionDefinitions::fiberStretch( C, Mith );

                // Definition of the fouth isochoric invariant : J^(-2.0/3.0) * I_4^i
                ExpressionDefinitions::isochoricStretch_Type IVithBar =
                    ExpressionDefinitions::isochoricFourthInvariant( Jel, IVith );

                *stretch[ j ] *= 0.0;
                *activationFunction[ j ] *= 0.0;
                *deformedVector[ j ] *= 0.0;

                evaluateNode( elements ( dScalarETFESpace->mesh() ),
                              fakeQuadratureRule,
                              dScalarETFESpace,
                              meas_K * IVithBar  * phi_i
                              ) >> stretch[ j ];
                stretch[ j ]->globalAssemble();
                *(stretch[ j ]) = *(stretch[ j ]) / *patchAreaVectorScalar;

                Real stretch = dataStructure->ithCharacteristicStretch( j ) * dataStructure->ithCharacteristicStretch( j );
                evaluateNode( elements ( dScalarETFESpace->mesh() ),
                              fakeQuadratureRule,
                              dScalarETFESpace,
                              meas_K * atan( IVithBar - value( stretch ) , dataStructure->smoothness() , ( 1 / 3.14159265359 ), ( 1.0/2.0 )  )  * phi_i
                              ) >> activationFunction[ j ];
                activationFunction[ j ]->globalAssemble();
                *(activationFunction[ j ]) = *(activationFunction[ j ]) / *patchAreaVectorScalar;

                evaluateNode( elements ( dETFESpace->mesh() ),
                              fakeQuadratureRule,
                              dETFESpace,
                              meas_K *  dot( ( F * value( dETFESpace ,*vectorInterpolated[ j ] ) )  , phi_i )
                              ) >> deformedVector[ j ];
                deformedVector[ j ]->globalAssemble();
                *deformedVector[ j ] = *deformedVector[ j ] / *patchAreaVector;

                evaluateNode( elements ( dScalarETFESpace->mesh() ),
                              fakeQuadratureRule,
                              dScalarETFESpace,
                              meas_K * dot( normalize( value( dETFESpace, *deformedVector[j] ) ) , value( dETFESpace, *referenceDirectionVector ) )  * phi_i
                              ) >> characteristicAngle[ j ];
                characteristicAngle[ j ]->globalAssemble();
                *(characteristicAngle[ j ]) = *(characteristicAngle[ j ]) / *patchAreaVectorScalar;

                evaluateNode( elements ( dScalarETFESpace->mesh() ),
                              fakeQuadratureRule,
                              dScalarETFESpace,
                              meas_K * dot( normalize( value( dETFESpace, *vectorInterpolated[j] ) ) , value( dETFESpace, *referenceDirectionVector ) )  * phi_i
                              ) >> referenceAngle[ j ];
                referenceAngle[ j ]->globalAssemble();
                *(referenceAngle[ j ]) = *(referenceAngle[ j ]) / *patchAreaVectorScalar;

            }

        }
        exporter->postProcess( dataStructure->dataTime()->initialTime() + k * dataStructure->dataTime()->timeStep() );
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

    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
