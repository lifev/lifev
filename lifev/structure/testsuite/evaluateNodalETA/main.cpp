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

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/structure/fem/ExpressionDefinitions.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

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

    void checkResults( const Real a, const Real b, const Real c,
		       const Real d, const Real e, const Real f,
		       const Real g, const Real h, const Real i,
		       const Real l, const Real m, const Real n,
		       const Real o);


private:
    struct Private;
    boost::shared_ptr<Private> parameters;
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

    // General typedefs
    typedef VectorEpetra                                          vector_Type;
    typedef boost::shared_ptr<vector_Type>                        vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >                       solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                  solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                     scalarETFESpacePtr_Type;

    // typedefs for fibers
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fiberFunction_Type;
    typedef boost::shared_ptr<fiberFunction_Type> fiberFunctionPtr_Type;

    typedef std::vector<fiberFunctionPtr_Type>                          vectorFiberFunction_Type;
    typedef boost::shared_ptr<vectorFiberFunction_Type>                 vectorFiberFunctionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 vectorInterpolatedFibers_Type;



    bool verbose = (parameters->comm->MyPID() == 0);

    //! dataElasticStructure for parameters
    GetPot dataFile ( parameters->data_file_name.c_str() );
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    //! Read and partition mesh
    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr (new RegionMesh<LinearTetra> ( ( parameters->comm ) ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, parameters->comm );

    //! Functional spaces - needed for the computations of the gradients
    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (meshPart, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );

    //! Scalar ETFEspace to evaluate scalar quantities
    solidFESpacePtr_Type dScalarFESpace ( new solidFESpace_Type (meshPart, dOrder, 1, parameters->comm) );
    scalarETFESpacePtr_Type dScalarETFESpace ( new scalarETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );


    if (verbose)
    {
        std::cout << std::endl;
    }

    // The vector where the solution will be stored
    vectorPtr_Type solidDisp (new vector_Type (dFESpace->map(), LifeV::Unique ) );
    *solidDisp *= 0.0;


    // Set up of the fiber field

    // Setting the fibers
    vectorFiberFunctionPtr_Type pointerToVectorOfFamilies( new vectorFiberFunction_Type( ) );
    (*pointerToVectorOfFamilies).resize( dataStructure->numberFibersFamilies( ) );

    if( verbose )
        std::cout << "Size of the number of families: " << (*pointerToVectorOfFamilies).size() << std::endl;

    fibersDirectionList setOfFiberFunctions;
    setOfFiberFunctions.setupFiberDefinitions( dataStructure->numberFibersFamilies( ) );

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


    vectorInterpolatedFibers_Type vectorInterpolated(0);

    // Interpolating fiber fields
    vectorInterpolated.resize( (*pointerToVectorOfFamilies).size() );

    for( UInt k(0); k < (*pointerToVectorOfFamilies).size(); k++ )
    {
        vectorInterpolated[ k ].reset( new vector_Type( dFESpace->map() ) );
        dFESpace->interpolate ( *( ( *(pointerToVectorOfFamilies) )[ k ] ) ,
                                * ( ( vectorInterpolated )[ k ] ),
                                0.0 );
    }



    //! 6. Post-processing setting
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "hdf5");
    std::string const nameExporter =  dataFile ( "exporter/name", "jacobian");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new hdf5Filter_Type ( dataFile, nameExporter ) );
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

    exporter->setMeshProcId (dFESpace->mesh(), dFESpace->map().comm().MyPID() );

    // Scalar vector to have scalar quantities
    vectorPtr_Type patchAreaVector ( new vector_Type ( dETFESpace->map(),  LifeV::Unique ) );
    vectorPtr_Type patchAreaVectorScalar ( new vector_Type ( dScalarETFESpace->map(), Unique ) );
    vectorPtr_Type JacobianZero( new vector_Type( dScalarETFESpace->map(), Unique ) );
    vectorPtr_Type JacobianZeroA( new vector_Type( dScalarETFESpace->map(), Unique ) );
    vectorPtr_Type JacobianA( new vector_Type( dScalarETFESpace->map(), Unique ) );

    // vectorPtr_Type jacobianVector ( new vector_Type (solid.displacement(),  LifeV::Unique ) );
    // vectorPtr_Type jacobianReference ( new vector_Type (solid.displacement(),  LifeV::Unique ) );


    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "J_0", dScalarFESpace, JacobianA, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "J_0(ta)", dScalarFESpace, JacobianZeroA, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "Ja", dScalarFESpace, JacobianA, UInt (0) );

    vectorInterpolatedFibers_Type stretchesVector(0);
    stretchesVector.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type atanStretchesVector(0);
    atanStretchesVector.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type scalarExpressionMultimechanism(0);
    scalarExpressionMultimechanism.resize( (*pointerToVectorOfFamilies).size() );

    // Deformation Gradient Fa
    vectorInterpolatedFibers_Type Fa_col1(0);
    Fa_col1.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type Fa_col2(0);
    Fa_col2.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type Fa_col3(0);
    Fa_col3.resize( (*pointerToVectorOfFamilies).size() );

    // Outer product Tensor
    vectorInterpolatedFibers_Type Mith_col1(0);
    Mith_col1.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type Mith_col2(0);
    Mith_col2.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type Mith_col3(0);
    Mith_col3.resize( (*pointerToVectorOfFamilies).size() );

    // FzeroAminusT
    vectorInterpolatedFibers_Type FzeroAminusT_col1(0);
    FzeroAminusT_col1.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type FzeroAminusT_col2(0);
    FzeroAminusT_col2.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type FzeroAminusT_col3(0);
    FzeroAminusT_col3.resize( (*pointerToVectorOfFamilies).size() );

    // Piola Kirchhoff
    vectorInterpolatedFibers_Type P_col1(0);
    P_col1.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type P_col2(0);
    P_col2.resize( (*pointerToVectorOfFamilies).size() );

    vectorInterpolatedFibers_Type P_col3(0);
    P_col3.resize( (*pointerToVectorOfFamilies).size() );


    // Adding the fibers vectors
    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
    {
        stretchesVector[ k - 1 ].reset( new vector_Type( dScalarFESpace->map() ) );
        atanStretchesVector[ k - 1 ].reset( new vector_Type( dScalarFESpace->map() ) );
        scalarExpressionMultimechanism[ k - 1 ].reset( new vector_Type( dScalarFESpace->map() ) );

	Fa_col1[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        Fa_col2[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        Fa_col3[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );

	Mith_col1[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        Mith_col2[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        Mith_col3[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );

	FzeroAminusT_col1[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        FzeroAminusT_col2[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        FzeroAminusT_col3[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );

        P_col1[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        P_col2[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );
        P_col3[ k - 1 ].reset( new vector_Type( dFESpace->map() ) );


        // Setting up the name of the function to define the family
        std::string stretchfamily="stretchFamily-";
        std::string familyAtan="atanStretchFamily-";
        std::string familyScalar="scalarQuantityFamily-";

        std::string Fa1="Fa1-";
        std::string Fa2="Fa2-";
        std::string Fa3="Fa3-";

        std::string Mith1="Mith1-";
        std::string Mith2="Mith2-";
        std::string Mith3="Mith3-";

        std::string FzeroAminusT1="FzeroAminusT1-";
        std::string FzeroAminusT2="FzeroAminusT2-";
        std::string FzeroAminusT3="FzeroAminusT3-";

        std::string P1="P1-";
        std::string P2="P2-";
        std::string P3="P3-";


        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = stretchfamily + familyNumber;
        std::string creationStringAtan = familyAtan + familyNumber;
        std::string creationStringScalar = familyScalar + familyNumber;

        std::string creationStringFa1 = Fa1 + familyNumber;
        std::string creationStringFa2 = Fa2 + familyNumber;
        std::string creationStringFa3 = Fa3 + familyNumber;

        std::string creationStringFzeroAminusT1 = FzeroAminusT1 + familyNumber;
        std::string creationStringFzeroAminusT2 = FzeroAminusT2 + familyNumber;
        std::string creationStringFzeroAminusT3 = FzeroAminusT3 + familyNumber;

        std::string creationStringMith1 = Mith1 + familyNumber;
        std::string creationStringMith2 = Mith2 + familyNumber;
        std::string creationStringMith3 = Mith3 + familyNumber;

        std::string creationStringP1 = P1 + familyNumber;
        std::string creationStringP2 = P2 + familyNumber;
        std::string creationStringP3 = P3 + familyNumber;

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, creationString, dScalarFESpace, stretchesVector[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, creationStringAtan, dScalarFESpace, atanStretchesVector[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, creationStringScalar, dScalarFESpace, scalarExpressionMultimechanism[ k-1 ], UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringFa1, dFESpace, Fa_col1[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringFa2, dFESpace, Fa_col2[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringFa3, dFESpace, Fa_col3[ k-1 ], UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringFzeroAminusT1, dFESpace, FzeroAminusT_col1[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringFzeroAminusT2, dFESpace, FzeroAminusT_col2[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringFzeroAminusT3, dFESpace, FzeroAminusT_col3[ k-1 ], UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringMith1, dFESpace, Mith_col1[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringMith2, dFESpace, Mith_col2[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringMith3, dFESpace, Mith_col3[ k-1 ], UInt (0) );

        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringP1, dFESpace, P_col1[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringP2, dFESpace, P_col2[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, creationStringP3, dFESpace, P_col3[ k-1 ], UInt (0) );

    }


    exporter->postProcess ( 0.0 );

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

    // //Compute the gradient along X of the displacement field
    // *grDisplX = dFESpace->gradientRecovery (*solidDisp, 0);

    // //Compute the gradient along Y of the displacement field
    // *grDisplY = dFESpace->gradientRecovery (*solidDisp, 1);

    // //Compute the gradient along Z of the displacement field
    // *grDisplZ = dFESpace->gradientRecovery (*solidDisp, 2);

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



    // Definition of F
    ExpressionDefinitions::deformationGradient_Type  F =
        ExpressionDefinitions::deformationGradient ( dETFESpace,  *solidDisp, 0, identity );

    // Definition of J
    ExpressionDefinitions::determinantTensorF_Type J = ExpressionDefinitions::determinantF( F );

    // Definition of F_0(ta)
    ExpressionDefinitions::deformationGradient_Type  FzeroA =
        ExpressionDefinitions::deformationGradient ( dETFESpace,  *solidDisp, 0, identity );

    // Definition of J
    ExpressionDefinitions::determinantTensorF_Type JzeroA = ExpressionDefinitions::determinantF( FzeroA );

    // Definition of J_0(ta) ^ -1
    ExpressionDefinitions::powerExpression_Type  JzeroAminus1 = ExpressionDefinitions::powerExpression( JzeroA , (-1.0) );

    // Definition of J_a
    ExpressionMultimechanism::activatedDeterminantF_Type Ja =
        ExpressionMultimechanism::activateDeterminantF( J, JzeroAminus1 );

    // Definition of J_a^(-2.0/3.0)
    ExpressionMultimechanism::activePowerExpression_Type  JactiveEl =
        ExpressionMultimechanism::activePowerExpression( Ja , (-2.0/3.0) );

    // Definition of C = F^T F
    ExpressionDefinitions::rightCauchyGreenTensor_Type C =
        ExpressionDefinitions::tensorC( transpose(F), F );

    LifeChrono chrono;
    chrono.start();

    evaluateNode( elements ( dScalarETFESpace->mesh() ),
                  fakeQuadratureRule,
                  dScalarETFESpace,
                  meas_K * phi_i
                  ) >> patchAreaVectorScalar;

    evaluateNode( elements ( dScalarETFESpace->mesh() ),
                  fakeQuadratureRule,
                  dScalarETFESpace,
                  meas_K * J  * phi_i
                  ) >> JacobianZero;

    *JacobianZero = *JacobianZero / *patchAreaVectorScalar;


    evaluateNode( elements ( dScalarETFESpace->mesh() ),
                  fakeQuadratureRule,
                  dScalarETFESpace,
                  meas_K * JzeroA  * phi_i
                  ) >> JacobianZeroA;

    *JacobianZeroA = *JacobianZeroA / *patchAreaVectorScalar;


    evaluateNode( elements ( dScalarETFESpace->mesh() ),
                  fakeQuadratureRule,
                  dScalarETFESpace,
                  meas_K * JactiveEl  * phi_i
                  ) >> JacobianA;

    *JacobianA = *JacobianA / *patchAreaVectorScalar;

    //Extracting the tensions

    // The patch area in vectorial form
    ExpressionVectorFromNonConstantScalar<ExpressionMeas, 3  > vMeas( meas_K );
    evaluateNode( elements ( dETFESpace->mesh() ),
		  fakeQuadratureRule,
		  dETFESpace,
		  dot( vMeas , phi_i )
		  ) >> patchAreaVector;

    for( UInt i(0); i < pointerToVectorOfFamilies->size( ); i++ )
    {

        // Fibers definitions
        ExpressionDefinitions::interpolatedValue_Type fiberIth =
            ExpressionDefinitions::interpolateFiber( dETFESpace, *(vectorInterpolated[ i ] ) );

        // FzeroA ( for zero displacement and activation at 1 it's equal to F )
        ExpressionDefinitions::deformationGradient_Type  ithFzeroA =
            ExpressionDefinitions::deformationGradient ( dETFESpace,  *solidDisp, 0, identity );

        // Definiton of FzeroA^{-T}
        ExpressionDefinitions::minusTransposedTensor_Type FzeroAminusT =
            ExpressionDefinitions::minusT( ithFzeroA );

        // Definition of FzeroA^{-1}
        ExpressionDefinitions::inverseTensor_Type FzeroAminus1 =
            ExpressionDefinitions::inv( ithFzeroA );

        // Definition of Ca = F_0^{-T}(ta) * C_0 * F_0^{-1}(ta)
        ExpressionMultimechanism::rightCauchyGreenMultiMechanism_Type Ca =
            ExpressionMultimechanism::activationRightCauchyGreen( FzeroAminusT, C, FzeroAminus1 );

        // Definition of F_0(ta) * f_0
        ExpressionMultimechanism::activatedFiber_Type activeIthFiber =
            ExpressionMultimechanism::activateFiberDirection( ithFzeroA, fiberIth );

        // Definition of M = f_a \otimes f_a
        ExpressionMultimechanism::activeOuterProduct_Type Mith =
            ExpressionMultimechanism::activeOuterProduct( activeIthFiber );

        // Definition of the stretch with respect the activation configuration
        ExpressionMultimechanism::activeStretch_Type IVith =
            ExpressionMultimechanism::activeFiberStretch( Ca, Mith );

        // Definition of the fouth isochoric invariant : J^(-2.0/3.0) * I_4^i
        ExpressionMultimechanism::activeIsochoricStretch_Type IVithBar =
            ExpressionMultimechanism::activeIsochoricFourthInvariant( JactiveEl, IVith );

        evaluateNode( elements ( dScalarETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dScalarETFESpace,
                      meas_K * IVithBar  * phi_i
                      ) >> stretchesVector[ i ];

        *( stretchesVector[ i ] ) = *( stretchesVector[ i ] ) / *patchAreaVectorScalar;

        evaluateNode( elements ( dScalarETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dScalarETFESpace,
                      meas_K * atan( IVithBar - value( dataStructure->ithCharacteristicStretch( i ) ),
                                     dataStructure->smoothness() , ( 1.0 / 3.14159265359 ), (1.0/2.0) ) * phi_i
                      ) >> atanStretchesVector[ i ];

        *( atanStretchesVector[ i ] ) = *( atanStretchesVector[ i ] ) / *patchAreaVectorScalar;

	Real stretch = dataStructure->ithCharacteristicStretch( i );
	Real pi = 3.14159265359;

        evaluateNode( elements ( dScalarETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dScalarETFESpace,
                      meas_K * 
		      JzeroA * atan( IVithBar - value( stretch ) , dataStructure->smoothness(), ( 1 / pi ), ( 1.0/2.0 )  )  * JactiveEl *
                      (value( 2.0 ) * value( dataStructure->ithStiffnessFibers( i ) ) * JactiveEl * ( IVithBar - value( stretch ) ) *
                       exp( value( dataStructure->ithNonlinearityFibers( i ) ) * ( IVithBar- value( stretch ) ) * ( IVithBar- value( stretch ) )  ) )  * phi_i
                      ) >> scalarExpressionMultimechanism[ i ];
        *( scalarExpressionMultimechanism[ i ] ) = *( scalarExpressionMultimechanism[ i ] ) / *patchAreaVectorScalar;


	// exporting the components of the Piola-Kirchhoff tensor
	// Definition of the expression definition the portion of the Piola-Kirchhoff 

	ExpressionMultimechanism::deformationActivatedTensor_Type Fa = 
	  ExpressionMultimechanism::createDeformationActivationTensor( F , FzeroAminus1);

	typedef ExpressionProduct<ExpressionMultimechanism::deformationActivatedTensor_Type,
				  ExpressionMultimechanism::activeOuterProduct_Type>  productFaMith_Type;

	typedef ExpressionProduct< productFaMith_Type,
				   ExpressionDefinitions::minusTransposedTensor_Type> firstPartPiolaMultimech_Type;

	productFaMith_Type FaMith( Fa, Mith );
	firstPartPiolaMultimech_Type firstPartPiola( FaMith, FzeroAminusT );

	//Extract columns
	ExpressionVectorFromNonConstantMatrix< ExpressionMultimechanism::deformationActivatedTensor_Type,3 ,3 > Fa_i1( Fa, 0 );
	ExpressionVectorFromNonConstantMatrix< ExpressionMultimechanism::deformationActivatedTensor_Type,3 ,3 > Fa_i2( Fa, 1 );
	ExpressionVectorFromNonConstantMatrix< ExpressionMultimechanism::deformationActivatedTensor_Type,3 ,3 > Fa_i3( Fa, 2 );

	ExpressionVectorFromNonConstantMatrix< ExpressionMultimechanism::activeOuterProduct_Type,3 ,3 > Mith_i1( Mith, 0 );
	ExpressionVectorFromNonConstantMatrix< ExpressionMultimechanism::activeOuterProduct_Type,3 ,3 > Mith_i2( Mith, 1 );
	ExpressionVectorFromNonConstantMatrix< ExpressionMultimechanism::activeOuterProduct_Type,3 ,3 > Mith_i3( Mith, 2 );

	ExpressionVectorFromNonConstantMatrix< ExpressionDefinitions::minusTransposedTensor_Type,3 ,3 > FzeroAminusT_i1( FzeroAminusT, 0 );
	ExpressionVectorFromNonConstantMatrix< ExpressionDefinitions::minusTransposedTensor_Type,3 ,3 > FzeroAminusT_i2( FzeroAminusT, 1 );
	ExpressionVectorFromNonConstantMatrix< ExpressionDefinitions::minusTransposedTensor_Type,3 ,3 > FzeroAminusT_i3( FzeroAminusT, 2 );

	ExpressionVectorFromNonConstantMatrix< firstPartPiolaMultimech_Type,3 ,3  > P_i1( firstPartPiola, 0 );
	ExpressionVectorFromNonConstantMatrix< firstPartPiolaMultimech_Type,3 ,3 > P_i2( firstPartPiola, 1 );
	ExpressionVectorFromNonConstantMatrix< firstPartPiolaMultimech_Type,3 ,3 > P_i3( firstPartPiola, 2 );

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( Fa_i1, phi_i) 
                      ) >> Fa_col1[ i ];
	*(Fa_col1[i]) = *(Fa_col1[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( Fa_i2, phi_i) 
                      ) >> Fa_col2[ i ];
	*(Fa_col2[i]) = *(Fa_col2[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( Fa_i3, phi_i) 
                      ) >> Fa_col3[ i ];
	*(Fa_col3[i]) = *(Fa_col3[i]) / *patchAreaVector;


        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( Mith_i1, phi_i) 
                      ) >> Mith_col1[ i ];
	*(Mith_col1[i]) = *(Mith_col1[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( Mith_i2, phi_i) 
                      ) >> Mith_col2[ i ];
	*(Mith_col2[i]) = *(Mith_col2[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( Mith_i3, phi_i) 
                      ) >> Mith_col3[ i ];
	*(Mith_col3[i]) = *(Mith_col3[i]) / *patchAreaVector;


        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( FzeroAminusT_i1, phi_i) 
                      ) >> FzeroAminusT_col1[ i ];
	*(FzeroAminusT_col1[i]) = *(FzeroAminusT_col1[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( FzeroAminusT_i2, phi_i) 
                      ) >> FzeroAminusT_col2[ i ];
	*(FzeroAminusT_col2[i]) = *(FzeroAminusT_col2[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( FzeroAminusT_i3, phi_i) 
                      ) >> FzeroAminusT_col3[ i ];
	*(FzeroAminusT_col3[i]) = *(FzeroAminusT_col3[i]) / *patchAreaVector;


        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( P_i1, phi_i) 
                      ) >> P_col1[ i ];
	*(P_col1[i]) = *(P_col1[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( P_i2, phi_i) 
                      ) >> P_col2[ i ];
	*(P_col2[i]) = *(P_col2[i]) / *patchAreaVector;

        evaluateNode( elements ( dETFESpace->mesh() ),
                      fakeQuadratureRule,
                      dETFESpace,
                      meas_K *  dot ( P_i3, phi_i) 
                      ) >> P_col3[ i ];
	*(P_col3[i]) = *(P_col3[i]) / *patchAreaVector;
		
    }

    exporter->postProcess ( 1.0 );

    if (verbose )
    {
        std::cout << "Analysis Completed!" << std::endl;
    }

    //Closing files
    exporter->closeFile();

    if (verbose )
    {
        std::cout << "finished" << std::endl;
    }

    checkResults( patchAreaVectorScalar->normInf(), patchAreaVectorScalar->norm2(),
		  patchAreaVector->normInf(), patchAreaVector->norm2(),
		  JacobianZero->normInf(), JacobianZeroA->normInf(), JacobianA->normInf(),
		  atanStretchesVector[0]->normInf(),atanStretchesVector[1]->normInf(),
		  scalarExpressionMultimechanism[0]->norm2(), scalarExpressionMultimechanism[1]->norm2(),
		  Mith_col1[0]->normInf(), Mith_col1[1]->normInf());

    MPI_Barrier (MPI_COMM_WORLD);
    //!---------------------------------------------.-----------------------------------------------------
}

void
Structure::checkResults(const Real a, const Real b, const Real c,
			const Real d, const Real e, const Real f,
			const Real g, const Real h, const Real i,
			const Real l, const Real m, const Real n,
			const Real o)
{

  if( std::fabs( a - 0.072916) < 1e-6 &&  
      std::fabs( b - 0.601692) < 1e-6  && 
      std::fabs( c - 0.072916) < 1e-6  && 
      std::fabs( d - 1.042161) < 1e-6  && 
      std::fabs( e - 1 ) < 1e-6 && 
      std::fabs( f - 1 ) < 1e-6 && 
      std::fabs( g - 1 ) < 1e-6 &&
      std::fabs( h - 3.18309e-05 ) < 1e-9 &&
      std::fabs( i - 3.18309e-05 ) < 1e-9 &&
      std::fabs( l - 6968.505143 ) < 1e-6 &&
      std::fabs( m - 6968.505143 ) < 1e-6 &&
      std::fabs( n - 0.433012 ) < 1e-6 &&
      std::fabs( o - 0.433012 ) < 1e-6 )
    {
      std::cout << "Correct result! "<< std::endl;
      returnValue = EXIT_SUCCESS;
    }

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
