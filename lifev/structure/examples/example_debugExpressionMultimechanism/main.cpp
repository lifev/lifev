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


    // Adding the fibers vectors
    // Setting the vector of fibers functions
    for( UInt k(1); k <= pointerToVectorOfFamilies->size( ); k++ )
    {
        stretchesVector[ k - 1 ].reset( new vector_Type( dScalarFESpace->map() ) );
        atanStretchesVector[ k - 1 ].reset( new vector_Type( dScalarFESpace->map() ) );

        // Setting up the name of the function to define the family
        std::string stretchfamily="stretchFamily-";
        std::string familyAtan="atanStretchFamily-";
        // adding the number of the family
        std::string familyNumber;
        std::ostringstream number;
        number << ( k );
        familyNumber = number.str();

        // Name of the function to create
        std::string creationString = stretchfamily + familyNumber;
        std::string creationStringAtan = familyAtan + familyNumber;
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, creationString, dScalarFESpace, stretchesVector[ k-1 ], UInt (0) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, creationStringAtan, dScalarFESpace, atanStretchesVector[ k-1 ], UInt (0) );
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
    std::cout << "Norm of the J       : " << JacobianZero->normInf() << std::endl;
    std::cout << "Norm of the J_0(ta) : " << JacobianZeroA->normInf() << std::endl;
    std::cout << "Norm of the J_a     : " << JacobianA->normInf() << std::endl;


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
