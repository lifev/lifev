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
   \file fefct.cpp
   \author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   \date 2011-05-02
 */

// ===================================================
//! Includes
// ===================================================

#include "fefct.hpp"
#include "user_fun.hpp"

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

// ===================================================
//!              Standard functions
// ===================================================

Real UOne( const Real& /* t */,
           const Real& /* x */,
           const Real& /* y */,
           const Real& /* z */,
           const ID&   /* icomp */)
{
    return 1.;
}

Real UZero( const Real& /* t */,
            const Real& /* x */,
            const Real& /* y */,
            const Real& /* z */,
            const ID&   /* icomp */)
{
    return 0.;
}


// ===================================================
//!                  Private Members
// ===================================================

struct fefct::Private
{
    Private() {}

    std::string    data_file_name;
    std::string    discretization_section;

    boost::shared_ptr<Epetra_Comm>   comm;

};

// ===================================================
//!                  Constructors
// ===================================================

fefct::fefct( int argc,
              char** argv )
        : Members( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "fefct";

#ifdef EPETRA_MPI
    std::cout << "Epetra Initialization" << std::endl;
    Members->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    Members->comm.reset( new Epetra_SerialComm() );
#endif

}

// ===================================================
//!                      Methods
// ===================================================

Real
fefct::run()
{
    LifeChrono chronoTotal;
    LifeChrono chronoReadAndPartitionMesh;
    LifeChrono chronoFiniteElementSpace;
    LifeChrono chronoFEFieldAndFEFct;
    LifeChrono chronoAddFields;
    LifeChrono chronoProblem;
    LifeChrono chronoProcess;
    LifeChrono chronoError;

    // Start chronoTotal for measure the total time for the computation
    chronoTotal.start();

    // Reading from data file
    GetPot dataFile( Members->data_file_name.c_str() );

    // Create the leader process, i.e. the process with MyPID equal to zero
    bool isLeader = ( Members->comm->MyPID() == 0 );

    if ( isLeader )
        std::cout << "The FEField and FEFunction test" << std::endl;

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes
    chronoReadAndPartitionMesh.start();

    // Create the mesh file handler
    MeshData meshData;

    // Set up the mesh file
    meshData.setup( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the the mesh
    regionMeshPtr_Type fullMeshPtr( new regionMesh_Type );

    // Set up the mesh
    readMesh( *fullMeshPtr, meshData );

    // Partition the mesh using ParMetis
    MeshPartitioner< regionMesh_Type >  meshPart( fullMeshPtr, Members->comm );

    // Stop chronoReadAndPartitionMesh
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    if ( isLeader )
        std::cout << "Time for read and partition the mesh " <<
                  chronoReadAndPartitionMesh.diff() << std::endl << std::flush;

    // Create the FEField and FEFunction spaces

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces
    chronoFiniteElementSpace.start();

    // First FEScalarField parameters
    const ReferenceFE*    refFE_scalarField1 ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_scalarField1    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_scalarField1  ( static_cast<QuadratureRule*>(NULL) );

    refFE_scalarField1 = &feTetraP0;
    qR_scalarField1    = &quadRuleTetra15pt;
    bdQr_scalarField1  = &quadRuleTria4pt;

    // Second FEScalarField parameters
    const ReferenceFE*    refFE_scalarField2 ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_scalarField2    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_scalarField2  ( static_cast<QuadratureRule*>(NULL) );

    refFE_scalarField2 = &feTetraP1;
    qR_scalarField2    = &quadRuleTetra15pt;
    bdQr_scalarField2  = &quadRuleTria4pt;

    // FEVectorField parameters
    const ReferenceFE*    refFE_vectorField ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_vectorField    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_vectorField  ( static_cast<QuadratureRule*>(NULL) );

    refFE_vectorField = &feTetraP0;
    qR_vectorField    = &quadRuleTetra15pt;
    bdQr_vectorField  = &quadRuleTria4pt;

    // Parameters for function visualization
    const ReferenceFE*    refFE_fct ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR_fct    ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr_fct  ( static_cast<QuadratureRule*>(NULL) );

    refFE_fct = &feTetraP0;
    qR_fct    = &quadRuleTetra15pt;
    bdQr_fct  = &quadRuleTria4pt;

    // Finite element space of the first scalar field
    FESpacePtr_Type scalarField1_FESpace ( new FESpace_Type ( meshPart, *refFE_scalarField1, 
                                                              *qR_scalarField1, *bdQr_scalarField1,
                                                              1, Members->comm ) );

    // Finite element space of the second scalar field
    FESpacePtr_Type scalarField2_FESpace ( new FESpace_Type ( meshPart, *refFE_scalarField2,
                                                              *qR_scalarField2, *bdQr_scalarField2,
                                                              1, Members->comm ) );

    // Finite element space of the vector field
    FESpacePtr_Type vectorField_FESpace ( new FESpace_Type ( meshPart, *refFE_vectorField,
                                                             *qR_vectorField, *bdQr_vectorField,
                                                             3, Members->comm ) );

    // Finite element space for the function visualization
    FESpacePtr_Type function_FESpace ( new FESpace_Type ( meshPart, *refFE_fct, *qR_fct,
                                                          *bdQr_fct, 1, Members->comm ) );

    // Stop chronoFiniteElementSpace
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace
    if ( isLeader )
        std::cout << "Time for create the finite element spaces " <<
                  chronoFiniteElementSpace.diff() << std::endl << std::flush;

    // Start chronoFEFieldAndFEFct for measure the total time for create the fields and the function
    chronoFEFieldAndFEFct.start();

    // First scalar field
    FEScalarFieldPtr_Type scalarField1 ( new FEScalarField_Type ( scalarField1_FESpace ) );
    
    // Create a dummy epetra vector
    scalarField1->getVector() = 1. ;

    // Second scalar field
    FEScalarFieldPtr_Type scalarField2 ( new FEScalarField_Type ( scalarField2_FESpace ) );
    scalarField2->getVector() = 2.;

    // Vector field
    FEVectorFieldPtr_Type vectorField ( new FEVectorField_Type ( vectorField_FESpace ) );
    vectorField->getVector() = 3.;

    // Function
    dataProblem::MyFun function;

    // Scalar field for visualize the function
    FEScalarFieldPtr_Type scalarFieldFunction ( new FEScalarField_Type ( function_FESpace ) );

    // Stop chronoFEFieldAndFEFct
    chronoFEFieldAndFEFct.stop();

    // The leader process print chronoFEFieldAndFEFct
    if ( isLeader )
        std::cout << "Time for create the fields and the function " <<
                  chronoFEFieldAndFEFct.diff() << std::endl;

    // Add the fields to the function

    // Start chronoAddFields
    chronoAddFields.start();

    // Add the first scalar field
    function.addScalarField ( scalarField1 );

    // Add the second scalar field
    function.addScalarField ( scalarField2 );

    // Add the vector field
    function.addVectorField ( vectorField );

    // Stop chronoAddFields
    chronoAddFields.stop();

    // The leader process print chronoAddFields
    if ( isLeader )
        std::cout << "Time for create the add the fields to the function " <<
                  chronoAddFields.diff() << std::endl;

    // Set the exporter for the fields
    exporterPtr_Type exporter;

    // Type of the exporter
    std::string const exporterType =  dataFile( "exporter/type", "hdf5");

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare("hdf5") == 0 )
    {
        exporter.reset( new ExporterHDF5< regionMesh_Type > ( dataFile, dataFile( "exporter/file_name", "FieldFunction" ) ) );

        // Set directory where to save the solution
        exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

        exporter->setMeshProcId( meshPart.meshPartition(), Members->comm->MyPID() );
    }
    else
#endif
    {
        if ( exporterType.compare("none") == 0 )
        {
            exporter.reset( new ExporterEmpty< regionMesh_Type > ( dataFile, dataFile( "exporter/file_name", "FieldFunction" ) ) );

            // Set directory where to save the solution
            exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPart.meshPartition(), Members->comm->MyPID() );
        }
        else
        {
            exporter.reset( new ExporterEnsight< regionMesh_Type > ( dataFile, dataFile( "exporter/file_name", "FieldFunction" ) ) );

            // Set directory where to save the solution
            exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPart.meshPartition(), Members->comm->MyPID() );
        }
    }

    const std::string fieldName = dataFile( "exporter/name_field", "Field" );

    // Add the first scalar field to the exporter
    exporter->addVariable( ExporterData< regionMesh_Type >::ScalarField,
                           fieldName + "1",
                           scalarField1_FESpace,
                           scalarField1->getVectorPtr(),
                           static_cast<UInt>( 0 ),
                           ExporterData< regionMesh_Type >::UnsteadyRegime,
                           ExporterData< regionMesh_Type >::Cell );

    // Add the second scalar field to the exporter
    exporter->addVariable( ExporterData< regionMesh_Type >::ScalarField,
                           fieldName + "2",
                           scalarField2_FESpace,
                           scalarField2->getVectorPtr(),
                           static_cast<UInt>( 0 ),
                           ExporterData< regionMesh_Type >::UnsteadyRegime );

    // Add the vector field to the exporter
    exporter->addVariable( ExporterData< regionMesh_Type >::VectorField,
                           fieldName + "3",
                           vectorField_FESpace,
                           vectorField->getVectorPtr(),
                           static_cast<UInt>( 0 ),
                           ExporterData< regionMesh_Type >::UnsteadyRegime,
                           ExporterData< regionMesh_Type >::Cell );

    const std::string fctName =  dataFile( "exporter/name_fct", "Function" );

    // Add the scalar field representing the function to the exporter
    exporter->addVariable( ExporterData< regionMesh_Type >::ScalarField,
                           fctName + "1",
                           function_FESpace,
                           scalarFieldFunction->getVectorPtr(),
                           static_cast<UInt>( 0 ),
                           ExporterData< regionMesh_Type >::UnsteadyRegime,                           
                           ExporterData< regionMesh_Type >::Cell );

    // Interpolate the value of the function
    function.interpolate( *scalarFieldFunction );

    // Export all the solutions
    exporter->postProcess(0);

    // Define the dummy element and coordinate for the evaluation of the function
    UInt iElem = 0;
    std::vector<Real> point(3,0);

    // Compute the error
    Real error = 0;
    error = std::fabs( ( std::sin(1.) + 4. ) / 3. - function.eval( iElem, point, 0. ) );

    // Return the error
    return error;

} // run
