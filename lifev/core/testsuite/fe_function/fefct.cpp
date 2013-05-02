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
/*
  This simple test show how to use FEField (scalar and vector) and FEFunction.
  After reading a mesh and creating the FESpaces the programs create three FEFields:
  the first is a Raviart-Thomas field, the second is a P1 vector field and the third is a
  P0 scalar field. Then the function is created, using the one written in user_fun.hpp,
  and the fields are added. After that there is a computation of an error between the FEFunction
  and its "free function" equivalence. Since we want to export the FEFunction we need to
  interpolate it on a finite element space, in this case P0, so we call the interpolate over a new
  field. The last part of the code contains the standard part to export the fields and the
  interpolated function.
*/

// ===================================================
//!                     Includes
// ===================================================

#include "fefct.hpp"
#include "user_fun.hpp"

// ===================================================
//!                    Namespaces
// ===================================================

using namespace LifeV;

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

fefct::fefct ( int argc,
               char** argv )
    : Members ( new Private )
{
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "fefct";

#ifdef EPETRA_MPI
    std::cout << "Epetra Initialization" << std::endl;
    Members->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size (MPI_COMM_WORLD, &ntasks);
#else
    Members->comm.reset ( new Epetra_SerialComm() );
#endif

}

// ===================================================
//!                      Methods
// ===================================================

Real
fefct::run()
{
    using namespace dataProblem;

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
    GetPot dataFile ( Members->data_file_name.c_str() );

    // Create the displayer
    Displayer displayer ( Members->comm );

    displayer.leaderPrint ( "The FEField and FEFunction test\n" );

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes
    chronoReadAndPartitionMesh.start();

    // Create the mesh file handler
    MeshData meshData;

    // Set up the mesh file
    meshData.setup ( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the the mesh
    regionMeshPtr_Type fullMeshPtr ( new regionMesh_Type ( Members->comm ) );

    // Set up the mesh
    readMesh ( *fullMeshPtr, meshData );

    // Create local mesh using ParMetis partitioner
    regionMeshPtr_Type meshPtr;
    {
        // local scope to properly delete the meshPart object
        MeshPartitioner < regionMesh_Type >  meshPart ( fullMeshPtr, Members->comm );
        meshPtr = meshPart.meshPartition();
    }

    // Stop chronoReadAndPartitionMesh
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    displayer.leaderPrintMax ( "Time for read and partition the mesh ",
                               chronoReadAndPartitionMesh.diff() );

    // Create the FEField and FEFunction spaces

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces
    chronoFiniteElementSpace.start();

    // Finite element space of the first scalar field - RT0
    FESpacePtr_Type scalarField1_FESpace ( new FESpace_Type ( meshPtr, feTetraP0, quadRuleTetra15pt,
                                                              quadRuleTria4pt, 1, Members->comm ) );

    // Finite element space of the second scalar field - P1
    FESpacePtr_Type scalarField2_FESpace ( new FESpace_Type ( meshPtr, feTetraP1, quadRuleTetra15pt,
                                                              quadRuleTria4pt, 1, Members->comm ) );

    // Finite element space of the vector field - P0
    FESpacePtr_Type vectorField_FESpace ( new FESpace_Type ( meshPtr, feTetraP0, quadRuleTetra15pt,
                                                             quadRuleTria4pt, 3, Members->comm ) );

    // Finite element space for the function visualization - P0
    FESpacePtr_Type function_FESpace ( new FESpace_Type ( meshPtr, feTetraP0, quadRuleTetra15pt,
                                                          quadRuleTria4pt, 1, Members->comm ) );

    // Stop chronoFiniteElementSpace
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace
    displayer.leaderPrintMax ( "Time for create the finite element spaces ",
                               chronoFiniteElementSpace.diff() );

    // Start chronoFEFieldAndFEFct for measure the total time for create the fields and the function
    chronoFEFieldAndFEFct.start();

    // First scalar field
    FEScalarFieldPtr_Type scalarField1 ( new FEScalarField_Type ( scalarField1_FESpace ) );
    scalarField1->getVector() = 1. ;

    // Second scalar field
    FEScalarFieldPtr_Type scalarField2 ( new FEScalarField_Type ( scalarField2_FESpace ) );
    scalarField2->getVector() = 2.;

    // Vector field
    FEVectorFieldPtr_Type vectorField ( new FEVectorField_Type ( vectorField_FESpace ) );
    vectorField->getVector() = 3.;

    // Function
    MyFun function;

    // Scalar field for visualize the function
    FEScalarFieldPtr_Type scalarFieldFunction ( new FEScalarField_Type ( function_FESpace ) );

    // Stop chronoFEFieldAndFEFct
    chronoFEFieldAndFEFct.stop();

    // The leader process print chronoFEFieldAndFEFct
    displayer.leaderPrintMax ( "Time for create the fields and the function ",
                               chronoFEFieldAndFEFct.diff() );

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
    displayer.leaderPrintMax ( "Time for create the add the fields to the function ",
                               chronoAddFields.diff() );

    // Define the dummy element and coordinate for the evaluation of the function
    UInt iElem = 0;
    Vector3D point;

    // Compute the error
    Real error = 0;
    error = std::fabs ( ( std::sin (1.) + 4. ) / 3. - function.eval ( iElem, point, 0. ) );

    // Interpolate the value of the function
    function.interpolate ( *scalarFieldFunction );

    // Set the exporter for the fields
    exporterPtr_Type exporter;

    // Type of the exporter
    std::string const exporterType =  dataFile ( "exporter/type", "hdf5");

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare ("hdf5") == 0 )
    {
        exporter.reset ( new ExporterHDF5< regionMesh_Type > ( dataFile, dataFile ( "exporter/file_name", "FieldFunction" ) ) );

        // Set directory where to save the solution
        exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );

        exporter->setMeshProcId ( meshPtr, Members->comm->MyPID() );
    }
    else
#endif
    {
        if ( exporterType.compare ("none") == 0 )
        {
            exporter.reset ( new ExporterEmpty< regionMesh_Type > ( dataFile, dataFile ( "exporter/file_name", "FieldFunction" ) ) );

            // Set directory where to save the solution
            exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );

            exporter->setMeshProcId ( meshPtr, Members->comm->MyPID() );
        }
        else
        {
            exporter.reset ( new ExporterEnsight< regionMesh_Type > ( dataFile, dataFile ( "exporter/file_name", "FieldFunction" ) ) );

            // Set directory where to save the solution
            exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );

            exporter->setMeshProcId ( meshPtr, Members->comm->MyPID() );
        }
    }

    const std::string fieldName = dataFile ( "exporter/name_field", "Field" );

    // Add the first scalar field to the exporter
    exporter->addVariable ( ExporterData< regionMesh_Type >::ScalarField,
                            fieldName + "1",
                            scalarField1_FESpace,
                            scalarField1->getVectorPtr(),
                            static_cast<UInt> ( 0 ),
                            ExporterData< regionMesh_Type >::UnsteadyRegime,
                            ExporterData< regionMesh_Type >::Cell );

    // Add the second scalar field to the exporter
    exporter->addVariable ( ExporterData< regionMesh_Type >::ScalarField,
                            fieldName + "2",
                            scalarField2_FESpace,
                            scalarField2->getVectorPtr(),
                            static_cast<UInt> ( 0 ),
                            ExporterData< regionMesh_Type >::UnsteadyRegime );

    // Add the vector field to the exporter
    exporter->addVariable ( ExporterData< regionMesh_Type >::VectorField,
                            fieldName + "3",
                            vectorField_FESpace,
                            vectorField->getVectorPtr(),
                            static_cast<UInt> ( 0 ),
                            ExporterData< regionMesh_Type >::UnsteadyRegime,
                            ExporterData< regionMesh_Type >::Cell );

    const std::string fctName =  dataFile ( "exporter/name_fct", "Function" );

    // Add the scalar field representing the function to the exporter
    exporter->addVariable ( ExporterData< regionMesh_Type >::ScalarField,
                            fctName + "1",
                            function_FESpace,
                            scalarFieldFunction->getVectorPtr(),
                            static_cast<UInt> ( 0 ),
                            ExporterData< regionMesh_Type >::UnsteadyRegime,
                            ExporterData< regionMesh_Type >::Cell );

    // Export all the solutions
    exporter->postProcess (0);

    // Return the error
    return error;

} // run
