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

/*!
    @file
    @brief Generation muscular fibers and sheets

    @author Simone Rossi <simone.rossi@epfl.ch>
    @maintainer Simone Palamara <palamara.simone@gmail.com>
    @date 31-01-2014


    Generation of the muscular fibers and sheets on a generic
    geometry representing the left or right ventricle, generated
    according to geometrical rules based on anatomical knowledge.
    For more details about the method see [S.Rossi et al,European Journal of
    Mechanics A/Solids (2013), http://dx.doi.org/10.1016/j.euromechsol.2013.10.009]

 */

//#include <Epetra_ConfigDefs.h>

// ------------------------------------------------------------------------------
//  Include MPI for parallel simulations
// ------------------------------------------------------------------------------
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// ------------------------------------------------------------------------------
//  To set up the linear solver we need a Teuchos parameter list
// ------------------------------------------------------------------------------
#include <Teuchos_ParameterList.hpp>

// ------------------------------------------------------------------------------
//  Needed to generate the ouput folder
// ------------------------------------------------------------------------------
#include <sys/stat.h>

// ------------------------------------------------------------------------------
// BCInterface is the interface that between the datafile and the
// boundary conditions. We will create a ummy physical solver in order to use it.
// ------------------------------------------------------------------------------
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/bc_interface/core/solver/DefaultPhysicalSolver.hpp>

// ------------------------------------------------------------------------------
//  Usefule utility to load the mesh in one line.
// Not working with partitioned meshes
// ------------------------------------------------------------------------------
#include <lifev/core/mesh/MeshLoadingUtility.hpp>

// ------------------------------------------------------------------------------
//  The sheets are defined as the gradient of a scalar potential, we therefore
// need to be able to compute the gradient in the nodes.
// ------------------------------------------------------------------------------
#include <lifev/core/fem/GradientRecovery.hpp>
#include <lifev/core/fem/BCManage.hpp>

// ------------------------------------------------------------------------------
// To solve the laplacian, use the linear solver (AztecOO or Belos) with ML
// ------------------------------------------------------------------------------
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

// ------------------------------------------------------------------------------
// The laplacian is assembled using Expression Template Assembly
// ------------------------------------------------------------------------------
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

// ------------------------------------------------------------------------------
//  Cannot save without HDF5 !!! Please make sure you have it
// ------------------------------------------------------------------------------
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

// ------------------------------------------------------------------------------
//  In this file there are a bunch of useful functions.
// Here we use it only to normalize a vector.
// ------------------------------------------------------------------------------
#include <lifev/electrophysiology/util/HeartUtility.hpp>

// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. Moreover,
// to make the code more readable, we also make typedefs for the mesh type,
// matrix type, vector type, boundary condition
// ---------------------------------------------------------------

using namespace LifeV;

// ---------------------------------------------------------------
// We typedef some common type we will frequently
// ---------------------------------------------------------------

typedef RegionMesh<LinearTetra>                        		mesh_Type;
typedef boost::shared_ptr< mesh_Type >                 		meshPtr_Type;

typedef MatrixEpetra<Real>                             		matrix_Type;
typedef boost::shared_ptr< matrix_Type >               		matrixPtr_Type;

typedef VectorEpetra                                   		vector_Type;
typedef boost::shared_ptr< vector_Type >               		vectorPtr_Type;



typedef FESpace< mesh_Type, MapEpetra >					    fespace_Type;
typedef boost::shared_ptr<fespace_Type >				    fespacePtr_Type;


// ---------------------------------------------------------------
// In order to keep the code more readble I created a couple
// of auxiliary functions.
// In this test we only have the datafile that will be
// a GetPot object. To set up the linear solver we need a Teuchos::ParameterList
// that typically reads  xml files. Therefore,
// the createListFromGetPot function create the requested Teuchos list
// from the datafile we have.
// In the end we will export three vector field in three different
// files. To avoid code repetition, I created this function
// that exports the requested vecotr fields.
// ---------------------------------------------------------------
void createListFromGetPot(Teuchos::ParameterList& solverList, const GetPot& dataFile);
void exportVectorField(boost::shared_ptr<Epetra_Comm> comm,
		               meshPtr_Type mesh,
		               fespacePtr_Type fespace,
		               vectorPtr_Type vector,
		               std::string postDir,
		               std::string outputName,
		               std::string hdf5name);


// ------------------------------------------------------------------------------
// Dummy class for to be used with the BCInterfac3D.
// We could have used the ElectroETAMonodomainSolver as physical solver,
// but in this way, this test is totally independent.
// ------------------------------------------------------------------------------

Real fzero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}
// Starting ...
int main ( int argc, char** argv )
{

// ---------------------------------------------------------------
//  In parallel? Yes, we dare!!!
// ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    typedef BCHandler                                      		bc_Type;
    typedef boost::shared_ptr< bc_Type >                   		bcPtr_Type;

    typedef DefaultPhysicalSolver<VectorEpetra>  				physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >   	bcInterface_Type;

     typedef boost::shared_ptr< bcInterface_Type >          		bcInterfacePtr_Type;
     typedef MeshUtility::MeshTransformer<mesh_Type>        		meshTransformer_Type;

    //*************************************************************//
    // We create, as usual, the output folder where
    // we are going to save the solutions of
    // the simulation.
    // It requires to append "-o OutputFolderName" flag
    // when you laucnh the executable, e.g.
    // mpirun -n 2 Electrophysiology_thisTestName -o SolutionFolder
    // By default the name of the folder will be "Output"
    //*************************************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( Comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    //*************************************************************//
    // We create the datafile. The datafile is passed through the flag
    // "-f datafileName" when launching the executable.
    // By default the name of the datafiler is "data".
    // Usually you should not bother about this flag, as the
    // datafile has already the default filename.
    // Remember about this, though, if you are going to change the filename
    // or if you want to add another datafile and run your simulation
    // with that one.
    //*************************************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //*************************************************************//
    // We specified in the datafile the name and the path of the
    // mesh we are going to create the fibers on.
    //*************************************************************//

    std::string meshName = dataFile( "problem/space_discretization/mesh_file", "" );
    std::string meshPath = dataFile( "problem/space_discretization/mesh_dir", "./" );

    //*************************************************************//
    // Here we create a pointer to the mesh and then we load it.
    // Notice how cool is to load the mesh in just one line!!!
    // Note that the pointer will point to the partitioned mesh
    // If you would like to keep informations about the full mesh
    // create another meshPtr_Type and call for example
    // MeshUtility::loadMesh (meshPart, meshFull, meshName, meshPath);
    //*************************************************************//
    meshPtr_Type meshPart (new mesh_Type ( Comm ) );
    MeshUtility::loadMesh (meshPart, meshName, meshPath);

    //*************************************************************//
    // Here we define the finite element spaces. In particular
    // - uSpace is the finite element space for Expression Template Assembly
    // - uFESpace is the usual finite element space required for
    //   boundary conditions and for the exporter
    // - vectorESpace is the space for the vectorial fields (sheets and fibers)
    //   we will define. Even if we don't need to do finite element
    //   operation on them, this fespace will be required to export
    //   the solution.
    //*************************************************************//

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );

    fespacePtr_Type uFESpace( new FESpace< mesh_Type, MapEpetra > (meshPart, "P1", 1, Comm) );

    fespacePtr_Type vectorFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPart, "P1", 3, Comm) );

    //*************************************************************//
    // We asseble the stiffness matrix using expression template.
    // For more details, look at the ETA tutorial.
    //*************************************************************//

    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uSpace->map() ) );

    *systemMatrix *= 0.0;
    {
        using namespace ExpressionAssembly;


        integrate (  elements (uSpace->mesh() ),
                     quadRuleTetra4pt,
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
        )
        >> systemMatrix;
    }

    systemMatrix->globalAssemble();

    //*************************************************************//
    // Setting up the boundary conditions is always an issue.
    // Fortunately Cristiano Malossi implemented a way to read
    // the boundary condition from the datafile. This is achieved
    // by using the BCInterface class.
    // The BC is template on a physicalSolver_Type.
    // We have not created a specific BCInterface for the
    // monodomain model (as usually one imposes homogeneous Neumann
    // conditions). Therefore here I created a default physical solver
    // in the bc interface module which will be used to solve this
    // simple laplacian problem. Check the typedef.
    // Then to set up the boundary conditions, we need to
    // 1 -  create the BCHandler
    // 2 - fill the handler with the boundary conditions specified
    //     in the datafile
    // 3 - update the boundary condtions
    // For more information on how to use the BCInterface refer
    // to that module.
    //*************************************************************//
    bcInterfacePtr_Type                     BC ( new bcInterface_Type() );
    BC->createHandler();
    BC->fillHandler ( data_file_name, "problem" );
    BC->handler()->bcUpdate( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

    //**********************************************************//
    //  We are going to solve the laplace equation with an
    // iterative solver. We precondition the system with Ifpack.
    // The parameters of the preconditioner are in the datafile.
    // If you want to use ML, change the precType to
    // LifeV::PreconditionerML
    //**********************************************************//

    typedef LifeV::Preconditioner             basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack           prec_Type;
    typedef boost::shared_ptr<prec_Type>      precPtr_Type;

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    //**********************************************************//
    // As detailed above to setup the linear system we need a
    // Teuchos::ParmeterList list. Usually this is read from
    // an xml file. Since the xml is not available in this test
    // the createListFromGetPot reads the required parameters
    // from the datafile and put the in the Teuchos list.
    // You can see the actual implementation at the end of the
    // test.
    //**********************************************************//
	Teuchos::ParameterList solverList;
	createListFromGetPot(solverList, dataFile);

	LinearSolver linearSolver;
	linearSolver.setCommunicator (Comm);
    linearSolver.setParameters ( solverList );
	linearSolver.setPreconditioner ( precPtr );

    //*************************************************************//
	// We are going to solve the laplace equation. The right hand
	// is zero!
	//*************************************************************//
    vectorPtr_Type rhs (new vector_Type ( uSpace -> map() ) );
    *rhs *= 0.0;
    rhs -> globalAssemble();

    //*************************************************************//
    // We impose the boundary conditions, on our system.
    //*************************************************************//

    bcManage ( *systemMatrix, *rhs, *uSpace->mesh(), uSpace->dof(), *BC -> handler(), uFESpace->feBd(), 1.0, 0.0 );

    //*************************************************************//
    // We declare the vector where we want to put the solution.
    // Here we want to solve the linear system Ax=b;
    // we tell the linear solver which A to use (systemMatrix)
    // and which right hand side to use (rhs).
    // Then we solve telling the solver to put the solution in the
    // vector x (solution).
    //*************************************************************//
    vectorPtr_Type solution ( new vector_Type ( uFESpace -> map() ) );


    linearSolver.setOperator (systemMatrix);
    linearSolver.setRightHandSide (rhs);
    linearSolver.solve (solution);

    //*************************************************************//
    // We save the potential field just computed on a file
    // using HDF5.
    //*************************************************************//
    ExporterHDF5< mesh_Type > exporter;
    exporter.setMeshProcId ( meshPart, Comm -> MyPID() );
    exporter.setPostDir (problemFolder);
    exporter.setPrefix ("Potential");
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential", uFESpace, solution, UInt (0) );
    exporter.postProcess (0);
    exporter.closeFile();

    //*************************************************************//
    // The sheets are defined as the gradient of the scalar
    // potential just computed.  We create three vector where we
    // store the components of the sheets direction.
    // We computed the gradient using the superconvergent
    // gradient recovery patch of ZZ (check that file for more infos)
    //*************************************************************//
    vectorPtr_Type sx (new vector_Type ( uSpace -> map() ) );
    vectorPtr_Type sy (new vector_Type ( uSpace -> map() ) );
    vectorPtr_Type sz (new vector_Type ( uSpace -> map() ) );

    *sx = GradientRecovery::ZZGradient (uSpace, *solution, 0);
    *sy = GradientRecovery::ZZGradient (uSpace, *solution, 1);
    *sz = GradientRecovery::ZZGradient (uSpace, *solution, 2);

    //*************************************************************//
    // We declare the rule-based sheets, fibers, and the projection
    // vectors. All these vectors have three components. Therefore
    // we use the vectorFESpace to create them.
    //*************************************************************//

    vectorPtr_Type rbSheet ( new vector_Type ( vectorFESpace -> map() ) );
    vectorPtr_Type rbFiber ( new vector_Type ( vectorFESpace -> map() ) );
    vectorPtr_Type projection ( new vector_Type ( vectorFESpace -> map() ) );

    //*************************************************************//
    // From now on all the operations (except for the export)
    // will not require communication between the processors.
    // Using more processor will surely speed up the computations.
    // Therefore the first step is to understand what is the length
    // of each of the above vectors on each processor.
    // Since we used the same map to declare the vectors, they
    // are distributed on the processors in the same way.
    // We check only the lenght one of them calling the
    // MyLength() method (this method is defined in Trilinos).
    // The length n is an integer, and since we have 3 components
    // each component will hav n/3 entries.
    //*************************************************************//


    int n = (*rbSheet).epetraVector().MyLength();
    int d = n / 3;

    //*************************************************************//
    // We loop over the number of entries for each component and
    // we fill the rule-based sheet vector field.
    // To access the components of the vector we use the global IDs.
    // Therefore, we first compute the GIDs and the we assign to
    // rbSheet the components of the sheets sx, sy and sz.
    // After that we normalize this vector field.
    // In this way we have computed the sheet field.
    //*************************************************************//
    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

        (*rbSheet) [i] = (*sx) [i];
        (*rbSheet) [j] = (*sy) [i];
        (*rbSheet) [k] = (*sz) [i];
    }

    ElectrophysiologyUtility::normalize (*rbSheet);

    //*************************************************************//
    // The algorithm requires to give as input the centerline of
    // the left ventricle. This can be read from datafile
    // so that if you want to change the mesh/geoemtry you don't have
    // to recompile
    //*************************************************************//
    Real cx = dataFile("problem/centerline_x", 0.0);
    Real cy = dataFile("problem/centerline_y", 0.0);
    Real cz = dataFile("problem/centerline_z", 1.0);


    //*************************************************************//
    // What I call the projection field is the projection of the
    // centerline to the plane orthogonal to the sheets.
    // Given the vector of the centerline c and the sheets s.
    // then this projection vector p is given by
    // p = c - (c,s) s,
    // where (c,s) is the scalar product of c with s.
    // After we normalize the projection vector.
    //*************************************************************//

    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

        Real cdot = cx * (*rbSheet) [i] + cy * (*rbSheet) [j] + cz * (*rbSheet) [k];

        (*projection) [i] = cx - cdot *  (*rbSheet) [i];
        (*projection) [j] = cy - cdot *  (*rbSheet) [j];
        (*projection) [k] = cz - cdot *  (*rbSheet) [k];
    }

    ElectrophysiologyUtility::normalize (*projection);

    //*************************************************************//
    // In each point we have defined to orthogonal direction s and p
    // We finish computing the local coordinate system computing the
    // third orthonormal component. This is a prototype of the fiber
    // field. I usually call it flat fiber field in the sense that
    // is a vector field which has zero component in the direction
    // of the left ventricle centerline and it can be rotated in
    // order to define the seeked fiber field.
    // This flat fiber vector f is the vector product of s and p:
    // f = s x p;
    // Although f should already be normalized, since the vector
    // product of two unit vector is a unit vector, we call again
    // the normalize function.
    //*************************************************************//

    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

        (*rbFiber) [i] = (*rbSheet) [j]  * (*projection) [k] - (*rbSheet) [k]  * (*projection) [j];
        (*rbFiber) [j] = (*rbSheet) [k]  * (*projection) [i] - (*rbSheet) [i]  * (*projection) [k];
        (*rbFiber) [k] = (*rbSheet) [i]  * (*projection) [j] - (*rbSheet) [j]  * (*projection) [i];
    }

    ElectrophysiologyUtility::normalize (*rbFiber);

    //*************************************************************//
    // The last step is to rotate the flat fiber field just computed
    // We rotate the fibers with respect to the sheet field in order
    // to keep the orthogonal.
    // We make use of Rodrigues formula (you can check wikipedia :) )
    // The angle of rotation changes through the wall thickness.
    // The hypothesis is that the angle of rotation can be written
    // as a function of the scalar potential we computed at the
    // beginning.
    // We first read from the data file the angle of rotation we
    // want to impose on the endocardium and on the epicardium.
    //*************************************************************//
    Real epi_angle = dataFile("problem/epi_angle", -60.0);
    Real endo_angle = dataFile("problem/endo_angle", 60.0);

    for ( int l (0); l < d; l++)
    {
        int i = (*rbSheet).blockMap().GID (l);
        int j = (*rbSheet).blockMap().GID (l + d);
        int k = (*rbSheet).blockMap().GID (l + 2 * d);

        //*************************************************************//
        // We define a linear relation between the scalar potential u
        // (the vector solution) and the angle teta, such that
        //  teta = m * u + q;
        // We need to compute m and q!
        // First: transform the angles in radiants.
        // Second: define m as the difference betwee the epi and the endo
        // angles.
        // Third: define q as the endo angle.
        // Fourth: compute m * u + q
        //*************************************************************//
        Real p = 3.14159265358979;
        Real teta1 = p * epi_angle / 180;
        Real teta2 = p * endo_angle / 180;
        Real m = (teta1 - teta2 );
        Real q = teta2;
        Real teta;

        teta = m * (*solution) [i] + q;

        //*************************************************************//
        // For easiness of implementation and to be more readable
        // we define the local components of the sheets and of the fibers.
        // Note that the component of the fibers are the one of the
        // flat fiber field.
        //*************************************************************//

        Real s01 = (*rbSheet) [i];
        Real s02 = (*rbSheet) [j];
        Real s03 = (*rbSheet) [k];
        Real f01 = (*rbFiber) [i];
        Real f02 = (*rbFiber) [j];
        Real f03 = (*rbFiber) [k];

        //*************************************************************//
        // The fiber field F is a rotation of the flat fiber field f
        // F = R f
        // where R is the rotation matrix.
        // To compute R we need the sin(teta) and
        // the sin(teta)^2 and the cross-product matrix W (check
        // rodrigues formula on wikipedia :) )
        //*************************************************************//
        Real sa = std::sin (teta);
        Real sa2 = 2.0 * std::sin (0.5 * teta) * std::sin (0.5 * teta);

        Real W11 = 0.0;
        Real W12 = -s03;
        Real W13 = s02;
        Real W21 = s03;
        Real W22 = 0.0;
        Real W23 = -s01;
        Real W31 = -s02;
        Real W32 = s01;
        Real W33 = 0.0;
        //
        Real R11 = 1.0 + sa * W11 + sa2 * ( s01 * s01 - 1.0 );
        Real R12 = 0.0 + sa * W12 + sa2 * ( s01 * s02 );
        Real R13 = 0.0 + sa * W13 + sa2 * ( s01 * s03 );
        Real R21 = 0.0 + sa * W21 + sa2 * ( s02 * s01 );
        Real R22 = 1.0 + sa * W22 + sa2 * ( s02 * s02 - 1.0 );
        Real R23 = 0.0 + sa * W23 + sa2 * ( s02 * s03 );
        Real R31 = 0.0 + sa * W31 + sa2 * ( s03 * s01 );
        Real R32 = 0.0 + sa * W32 + sa2 * ( s03 * s02 );
        Real R33 = 1.0 + sa * W33 + sa2  * ( s03 * s03 - 1.0 );

        //*************************************************************//
        // Finally we perform the rotation of the flat fiber field
        //*************************************************************//

        (*rbFiber) [i] = R11 * f01 + R12 * f02 + R13 * f03;
        (*rbFiber) [j] = R21 * f01 + R22 * f02 + R23 * f03;
        (*rbFiber) [k] = R31 * f01 + R32 * f02 + R33 * f03;
    }


    //*************************************************************//
    // Now that we have computed all the desired vector fields
    // we export them using HDF5 format.
    // From datafile we can choose the name of the files of the
    // fiber and sheet vectors.
    // Also the field are saved on the h5 file with a identification
    // name (that is the one you can select on paraview). I called
    // this fiberHDF5Name and sheetHDF5Name. This name is important
    // to know if we want to import the computed fields in other
    // simulations
    //*************************************************************//
    std::string outputFiberFileName = dataFile("problem/output_fiber_filename", "FiberDirection");
    std::string fiberHDF5Name = dataFile("problem/hdf5_fiber_name", "fibers");

    std::string outputSheetsFileName = dataFile ("problem/output_sheets_filename", "SheetsDirection");
    std::string sheetsHDF5Name = dataFile("problem/hdf5_sheets_name", "sheets");

    exportVectorField(Comm, meshPart, vectorFESpace, rbFiber, problemFolder, outputFiberFileName, fiberHDF5Name );
    exportVectorField(Comm, meshPart, vectorFESpace, rbSheet, problemFolder, outputSheetsFileName, sheetsHDF5Name );
    exportVectorField(Comm, meshPart, vectorFESpace, projection, problemFolder, "Projection", "projection" );


    //*************************************************************//
    // We test if the norm of the fuber field and the norm of the
    // sheet field are the same
    //*************************************************************//

    Real normS = rbSheet-> norm2();
    Real normF = rbFiber-> norm2();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    Real err = std::abs (normF - normS) / std::abs(normS);
	std::cout << std::setprecision(20) << "\nError: " <<  err << "\nFiber Norm: " <<  normF << "\n";
	std::cout << std::setprecision(20) << "Sheet Norm: " <<  normS << "\n";
    if ( err > 1e-13 )
    {
        return EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        return EXIT_SUCCESS;
    }
}

//Implementation of the predeclared functions:
//
//Read the parameters from a datafile and
// put them in a Teuchos:ParameterList
void createListFromGetPot(Teuchos::ParameterList& solverList, const GetPot& dataFile)
{
	std::string solverName   = dataFile( "problem/solver/solver_name", "AztecOO");
	std::string solver       = dataFile( "problem/solver/solver", "gmres");
	std::string conv         = dataFile( "problem/solver/conv", "rhs");
	std::string scaling      = dataFile( "problem/solver/scaling", "none");
	std::string output       = dataFile( "problem/solver/output", "all");
	Int maxIter              = dataFile( "problem/solver/max_iter", 200);
	Int maxIterForReuse      = dataFile( "problem/solver/max_iter_reuse", 250);
	Int kspace               = dataFile( "problem/solver/kspace", 100);
	Int orthog               = dataFile( "problem/solver/orthog", 0);
	Int auxvec               = dataFile( "problem/solver/aux_vec", 0);
	double tol               = dataFile( "problem/solver/tol", 1e-10);
	bool reusePreconditioner = dataFile( "problem/solver/reuse", true);
	bool quitOnFailure       = dataFile( "problem/solver/quit", false);
	bool silent              = dataFile( "problem/solver/silent", false);

	solverList.set("Solver Type", solverName);
	solverList.set("Maximum Iterations", maxIter);
	solverList.set("Max Iterations For Reuse", maxIterForReuse);
	solverList.set("Reuse Preconditioner", reusePreconditioner);
	solverList.set("Quit On Failure", quitOnFailure);
	solverList.set("Silent", silent);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("solver", solver);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("conv", conv);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("scaling", scaling);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("output", output);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("tol", tol);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("max_iter", maxIter);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("kspace", kspace);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("orthog", orthog);
	solverList.sublist("Solver: Operator List").sublist("Trilinos: AztecOO List").set("aux_vec", auxvec);;
}

//Export vector to file using HDF5 exporter
void exportVectorField(boost::shared_ptr<Epetra_Comm> comm,
		               meshPtr_Type mesh,
		               fespacePtr_Type fespace,
		               vectorPtr_Type vector,
		               std::string postDir,
		               std::string outputName,
		               std::string hdf5name)
{
    ExporterHDF5< mesh_Type > exporter;
    exporter.setMeshProcId ( mesh, comm -> MyPID() );
    exporter.setPostDir (postDir);
    exporter.setPrefix (outputName);
    exporter.addVariable ( ExporterData<mesh_Type>::VectorField,  hdf5name, fespace, vector, UInt (0) );
    exporter.postProcess (0);
    exporter.closeFile();
}


