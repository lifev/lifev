/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
             Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2011-03-09

  Copyright (C) 2010 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/*!
    @file ethiersteiman.hpp
    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 2011-03-08
 */

#ifndef NAVIERSTOKES_H
#define NAVIERSTOKES_H 1

#include <lifev/navier_stokes/solver/OseenSolver.hpp>
#include <lifev/core/mesh/ElementShapes.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

/*! @enum TimeScheme
    Order of the BDF
 */
enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
    @class NavierStokes
    @brief NavierStokes Simulation class

    @author Christophe Prud'homme
    @author Gwenol Grandperrin
    @see C.R. Ethier and D.A. Steinman. Exact fully 3D Navier-Stokes solutions for benchmarking. Int. J. Numer. Methods Fluids, 19(5):369-375, 1994

    This class tests many settings, we propose here a list of the available options.
    <ul>
        <li> NavierStokes/test (none/accuracy/space_convergence)
        <li> NavierStokes/initialization (interpolation/projection)
        <li> NavierStokes/export_norms
        <li> NavierStokes/export_exact_solutions
        <li> NavierStokes/mesh_source (regular_mesh/file)
        <li> exporter/type (ensight/hdf5)
        <li> fluid/problem/Re
        <li> fluid/physics/viscosity
        <li> fluid/physics/density
    </ul>

    In the case of a space convergence test the following options are available
    <ul>
        <li> NavierStokes/space_convergence_tolerance
        <li> fluid/space_discretization/mesh_number
        <li> fluid/space_discretization/mesh_size
        <li> fluid/space_discretization/FE_number
        <li> fluid/space_discretization/vel_order
        <li> fluid/space_discretization/vel_conv_order_order
        <li> fluid/space_discretization/press_order
        <li> fluid/space_discretization/press_conv_order
    </ul>

    In the case of an accuracy test the following options are available
    <ul>
        <li> NavierStokes/accuracy_tolerance
    </ul>

    If the mesh source is a file, the file must be specified in
    <ul>
        <li> fluid/space_discretization
    </ul>
 */

template<typename MeshType, typename Problem>
class NavierStokes
{
public:
    typedef MeshType                                      mesh_Type;
    typedef Problem                                       problem_Type;
    typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;
    typedef LifeV::OseenSolver< mesh_Type >               fluid_Type;
    typedef typename fluid_Type::vector_Type              vector_Type;
    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;
    typedef typename fluid_Type::matrix_Type              matrix_Type;

    /** @name Constructors, destructor
     */
    //@{

    //! Constructor
    /*!
        @param argc number of parameter passed through the command line
        @param argv char passed through the command line
     */
    NavierStokes ( int argc,
                   char** argv,
                   const std::string defaultDataName = "data",
                   const std::string outputName = "navierStokes");

    //! Destructor
    ~NavierStokes()
    {}

    //@}

    /** @name  Methods
     */
    //@{

    //! Launches the simulation
    void run();

    //@}


private:
    /*! @enum TestType
        Order of the BDF
     */
    enum TestType {None, Accuracy, SpaceConvergence};

    /*! @enum InitializationType
        Type of initialization. "Interpolation" just interpolates the value of the exact solution to the DoFs.
        "Projection" solves an Oseen problem where alpha=0, the convective term is linearized by using the exact solution for beta,
        and the time derivative is passed to the right hand side and computed from the exact solution.
     */
    enum InitializationType {Projection, Interpolation};

    /*! @enum MeshSourceType
        Type of mesh source. It can be a "File" or a "RegularMesh" which is generated during the simulation
     */
    enum MeshSourceType {File, RegularMesh};

    struct RESULT_CHANGED_EXCEPTION
    {
    public:
        RESULT_CHANGED_EXCEPTION()
        {
            std::cout << "Some modifications led to changes in the l2 norm of the solution" << std::endl;
        }
    };

    //! Computes the L2 error and the relative L2 error with the exact solution
    /*!
        @param velocityAndPressureSolution number of parameter passed through the command line
        @param uL2Error Variable to store the L2 error for the veloctity
        @param uRelError Variable to store the Relative L2 error for the velocity
        @param uFESpace Variable to store the FE space for the velocity
        @param pL2Error Variable to store the L2 error the pressure
        @param pRelError Variable to store the Relative L2 error for the pressure
        @param pFESpace Variable to store the FE space for the pressure
        @param time Actual timestep of the simulation
     */
    void computeErrors (const vector_Type& velocityAndPressureSolution,
                        LifeV::Real& uL2Error, LifeV::Real& uRelError, feSpacePtr_Type& uFESpace,
                        LifeV::Real& pL2Error, LifeV::Real& pRelError, feSpacePtr_Type& pFESpace,
                        LifeV::Real time);

    //! Method to check the convergence rate of the solution
    /*!
        For each type of finite elements the method uses the computed errors obtained with the
        different meshes to check if the order of convergence follows the theory predictions
        @param uFELabels Vector containing the FE names for the velocity (e.g. P2, P1Bubble, P1)
        @param uL2Error Vector containing the computed errors for the velocity
        @param uConvergenceOrder Vector containing the convergence order corresponding to uFELabel
        @param pFELabels Vector containing the FE names for the pressure (e.g. P2, P1Bubble, P1)
        @param pL2Error Vector containing the computed errors for the pressure
        @param pConvergenceOrder Vector containing the convergence order corresponding to pFELabel
        @param meshDiscretization Vector containing the subdivisions values used to generate the meshes
        @param convTolerance Tolerance for the test. The test is passed if (observed convergence)>convTolerance*(theory error prediction)
     */
    bool checkConvergenceRate (const std::vector<std::string>& uFELabels,
                               const std::vector<std::vector<LifeV::Real> >& uL2Error,
                               const std::vector<LifeV::UInt>& uConvergenceOrder,
                               const std::vector<std::string>& pFELabels,
                               const std::vector<std::vector<LifeV::Real> > pL2Error,
                               const std::vector<LifeV::UInt>& pConvergenceOrder,
                               const std::vector<LifeV::UInt>& meshDiscretizations,
                               LifeV::Real convTolerance);


    struct Private;
    boost::shared_ptr<Private> M_data;

    std::vector<LifeV::UInt>   M_meshDiscretization;
    std::vector<std::string>   M_uFELabels;
    std::vector<std::string>   M_pFELabels;
    std::vector<LifeV::UInt>   M_uConvergenceOrder;
    std::vector<LifeV::UInt>   M_pConvergenceOrder;

    // Test to be performed (accuracy or convergence in space)
    TestType                   M_test;
    LifeV::Real                M_convTol; // Tolerance of the test (should be <1)
    // Actually for convTol=1, the test failed
    // if the improvement of accuracy is less
    // than predicted by the theory.
    // convTol lower down the theoretical bounds
    LifeV::Real                M_accuracyTol;

    // Data related to norm export
    bool                       M_exportNorms;
    std::ofstream              M_outNorm;

    // Data related to solution export
    bool                       M_exportExactSolutions;

    // Initialization method
    InitializationType         M_initMethod;

    // Mesh source
    MeshSourceType             M_meshSource;

    // output file name
    std::string                M_outputName;
};


using namespace LifeV;

//! Null function used to define the boundary conditions
/*!
    @param t Current time
    @param x x-position
    @param y y-position
    @param z z-position
    @param i ith component of the returned value of the function
 */
Real zero_scalar ( const Real& /* t */,
                   const Real& /* x */,
                   const Real& /* y */,
                   const Real& /* z */,
                   const ID& /* i */ )
{
    return 0.;
}


//! Parse a string using "," to separate values and return a set of value
/*!
    @param list String containing the list of desired value separated by ","
 */
std::set<UInt> parseList ( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    std::string::size_type commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find ( "," );
        setList.insert ( atoi ( stringList.substr ( 0, commaPos ).c_str() ) );
        stringList = stringList.substr ( commaPos + 1 );
    }
    setList.insert ( atoi ( stringList.c_str() ) );
    return setList;
}


template<typename MeshType, typename Problem>
struct NavierStokes<MeshType, Problem>::Private
{
    Private() :
        nu    (1),
        steady (0)
    {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_Type;

    double         Re;

    std::string    data_file_name;

    double         nu;  /* < viscosity (in m^2/s) */
    //const double rho; /* < density is constant (in kg/m^3) */

    bool                             steady;
    boost::shared_ptr<Epetra_Comm>   comm;
};

template<typename MeshType, typename Problem>
NavierStokes<MeshType, Problem>::NavierStokes ( int argc,
                                                char** argv,
                                                const std::string defaultDataName,
                                                const std::string outputName)
    :
    M_data ( new Private ),
    M_outputName (outputName)
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow (defaultDataName.c_str(), 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    M_data->data_file_name = data_file_name;

    M_data->Re = dataFile ( "fluid/problem/Re", 1. );
    M_data->nu = dataFile ( "fluid/physics/viscosity", 1. ) /
                 dataFile ( "fluid/physics/density", 1. );

    // Test type
    string testType = dataFile ("NavierStokes/test", "none");
    if (testType == "none")
    {
        M_test = None;
    }
    else if (testType == "accuracy")
    {
        M_test = Accuracy;
    }
    else if (testType == "space_convergence")
    {
        M_test = SpaceConvergence;
    }
    else
    {
        std::cout << "[Error] Unknown test method" << std::endl;
        exit (1);
    }

    M_convTol     = dataFile ("NavierStokes/space_convergence_tolerance", 1.0);
    M_accuracyTol = dataFile ("NavierStokes/accuracy_tolerance", 1.0);

    // Method of initialization
    string initType = dataFile ("NavierStokes/initialization", "projection");
    if (initType == "projection")
    {
        M_initMethod = Projection;
    }
    else if (initType == "interpolation")
    {
        M_initMethod = Interpolation;
    }
    else
    {
        std::cout << "[Error] Unknown initialization method" << std::endl;
        exit (1);
    }

    M_exportNorms = dataFile ("NavierStokes/export_norms", false);
    M_exportExactSolutions = dataFile ("NavierStokes/export_exact_solutions", false);

    std::string meshSource =  dataFile ( "NavierStokes/mesh_source", "regular_mesh");
    if (meshSource == "regular_mesh")
    {
        M_meshSource = RegularMesh;
    }
    else if (meshSource == "file")
    {
        M_meshSource = File;
    }
    else
    {
        std::cout << "[Error] Unknown mesh source" << std::endl;
        exit (1);
    }

    //Checking the consistency of the data
    if (M_meshSource == File && M_test == SpaceConvergence)
    {
        std::cout << "[Error] You cannot use mesh files to test the space convergence." << std::endl;
        exit (1);
    }


#ifdef EPETRA_MPI

    //    MPI_Init(&argc,&argv);

    M_data->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size (MPI_COMM_WORLD, &ntasks);
#else
    M_data->comm.reset ( new Epetra_SerialComm() );
#endif

}

template<typename MeshType, typename Problem>
void
NavierStokes<MeshType, Problem>::computeErrors (const vector_Type& velocityAndPressureSolution,
                                                LifeV::Real& uL2Error, LifeV::Real& uRelError, feSpacePtr_Type& uFESpace,
                                                LifeV::Real& pL2Error, LifeV::Real& pRelError, feSpacePtr_Type& pFESpace,
                                                LifeV::Real time)
{
    // Computation of the error
    vector_Type vel  (uFESpace->map(), Repeated);
    vector_Type press (pFESpace->map(), Repeated);
    vector_Type velpressure ( velocityAndPressureSolution, Repeated );

    velpressure = velocityAndPressureSolution;
    vel.subset (velpressure);
    press.subset (velpressure, uFESpace->dim() *uFESpace->fieldDim() );

    uL2Error = uFESpace->l2Error (problem_Type::uexact, vel  , time, &uRelError );
    pL2Error = pFESpace->l20Error (problem_Type::pexact, press, time, &pRelError );
}

template<typename MeshType, typename Problem>
bool
NavierStokes<MeshType, Problem>::checkConvergenceRate (const std::vector<std::string>& uFELabels,
                                                       const std::vector<std::vector<LifeV::Real> >& uL2Error,
                                                       const std::vector<UInt>& uConvergenceOrder,
                                                       const std::vector<std::string>& pFELabels,
                                                       const std::vector<std::vector<LifeV::Real> > pL2Error,
                                                       const std::vector<UInt>& pConvergenceOrder,
                                                       const std::vector<UInt>& meshDiscretizations,
                                                       LifeV::Real convTolerance)
{
    // We want to check the convergence of the error and
    // see if it matches the theory.
    std::cout << "Checking the convergence:" << std::endl;

    // Test variable
    bool success (true); // Variable to keep trace of a previous error
    Real h1 (0.0), h2 (0.0); // Space discretization step
    Real uBound (0.0), pBound (0.0); // Velocity and pressure bounds
    Real uErrRatio (0.0), pErrRatio (0.0); // Ratio of the error E1/E2
    std::string status (""); // Information string

    UInt numFELabels (uFELabels.size() );
    UInt numDiscretizations (meshDiscretizations.size() );

    for (UInt iFELabel (0); iFELabel < numFELabels; ++iFELabel)
    {
        std::cout << "    - " << uFELabels[iFELabel] << "-" << pFELabels[iFELabel] << " ... " << std::endl;

        // Everything is OK a priori
        status = "OK";

        for (UInt jDiscretization (0); jDiscretization < numDiscretizations - 1; ++jDiscretization)
        {
            h1 = 1.0 / meshDiscretizations[jDiscretization];
            h2 = 1.0 / meshDiscretizations[jDiscretization + 1];

            uBound = convTolerance * pow (h1 / h2, int (uConvergenceOrder[iFELabel]) );
            pBound = convTolerance * pow (h1 / h2, int (pConvergenceOrder[iFELabel]) );

            uErrRatio = uL2Error[iFELabel][jDiscretization] / uL2Error[iFELabel][jDiscretization + 1]; // E1/E2
            pErrRatio = pL2Error[iFELabel][jDiscretization] / pL2Error[iFELabel][jDiscretization + 1];

            if (uErrRatio < uBound)
            {
                status = "FAILED";
                success = false;
            }
            if (pErrRatio < pBound)
            {
                status = "FAILED";
                success = false;
            }
            std::cout << "      " << " (velocity: " << uErrRatio << ">=?" << uBound
                      << ", pressure: " << pErrRatio << ">=?" << pBound << std::endl;
        }
        std::cout << "      Status: " << status << std::endl;

    }

    return success;
}

template<typename MeshType, typename Problem>
void
NavierStokes<MeshType, Problem>::run()
{
    bool verbose = (M_data->comm->MyPID() == 0);
    int nproc;
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    if (verbose)
    {
        std::cout << "[[BEGIN_SIMULATION]]" << std::endl << std::endl;

        std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
        std::cout << "Using MPI (" << nproc << " proc.)" << std::endl;
#else
        std::cout << "Using serial version" << std::endl;
#endif
    }

    // +-----------------------------------------------+
    // |             Begining of the test              |
    // +-----------------------------------------------+
    LifeChrono globalChrono;
    LifeChrono runChrono;
    LifeChrono initChrono;
    LifeChrono iterChrono;

    globalChrono.start();
    initChrono.start();

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Loading the data]" << std::endl;
    }
    GetPot dataFile ( M_data->data_file_name.c_str() );
    if (verbose)
    {
        switch (M_test)
        {
            case Accuracy:
                std::cout << "Test : checks the accuracy of the solution" << std::endl;
                break;
            case SpaceConvergence:
                std::cout << "Test : checks the convergence in space of the solution" << std::endl;
                break;
            case None:
                break;
        }
    }
    problem_Type::setParamsFromGetPot ( dataFile );

    UInt numDiscretizations;

    // Loading the discretization to be tested
    numDiscretizations = dataFile ( "fluid/space_discretization/mesh_number", 1 );
    for ( UInt i ( 0 ); i < numDiscretizations; ++i )
    {
        M_meshDiscretization.push_back (dataFile ( "fluid/space_discretization/mesh_size", 8, i ) );
    }

    UInt numFELabels = dataFile ( "fluid/space_discretization/FE_number", 1 );
    if (M_test == SpaceConvergence)
    {
        // Loading the convergence rate for the finite elements tested
        for ( UInt i ( 0 ); i < numFELabels; ++i )
        {
            M_uConvergenceOrder.push_back (dataFile ( "fluid/space_discretization/vel_conv_order", 2, i ) );
        }
        for ( UInt i ( 0 ); i < numFELabels; ++i )
        {
            M_pConvergenceOrder.push_back (dataFile ( "fluid/space_discretization/press_conv_order", 2, i ) );
        }
    }

    // Loading the Finite element to be tested
    for ( UInt i ( 0 ); i < numFELabels; ++i )
    {
        M_uFELabels.push_back (dataFile ( "fluid/space_discretization/vel_order", "P1", i ) );
    }
    for ( UInt i ( 0 ); i < numFELabels; ++i )
    {
        M_pFELabels.push_back (dataFile ( "fluid/space_discretization/press_order", "P1", i ) );
    }

    // Initialization of the errors array
    std::vector<std::vector<LifeV::Real> > uL2Error;
    std::vector<std::vector<LifeV::Real> > pL2Error;
    uL2Error.clear();
    pL2Error.clear();
    std::vector<LifeV::Real> tmpVec (numDiscretizations, 0.0);
    for (UInt iFELabel (0); iFELabel < numFELabels; ++iFELabel)
    {
        uL2Error.push_back (tmpVec);
        pL2Error.push_back (tmpVec);
    }

    initChrono.stop();
    if (verbose)
    {
        std::cout << "Initialization time (pre-run): " << initChrono.diff() << " s." << std::endl;
    }

    // Loop on the mesh refinement
    for (UInt jDiscretization (0); jDiscretization < numDiscretizations; ++jDiscretization)
    {
        UInt mElem = M_meshDiscretization[jDiscretization];

        // Loop on the finite element
        for (UInt iFELabel (0); iFELabel < numFELabels; ++iFELabel)
        {
            if (verbose)
            {
                std::cout << std::endl << "[[BEGIN_RUN_" << jDiscretization* numFELabels + iFELabel << "]]" << std::endl;
            }
            runChrono.reset();
            runChrono.start();
            initChrono.reset();
            initChrono.start();

            if (verbose && M_exportNorms)
            {
                std::string fileName ("norm_");
                std::ostringstream oss;
                oss << mElem;
                fileName.append (oss.str() );
                fileName.append ("_");
                fileName.append (M_uFELabels[iFELabel]);
                fileName.append (M_pFELabels[iFELabel]);
                fileName.append (".txt");
                M_outNorm.open (fileName.c_str() );
                M_outNorm << "% time / u L2 error / u H1 error / p L2 error \n" << std::flush;
            }

            // +-----------------------------------------------+
            // |               Loading the mesh                |
            // +-----------------------------------------------+
            if (verbose)
            {
                std::cout << "[Loading the mesh]" << std::endl;
            }

            boost::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( M_data->comm ) );

            Int geoDimensions = mesh_Type::S_geoDimensions;
            // Building the mesh from the source
            if(M_meshSource == RegularMesh) 
            {
                    regularMesh3D( *fullMeshPtr, 2, mElem, mElem, mElem, false, 1.0, 1.0, 1.0, -1.0,  -1.0,  -1.0);

                    if (verbose) std::cout << "Mesh source: regular mesh("
                                           << mElem << "x" << mElem << "x" << mElem << ")" << std::endl;
                
                if (verbose){
                        std::cout << "\n[Mesh size max : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).maxH << " ]\n";
                }
            }
            else if(M_meshSource == File)
            {
                MeshData meshData;
                meshData.setup (dataFile, "fluid/space_discretization");
                readMesh (*fullMeshPtr, meshData);

                if (verbose) std::cout << "Mesh source: file("
                                           << meshData.meshDir() << meshData.meshFile() << ")" << std::endl;
            }
            else
            {
                if (verbose)
                {
                    std::cout << std::endl << "Error: Unknown source type for the mesh" << std::endl;
                }
                exit (1);
            }
            
            if (verbose)
                std::cout << "Partitioning the mesh ... " << std::flush;
            
            
            boost::shared_ptr<mesh_Type > localMeshPtr;
        
            MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_data->comm);
            localMeshPtr = meshPart.meshPartition();
            
            
            fullMeshPtr.reset(); //Freeing the global mesh to save memory

            // +-----------------------------------------------+
            // |            Creating the FE spaces             |
            // +-----------------------------------------------+
            if (verbose)
            {
                std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
            }
            std::string uOrder =  M_uFELabels[iFELabel];
            std::string pOrder =  M_pFELabels[iFELabel];

            if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                                       << "FE for the pressure: " << pOrder << std::endl;

            if (verbose)
            {
                std::cout << "Building the velocity FE space ... " << std::flush;
            }
            feSpacePtr_Type uFESpace;
            uFESpace.reset (new feSpace_Type (localMeshPtr, uOrder, geoDimensions, M_data->comm) );
            if (verbose)
            {
                std::cout << "ok." << std::endl;
            }

            if (verbose)
            {
                std::cout << "Building the pressure FE space ... " << std::flush;
            }
            feSpacePtr_Type pFESpace;
            pFESpace.reset (new feSpace_Type (localMeshPtr, pOrder, 1, M_data->comm) );
            if (verbose)
            {
                std::cout << "ok." << std::endl;
            }

            UInt totalVelDof   = uFESpace->dof().numTotalDof();
            UInt totalPressDof = pFESpace->dof().numTotalDof();

            // Pressure offset in the vector
            UInt pressureOffset =  uFESpace->fieldDim() * uFESpace->dof().numTotalDof();

            if (verbose)
            {
                std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
            }
            if (verbose)
            {
                std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;
            }

            // +-----------------------------------------------+
            // |             Boundary conditions               |
            // +-----------------------------------------------+
            if (verbose)
            {
                std::cout << std::endl << "[Boundary conditions]" << std::endl;
            }
            std::string dirichletList = dataFile ( "fluid/problem/dirichletList", "" );
            std::set<UInt> dirichletMarkers = parseList ( dirichletList );
            std::string neumannList = dataFile ( "fluid/problem/neumannList", "" );
            std::set<UInt> neumannMarkers = parseList ( neumannList );

            BCHandler bcH;
            BCFunctionBase uWall ( problem_Type::uexact );
            BCFunctionBase uNeumann ( problem_Type::fNeumann );

            for (std::set<UInt>::const_iterator it = dirichletMarkers.begin();
                    it != dirichletMarkers.end(); ++it)
            {
                bcH.addBC ( "Wall", *it, Essential, Full, uWall, geoDimensions );
            }
            for (std::set<UInt>::const_iterator it = neumannMarkers.begin();
                    it != neumannMarkers.end(); ++it)
            {
                bcH.addBC ( "Flux", *it, Natural, Full, uNeumann, geoDimensions );
            }

            // If we change the FE we have to update the BCHandler (internal data)
            bcH.bcUpdate ( *localMeshPtr, uFESpace->feBd(), uFESpace->dof() );

            // +-----------------------------------------------+
            // |             Creating the problem              |
            // +-----------------------------------------------+
            if (verbose)
            {
                std::cout << std::endl << "[Creating the problem]" << std::endl;
            }
            boost::shared_ptr<OseenData> oseenData (new OseenData() );
            oseenData->setup ( dataFile );

            if (verbose)
            {
                std::cout << "Time discretization order " << oseenData->dataTimeAdvance()->orderBDF() << std::endl;
            }

            OseenSolver< mesh_Type > fluid (oseenData,
                                            *uFESpace,
                                            *pFESpace,
                                            M_data->comm);

            MapEpetra fullMap (fluid.getMap() );

            fluid.setUp (dataFile);

            fluid.buildSystem();

            MPI_Barrier (MPI_COMM_WORLD);

            // +-----------------------------------------------+
            // |       Initialization of the simulation        |
            // +-----------------------------------------------+
            if (verbose)
            {
                std::cout << std::endl << "[Initialization of the simulation]" << std::endl;
            }
            Real dt     = oseenData->dataTime()->timeStep();
            Real t0     = oseenData->dataTime()->initialTime();
            Real tFinal = oseenData->dataTime()->endTime();


            // bdf object to store the previous solutions
            TimeAdvanceBDFNavierStokes<vector_Type> bdf;
            bdf.setup (oseenData->dataTimeAdvance()->orderBDF() );

            /*
                Initialization with exact solution: either interpolation or "L2-NS"-projection
                Depending on which order scheme we want for the time derivative, we need a to
                setup a fixed number of solution required by the scheme. Therefore we need to
                compute the solution for some timestep before t0.
             */
            t0 -= dt * bdf.bdfVelocity().order();

            if (verbose)
            {
                std::cout << "Computing the initial solution ... " << std::endl;
            }

            vector_Type beta ( fullMap );
            vector_Type rhs ( fullMap );

            MPI_Barrier (MPI_COMM_WORLD);

            oseenData->dataTime()->setTime (t0);
            fluid.initialize ( problem_Type::uexact, problem_Type::pexact );

            bdf.bdfVelocity().setInitialCondition ( *fluid.solution() );

            /*
                Initial solution loading (interpolation or projection)
                This part of the code take advantage of the fact that the Projection
                method adds only a few lines to the Interpolation initialization method.
                First notice that the loop ensure that enough solutions are computed in
                order to use the BDF scheme chose previously.

                Interpolation case:
                We start by setting the current time then we initialize the OseenSolver
                using fluid.initialize. Therefore the exact solution is interpolated to
                obtain a solution. Then we store this solution for the BDF scheme using
                bdf.bdfVelocity().shiftRight(...).

                Projection case:
                In that case the solution obtained in fluid.initialize is used as a starting
                point to solve the steady-state problem:
                \Delta u + u^*\nabla u +\nabla p = -(\frace{\partial u}{\partial t})^*
                where the * means that the value is obtained by interpolating the quantity
                using the exact solution.
             */
            
            Real time = t0 + dt;
            for (  ; time <=  oseenData->dataTime()->initialTime() + dt / 2.; time += dt)
            {

                oseenData->dataTime()->setTime (time);

                beta *= 0.;
                rhs  *= 0.;

                fluid.initialize ( problem_Type::uexact, problem_Type::pexact );

                beta = *fluid.solution();

                if (M_initMethod == Projection)
                {
                    uFESpace->interpolate ( static_cast<typename feSpace_Type::function_Type> ( problem_Type::uderexact ), rhs, time);
                    rhs *= -1.;
                    rhs = fluid.matrixMass() * rhs;
                    fluid.setVelocityRhs(bdf.bdfVelocity().rhsContributionFirstDerivative());
                    fluid.updateSystem ( 0., beta, rhs );
                    fluid.iterate (bcH);
                }

                // Computation of the error
                LifeV::Real urelerr, prelerr, ul2error, pl2error;

                Real uh1error = 0.0;
                fluid.h1normVelocity(uh1error);

                computeErrors (*fluid.solution(),
                               ul2error, urelerr, uFESpace,
                               pl2error, prelerr, pFESpace,
                               time);

                if (verbose)
                {
                	std::cout << "\n[ERRORS]:\n";
                	std::cout << "Time: " << time  << ", "
                			<< "L2 velocity error: " << ul2error << ", "
                			<< "H1 velocity error: " << uh1error << ", "
                			<< "L2 pressure error: " << pl2error << "\n" << std::flush;
                }

                if (verbose && M_exportNorms)
                {
                    M_outNorm << time  << " "
                              << ul2error << " "
                              << uh1error << " "
                              << pl2error << "\n" << std::flush;
                }

                // Updating bdf
                bdf.bdfVelocity().shiftRight ( *fluid.solution() );

            }

            fluid.resetPreconditioner();

            boost::shared_ptr< Exporter<mesh_Type > > exporter;

            vectorPtr_Type velAndPressure;
            // only for export -->
            vectorPtr_Type exactPressPtr;
            vector_Type exactPress (pFESpace->map(), Repeated);
            vectorPtr_Type exactVelPtr;
            vector_Type exactVel (uFESpace->map(), Repeated);
            // <--

            std::string const exporterType =  dataFile ( "exporter/type", "ensight");

#ifdef HAVE_HDF5
            if (exporterType.compare ("hdf5") == 0)
            {
                exporter.reset ( new ExporterHDF5<mesh_Type > ( dataFile, M_outputName ) );
                exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
                exporter->setMeshProcId ( localMeshPtr, M_data->comm->MyPID() );
            }
#endif
            else if(exporterType.compare ("vtk") == 0)
            {
                exporter.reset ( new ExporterVTK<mesh_Type > ( dataFile, M_outputName ) );
                exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
                exporter->setMeshProcId ( localMeshPtr, M_data->comm->MyPID() );
            }
            else
            {
                if (exporterType.compare ("none") == 0)
                {
                    exporter.reset ( new ExporterEmpty<mesh_Type > ( dataFile, localMeshPtr, M_outputName, M_data->comm->MyPID() ) );
                }
                else
                {
                    exporter.reset ( new ExporterEnsight<mesh_Type > ( dataFile, localMeshPtr, M_outputName, M_data->comm->MyPID() ) );
                }
            }

            velAndPressure.reset ( new vector_Type (*fluid.solution(), exporter->mapType() ) );
            if (M_exportExactSolutions)
            {
                exactPressPtr.reset ( new vector_Type (exactPress, exporter->mapType() ) );
                pFESpace->interpolate ( static_cast<typename feSpace_Type::function_Type> ( problem_Type::pexact ), *exactPressPtr, 0 );
                exactVelPtr.reset ( new vector_Type (exactVel, exporter->mapType() ) );
                uFESpace->interpolate ( static_cast<typename feSpace_Type::function_Type> ( problem_Type::uexact ), *exactVelPtr, 0 );
            }

            exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace,
                                    velAndPressure, UInt (0) );
            exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace,
                                    velAndPressure, pressureOffset );
            if (M_exportExactSolutions)
            {
                exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "exactVelocity", uFESpace,
                                        exactVelPtr, UInt (0) );
                exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "exactPressure", pFESpace,
                                        exactPressPtr, UInt (0) );
            }
            exporter->postProcess ( 0 );

            initChrono.stop();
            if (verbose)
            {
                std::cout << "Initialization time: " << initChrono.diff() << " s." << std::endl;
            }

            //vector_Type oldVel(uFESpace->map(), Unique);
            //oldVel.subset(*velAndPressure, uFESpace->map(),0,0);


            // +-----------------------------------------------+
            // |             Solving the problem               |
            // +-----------------------------------------------+
            if (verbose)
            {
                std::cout << std::endl << "[Solving the problem]" << std::endl;
            }
            int iter = 1;

            for ( ; time <= tFinal + dt / 2.; time += dt, iter++)
            {
                iterChrono.reset();
                iterChrono.start();

                oseenData->dataTime()->setTime (time);

                if (verbose)
                {
                    std::cout << "[t = " << oseenData->dataTime()->time() << " s.]" << std::endl;
                }

                double alpha = bdf.bdfVelocity().coefficientFirstDerivative ( 0 ) / oseenData->dataTime()->timeStep();

                bdf.bdfVelocity().extrapolation (beta); // Extrapolation for the convective term

                bdf.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );

                fluid.setVelocityRhs(bdf.bdfVelocity().rhsContributionFirstDerivative());

                fluid.getDisplayer().leaderPrint ("alpha ", alpha);
                fluid.getDisplayer().leaderPrint ("\n");
                fluid.getDisplayer().leaderPrint ("norm beta ", beta.norm2() );
                fluid.getDisplayer().leaderPrint ("\n");
                fluid.getDisplayer().leaderPrint ("norm rhs  ", rhs.norm2() );
                fluid.getDisplayer().leaderPrint ("\n");

//                if (oseenData->conservativeFormulation() )
//                {
                    rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();
//                }

                fluid.updateSystem ( alpha, beta, rhs );

//                if (!oseenData->conservativeFormulation() )
//                {
//                    rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();
//                }

                fluid.iterate ( bcH );

                Real uh1error = 0.0;
                fluid.h1normVelocity(uh1error);

                bdf.bdfVelocity().shiftRight ( *fluid.solution() );

                // Computation of the error
                LifeV::Real urelerr, prelerr, ul2error, pl2error;

                computeErrors (*fluid.solution(),
                               ul2error, urelerr, uFESpace,
                               pl2error, prelerr, pFESpace,
                               time);
                if (verbose)
                {
                	std::cout << "\n[ERRORS]:\n";
                	std::cout << "Time: " << time  << ", "
                			<< "L2 velocity error: " << ul2error << ", "
                			<< "H1 velocity error: " << uh1error << ", "
                			<< "L2 pressure error: " << pl2error << "\n" << std::flush;
                }

                if (verbose && M_exportNorms)
                {
                	 M_outNorm << time  << " "
                			 << ul2error << " "
                			 << uh1error << " "
                			 << pl2error << "\n" << std::flush;
                }

                // Saving the errors for the final test
                uL2Error[iFELabel][jDiscretization] = ul2error;
                pL2Error[iFELabel][jDiscretization] = pl2error;

                // Exporting the solution
                *velAndPressure = *fluid.solution();

                //oldVel *= 0;
                //oldVel.subset(*velAndPressure, uFESpace->map(),0,0);

                if (M_exportExactSolutions)
                {
                    pFESpace->interpolate ( static_cast<typename feSpace_Type::function_Type> ( problem_Type::pexact ), *exactPressPtr, time );
                    uFESpace->interpolate ( static_cast<typename feSpace_Type::function_Type> ( problem_Type::uexact ), *exactVelPtr, time );
                }
                exporter->postProcess ( time );


                MPI_Barrier (MPI_COMM_WORLD);

                iterChrono.stop();
                if (verbose)
                {
                    std::cout << "Iteration time: " << initChrono.diff() << " s." << std::endl << std::endl;
                }
            }

            if (verbose && M_exportNorms)
            {
                M_outNorm.close();
            }

            // ** BEGIN Accuracy test **
            if (M_test == Accuracy)
            {
                // Computation of the error
                LifeV::Real urelerr, prelerr, ul2error, pl2error;

                computeErrors (*fluid.solution(),
                               ul2error, urelerr, uFESpace,
                               pl2error, prelerr, pFESpace,
                               time);

                if (verbose) std::cout << "Relative error: E(u)=" << urelerr << ", E(p)=" << prelerr << std::endl
                                           << "Tolerance=" << M_accuracyTol << std::endl;

                /*
                if (urelerr > M_accuracyTol || prelerr > M_accuracyTol)
                {
                    if (verbose)
                    {
                        std::cout << "TEST_NAVIERSTOKES STATUS: ECHEC" << std::endl;
                    }
                    throw typename NavierStokes::RESULT_CHANGED_EXCEPTION();
                }
                */
            }
            // ** END Accuracy test **

            runChrono.stop();
            if (verbose)
            {
                std::cout << "Total run time: " << runChrono.diff() << " s." << std::endl;
            }
            if (verbose)
            {
                std::cout << "[[END_RUN_" << jDiscretization* numFELabels + iFELabel << "]]" << std::endl;
            }

        } // End of loop on the finite elements
    } // End of loop on the mesh refinement

    // ** BEGIN Space convergence test **
    /*
    if (verbose && (M_test == SpaceConvergence) )
    {
        bool success;
        success = checkConvergenceRate (M_uFELabels, uL2Error, M_uConvergenceOrder,
                                        M_pFELabels, pL2Error, M_pConvergenceOrder,
                                        M_meshDiscretization,
                                        M_convTol);

        if (!success)
        {
            if (verbose)
            {
                std::cout << "TEST_NAVIERSTOKES STATUS: ECHEC" << std::endl;
            }
            throw typename NavierStokes::RESULT_CHANGED_EXCEPTION();
        }
    }
    // ** END Space convergence test **
    */
    globalChrono.stop();
    if (verbose)
    {
        std::cout << std::endl << "Total simulation time:" << globalChrono.diff() << " s." << std::endl;
    }
    if (verbose && (M_test != None) )
    {
        std::cout << "TEST_NAVIERSTOKES_STATUS: SUCCESS" << std::endl;
    }
    if (verbose)
    {
        std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;
    }
}

#endif /* NAVIERSTOKES_H */
