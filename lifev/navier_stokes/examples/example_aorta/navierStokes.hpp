/* -*- mode: c++ -*-

  This file is part of the LifeV library.
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
    @file navierStokes.hpp
    @author Davide Forti <davide.forti@epfl.ch>
    @date 2014-02-06
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
#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <algorithm>    // std::reverse

using namespace LifeV;

typedef VectorEpetra vector_Type;
typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

typedef FESpace< RegionMesh<LinearTetra>, MapEpetra > feSpace_Type;
typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

class NavierStokes
{
public:
    typedef RegionMesh<LinearTetra>                       mesh_Type;
    typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;
    typedef LifeV::OseenSolver< mesh_Type >               fluid_Type;
    typedef typename fluid_Type::vector_Type              vector_Type;
    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;
    typedef typename fluid_Type::matrix_Type              matrix_Type;

    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > > filterPtr_Type;

#ifdef HAVE_HDF5
    typedef LifeV::ExporterHDF5<mesh_Type>      hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>  hdf5FilterPtr_Type;
#endif

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
                   const std::string defaultDataName = "data");

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

    /*! @enum InitializationType
        Type of initialization. "Interpolation" just interpolates the value of the exact solution to the DoFs.
        "Projection" solves an Oseen problem where alpha=0, the convective term is linearized by using the exact solution for beta,
        and the time derivative is passed to the right hand side and computed from the exact solution.
     */
    enum InitializationType {Projection, Interpolation};


    struct Private;
    
    boost::shared_ptr<Private> M_data;
    
    // Initialization method
    InitializationType         M_initMethod;

    // output file name
    std::string                M_outputName;
    
    std::ofstream              M_out;
    
    bool                       M_exportCoeff;
};


using namespace LifeV;

Real zeroFunction(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0;
}

void preprocessBoundary(const Real& nx, const Real& ny, const Real& nz, BCHandler& bc, Real& Q_hat, const vectorPtr_Type& Phi_h, const UInt flag,
					    const feSpacePtr_Type& uFESpace_scalar,
						const boost::shared_ptr< ETFESpace<RegionMesh<LinearTetra>, MapEpetra, 3, 1 > >& uFESpace_ETA,
						vectorPtr_Type& V_hat_x, vectorPtr_Type& V_hat_y, vectorPtr_Type& V_hat_z)
{
	vectorPtr_Type Phi_h_inflow;
	Phi_h_inflow.reset ( new vector_Type ( uFESpace_scalar->map() ) );
	Phi_h_inflow->zero();

	bcManageRhs ( *Phi_h_inflow, *uFESpace_scalar->mesh(), uFESpace_scalar->dof(),  bc, uFESpace_scalar->feBd(), 1., 0.);

	*Phi_h_inflow *= *Phi_h;

	Q_hat = 0.0;

	// Computing the flowrate associated to Phi_h_inflow

	vectorPtr_Type Phi_h_inflow_rep;
	Phi_h_inflow_rep.reset ( new vector_Type ( *Phi_h_inflow, Repeated ) );

	vectorPtr_Type Q_hat_vec;
	Q_hat_vec.reset ( new vector_Type ( uFESpace_scalar->map() ) );
	Q_hat_vec->zero();

	vectorPtr_Type ones_vec;
	ones_vec.reset ( new vector_Type ( uFESpace_scalar->map() ) );
	ones_vec->zero();
	*ones_vec += 1.0;

	{
		using namespace ExpressionAssembly;
		QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria6pt) );
		integrate (
				boundary (uFESpace_ETA->mesh(), flag),
				myBDQR,
				uFESpace_ETA,
				value(uFESpace_ETA, *Phi_h_inflow_rep)*phi_i
		)
		>> Q_hat_vec;
	}

	Q_hat = Q_hat_vec->dot(*ones_vec);

	V_hat_x.reset ( new vector_Type ( *Phi_h_inflow ) );
	V_hat_y.reset ( new vector_Type ( *Phi_h_inflow ) );
	V_hat_z.reset ( new vector_Type ( *Phi_h_inflow ) );

	*V_hat_x *= nx;
	*V_hat_y *= ny;
	*V_hat_z *= nz;
}

Real oneFunction(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1.0;
}

Real inflowFunction(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
	Real Velocity = 1.0;

	switch (i)
	{
	case 0:
		return Velocity*0.07780;
		break;
	case 2:
		return Velocity*0.0;
		break;
	case 1:
		return Velocity*0.99696;
		break;
	}
	return 0;
}

struct NavierStokes::Private
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

 
NavierStokes::NavierStokes ( int argc,
                                                char** argv,
                                                const std::string defaultDataName)
    :
    M_data ( new Private )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow (defaultDataName.c_str(), 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    M_data->data_file_name = data_file_name;

    M_data->Re = dataFile ( "fluid/problem/Re", 1. );
    M_data->nu = dataFile ( "fluid/physics/viscosity", 1. ) /
                 dataFile ( "fluid/physics/density", 1. );



#ifdef EPETRA_MPI

    //    MPI_Init(&argc,&argv);

    M_data->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size (MPI_COMM_WORLD, &ntasks);
#else
    M_data->comm.reset ( new Epetra_SerialComm() );
#endif

}

 
void
NavierStokes::run()
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
    
    if (verbose)
        std::cout << std::endl << "[Loading the data]" << std::endl;
    
    GetPot dataFile ( M_data->data_file_name.c_str() );
    initChrono.stop();
    
    if (verbose)
        std::cout << "Initialization time (pre-run): " << initChrono.diff() << " s." << std::endl;
    
    if (verbose)
        std::cout << std::endl << "[[BEGIN_RUN]]" << std::endl;

    runChrono.reset();
    runChrono.start();
    initChrono.reset();
    initChrono.start();
    
    if (verbose)
        std::cout << "[Loading the mesh]" << std::endl;

    Int geoDimensions = mesh_Type::S_geoDimensions;
    
    /*
     * 	Handling offline/online mesh partitioning - BEGIN -
     */

    bool offlinePartio = dataFile ("offline_partitioner/useOfflinePartitionedMesh", false);

    boost::shared_ptr<mesh_Type > localMeshPtr;

    if ( offlinePartio )
    {
    	const std::string partsFileName (dataFile ("offline_partitioner/hdf5_file_name", "name.h5") );

    	if(verbose)
    		std::cout<< "   -- Loading a mesh which was already partitioned" << std::endl;

    	if (verbose)
    		std::cout << "   -- Partitioned mesh file: " << partsFileName << std::endl;

    	boost::shared_ptr<Epetra_MpiComm> comm = boost::dynamic_pointer_cast<Epetra_MpiComm>(M_data->comm);

    	PartitionIO<mesh_Type > partitionIO (partsFileName, comm);
    	partitionIO.read (localMeshPtr);
    }
    else
    {
        boost::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( M_data->comm ) );

        MeshData meshData;
        meshData.setup (dataFile, "fluid/space_discretization");
        readMesh (*fullMeshPtr, meshData);

        if (verbose) std::cout << "Mesh source: file("
            << meshData.meshDir() << meshData.meshFile() << ")" << std::endl;

        if ( verbose )
        {
            std::cout << "Mesh size max : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).maxH << std::endl;
        }

        if ( verbose )
        {
            std::cout << "Mesh size mean : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).meanH << std::endl;
        }


        if ( verbose )
        {
            std::cout << "Mesh size min : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).minH << std::endl;
        }

        if (verbose)
            std::cout << "Partitioning the mesh ... " << std::flush;

        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_data->comm);
        localMeshPtr = meshPart.meshPartition();

        fullMeshPtr.reset(); //Freeing the global mesh to save memory
    }
    
    /*
     * 	Handling offline/online mesh partitioning - END -
     */

    /*
     *	Reading if we need to store at each timestep or not - BEGIN
     */

    int saveEvery = dataFile("fluid/save_every",1);

    /*
     *	Reading if we need to store at each timestep or not - END
     */


    if (verbose)
        std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    
    std::string uOrder = dataFile("fluid/space_discretization/vel_order","P1");
    std::string pOrder = dataFile("fluid/space_discretization/pres_order","P1");;
    
    if (verbose)
        std::cout << "\tFE for the velocity: " << uOrder << std::endl
        << "\tFE for the pressure: " << pOrder << std::endl;
    
    if (verbose)
        std::cout << "\tBuilding the velocity FE space ... " << std::flush;
    
    feSpacePtr_Type uFESpace;
    uFESpace.reset (new feSpace_Type (localMeshPtr, uOrder, geoDimensions, M_data->comm) );
    
    if (verbose)
        std::cout << "ok." << std::endl;
    
    if (verbose)
        std::cout << "\tBuilding the pressure FE space ... " << std::flush;
    
    feSpacePtr_Type pFESpace;
    pFESpace.reset (new feSpace_Type (localMeshPtr, pOrder, 1, M_data->comm) );
    
    if (verbose)
        std::cout << "ok." << std::endl;
    
    UInt totalVelDof   = uFESpace->dof().numTotalDof();
    UInt totalPressDof = pFESpace->dof().numTotalDof();
    
    // Pressure offset in the vector
    UInt pressureOffset =  uFESpace->fieldDim() * uFESpace->dof().numTotalDof();
    
    if (verbose)
        std::cout << "\tTotal Velocity Dof = " << totalVelDof*3 << std::endl;
    
    if (verbose)
        std::cout << "\tTotal Pressure Dof = " << totalPressDof << std::endl;
    
    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose)
        std::cout << std::endl << "[Boundary conditions]" << std::endl;

    BCFunctionBase zero( zeroFunction );
    BCFunctionBase one ( oneFunction  );
	BCFunctionBase uInflow( inflowFunction );
    
	// BCs for the NS problem

    BCHandler bcH_aorta;

    bcH_aorta.addBC( "Inflow",         3, Essential,      Full,        zero, 3 );
    bcH_aorta.addBC( "Walls",        200, Essential, 	  Full, 	   zero, 3 );
    bcH_aorta.addBC( "Outflow2",       2,   Natural, 	Normal,    	   zero    );
    bcH_aorta.addBC( "Outflow4",       4, 	Natural, 	Normal,    	   zero    );
    bcH_aorta.addBC( "Outflow5",       5, 	Natural, 	Normal,    	   zero    );
    bcH_aorta.addBC( "Outflow6",       6, 	Natural, 	Normal,    	   zero    );
    bcH_aorta.addBC( "Outflow7",       7, 	Natural, 	Normal,        zero    );
    bcH_aorta.addBC( "Outflow8",       8, 	Natural, 	Normal,    	   zero    );
    bcH_aorta.addBC( "Outflow9",       9, 	Natural, 	Normal,    	   zero    );

    bcH_aorta.bcUpdate ( *localMeshPtr, uFESpace->feBd(), uFESpace->dof() );

    // +-----------------------------------------------+
    // |             Solving the laplacian             |
    // +-----------------------------------------------+

    if (verbose)
            std::cout << std::endl << "[Solving the laplacian for the inflow profile]" << std::endl;

    feSpacePtr_Type uFESpace_scalar;
    uFESpace_scalar.reset (new feSpace_Type (localMeshPtr, uOrder, 1, M_data->comm) );

    // BCs for the laplacian problem

    BCHandler bcH_laplacian;
    bcH_laplacian.addBC( "Walls", 200, Essential, Full, zero, 1 );
    bcH_laplacian.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    vectorPtr_Type Phi_h;
    Phi_h.reset ( new vector_Type ( uFESpace_scalar->map() ) );
    Phi_h->zero();

    vectorPtr_Type rhs_laplacian;
    rhs_laplacian.reset ( new vector_Type ( uFESpace_scalar->map() ) );
    rhs_laplacian->zero();
    *rhs_laplacian += 1.0;

    boost::shared_ptr<MatrixEpetra<Real> > Laplacian;
    Laplacian.reset ( new MatrixEpetra<Real> ( uFESpace_scalar->map() ) );

    boost::shared_ptr< ETFESpace<mesh_Type, MapEpetra, 3, 1 > > uFESpace_ETA;
    uFESpace_ETA.reset( new ETFESpace<mesh_Type, MapEpetra, 3, 1 > ( uFESpace_scalar->mesh(), &(uFESpace_scalar->refFE()), M_data->comm) );

    {
    	using namespace ExpressionAssembly;

    	integrate(
    				elements(uFESpace_ETA->mesh()),
    				uFESpace_scalar->qr(),
    				uFESpace_ETA,
    				uFESpace_ETA,
    				dot( grad(phi_i) , grad(phi_j) )
    			 ) >> Laplacian;
    }

    bcManage ( *Laplacian, *rhs_laplacian, *uFESpace_scalar->mesh(), uFESpace_scalar->dof(), bcH_laplacian, uFESpace_scalar->feBd(), 1.0, 0.0 );
    Laplacian->globalAssemble();

    boost::shared_ptr<SolverAztecOO> linearSolver_laplacian;
    linearSolver_laplacian.reset( new SolverAztecOO (uFESpace_scalar->map().commPtr()) );

    linearSolver_laplacian->setupPreconditioner ( dataFile, "laplacian/prec" );
    linearSolver_laplacian->setDataFromGetPot ( dataFile, "laplacian/solver" );

    linearSolver_laplacian->setMatrix ( *Laplacian );

    boost::shared_ptr<MatrixEpetra<Real> > staticCast_laplacian = boost::static_pointer_cast<MatrixEpetra<Real> > (Laplacian);
    Int numIter_laplacian = linearSolver_laplacian->solveSystem ( *rhs_laplacian, *Phi_h, staticCast_laplacian );

    // Inflow

    Real nx_inflow = 0.07780;
    Real ny_inflow = 0.0;
    Real nz_inflow = 0.99696;
    Real Q_hat_inflow = 0.0;
    vectorPtr_Type V_hat_x_inflow;
    vectorPtr_Type V_hat_y_inflow;
    vectorPtr_Type V_hat_z_inflow;
    BCHandler bcH_laplacian_inflow;
    bcH_laplacian_inflow.addBC( "Inflow", 3, Essential, Full, one, 1 );
    bcH_laplacian_inflow.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    preprocessBoundary(nx_inflow, ny_inflow, nz_inflow, bcH_laplacian_inflow, Q_hat_inflow, Phi_h, 3, uFESpace_scalar, uFESpace_ETA, V_hat_x_inflow, V_hat_y_inflow, V_hat_z_inflow);

    if (verbose)
        	std::cout << "\tValue of inflow, Q_hat = " << Q_hat_inflow << std::endl;

    // Outflow 4
    Real nx_flag4 = 0.0;
    Real ny_flag4 = 0.206;
    Real nz_flag4 = 0.978;
    Real Q_hat_flag4 = 0.0;
    vectorPtr_Type V_hat_x_flag4;
    vectorPtr_Type V_hat_y_flag4;
    vectorPtr_Type V_hat_z_flag4;
    BCHandler bcH_laplacian_flag4;
    bcH_laplacian_flag4.addBC( "Outflow4", 4, Essential, Full, one, 1 );
    bcH_laplacian_flag4.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    preprocessBoundary(nx_flag4, ny_flag4, nz_flag4, bcH_laplacian_flag4, Q_hat_flag4, Phi_h, 4, uFESpace_scalar, uFESpace_ETA, V_hat_x_flag4, V_hat_y_flag4, V_hat_z_flag4);

    if (verbose)
            	std::cout << "\tValue of outflow 4, Q_hat = " << Q_hat_flag4 << std::endl;

    // Outflow 5
    Real nx_flag5 = 0.0;
    Real ny_flag5 = 0.0;
    Real nz_flag5 = 1.0;
    Real Q_hat_flag5 = 0.0;
    vectorPtr_Type V_hat_x_flag5;
    vectorPtr_Type V_hat_y_flag5;
    vectorPtr_Type V_hat_z_flag5;
    BCHandler bcH_laplacian_flag5;
    bcH_laplacian_flag5.addBC( "Outflow5", 5, Essential, Full, one, 1 );
    bcH_laplacian_flag5.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    preprocessBoundary(nx_flag5, ny_flag5, nz_flag5, bcH_laplacian_flag5, Q_hat_flag5, Phi_h, 5, uFESpace_scalar, uFESpace_ETA, V_hat_x_flag5, V_hat_y_flag5, V_hat_z_flag5);

    if (verbose)
    	std::cout << "\tValue of outflow 5, Q_hat = " << Q_hat_flag5 << std::endl;

    // Outflow 6
    Real nx_flag6 = 0.0;
    Real ny_flag6 = 0.51;
    Real nz_flag6 = 0.86;
    Real Q_hat_flag6 = 0.0;
    vectorPtr_Type V_hat_x_flag6;
    vectorPtr_Type V_hat_y_flag6;
    vectorPtr_Type V_hat_z_flag6;
    BCHandler bcH_laplacian_flag6;
    bcH_laplacian_flag6.addBC( "Outflow6", 6, Essential, Full, one, 1 );
    bcH_laplacian_flag6.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    preprocessBoundary(nx_flag6, ny_flag6, nz_flag6, bcH_laplacian_flag6, Q_hat_flag6, Phi_h, 6, uFESpace_scalar, uFESpace_ETA, V_hat_x_flag6, V_hat_y_flag6, V_hat_z_flag6);

    if (verbose)
    	std::cout << "\tValue of outflow 6, Q_hat = " << Q_hat_flag6 << std::endl;

    // Outflow 7
    Real nx_flag7 = 0.0;
    Real ny_flag7 = 0.807;
    Real nz_flag7 = 0.589;
    Real Q_hat_flag7 = 0.0;
    vectorPtr_Type V_hat_x_flag7;
    vectorPtr_Type V_hat_y_flag7;
    vectorPtr_Type V_hat_z_flag7;
    BCHandler bcH_laplacian_flag7;
    bcH_laplacian_flag7.addBC( "Outflow7", 7, Essential, Full, one, 1 );
    bcH_laplacian_flag7.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    preprocessBoundary(nx_flag7, ny_flag7, nz_flag7, bcH_laplacian_flag7, Q_hat_flag7, Phi_h, 7, uFESpace_scalar, uFESpace_ETA, V_hat_x_flag7, V_hat_y_flag7, V_hat_z_flag7);

    if (verbose)
    	std::cout << "\tValue of outflow 7, Q_hat = " << Q_hat_flag7 << std::endl;

    // Outflow 8
    Real nx_flag8 = 0.0;
    Real ny_flag8 = 0.185;
    Real nz_flag8 = 0.983;
    Real Q_hat_flag8 = 0.0;
    vectorPtr_Type V_hat_x_flag8;
    vectorPtr_Type V_hat_y_flag8;
    vectorPtr_Type V_hat_z_flag8;
    BCHandler bcH_laplacian_flag8;
    bcH_laplacian_flag8.addBC( "Outflow5", 8, Essential, Full, one, 1 );
    bcH_laplacian_flag8.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    preprocessBoundary(nx_flag8, ny_flag8, nz_flag8, bcH_laplacian_flag8, Q_hat_flag8, Phi_h, 8, uFESpace_scalar, uFESpace_ETA, V_hat_x_flag8, V_hat_y_flag8, V_hat_z_flag8);

    if (verbose)
    	std::cout << "\tValue of outflow 8, Q_hat = " << Q_hat_flag8 << std::endl;

    // Outflow 9
    Real nx_flag9 = 0.0;
    Real ny_flag9 = 0.851;
    Real nz_flag9 = -0.525;
    Real Q_hat_flag9 = 0.0;
    vectorPtr_Type V_hat_x_flag9;
    vectorPtr_Type V_hat_y_flag9;
    vectorPtr_Type V_hat_z_flag9;
    BCHandler bcH_laplacian_flag9;
    bcH_laplacian_flag9.addBC( "Outflow5", 9, Essential, Full, one, 1 );
    bcH_laplacian_flag9.bcUpdate ( *localMeshPtr, uFESpace_scalar->feBd(), uFESpace_scalar->dof() );

    preprocessBoundary(nx_flag9, ny_flag9, nz_flag9, bcH_laplacian_flag9, Q_hat_flag9, Phi_h, 9, uFESpace_scalar, uFESpace_ETA, V_hat_x_flag9, V_hat_y_flag9, V_hat_z_flag9);

    if (verbose)
    	std::cout << "\tValue of outflow 9, Q_hat = " << Q_hat_flag9 << std::endl;

    // Cleaning useless things
    linearSolver_laplacian.reset();
    Phi_h.reset();
    Laplacian.reset();
    staticCast_laplacian.reset();
    rhs_laplacian.reset();
    Phi_h_inflow_rep.reset();

    /*
    // +-----------------------------------------------+
    // |             Creating the problem              |
    // +-----------------------------------------------+

    if (verbose)
        std::cout << std::endl << "[Creating the problem]" << std::endl;
    
    boost::shared_ptr<OseenData> oseenData (new OseenData() );
    oseenData->setup ( dataFile );
    
    if (verbose)
        std::cout << "Time discretization order " << oseenData->dataTimeAdvance()->orderBDF() << std::endl;
    
    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpace,
                                    *pFESpace,
                                    M_data->comm);
    
    MapEpetra fullMap (fluid.getMap() );
    
    fluid.setUp (dataFile);
    
    fluid.buildSystem();
    
    // +-----------------------------------------------+
    // |       Initialization of the simulation        |
    // +-----------------------------------------------+
    if (verbose)
        std::cout << std::endl << "[Initialization of the simulation]" << std::endl;
    
    Real dt       = dataFile("fluid/time_discretization/timestep",0.0);
    Real t0       = dataFile("fluid/time_discretization/initialtime",0.0);
    Real tFinal   = dataFile("fluid/time_discretization/endtime",0.0);
    UInt orderBDF = dataFile("fluid/time_discretization/BDF_order",2);
    
    // Time handler objects to deal with time advancing and extrapolation
    TimeAndExtrapolationHandler timeVelocity;
    TimeAndExtrapolationHandler timePressure;

    // Order of BDF and extrapolation for the velocity
    timeVelocity.setBDForder(orderBDF);
    timeVelocity.setMaximumExtrapolationOrder(orderBDF);
    timeVelocity.setTimeStep(dt);

    // Order of BDF and extrapolation for the pressure
    timePressure.setBDForder(orderBDF);
    timePressure.setMaximumExtrapolationOrder(orderBDF);
    timePressure.setTimeStep(dt);

    if (verbose)
        std::cout << "Computing the initial solution ... " << std::endl;
    
    vector_Type beta ( uFESpace->map() );
    
    M_outputName = dataFile ( "exporter/filename", "result");
    boost::shared_ptr< Exporter<mesh_Type > > exporter;

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

    vectorPtr_Type velAndPressure;
    velAndPressure.reset ( new vector_Type (*fluid.solution(), exporter->mapType() ) );

    vectorPtr_Type velAndPressure_inflow_reference;
    velAndPressure_inflow_reference.reset ( new vector_Type (*fluid.solution() ) );

    vectorPtr_Type velAndPressure_inflow;
    velAndPressure_inflow.reset ( new vector_Type (*fluid.solution() ) );

    vectorPtr_Type betaFull;
    betaFull.reset ( new vector_Type (*fluid.solution(), exporter->mapType() ) );

    vectorPtr_Type pressureFull;
    pressureFull.reset ( new vector_Type (*fluid.solution(), exporter->mapType() ) );

    vectorPtr_Type rhs;
    rhs.reset ( new vector_Type (*fluid.solution(), exporter->mapType() ) );

    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace, velAndPressure, UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace, velAndPressure, pressureOffset );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "Laplacian", uFESpace_scalar, Phi_h_inflow, UInt (0) );
	*/

    /*
     *  Starting from scratch or restarting? -BEGIN-
     */

    /*
    bool doRestart = dataFile("importer/restart", false);

    Real time = t0;

    oseenData->dataTime()->setTime (t0);

    // initialize stencils
	vector_Type velocityInitial ( uFESpace->map() );
	vector_Type pressureInitial ( pFESpace->map() );

	std::vector<vector_Type> initialStateVelocity;
	std::vector<vector_Type> initialStatePressure;

	vector_Type vectorForInitializationVelocitySolution(*fluid.solution(), Unique);
	vector_Type vectorForInitializationPressureSolution(*fluid.solution(), Unique);

    if (!doRestart)
    {
    	// if you start from scratch
    	velocityInitial *= 0 ;
    	pressureInitial *= 0;

    	if(orderBDF==1)
    	{
    		initialStateVelocity.push_back(velocityInitial);
    		initialStatePressure.push_back(pressureInitial);
    	}
    	else if (orderBDF==2)
    	{
    		initialStateVelocity.push_back(velocityInitial);
    		initialStateVelocity.push_back(velocityInitial);

    		initialStatePressure.push_back(pressureInitial);
    		initialStatePressure.push_back(pressureInitial);
    	}

    	timeVelocity.initialize(initialStateVelocity);
    	timePressure.initialize(initialStatePressure);
    }
    else
    {
    	std::string const importerType  =  dataFile ( "importer/type", "hdf5");
    	std::string const fileName      =  dataFile ( "importer/filename", "SolutionRestarted");
    	std::string const initialLoaded =  dataFile ( "importer/initSol", "NO_DEFAULT_VALUE");

    	boost::shared_ptr<hdf5Filter_Type> importer;
    	importer.reset ( new  hdf5Filter_Type ( dataFile, fileName) );

    	importer->setMeshProcId (uFESpace->mesh(), M_data->comm->MyPID() );

    	std::string iterationString;
    	iterationString = initialLoaded;

    	for (UInt iterInit = 0; iterInit < orderBDF; iterInit++ )
    	{
    		vectorPtr_Type velocityRestart;
    		velocityRestart.reset ( new vector_Type (uFESpace->map(),  Unique ) );
    		*velocityRestart *= 0.0;

    		vectorPtr_Type pressureRestart;
    		pressureRestart.reset ( new vector_Type (pFESpace->map(),  Unique ) );
    		*pressureRestart *= 0.0;

    		LifeV::ExporterData<mesh_Type> velocityReader (LifeV::ExporterData<mesh_Type>::VectorField,
    													   std::string ("velocity." + iterationString),
    													   uFESpace,
    													   velocityRestart,
    													   UInt (0),
    													   LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    		LifeV::ExporterData<mesh_Type> pressureReader (LifeV::ExporterData<mesh_Type>::ScalarField,
    													   std::string ("pressure." + iterationString),
    													   pFESpace,
    													   pressureRestart,
    													   UInt (0),
    													   LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    		importer->readVariable (velocityReader);
    		importer->readVariable (pressureReader);

    		int iterations = std::atoi (iterationString.c_str() );
    		Real iterationsReal = iterations;

    		velocityInitial = *velocityRestart;
    		pressureInitial = *pressureRestart;
    		velAndPressure->subset(velocityInitial,  uFESpace->map(), 0, 0);
    		velAndPressure->subset(pressureInitial,  pFESpace->map(), 0, pressureOffset);
    		//exporter->postProcess ( iterationsReal*0.0025 + 0.0025);

    		if ( iterInit == 0 )
    		{
    			vectorForInitializationVelocitySolution = velocityInitial;
    			vectorForInitializationPressureSolution = pressureInitial;
    		}


    		iterations--;

    		std::ostringstream iter;
    		iter.fill ( '0' );
    		iter << std::setw (5) << ( iterations );
    		iterationString = iter.str();

    		initialStateVelocity.push_back(*velocityRestart);
    		initialStatePressure.push_back(*pressureRestart);
    	}

    	// For BDF 1 it does not change anything, for BDF2 it is necessary
    	std::reverse(initialStateVelocity.begin(),initialStateVelocity.end());
    	std::reverse(initialStatePressure.begin(),initialStatePressure.end());

    	timeVelocity.initialize(initialStateVelocity);
    	timePressure.initialize(initialStatePressure);
    	importer->closeFile();
    }
	*/

    /*
     *  Starting from scratch or restarting? -END-
     */

    /*
    velAndPressure_inflow_reference->zero();
    velAndPressure_inflow_reference->subset(*V_hat_x,uFESpace_scalar->map(), 0, 0);
    velAndPressure_inflow_reference->subset(*V_hat_y,uFESpace_scalar->map(), 0, uFESpace_scalar->dof().numTotalDof() );
    velAndPressure_inflow_reference->subset(*V_hat_z,uFESpace_scalar->map(), 0, 2*uFESpace_scalar->dof().numTotalDof() );

    initChrono.stop();
    
    if (verbose)
        std::cout << "Initialization time: " << initChrono.diff() << " s." << std::endl;
    
    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    if (verbose)
        std::cout << std::endl << "[Solving the problem]" << std::endl;
    
    int iter = 1;
    time = t0 + dt;
    
    int n_iter = 0;

    // initialize stencils
    vector_Type pressure ( pFESpace->map() );
    vector_Type rhsVelocity ( uFESpace->map() );
    vector_Type rhsVelocityFull ( fullMap );
	*/

    /*
     * 	Compute the vector for computing Loads - END
     */

    /*

    if (doRestart)
    {
    	fluid.initializeVelocitySolution(vectorForInitializationVelocitySolution);
    	fluid.initializePressureSolution(vectorForInitializationPressureSolution);
    }

    for ( ; time <= tFinal + dt / 2.; time += dt, iter++)
    {
        n_iter++;
        iterChrono.reset();
        iterChrono.start();
        
        oseenData->dataTime()->setTime (time);
        
        if (verbose)
            std::cout << "[t = " << oseenData->dataTime()->time() << " s.]" << std::endl;
        
        double alpha = timeVelocity.alpha() / dt;
        
        beta *= 0;
        timeVelocity.extrapolate (orderBDF, beta); // Extrapolation for the convective term
        *betaFull *= 0;
        betaFull->subset(beta, uFESpace->map(), 0, 0);

        pressure *= 0;
        timePressure.extrapolate (orderBDF, pressure); // Extrapolation for the LES terms
        *pressureFull *= 0;
        pressureFull->subset(pressure, pFESpace->map(), 0, pressureOffset);

        rhsVelocity *= 0;
        timeVelocity.rhsContribution (rhsVelocity);
        rhsVelocityFull *= 0;
        rhsVelocityFull.subset(rhsVelocity, uFESpace->map(), 0, 0);

        fluid.setVelocityRhs(rhsVelocityFull);
        fluid.setPressureExtrapolated(*pressureFull);

        *rhs  = fluid.matrixMass() * rhsVelocityFull;
        
        fluid.updateSystem ( alpha, *betaFull, *rhs );
        

        // ---------------------------------
        // Evaluation of the inflow velocity
        // ---------------------------------

        Real Q_in = 0;
        Real T_heartbeat = 0.8;

        if ( (time >= 0.05 && time <= 0.42) || (time >= (0.05+T_heartbeat) && time <= (0.42+T_heartbeat) ) || (time >= (0.05+2*T_heartbeat) && time <= (0.42+2*T_heartbeat) ) || (time >= (0.05+3*T_heartbeat) && time <= (0.42+3*T_heartbeat) ) )
        	Q_in = 2.422818092859456e+8*std::pow(time,8)-4.764207344433996e+8*std::pow(time,7) + 3.993883831476327e+8*std::pow(time,6) -1.867066900011057e+8*std::pow(time,5) +0.533079809563519e+8*std::pow(time,4) -0.094581323616832e+8*std::pow(time,3) +0.009804512311267e+8*std::pow(time,2) -0.000482942399225e+8*time+0.000008651437192e+8;
        else
        	Q_in = 0.0;

        Real alpha_flowrate = Q_in/Q_hat;

        if (verbose)
        {
        	std::cout << "Q_in: " << Q_in << std::endl << std::endl;
        	std::cout << "alpha_flowrate: " << alpha_flowrate << std::endl << std::endl;
        }

        velAndPressure_inflow->zero();
        *velAndPressure_inflow += *velAndPressure_inflow_reference;
        *velAndPressure_inflow *= alpha_flowrate;

        fluid.iterate ( bcH_aorta, velAndPressure_inflow );
        
        vector_Type computedVelocity(uFESpace->map());
        computedVelocity.subset(*fluid.solution(), uFESpace->map(), 0, 0);

        vector_Type computedPressure(pFESpace->map());
        computedPressure.subset(*fluid.solution(), pressureOffset );

        timeVelocity.shift(computedVelocity);
        timePressure.shift(computedPressure);
        
        // Exporting the solution
        *velAndPressure = *fluid.solution();
        
        if( n_iter%saveEvery == 0 || (n_iter+orderBDF-1)%saveEvery == 0)
        {
        	exporter->postProcess ( time );
        }
        
        iterChrono.stop();
        
        if (verbose)
            std::cout << "Iteration time: " << initChrono.diff() << " s." << std::endl << std::endl;
    }
    
    exporter->closeFile();
	*/

    runChrono.stop();
    
    if (verbose)
        std::cout << "Total run time: " << runChrono.diff() << " s." << std::endl;
    
    if (verbose)
        std::cout << "[[END_RUN]]" << std::endl;

    globalChrono.stop();

    if (verbose)
        std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;
}

#endif /* NAVIERSTOKES_H */
