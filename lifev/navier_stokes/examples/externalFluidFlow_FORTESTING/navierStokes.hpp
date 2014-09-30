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
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

using namespace LifeV;


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

Real oneFunction(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1;
}

Real inflowFunction(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if (i == 0)
    {
        Real ux = 22.0;              // flat velocity profile
        return (1 * ux * t/0.15 * (t < 0.15) + 1 * ux * (t >= 0.15));
    }
    else
    {
        return 0;
    }
    
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
    
    M_exportCoeff = dataFile ("fluid/export_coefficients", false);
    
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
    int orderBDF  = dataFile("fluid/time_discretization/BDF_order",1);

    /*
     *	Reading if we need to store at each timestep or not - END
     */

    if (verbose)
        std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    
    std::string uOrder = dataFile("fluid/space_discretization/vel_order","P1");
    std::string pOrder = dataFile("fluid/space_discretization/pres_order","P1");;
    
    if (verbose)
        std::cout << "FE for the velocity: " << uOrder << std::endl
        << "FE for the pressure: " << pOrder << std::endl;
    
    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;
    
    feSpacePtr_Type uFESpace;
    uFESpace.reset (new feSpace_Type (localMeshPtr, uOrder, geoDimensions, M_data->comm) );
    
    if (verbose)
        std::cout << "ok." << std::endl;
    
    
    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;
    
    feSpacePtr_Type pFESpace;
    pFESpace.reset (new feSpace_Type (localMeshPtr, pOrder, 1, M_data->comm) );
    
    if (verbose)
        std::cout << "ok." << std::endl;
    
    
    UInt totalVelDof   = uFESpace->dof().numTotalDof();
    UInt totalPressDof = pFESpace->dof().numTotalDof();
    
    // Pressure offset in the vector
    UInt pressureOffset =  uFESpace->fieldDim() * uFESpace->dof().numTotalDof();
    
    if (verbose)
        std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
    
    if (verbose)
        std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;
    
    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose)
        std::cout << std::endl << "[Boundary conditions]" << std::endl;
    
    
    BCFunctionBase uZero( zeroFunction );
	BCFunctionBase uInflow( inflowFunction );
    
	std::vector<LifeV::ID> zComp(1), yComp(1), xComp(1);
    xComp[0] = 0;
	yComp[0] = 1;
	zComp[0] = 2;
    
    BCHandler bcH;

    bcH.addBC( "Outflow",        3, Natural,   Full,      uZero,   3 );
    bcH.addBC( "Inflow",         1, Essential, Full,      uInflow, 3 );
    bcH.addBC( "WallUpDown",     2, Essential, Component, uZero,   yComp );
    bcH.addBC( "Cylinder",       4, Essential, Full,      uZero,   3 );
    bcH.addBC( "WallLeftRight",  5, Essential, Component, uZero,   zComp );

    // If we change the FE we have to update the BCHandler (internal data)
    bcH.bcUpdate ( *localMeshPtr, uFESpace->feBd(), uFESpace->dof() );
    
    BCFunctionBase uOne( oneFunction );
    
    BCHandler bcHDrag;
    bcHDrag.addBC( "Cylinderr",   4, Essential, Component, uOne,   xComp ); // ATTENTO CAMBIA A SECONDA DEL FLAG DEL CILINDRO
    bcHDrag.bcUpdate ( *localMeshPtr, uFESpace->feBd(), uFESpace->dof() );
    
    BCHandler bcHLift;
    bcHLift.addBC( "Cylinderr",   4, Essential, Component, uOne,   yComp ); // ATTENTO CAMBIA A SECONDA DEL FLAG DEL CILINDRO
    bcHLift.bcUpdate ( *localMeshPtr, uFESpace->feBd(), uFESpace->dof() );
    
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
    
    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();
    
    // bdf object to store the previous solutions
    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup (oseenData->dataTimeAdvance()->orderBDF() );
    
    if (verbose)
        std::cout << "Computing the initial solution ... " << std::endl;
    
    vector_Type beta ( fullMap );
    vector_Type rhs ( fullMap );
    
    /*
     *  Starting from scratch or restarting? -BEGIN-
     */

    bool doRestart = dataFile("importer/restart", false);

    Real time = t0;
    oseenData->dataTime()->setTime (t0);

    // initialize stencils
    std::vector<vectorPtr_Type> solutionStencilVelocity;
    vectorPtr_Type velocity ( new vector_Type(uFESpace->map() ) );
    solutionStencilVelocity.resize ( bdf.bdfVelocity().size() );

    std::vector<vectorPtr_Type> solutionStencilPressure;
    vectorPtr_Type pressure ( new vector_Type(pFESpace->map() ) );
	solutionStencilPressure.resize ( bdf.bdfPressure().size() );

	// initialize vectors for printing the solution when we do the restart
	vectorPtr_Type velocityInitial ( new vector_Type(uFESpace->map() ) );
	vectorPtr_Type pressureInitial ( new vector_Type(pFESpace->map() ) );

    if (!doRestart)
    {
    	*velocity *= 0;
    	for ( UInt i(0); i < bdf.bdfVelocity().size() ; ++i)
    		solutionStencilVelocity[ i ] = fluid.solution();

    	*pressure *= 0;
    	for ( UInt i(0); i < bdf.bdfPressure().size() ; ++i)
    		solutionStencilPressure[ i ] = pressure;

    	bdf.bdfVelocity().setInitialCondition (solutionStencilVelocity);
    	bdf.bdfPressure().setInitialCondition (solutionStencilPressure);

    	*velocityInitial = *velocity;
    	*pressureInitial = *pressure;

    }
    else
    {
    	std::string const importerType  =  dataFile ( "importer/type", "hdf5");
    	std::string const fileName      =  dataFile ( "importer/filename", "SolutionRestarted");
    	std::string const initialLoaded =  dataFile ( "importer/initSol", "NO_DEFAULT_VALUE");

    	boost::shared_ptr<hdf5Filter_Type> importer;
    	importer.reset ( new  hdf5Filter_Type ( dataFile, fileName) );

    	importer->setMeshProcId (uFESpace->mesh(), M_data->comm->MyPID() );
    	vectorPtr_Type velocity (new vector_Type (uFESpace->map(), importer->mapType() ) );
    	vectorPtr_Type pressure (new vector_Type (pFESpace->map(), importer->mapType() ) );

    	std::string iterationString;
    	iterationString = initialLoaded;


    	for (UInt iterInit = 0; iterInit < bdf.bdfVelocity().size(); iterInit++ )
    	{
    		velocity.reset ( new vector_Type (uFESpace->map(),  Unique ) );
    		*velocity *= 0.0;


    		pressure.reset ( new vector_Type (pFESpace->map(),  Unique ) );
    		*pressure *= 0.0;

    		LifeV::ExporterData<mesh_Type> velocityReader (LifeV::ExporterData<mesh_Type>::VectorField,
    													   std::string ("velocity." + iterationString),
    													   uFESpace,
    													   velocity,
    													   UInt (0),
    													   LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    		LifeV::ExporterData<mesh_Type> pressureReader (LifeV::ExporterData<mesh_Type>::ScalarField,
    													   std::string ("pressure." + iterationString),
    													   pFESpace,
    													   pressure,
    													   UInt (0),
    													   LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    		importer->readVariable (velocityReader);
    		importer->readVariable (pressureReader);

    		vectorPtr_Type solutionVelocityPart (new vector_Type (fluid.solution()->map(), Unique, Zero));
    		*solutionVelocityPart *= 0;
    		solutionVelocityPart->subset(*velocity, velocity->map(), UInt(0), UInt(0) );

    		solutionStencilVelocity[ iterInit ] = solutionVelocityPart;
    		solutionStencilPressure[ iterInit ] = pressure;

    		int iterations = std::atoi (iterationString.c_str() );
    		iterations--;

    		std::ostringstream iter;
    		iter.fill ( '0' );
    		iter << std::setw (5) << ( iterations );
    		iterationString = iter.str();

    		if (iterInit == 0)
    		{
    			*velocityInitial = *velocity;
    			*pressureInitial = *pressure;
    		}

    	}

    	bdf.bdfVelocity().setInitialCondition (solutionStencilVelocity);
    	bdf.bdfPressure().setInitialCondition (solutionStencilPressure);
    	importer->closeFile();
    }

    bdf.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );

    bdf.bdfPressure().updateRHSContribution ( oseenData->dataTime()->timeStep() );

    /*
     *  Starting from scratch or restarting? -END-
     */

    M_outputName = dataFile ( "exporter/filename", "result");
    boost::shared_ptr< Exporter<mesh_Type > > exporter;
    
    vectorPtr_Type velAndPressure;
    
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

    velAndPressure->subset(*velocityInitial,  uFESpace->map(), 0, 0);
    velAndPressure->subset(*pressureInitial,  pFESpace->map(), 0, pressureOffset);

    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace, velAndPressure, UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace, velAndPressure, pressureOffset );
    
    exporter->postProcess ( time );
    
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
    Real drag(0.0);
    VectorSmall<2> AerodynamicCoefficients;
    
    if(verbose && M_exportCoeff)
    {
        M_out.open ("Coefficients.txt" );
        M_out << "% time / drag / lift / energy \n" << std::flush;
    }
    
    int n_iter = 0;

    for ( ; time <= tFinal + dt / 2.; time += dt, iter++)
    {
        n_iter++;
        iterChrono.reset();
        iterChrono.start();
        
        oseenData->dataTime()->setTime (time);
        
        if (verbose)
            std::cout << "[t = " << oseenData->dataTime()->time() << " s.]" << std::endl;
        
        double alpha = bdf.bdfVelocity().coefficientFirstDerivative ( 0 ) / oseenData->dataTime()->timeStep();
        

        bdf.bdfVelocity().extrapolation (beta); // Extrapolation for the convective term
        
        *pressure *= 0;
        bdf.bdfPressure().extrapolation (*pressure); // Extrapolation for the LES terms

        bdf.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );
        
        bdf.bdfPressure().updateRHSContribution ( oseenData->dataTime()->timeStep() );

        fluid.setVelocityRhs(bdf.bdfVelocity().rhsContributionFirstDerivative());
        
        fluid.setPressureExtrapolated(*pressure);

        rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();
        
        fluid.updateSystem ( alpha, beta, rhs );
        
        fluid.iterate ( bcH );
        
        Real velInfty = (1 * 22 * time/0.15 * (time < 0.15) + 1 * 22 * (time >= 0.15));
        
        AerodynamicCoefficients = fluid.computeDrag(1, bcHDrag, bcHLift, velInfty, 0.25*fluid.area(4) );
        
        Real energy = fluid.energy();

        if ( verbose && M_exportCoeff )
        {
            M_out << time  << " "
            << AerodynamicCoefficients[0] << " "
            << AerodynamicCoefficients[1] << " "
            << energy << "\n" << std::flush;
        }
        

        vector_Type computedPressure(pFESpace->map());
        computedPressure.subset(*fluid.solution(), pressureOffset );

        bdf.bdfVelocity().shiftRight ( *fluid.solution() );
        bdf.bdfPressure().shiftRight ( computedPressure );
        
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
    
    runChrono.stop();
    
    if (verbose)
        std::cout << "Total run time: " << runChrono.diff() << " s." << std::endl;
    
    if (verbose)
        std::cout << "[[END_RUN]]" << std::endl;
    
    if (verbose && M_exportCoeff)
    {
        M_out.close();
    }
    
    globalChrono.stop();
    
    if (verbose)
        std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;
}

#endif /* NAVIERSTOKES_H */
