/* davide forti*/

#undef HAVE_HDF5

#include <cassert>
#include <cstdlib>

#include <boost/timer.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <sys/stat.h>

// Exporters
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

// Files
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>

using namespace LifeV;

typedef VectorEpetra vector_Type;
typedef RegionMesh<LinearTetra> mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
typedef FESpace<mesh_Type, MapEpetra >  FESpace_Type;
typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;

int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    bool verbose = Comm->MyPID() == 0;
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm() );
    bool verbose = true;
#endif

    GetPot command_line (argc, argv);
    const bool check = command_line.search (2, "-c", "--check");

    const std::string dataFile = command_line.follow ("data_TimeAndExtrapolationHandler", 2, "-d", "--data_TimeAndExtrapolationHandler");
    GetPot data_file(dataFile);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  LOADING THE MESH
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    MeshData meshDataFluid;
    meshDataFluid.setup(data_file, "mesh_fluid");

    meshPtr_Type fluidMeshFull;

    if(Comm->MyPID()==0)
        std::cout << "\n[Loading the mesh ] .." << std::flush;

    fluidMeshFull.reset(new mesh_Type());
    readMesh(*fluidMeshFull, meshDataFluid);

    if(Comm->MyPID()==0)
        std::cout << " done! \n\n" << std::flush;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  PARTITIONING THE MESH
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(Comm->MyPID()==0)
        std::cout << "[Partitioning the fluid mesh ] \n\n";

    MeshPartitioner<mesh_Type > fluidPartio(fluidMeshFull, Comm);

    meshPtr_Type fluidMeshLocal;

    fluidMeshLocal.reset(new mesh_Type(*fluidPartio.meshPartition()));

    if(Comm->MyPID()==0)
        std::cout << "\n";

    FESpacePtr_Type velocityFESpace;
    velocityFESpace.reset (new FESpace_Type (fluidMeshLocal, "P1", 1, Comm) );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  TESTING THE TIME HANDLER OBJECT
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(Comm->MyPID()==0)
    	std::cout << "[Test with empty constructor ]";

    UInt orderBDF = data_file("dataTimeHandler/orderBDF",2);
    UInt orderExtrapolation = data_file("dataTimeHandler/orderExtrapolation",2);

    // Initialize the time handler object, set the order of the BDF and for the extrapolation
    TimeAndExtrapolationHandler timeVelocity;
    timeVelocity.setBDForder(orderBDF);
    timeVelocity.setMaximumExtrapolationOrder(orderExtrapolation);

    // Set the timestep
    timeVelocity.setTimeStep(0.2);

    // Initialize 3 vectors
    vector_Type firstVector(velocityFESpace->map());
    vector_Type secondVector(velocityFESpace->map());
    vector_Type thirdVector(velocityFESpace->map());
    firstVector  += 1;
    secondVector += 2;
    thirdVector  += 3;

    // Initialize the state with two vectors (we suppose that in the datafile orderBDF=2)
    std::vector<vector_Type> initialState;
    initialState.push_back(firstVector);
    initialState.push_back(secondVector);
    timeVelocity.initialize(initialState);

    if(Comm->MyPID()==0)
        std::cout << "\n";

    // empty vector to check the extrapolation
    vector_Type extrapolation(velocityFESpace->map());

    // Do an extrapolation
    timeVelocity.extrapolate(orderExtrapolation, extrapolation);
    extrapolation.spy("FirstExtrapolation");

    // Put in the stencil a new vector, and do a new extrapolation
    timeVelocity.shift(thirdVector);
    timeVelocity.extrapolate(orderExtrapolation, extrapolation);
    extrapolation.spy("SecondExtrapolation");

    // Compute the term that should go to the right hand side
    vector_Type rhsTerm(velocityFESpace->map());
    timeVelocity.rhsContribution(rhsTerm);
    rhsTerm.spy("rhsTerm");

    timeVelocity.state()[0].spy("StatoZero");
    timeVelocity.state()[1].spy("StatoUno");

    // Print the value of alpha
    std::cout << "\nAlpha = " << timeVelocity.alpha() << "\n\n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return 0;

}
