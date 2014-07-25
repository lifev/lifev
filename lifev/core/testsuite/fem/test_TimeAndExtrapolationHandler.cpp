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
    	std::cout << "[Test with empty constructor ] \n\n";

    TimeAndExtrapolationHandler timeVelocity;
    timeVelocity.setBDForder(1);
    timeVelocity.setMaximumExtrapolationOrder(3);

    vector_Type firstVector(velocityFESpace->map());
    vector_Type secondVector(velocityFESpace->map());
    vector_Type thirdVector(velocityFESpace->map());
    firstVector  += 1;
    secondVector += 2;
    thirdVector  += 3;

    std::vector<vector_Type> initialState;
    initialState.push_back(firstVector);
    initialState.push_back(secondVector);
    initialState.push_back(thirdVector);

    if(Comm->MyPID()==0)
        std::cout << "\n";

    timeVelocity.initialize(initialState);

    timeVelocity.shift(thirdVector);

    timeVelocity.setTimeStep(0.1);

    vector_Type extrapolation(velocityFESpace->map());
    timeVelocity.extrapolate(1,extrapolation);

    vector_Type rhsTerm(velocityFESpace->map());
    timeVelocity.rhsContribution(rhsTerm);

    std::cout << "\n\nAlpha = " << timeVelocity.alpha() << "\n\n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return 0;

}
