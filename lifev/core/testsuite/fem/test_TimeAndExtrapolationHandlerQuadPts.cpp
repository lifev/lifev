/* davide forti*/

#undef HAVE_HDF5

#include <cassert>
#include <cstdlib>

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
#include <lifev/core/fem/TimeAndExtrapolationHandlerQuadPts.hpp>

using namespace LifeV;

typedef VectorEpetra vector_Type;
typedef RegionMesh<LinearTetra> mesh_Type;
typedef std::shared_ptr<mesh_Type> meshPtr_Type;
typedef FESpace<mesh_Type, MapEpetra >  FESpace_Type;
typedef std::shared_ptr<FESpace_Type> FESpacePtr_Type;

int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    bool verbose = Comm->MyPID() == 0;
#else
    std::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm() );
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

    const UInt sizeField = 3;
    FESpacePtr_Type velocityFESpace;
    velocityFESpace.reset (new FESpace_Type (fluidMeshLocal, data_file("mesh_fluid/fe_order","P1"), sizeField, Comm) );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  TESTING THE TIME HANDLER OBJECT
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(Comm->MyPID()==0)
    	std::cout << "[Test with empty constructor ]";

    UInt orderBDF = data_file("dataTimeHandler/orderBDF",2);
    UInt orderExtrapolation = data_file("dataTimeHandler/orderExtrapolation",2);

    // Initialize the time handler object, set the order of the BDF and for the extrapolation
    TimeAndExtrapolationHandlerQuadPts<sizeField> timeVelocity;
    timeVelocity.setBDForder(orderBDF);

    // Set the timestep
    timeVelocity.setTimeStep(0.2);

    std::vector<std::vector<VectorSmall<sizeField>>> firstVector;
    std::vector<std::vector<VectorSmall<sizeField>>> secondVector;
    std::vector<std::vector<VectorSmall<sizeField>>> thirdVector;
    
    firstVector.resize(fluidMeshLocal->numElements());
    secondVector.resize(fluidMeshLocal->numElements());
    thirdVector.resize(fluidMeshLocal->numElements());
    
    // Initialize 3 vectors
    for (int i = 0 ; i < fluidMeshLocal->numElements() ; ++i ) // loop elements
    {
        firstVector[i].resize(velocityFESpace->qr().nbQuadPt());
        secondVector[i].resize(velocityFESpace->qr().nbQuadPt());
        thirdVector[i].resize(velocityFESpace->qr().nbQuadPt());
        
        for (int j = 0 ; j < velocityFESpace->qr().nbQuadPt(); ++j ) // loop quadrature points
        {
            for (int k = 0 ; k < sizeField; ++k ) // loop dimensions
            {
                firstVector[i][j](k)  = 1.0;
                secondVector[i][j](k) = 2.0;
                thirdVector[i][j](k)  = 3.0;
            }
        }
    }

    // Initialize the state with two vectors (we suppose that in the datafile orderBDF=2)
    std::vector<std::vector<std::vector<VectorSmall<sizeField>>>> initialState;
    initialState.push_back(firstVector);
    initialState.push_back(secondVector);
    timeVelocity.initialize(initialState);

    if(Comm->MyPID()==0)
        std::cout << "\n";

    // Put in the stencil a new vector, and do a new extrapolation
    timeVelocity.shift(thirdVector);
    
    // Compute the term that should go to the right hand side
    std::vector<std::vector<VectorSmall<sizeField>>> rhsTerm;
    rhsTerm.resize(fluidMeshLocal->numElements());
    for (int i = 0 ; i < fluidMeshLocal->numElements() ; ++i ) // loop elements
    {
        rhsTerm[i].resize(velocityFESpace->qr().nbQuadPt());
    }
    
    timeVelocity.rhsContribution(rhsTerm);

    for ( int i = 0; i < fluidMeshLocal->numElements(); ++i )
    {
        for ( int j = 0; j < velocityFESpace->qr().nbQuadPt(); ++j )
        {
            std::cout << " (";
            for ( int k = 0; k < sizeField; ++k )
            {
                std::cout << " " << rhsTerm[i][j](k);
            }
            std::cout << "),";
        }
        
        std::cout << "\n";
    }
    
    
    // Print the value of alpha
    std::cout << "\nAlpha = " << timeVelocity.alpha() << "\n\n";
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return 0;

}
