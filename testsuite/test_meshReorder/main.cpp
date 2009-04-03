#include "Epetra_config.h"
#include <life/lifefilters/medit_wrtrs.hpp>
#include<life/lifemesh/dataMesh.hpp>
#include<life/lifecore/GetPot.hpp>
#include <string>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    using namespace LifeV;
    GetPot command_line(argc,argv);
    std::string data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);
    DataMesh<RegionMesh3D<LinearTetra> > mesh_data(data_file, "discretization");
    //DataMesh<RegionMesh3D<LinearTetra> > solidData(data_file, "solid/discretization");
    //const char* mesh_input = command_line.follow(data_file("fluid/discretization/mesh_file", "mesh", 0), 2, "-i","--input");
    std::string mesh_output = command_line.follow(data_file("fluid/discretization/output_mesh_file", "mesh", 0).c_str(), 2, "-o","--output");
    int numtasks;
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm MPIcomm;
    MPI_Group  originGroup, newGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &originGroup);
    int members[numtasks];
    for(ID i=0; i<numtasks; ++i)
        members[i]=0;
    int ierr = MPI_Group_incl(originGroup, 1, members, &newGroup);

    MPI_Comm_create(MPI_COMM_WORLD, newGroup, &MPIcomm);
    if(me==0)
        {
            mesh_data.mesh()->orderMesh( MPIcomm);
            //solidData.mesh()->orderMesh( MPIcomm);
            writeMesh( mesh_output , *mesh_data.mesh());
            //writeMesh( "solid_ord.mesh", *solidData.mesh());
        }
    MPI_Finalize();
#endif
}
