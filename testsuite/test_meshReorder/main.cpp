#ifdef TWODIM
#error test_reorder cannot be compiled in 2D
#endif


#include "Epetra_config.h"

#include <life/lifefilters/medit_wrtrs.hpp>

#include<life/lifemesh/dataMesh.hpp>

#include<life/lifecore/GetPot.hpp>
#include <string>

#include <life/lifefem/dofInterface3Dto3D.hpp>
#include <life/lifefem/FESpace.hpp>

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
    DataMesh<RegionMesh3D<LinearTetra> > mesh_data(data_file, "space_discretization");
    //DataMesh<RegionMesh3D<LinearTetra> > solidData(data_file, "solid/space_discretization");
    //const char* mesh_input = command_line.follow(data_file("fluid/space_discretization/mesh_file", "mesh", 0), 2, "-i","--input");

    const std::string mesh_output = command_line.follow((data_file("space_discretization/output_mesh_file", "mesh").c_str()), 2, "-o", "--output");
    bool ordering = data_file("space_discretization/ordering", false);
    bool create_edge = data_file("interface/create_edge", true);
    UInt newMarker = data_file("interface/edgeMarker", 20);

    if(ordering)
        {
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
        }
    else
        if(create_edge)
            {
                UInt FluidInterfaceFlag           = data_file("interface/fluidInterfaceFlag",      2 );
                UInt SolidInterfaceFlag           = data_file("interface/solidInterfaceFlag",      2 );
                EntityFlag edgeFlag               = data_file("interface/edgeFlag",      2 );

                boost::shared_ptr<Epetra_Comm> uselessComm(new Epetra_MpiComm(MPI_COMM_WORLD));
                DataMesh<RegionMesh3D<LinearTetra> > mesh_data2(data_file, "second_mesh/space_discretization");

                boost::shared_ptr<FESpace<RegionMesh3D<LinearTetra>, EpetraMap> > firstFESpace;
                firstFESpace.reset(new FESpace<RegionMesh3D<LinearTetra>, EpetraMap>(mesh_data.mesh(),"P1",3,*uselessComm));

                boost::shared_ptr<FESpace<RegionMesh3D<LinearTetra>, EpetraMap> > secondFESpace;
                secondFESpace.reset(new FESpace<RegionMesh3D<LinearTetra>, EpetraMap>(mesh_data2.mesh(),"P1",3,*uselessComm));


                boost::shared_ptr<DofInterface3Dto3D>  dofEdgeFluidToEdgeSolid( new DofInterface3Dto3D );

                dofEdgeFluidToEdgeSolid->setup(firstFESpace->refFE(), firstFESpace->dof(),
                                               secondFESpace->refFE(), secondFESpace->dof()
                                               );

                dofEdgeFluidToEdgeSolid->update(
                                                *mesh_data.mesh(), FluidInterfaceFlag,
                                                *mesh_data2.mesh(), SolidInterfaceFlag,
                                                0., &edgeFlag);
                mesh_data2.mesh()->edgeMarkers(dofEdgeFluidToEdgeSolid->locDofMap(), newMarker);
                writeMesh( mesh_output , *mesh_data2.mesh());
            }
    MPI_Finalize();
#endif
}
