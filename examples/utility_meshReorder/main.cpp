#ifdef TWODIM
#error test_reorder cannot be compiled in 2D
#endif


// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include<life/lifemesh/MeshData.hpp>

#include<life/lifefilters/GetPot.hpp>
#include <life/lifefilters/MeshWriter.hpp>
#include <string>

#include <life/lifefem/DOFInterface3Dto3D.hpp>
#include <life/lifefem/FESpace.hpp>


int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    using namespace LifeV;
    GetPot command_line(argc,argv);
    std::string data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);

    MeshData mesh_data;
    mesh_data.setup(data_file, "space_discretization");

    boost::shared_ptr<RegionMesh3D<LinearTetra> > mesh;
    mesh.reset(new RegionMesh3D<LinearTetra>);

    readMesh(*mesh, mesh_data);

    //MeshData<RegionMesh3D<LinearTetra> > solidData(data_file, "solid/space_discretization");
    //const char* mesh_input = command_line.follow(data_file("fluid/space_discretization/mesh_file", "mesh", 0), 2, "-i","--input");

    const std::string mesh_output = command_line.follow((data_file("space_discretization/output_mesh_file", "mesh").c_str()), 2, "-o", "--output");
    bool ordering = data_file("space_discretization/ordering", false);
    bool create_edge = data_file("interface/create_edge", true);
    UInt TimeAdvanceNewmarker = data_file("interface/edgeMarker", 20);

    if (ordering)
    {
        int numtasks;
        int me;
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
        MPI_Comm MPIcomm;
        MPI_Group  originGroup, newGroup;
        MPI_Comm_group(MPI_COMM_WORLD, &originGroup);
        int members[numtasks];
        for (Int i=0; i<numtasks; ++i)
            members[i]=0;
        MPI_Group_incl(originGroup, 1, members, &newGroup);

        MPI_Comm_create(MPI_COMM_WORLD, newGroup, &MPIcomm);
        if(me==0)
        {
            mesh->orderMesh( MPIcomm);
            //solidData.mesh()->orderMesh( MPIcomm);

            MeshWriter::writeMeshMedit<RegionMesh3D<LinearTetra> >( mesh_output , *mesh);
            //writeMesh( "solid_ord.mesh", *solidData.mesh());
        }

    }
    else if (create_edge)
    {
        UInt FluidInterfaceFlag           = data_file("interface/fluidInterfaceFlag",      2 );
        UInt SolidInterfaceFlag           = data_file("interface/solidInterfaceFlag",      2 );
        int const edgeFlag                 (data_file("interface/edgeFlag",      2 ) );

        boost::shared_ptr<Epetra_Comm> uselessComm(new Epetra_MpiComm(MPI_COMM_WORLD));
        MeshData mesh_data2;
        mesh_data2.setup(data_file, "second_mesh/space_discretization");

        boost::shared_ptr<RegionMesh3D<LinearTetra> > mesh2;
        mesh2.reset(new RegionMesh3D<LinearTetra>);

        readMesh(*mesh2, mesh_data2);

        boost::shared_ptr<FESpace<RegionMesh3D<LinearTetra>, MapEpetra> > firstFESpace;
        firstFESpace.reset(new FESpace<RegionMesh3D<LinearTetra>, MapEpetra> (mesh,  "P1", 3, uselessComm));

        boost::shared_ptr<FESpace<RegionMesh3D<LinearTetra>, MapEpetra> > secondFESpace;
        secondFESpace.reset(new FESpace<RegionMesh3D<LinearTetra>, MapEpetra>(mesh2, "P1", 3, uselessComm));


        boost::shared_ptr<DOFInterface3Dto3D>  dofEdgeFluidToEdgeSolid( new DOFInterface3Dto3D );

        dofEdgeFluidToEdgeSolid->setup(firstFESpace->refFE(), firstFESpace->dof(),
                                       secondFESpace->refFE(), secondFESpace->dof()
                                      );

        dofEdgeFluidToEdgeSolid->update(
            *mesh, FluidInterfaceFlag,
            *mesh, SolidInterfaceFlag,
            0., &edgeFlag);
        mesh2->edgeMarkers(dofEdgeFluidToEdgeSolid->localDofMap(), TimeAdvanceNewmarker);
        MeshWriter::writeMeshMedit<RegionMesh3D<LinearTetra> >( mesh_output , *mesh2);

    }
    MPI_Finalize();
#endif
}
