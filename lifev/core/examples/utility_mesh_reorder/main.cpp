#ifdef TWODIM
#error test_reorder cannot be compiled in 2D
#endif


#include <lifev/core/mesh/InternalEntitySelector.hpp>
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/MeshWriter.hpp>
#include <string>

#include <lifev/core/fem/DOFInterface3Dto3D.hpp>
#include <lifev/core/fem/FESpace.hpp>

/*!
 * @file
 * @brief File containing a utility to reorder a mesh
 *
 * This utility contains a method to reorder the meshes reducing the fill-in of the matrix,
 * and a method to create a different marker on a region of a mesh.
 * In particular this region is identified by a marker of another mesh that shares an interface
 * (e.g. solid interface and fluid edges in FSI).
 *
 * For instance suppose that you have marker=20 on the fluid 'edges'
 * and on the solid mesh you have just different markers for the internal and external walls,
 * then using this test you can create a marker '30' on the solid mesh whenever the corresponding marker on the fluid
 * mesh is '20', and rewrite the mesh.
 *
 * usage: test_meshReorder -o mesh_output.mesh
 *
 * without the options -o or --output the output mesh can be specified from this data file.
 * Otherwise the mesh will be written in a file 'mesh'.
 *
 * the first mesh is the one reordered. The second mesh is the one where the new marker is added
 * (not read when create_edge=false).
 *
 * @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *
 * @mantainer Paolo Crosetto <paolo.crosetto@epfl.ch>
 */


int main (int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    using namespace LifeV;
    GetPot command_line (argc, argv);
    std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot data_file (data_file_name);

    MeshData mesh_data;
    mesh_data.setup (data_file, "space_discretization");

    std::shared_ptr<Epetra_Comm> uselessComm (new Epetra_MpiComm (MPI_COMM_WORLD) );
    std::shared_ptr<RegionMesh<LinearTetra> > mesh (new RegionMesh<LinearTetra> ( uselessComm) );

    readMesh (*mesh, mesh_data);

    //MeshData<RegionMesh<LinearTetra> > solidData(data_file, "solid/space_discretization");
    //const char* mesh_input = command_line.follow(data_file("fluid/space_discretization/mesh_file", "mesh", 0), 2, "-i","--input");

    const std::string mesh_output = command_line.follow ( (data_file ("space_discretization/output_mesh_file", "mesh").c_str() ), 2, "-o", "--output");
    bool ordering = data_file ("space_discretization/ordering", false);
    bool create_edge = data_file ("interface/create_edge", true);
    UInt TimeAdvanceNewmarker = data_file ("interface/edgeMarker", 20);

    if (ordering)
    {
        Int numtasks;
        Int me;
        MPI_Comm_rank (MPI_COMM_WORLD, &me);
        MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
        MPI_Comm MPIcomm;
        MPI_Group  originGroup, newGroup;
        MPI_Comm_group (MPI_COMM_WORLD, &originGroup);
        std::vector<Int> members (numtasks);
        for (Int i = 0; i < numtasks; ++i)
        {
            members[i] = 0;
        }
        MPI_Group_incl (originGroup, 1, &members[0], &newGroup);

        MPI_Comm_create (MPI_COMM_WORLD, newGroup, &MPIcomm);
        if (me == 0)
        {
            // LF TAKEN AWAY. IT IS BROKEN
            // mesh->orderMesh( MPIcomm);
            //solidData.mesh()->orderMesh( MPIcomm);

            MeshWriter::writeMeshMedit<RegionMesh<LinearTetra> > ( mesh_output , *mesh);
            //writeMesh( "solid_ord.mesh", *solidData.mesh());
        }

    }
    else if (create_edge)
    {
        UInt FluidInterfaceFlag           = data_file ("interface/fluidInterfaceFlag",      2 );
        UInt SolidInterfaceFlag           = data_file ("interface/solidInterfaceFlag",      2 );
        Int const edgeFlag                 (data_file ("interface/edgeFlag",      2 ) );

        MeshData mesh_data2;
        mesh_data2.setup (data_file, "second_mesh/space_discretization");

        std::shared_ptr<RegionMesh<LinearTetra> > mesh2;
        mesh2.reset (new RegionMesh<LinearTetra> (uselessComm) );

        readMesh (*mesh2, mesh_data2);

        std::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra> > firstFESpace;
        firstFESpace.reset (new FESpace<RegionMesh<LinearTetra>, MapEpetra> (mesh,  "P1", 3, uselessComm) );

        std::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra> > secondFESpace;
        secondFESpace.reset (new FESpace<RegionMesh<LinearTetra>, MapEpetra> (mesh2, "P1", 3, uselessComm) );


        std::shared_ptr<DOFInterface3Dto3D>  dofEdgeFluidToEdgeSolid ( new DOFInterface3Dto3D );

        dofEdgeFluidToEdgeSolid->setup (firstFESpace->refFE(), firstFESpace->dof(),
                                        secondFESpace->refFE(), secondFESpace->dof()
                                       );

        dofEdgeFluidToEdgeSolid->update (
            *mesh, FluidInterfaceFlag,
            *mesh, SolidInterfaceFlag,
            0., &edgeFlag);

        ChangeMarkersAccordingToMap (mesh2->pointList,
                                     dofEdgeFluidToEdgeSolid->localDofMap(),
                                     TimeAdvanceNewmarker);
        //mesh2->edgeMarkers(dofEdgeFluidToEdgeSolid->localDofMap(), TimeAdvanceNewmarker);
        MeshWriter::writeMeshMedit<RegionMesh<LinearTetra> > ( mesh_output , *mesh2);

    }
    MPI_Finalize();
#endif
}
