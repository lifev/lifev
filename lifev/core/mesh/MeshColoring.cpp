#include <lifev/core/mesh/MeshColoring.hpp>


namespace LifeV
{

MeshColoring::MeshColoring(const commPtr_Type& communicator):
    M_comm(communicator)
{
}
    
MeshColoring::~MeshColoring()
{
}

void MeshColoring::setMesh ( const meshPtr_Type& local_mesh )
{
    M_mesh = local_mesh;
}

void MeshColoring::setup ( )
{
    int numVertices = M_mesh->numGlobalVertices();
    int numVolumes  = M_mesh->numGlobalElements();
    UInt numberLocalDof = 4;
 
    M_nodesToVolumes.resize(numVertices);
    M_volumesToNodes.resize(numVolumes);
    M_numNeighbors.resize(numVertices);
    
    for ( int i = 0; i < M_mesh->numElements(); ++i ) // local number of elements
    {
        M_volumesToNodes[M_mesh->element(i).id()].resize(4);
        for ( int j = 0; j < 4; ++j )
        {
            M_volumesToNodes[M_mesh->element(i).id()][j] = M_mesh->element(i).point(j).id();
            M_nodesToVolumes[M_mesh->element(i).point(j).id()].push_back(M_mesh->element(i).id());
        }
    }
    
    for ( int i = 0; i < numVertices; ++i )
    {
        M_numNeighbors[i] = M_nodesToVolumes[i].size();
    }
    
    std::vector<UInt>::iterator indexMaxNumNeighbors = std::max_element(M_numNeighbors.begin(), M_numNeighbors.end() );
    M_indexNodeMaxNumNeighbors = std::distance(M_numNeighbors.begin(), indexMaxNumNeighbors);
    M_maxNumNeighbors = M_numNeighbors[M_indexNodeMaxNumNeighbors];
    M_numColors = M_maxNumNeighbors;
    
    // debug info
    // std::cout << "\nNode with most elements connected is " << M_indexMaxNumNeighbors << " with " << M_maxNumNeighbors << " elements\n";
    
    // building element - neighbors map (element i has neighbors element j k z for example)
    M_volumesToVolumes.resize(numVolumes);
    for ( int i = 0; i < M_mesh->numElements(); ++i ) // local number of elements
    {
        for ( int j = 0; j < 4; ++j )
        {
            std::vector<UInt> volumes_tmp = M_nodesToVolumes[M_mesh->element(i).point(j).id()];
            for ( int k = 0; k < volumes_tmp.size(); ++k )
            {
                M_volumesToVolumes[M_mesh->element(i).id()].push_back( volumes_tmp[k] );
            }
        }
        std::sort( M_volumesToVolumes[M_mesh->element(i).id()].begin(), M_volumesToVolumes[M_mesh->element(i).id()].end() );
        M_volumesToVolumes[M_mesh->element(i).id()].erase( unique( M_volumesToVolumes[M_mesh->element(i).id()].begin(),
                                                                   M_volumesToVolumes[M_mesh->element(i).id()].end() ),
                                                                   M_volumesToVolumes[M_mesh->element(i).id()].end() );
        M_volumesToVolumes[M_mesh->element(i).id()].erase(std::remove(M_volumesToVolumes[M_mesh->element(i).id()].begin(),
                                                                      M_volumesToVolumes[M_mesh->element(i).id()].end(),
                                                                      M_mesh->element(i).id()),
                                                          M_volumesToVolumes[M_mesh->element(i).id()].end());
    }
    
    // build distributed map for exporter
    std::vector<int> id_elem_scalar;
    
    for ( int i = 0; i < M_mesh->numVolumes(); ++i )
    {
        id_elem_scalar.push_back ( M_mesh->element(i).id() );
    }
    
    int* pointerToDofs_scalar (0);
    pointerToDofs_scalar = &id_elem_scalar[0];
    boost::shared_ptr<MapEpetra> map_scalar ( new MapEpetra ( -1, static_cast<int> (id_elem_scalar.size() ), pointerToDofs_scalar, M_comm ) );
    
    M_vectorColors.reset ( new vector_Type ( *map_scalar,  Unique ) );
    M_vectorColors->zero();

}
    
void MeshColoring::colorMesh()
{
    M_colors.resize ( M_mesh->numGlobalElements() );
    std::for_each ( M_colors.begin(), M_colors.end(), [] (int&d) { d = -1.0;} );
    
    for ( int j = 0; j < M_numColors; ++j )
    {
        M_colors[ M_nodesToVolumes[M_indexNodeMaxNumNeighbors][j] ] = j;
    }
    
    std::vector<int> all_colors;
    all_colors.resize(M_numColors);
    std::iota (std::begin(all_colors), std::end(all_colors), 0);
    
    for ( int k = 0; k < M_mesh->numElements(); ++k ) // local number of elements
    {
        if ( M_colors[ M_mesh->element(k).id() ] == -1.0 )
        {
            std::vector<UInt> neigbors = M_volumesToVolumes[M_mesh->element(k).id()];
            std::vector<int> my_colors;
            my_colors.resize(neigbors.size());
            for ( int z = 0; z < neigbors.size(); ++z )
            {
                my_colors[z] = M_colors[neigbors[z]];
            }
            
            std::vector<int> vector_diff;
            std::sort (my_colors.begin(),my_colors.end());
            std::set_difference (all_colors.begin(), all_colors.end(),
                                 my_colors.begin(), my_colors.end(),
                                 std::inserter(vector_diff, vector_diff.begin()) );
            
            if ( vector_diff.size() == 0 )
            {
                M_numColors = M_numColors + 1;
                M_colors[M_mesh->element(k).id()] = M_numColors-1;
                all_colors.resize(M_numColors);
                std::iota (std::begin(all_colors), std::end(all_colors), 0);
            }
            else
            {
                std::vector<int>::iterator result = std::min_element(std::begin(vector_diff), std::end(vector_diff));
                M_colors[M_mesh->element(k).id()] = vector_diff[std::distance(std::begin(vector_diff), result)];
            }
        }
    }
    
    // fill epetra vector with colors
    for ( int i = 0; i < M_vectorColors->epetraVector().MyLength(); ++i )
    {
        (*M_vectorColors)[M_vectorColors->blockMap().GID(i)] = M_colors[M_vectorColors->blockMap().GID(i)];
    }
    colorsForAssembly ( );
}
 
void MeshColoring::colorsForAssembly ( )
{
    M_colorsForAssembly.resize(M_numColors);
    
    for ( int i = 0; i < M_vectorColors->epetraVector().MyLength(); ++i )
    {
        M_colorsForAssembly[ (*M_vectorColors)[M_vectorColors->blockMap().GID(i)] ].push_back(i);
    }
}
    
    
void MeshColoring::printVolumesToVolumes ( )
{
    for ( int i = 0; i < M_volumesToVolumes.size(); ++i )
    {
        std::cout << "\nElement " << i << " attached to elements " << std::flush;
        for ( int j = 0; j < M_volumesToVolumes[i].size(); ++j )
        {
            std::cout << M_volumesToVolumes[i][j] << " ";
        }
    }
    std::cout << "\n";
}
    
void MeshColoring::printNodesToVolumes ( )
{
    for ( int i = 0; i < M_nodesToVolumes.size(); ++i )
    {
        std::cout << "\nNode " << i << " attached to elements " << std::flush;
        for ( int j = 0; j < M_nodesToVolumes[i].size(); ++j )
        {
            std::cout << M_nodesToVolumes[i][j] << " ";
        }
    }
    std::cout << "\n";
}

void MeshColoring::printVolumesToNodes ( )
{
    for ( int i = 0; i < M_volumesToNodes.size(); ++i )
    {
        std::cout << "\nElement " << i << " has vertices " << std::flush;
        for ( int j = 0; j < M_volumesToNodes[i].size(); ++j )
        {
            std::cout << M_volumesToNodes[i][j] << " ";
        }
    }
    std::cout << "\n";
}
    
void MeshColoring::printNumNeighbors ( )
{
    for ( int i = 0; i < M_numNeighbors.size(); ++i )
    {
        std::cout << "\nVertex " << i << " touches " << M_numNeighbors[i] << " elements " << std::endl;
    }
}
    
void MeshColoring::printColors ( )
{
    for ( int i = 0; i < M_vectorColors->epetraVector().MyLength(); ++i )
    {
        std::cout << "\nElement " << M_vectorColors->blockMap().GID(i)
                  << " on proc ID " << M_comm->MyPID()
                  << " has color " << M_colors[M_vectorColors->blockMap().GID(i)]
                  << std::endl;
    }
    std::cout << "\n";
}
    
int MeshColoring::performCheck ( )
{
    int check = 0;
    
    for ( int k = 0; k < M_mesh->numVertices(); ++k ) // local number of vertices
    {
        std::vector<UInt> neigbors = M_nodesToVolumes[M_mesh->point(k).id()];
        std::vector<int> colors_neighbors;
        for ( int z = 0; z < neigbors.size(); ++z )
        {
            colors_neighbors.push_back(M_colors[neigbors[z]]);
        }
        // put in the colors_neighbors also the color of element k. Then check if size(unique(colors_neighbors))) is
        // equal to size colors_neighbors. If not there is an error.
        int size_initial = colors_neighbors.size();
        std::sort( colors_neighbors.begin(), colors_neighbors.end() );
        colors_neighbors.erase( unique( colors_neighbors.begin(), colors_neighbors.end() ), colors_neighbors.end() );
        int size_after = colors_neighbors.size();
        
        if ( size_after != size_initial )
        {
            check++;
            std::cout << "\nInitial = " << size_initial << ", after = " << size_after << std::endl;
            std::cout << "\nError in coloring around vertex " << M_mesh->point(k).id() << std::endl;
        }
    }
    return check;
}

void MeshColoring::printColorsForAssembly ( )
{
    for ( int i = 0; i < M_colorsForAssembly.size(); ++i )
    {
        if ( M_colorsForAssembly[i].size() > 0 )
        {
            std::cout << "\nColor " << i << " has elements with LID " << std::flush;
            for ( int j = 0; j < M_colorsForAssembly[i].size(); ++j )
            {
                std::cout << M_colorsForAssembly[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << "\n";
}
    
void MeshColoring::printInfo ( )
{
    // Number of colors per process
    std::cout << "\nOn proc " << M_comm->MyPID() << " there are " << M_numColors << " colors\n";
    
    M_comm->Barrier();
    
    for ( int k = 0; k < M_comm->NumProc(); ++k )
    {
        if ( M_comm->MyPID() == k )
        {
            for ( int i = 0; i < M_colorsForAssembly.size(); ++i )
            {
                if ( M_colorsForAssembly[i].size() > 0 )
                {
                    std::cout << "\nColor " << i << " has elements " << M_colorsForAssembly[i].size() << std::endl;
                }
            }
        }
    }
    
}
    
    
}