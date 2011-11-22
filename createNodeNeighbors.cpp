
//! this routine generates node neighbors for the given mesh
/*! the routine assumes that the mesh is not yet partitioned or reordered 
 *  (i.e. the local id and the global id are the same).
 *  if this is not true the method should be changed to use a less 
 *  cheap STL find on the mesh points to get the correct point that has 
 *  the given global id. 
 */
template <typename MeshType>
void createNodeNeighbors( MeshType const & mesh )
{
    // TODO: ASSERT_COMPILE_TIME that MeshType::pointMarker == NeighborMarker
    // this guarantees that the nodeNeighbors structure is available. 
    
    // generate node neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < mesh.numEdges(); ie++ )
    {
        ID id0 = mesh.edge( ie ).point( 0 ).id();
        ID id1 = mesh.edge( ie ).point( 1 ).id();

        ASSERT ( mesh.point( id0 ).id() == id0 , "the mesh has been reordered, the point must be found" );
        ASSERT ( mesh.point( id1 ).id() == id1 , "the mesh has been reordered, the point must be found" );

        mesh.point( id0 ).nodeNeighbors().insert( id1 );
        mesh.point( id1 ).nodeNeighbors().insert( id0 );
    }
}   

