
    // get list of ghost gid from partitioner
    typename MeshPartitioner<RegionMesh>::GhostEntityDataMap_Type ghostData = meshPart.ghostDataMap();
    typename MeshPartitioner<RegionMesh>::GhostEntityDataMap_Type::const_iterator pIt  = ghostData.begin();
    typename MeshPartitioner<RegionMesh>::GhostEntityDataMap_Type::const_iterator pEnd = ghostData.end();
    typename MeshPartitioner<RegionMesh>::GhostEntityDataContainer_Type::const_iterator dIt;

    // create unique list of element gid to search
    std::set<int> uniqueGhostList;
    std::vector<ID> multiGhostVector;
    for ( ; pIt != pEnd; ++pIt )
    {
        for ( dIt = pIt->second.begin(); dIt != pIt->second.end(); ++dIt )
        {
            uniqueGhostList.insert ( dIt->ghostElementGlobalId );
            multiGhostVector.push_back( dIt->ghostElementGlobalId );
        }
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<int> uniqueGhostVector ( uniqueGhostList.begin(), uniqueGhostList.end() );

    // get map and create directory
    const Epetra_BlockMap & map ( *( feSpacePtr->map().map( Unique ) ) );
    Epetra_Directory* directory = feSpacePtr->map().comm().CreateDirectory( map );

    // ask directory for proc and lid of ghost elements
    std::vector<int> uniqueGhostProc ( uniqueGhostVector.size() );
    std::vector<int> lIDList ( ghostList.size() );
    directory->GetDirectoryEntries( map,
                                    uniqueGhostVector.size(),
                                    &uniqueGhostVector[ 0 ],
                                    &uniqueGhostProc[ 0 ],
                                    &lIDlist[ 0 ],
                                    NULL );

    std::ofstream fout ( ( "ghostfaces." + boost::lexical_cast<std::string>( Members->comm->MyPID() ) + ".out").c_str() );

    // associate ghost gid with proc holding them
    std::map<ID,UInt> ghostProc;
    for ( UInt i = 0; i < uniqueGhostVector.size(); i++ )
    {
        ghostProc[ uniqueGhostVector[ i ] ] = uniqueGhostProc[ i ];
    }

    // DEBUG
//    for ( UInt i = 0; i < multiGhostVector.size(); i++ )
//    {
//        fout << multiGhostVector[ i ] << " " << ghostProc[ multiGhostVector[ i ] ] << std::endl;
//    }

    // reverse map to associate to each proc the list of ghost gid that it is holding
    std::map<UInt, std::vector<ID> > ghost;
    for ( std::map<ID,UInt>::const_iterator it = ghostProc.begin(); it != ghostProc.end(); ++it )
    {
        ghost[ it->second ].push_back( it->first );
    }

    Members->comm->Barrier();
    abort();

