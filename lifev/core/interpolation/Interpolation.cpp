#include <lifev/core/interpolation/Interpolation.hpp>

namespace LifeV
{

Interpolation::Interpolation():
		M_precBuilt ( false )
{
}

Interpolation::~Interpolation()
{}

void
Interpolation::setup( const GetPot& datafile, parameterList_Type belosList )
{
	M_datafile = datafile;
	M_belosList = belosList;
}

double
Interpolation::rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius)
{
    double distance = std::sqrt ( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2) );
    return (1 - distance / radius) * (1 - distance / radius) * (1 - distance / radius) * (1 - distance / radius) * (4 * distance / radius + 1);
}

void
Interpolation::setVectors (const vectorPtr_Type& KnownField, const vectorPtr_Type& UnknownField )
{
	M_knownField.reset( new vector_Type ( *KnownField ) );
	M_unknownField.reset( new vector_Type ( *UnknownField ) );
}

void
Interpolation::buildOperators()
{
	interpolationOperator();
	projectionOperator();
	buildRhs();
	interpolateCostantField();
	free_space ();
}

void
Interpolation::projectionOperator()
{
	buildProjectionOperatorMap();

	// Total number of dofs taken into account
	int LocalNodesNumber = M_GIdsUnknownMesh.size();

	std::vector<std::vector<int>> MatrixGraph (LocalNodesNumber);
	int* GlobalID = new int[LocalNodesNumber];
	int k = 0;

	// I need to find the closest point in the "known mesh" to use its radius
	double d;
	double d_min;
	int nearestPoint;

	for (std::set<ID>::iterator it = M_GIdsUnknownMesh.begin(); it != M_GIdsUnknownMesh.end(); ++it)
	{
		GlobalID[k] = *it;
		d_min = 100;
		for (int j = 0; j <  M_xcoord_known.size(); ++j)
		{
			if ( M_marker_known[j] == M_flag )
			{
				d = std::sqrt ( ( M_xcoord_known[j] - M_xcoord_unknown[*it] ) * ( M_xcoord_known[j] - M_xcoord_unknown[*it] ) +
								( M_ycoord_known[j] - M_ycoord_unknown[*it] ) * ( M_ycoord_known[j] - M_ycoord_unknown[*it] ) +
								( M_zcoord_known[j] - M_zcoord_unknown[*it] ) * ( M_zcoord_known[j] - M_zcoord_unknown[*it] ) );
				if (d < d_min)
				{
					d_min = d;
					nearestPoint = j;
				}
			}
		}

		// For each of them, identify the neighbors on the other mesh within a certain number of circles M_links
		MatrixGraph[k] = M_dof_connectivity_known[nearestPoint];
		// RBF_radius[k] = computeRBFradius_unknown ( *it, MatrixGraph[k]);
		++k;
	}
	M_projectionOperator.reset (new matrix_Type (*M_projectionOperatorMap, 100) );

	Real Values;

	for ( int i = 0 ; i < LocalNodesNumber; ++i )
	{
		for ( int it = 0; it < MatrixGraph[i].size(); ++it )
		{
			Values  = rbf ( M_xcoord_unknown[GlobalID[i]],
							M_ycoord_unknown[GlobalID[i]],
							M_zcoord_unknown[GlobalID[i]],
							M_xcoord_known[MatrixGraph[i][it]],
							M_ycoord_known[MatrixGraph[i][it]],
							M_zcoord_known[MatrixGraph[i][it]],
							M_links /*RBF_radius[i]*/ );

			M_projectionOperator->addToCoefficient( (*M_numerationInterfaceUnknown)( GlobalID[i] ),
													(*(M_numerationInterfaceKnownColumns[M_pid]))(MatrixGraph[i][it]),
													Values );
		}
	}
	M_projectionOperator->globalAssemble (M_interpolationOperatorMap, M_projectionOperatorMap);
	delete[] GlobalID;
}

double
Interpolation::computeRBFradius_known ( const ID& index, std::vector<int> indexes )
{
	double r = 0;
	double r_max = 0;
	for ( int i = 0; i < indexes.size(); ++i )
	{
		r = std::sqrt (  ( M_xcoord_known[index] - M_xcoord_known[indexes[i]] ) * ( M_xcoord_known[index] - M_xcoord_known[indexes[i]] ) +
						 ( M_ycoord_known[index] - M_ycoord_known[indexes[i]] ) * ( M_ycoord_known[index] - M_ycoord_known[indexes[i]] ) +
						 ( M_zcoord_known[index] - M_zcoord_known[indexes[i]] ) * ( M_zcoord_known[index] - M_zcoord_known[indexes[i]] ) );
		r_max = ( r > r_max ) ? r : r_max;
	}
	return r_max;
}

double
Interpolation::computeRBFradius_unknown ( const ID& index, std::vector<int> indexes )
{
	double r = 0;
	double r_max = 0;
	for ( int i = 0; i < indexes.size(); ++i )
	{
		r = std::sqrt (  ( M_xcoord_unknown[index] - M_xcoord_known[indexes[i]] ) * ( M_xcoord_unknown[index] - M_xcoord_known[indexes[i]] ) +
						 ( M_ycoord_unknown[index] - M_ycoord_known[indexes[i]] ) * ( M_ycoord_unknown[index] - M_ycoord_known[indexes[i]] ) +
						 ( M_zcoord_unknown[index] - M_zcoord_known[indexes[i]] ) * ( M_zcoord_unknown[index] - M_zcoord_known[indexes[i]] ) );
		r_max = ( r > r_max ) ? r : r_max;
	}
	return r_max;
}

void
Interpolation::interpolationOperator()
{
	buildInterpolationOperatorMap();

	int LocalNodesNumber = M_GIdsKnownMesh.size();

	std::vector<double> RBF_radius (LocalNodesNumber);
	std::vector<std::vector<int>> MatrixGraph (LocalNodesNumber);
	int* GlobalID = new int[LocalNodesNumber];

	int k = 0;
	int Max_entries = 0;

	for (std::set<ID>::iterator it = M_GIdsKnownMesh.begin(); it != M_GIdsKnownMesh.end(); ++it)
	{
		GlobalID[k] = *it;
		MatrixGraph[k] = M_dof_connectivity_known[*it];
		//RBF_radius[k] = computeRBFradius_known ( *it, MatrixGraph[k] );
		// std::cout << "DOF " << *it << " has radius " << RBF_radius[k] << std::endl;
		++k;
	}

	// Matrix for interpolation operator
	M_interpolationOperator.reset (new matrix_Type (*M_interpolationOperatorMap, 100) );

	Real Values;

	M_interpolationOperator->zero();

	// Loop over all the dofs taken into account by the interpolation
	for ( int i = 0 ; i < LocalNodesNumber; ++i )
	{
		// For each i-th dof, evaluate the rbf between the dof and its neighbors
		for ( int it = 0; it < MatrixGraph[i].size(); ++it)
		{
			Values = rbf (  M_xcoord_known[GlobalID[i]],
							M_ycoord_known[GlobalID[i]],
							M_zcoord_known[GlobalID[i]],
							M_xcoord_known[MatrixGraph[i][it]],
							M_ycoord_known[MatrixGraph[i][it]],
							M_zcoord_known[MatrixGraph[i][it]],
							M_links /*RBF_radius[i]*/ );

			M_interpolationOperator->addToCoefficient( (*M_numerationInterfaceKnown)( GlobalID[i] ),
													   (*(M_numerationInterfaceKnownColumns[M_pid]))(MatrixGraph[i][it]),
													   Values );
		}
	}

	M_interpolationOperator->globalAssemble();

	delete[] GlobalID;
}


void
Interpolation::buildTableDofs_known ( const FESpacePtr_Type& fespace )
{
	CurrentFEManifold feBd1 ( fespace->refFE().boundaryFE(), getGeometricMap ( *fespace->mesh() ).boundaryMap() );

	int numBFaces  = fespace->mesh()->numBFaces();
	int numTotalDof = fespace->dof().numTotalDof();

	std::vector< std::vector<int> > dof_element_connectivity(numTotalDof);

	M_marker_known.resize(numTotalDof);
	M_xcoord_known.resize(numTotalDof);
	M_ycoord_known.resize(numTotalDof);
	M_zcoord_known.resize(numTotalDof);

	for ( int i = 0; i < numBFaces; ++i )
	{
		if ( fespace->mesh()->boundaryFacet ( i ).markerID() == M_flag )
		{
			feBd1.update ( fespace->mesh()->boundaryFacet ( i ), UPDATE_ONLY_CELL_NODES );  // Updating facet information on mesh1
			std::vector<ID> localToGlobalMapOnBFacet1 = fespace->dof().localToGlobalMapOnBdFacet (i);

			// devo riempire vettore con marker di ogni dof
			VectorSmall<6> FaceDof;
			FaceDof *= 0;

			for (ID lDof1 = 0; lDof1 < localToGlobalMapOnBFacet1.size(); lDof1++)
			{
				dof_element_connectivity[localToGlobalMapOnBFacet1[lDof1]].push_back(i);

				feBd1.coorMap ( M_xcoord_known[localToGlobalMapOnBFacet1[lDof1]],
								M_ycoord_known[localToGlobalMapOnBFacet1[lDof1]],
								M_zcoord_known[localToGlobalMapOnBFacet1[lDof1]],
								feBd1.refFE().xi ( lDof1 ),
								feBd1.refFE().eta ( lDof1 ) );

				if ( lDof1 < 3 )
					FaceDof[lDof1] = fespace->mesh()->boundaryFacet ( i ).point ( lDof1 ).markerID ();

			}

			if ( localToGlobalMapOnBFacet1.size() > 3 )
			{
				// a p2 dof that was not in the mesh, so a dof the we inserted, has flag 1439 (random number) or it inherits the flag of the edge
				FaceDof[3] = ( FaceDof[0] == FaceDof[1] ) ? FaceDof[0] : M_flag;
				FaceDof[4] = ( FaceDof[1] == FaceDof[2] ) ? FaceDof[1] : M_flag;
				FaceDof[5] = ( FaceDof[0] == FaceDof[2] ) ? FaceDof[0] : M_flag;
			}

			for (ID lDof1 = 0; lDof1 < localToGlobalMapOnBFacet1.size(); lDof1++)
			{
				M_marker_known[localToGlobalMapOnBFacet1[lDof1]] = FaceDof[lDof1];
				if ( FaceDof[lDof1] == M_flag )
					M_GIdsKnownMesh_common.insert(localToGlobalMapOnBFacet1[lDof1]);
			}
		}
	}

	M_dof_connectivity_known.resize(numTotalDof);

    for ( int i = 0; i < numTotalDof; ++i )
    {
        if ( M_marker_known[i] == M_flag)
        {
            for ( int j = 0; j < numTotalDof; ++j )
            {
                if ( M_marker_known[j] == M_flag)
                {
                    Real distance = std::sqrt( (M_xcoord_known[i] - M_xcoord_known[j])*(M_xcoord_known[i] - M_xcoord_known[j]) +
                                               (M_ycoord_known[i] - M_ycoord_known[j])*(M_ycoord_known[i] - M_ycoord_known[j]) +
                                               (M_zcoord_known[i] - M_zcoord_known[j])*(M_zcoord_known[i] - M_zcoord_known[j]) );
                    if ( distance < M_links )
                    {
                        M_dof_connectivity_known[i].push_back(j);
                    }
                }
            }
        }
    }

    /*
    for ( int i = 0; i < numTotalDof; ++i )
    {
        if ( M_marker_known[i] == M_flag)
        {
            std::cout << "GID " << i;
            for ( int k = 0; k < M_dof_connectivity_known[i].size(); ++k )
            {
                std::cout << "  " << M_dof_connectivity_known[i][k];
            }
            std::cout << "\n";
        }
    }
    */

    /*
	for ( int i = 0; i < numTotalDof; ++i )
	{
		for ( int k = 0; k < dof_element_connectivity[i].size(); ++k )
		{
			feBd1.update ( fespace->mesh()->boundaryFacet ( dof_element_connectivity[i][k] ), UPDATE_ONLY_CELL_NODES );
			std::vector<ID> localToGlobalMapOnBFacet1 = fespace->dof().localToGlobalMapOnBdFacet ( dof_element_connectivity[i][k] );
			for (ID lDof1 = 0; lDof1 < localToGlobalMapOnBFacet1.size(); lDof1++)
			{
				M_dof_connectivity_known[i].push_back(localToGlobalMapOnBFacet1[lDof1]);
			}
		}
		std::sort( M_dof_connectivity_known[i].begin(), M_dof_connectivity_known[i].end() );
		M_dof_connectivity_known[i].erase( unique( M_dof_connectivity_known[i].begin(), M_dof_connectivity_known[i].end() ), M_dof_connectivity_known[i].end() );
	}
    */


	/*
	// stampa coordinate
	if ( fespace->map().commPtr()->MyPID() == 0 )
		for ( int i = 0 ; i < M_xcoord_known.size(); ++i )
			std::cout << "GID " << i << ", x = " << M_xcoord_known[i] << ", y = " << M_ycoord_known[i] << ", z = " << M_zcoord_known[i] << std::endl;
	*/

	/*
	// stampa connettivita tra dof
	if ( fespace->map().commPtr()->MyPID() == 0 )
	{
		for ( int i = 0 ; i < M_dof_connectivity_known.size(); ++i )
		{
			std::cout << "GID " << i << " ha " << M_dof_connectivity_known[i].size() << " vicini\n";
			for ( int j = 0 ; j < M_dof_connectivity_known[i].size(); ++j )
			{
				std::cout << " " << M_dof_connectivity_known[i][j];
			}
			std::cout << "\n";
		}
	}
	*/

	/*
	// stampa flag di ogni dof
	if ( fespace->map().commPtr()->MyPID() == 0 )
	{
		for ( int i = 0 ; i < M_marker_known.size(); ++i )
			std::cout << "DOF " << i << " has flag " << M_marker_known[i] << std::endl;
	}
	*/
}

void
Interpolation::free_space ( )
{
	M_marker_unknown.clear();
	M_xcoord_unknown.clear();
	M_ycoord_unknown.clear();
	M_zcoord_unknown.clear();
	M_marker_known.clear();
	M_xcoord_known.clear();
	M_ycoord_known.clear();
	M_zcoord_known.clear();
	M_GIdsUnknownMesh_common.clear();
	M_GIdsKnownMesh_common.clear();
	M_GIdsKnownMesh.clear();
	M_GIdsUnknownMesh.clear();
}

void
Interpolation::buildTableDofs_unknown ( const FESpacePtr_Type& fespace )
{
	CurrentFEManifold feBd1 ( fespace->refFE().boundaryFE(), getGeometricMap ( *fespace->mesh() ).boundaryMap() );

	int numBFaces  = fespace->mesh()->numBFaces();
	int numTotalDof = fespace->dof().numTotalDof();

	std::vector< std::vector<int> > dof_element_connectivity(numTotalDof);

	M_marker_unknown.resize(numTotalDof);
	M_xcoord_unknown.resize(numTotalDof);
	M_ycoord_unknown.resize(numTotalDof);
	M_zcoord_unknown.resize(numTotalDof);

	for ( int i = 0; i < numBFaces; ++i )
	{
		if ( fespace->mesh()->boundaryFacet ( i ).markerID() == M_flag )
		{
			feBd1.update ( fespace->mesh()->boundaryFacet ( i ), UPDATE_ONLY_CELL_NODES );  // Updating facet information on mesh1
			std::vector<ID> localToGlobalMapOnBFacet1 = fespace->dof().localToGlobalMapOnBdFacet (i);

			// devo riempire vettore con marker di ogni dof
			VectorSmall<6> FaceDof;
			FaceDof *= 0;

			for (ID lDof1 = 0; lDof1 < localToGlobalMapOnBFacet1.size(); lDof1++)
			{
				dof_element_connectivity[localToGlobalMapOnBFacet1[lDof1]].push_back(i);

				feBd1.coorMap ( M_xcoord_unknown[localToGlobalMapOnBFacet1[lDof1]],
								M_ycoord_unknown[localToGlobalMapOnBFacet1[lDof1]],
								M_zcoord_unknown[localToGlobalMapOnBFacet1[lDof1]],
								feBd1.refFE().xi ( lDof1 ),
								feBd1.refFE().eta ( lDof1 ) );

				if ( lDof1 < 3 )
					FaceDof[lDof1] = fespace->mesh()->boundaryFacet ( i ).point ( lDof1 ).markerID ();
			}

			if ( localToGlobalMapOnBFacet1.size() > 3 )
			{
				// a p2 dof that was not in the mesh, so a dof the we inserted, has flag 1439 (random number) or it inherits the flag of the edge
				FaceDof[3] = ( FaceDof[0] == FaceDof[1] ) ? FaceDof[0] : M_flag;
				FaceDof[4] = ( FaceDof[1] == FaceDof[2] ) ? FaceDof[1] : M_flag;
				FaceDof[5] = ( FaceDof[0] == FaceDof[2] ) ? FaceDof[0] : M_flag;
			}

			for (ID lDof1 = 0; lDof1 < localToGlobalMapOnBFacet1.size(); lDof1++)
			{
				M_marker_unknown[localToGlobalMapOnBFacet1[lDof1]] = FaceDof[lDof1];
				if ( FaceDof[lDof1] == M_flag )
					M_GIdsUnknownMesh_common.insert(localToGlobalMapOnBFacet1[lDof1]);
			}
		}
	}

	/*
	if ( fespace->map().commPtr()->MyPID() == 0 )
		for ( int i = 0 ; i < M_xcoord_unknown.size(); ++i )
			std::cout << "GID " << i << ", x = " << M_xcoord_unknown[i] << ", y = " << M_ycoord_unknown[i] << ", z = " << M_zcoord_unknown[i] << std::endl;
	*/

	/*
	// stampa connettivita tra dof
	if ( fespace->map().commPtr()->MyPID() == 0 )
	{
		for ( int i = 0 ; i < M_dof_connectivity_unknown.size(); ++i )
		{
			std::cout << "GID " << i << " ha " << M_dof_connectivity_unknown[i].size() << " vicini\n";
			for ( int j = 0 ; j < M_dof_connectivity_unknown[i].size(); ++j )
			{
				std::cout << " " << M_dof_connectivity_unknown[i][j];
			}
			std::cout << "\n";
		}
	}
	*/
	/*
	// stampa flag di ogni dof
	if ( fespace->map().commPtr()->MyPID() == 0 )
	{
		for ( int i = 0 ; i < M_marker_unknown.size(); ++i )
			std::cout << "DOF " << i << " has flag " << M_marker_unknown[i] << std::endl;
	}
	*/
}

void
Interpolation::identifyNodes_known ( )
{
	for (std::set<ID>::iterator it = M_GIdsKnownMesh_common.begin(); it != M_GIdsKnownMesh_common.end(); ++it)
	{
		if ( M_knownField->mapPtr()->map(Unique)->LID( static_cast<EpetraInt_Type> ( *it ) ) >= 0 )
		{
			M_GIdsKnownMesh.insert(*it);
		}
	}
}

void
Interpolation::identifyNodes_unknown ( )
{
	for (std::set<ID>::iterator it = M_GIdsUnknownMesh_common.begin(); it != M_GIdsUnknownMesh_common.end(); ++it)
	{
		if ( M_unknownField->mapPtr()->map(Unique)->LID( static_cast<EpetraInt_Type> ( *it ) ) >= 0 )
		{
			M_GIdsUnknownMesh.insert(*it);
		}
	}
}

void
Interpolation::buildKnownInterfaceMap()
{
	std::vector<int> dofKnown;
	dofKnown.reserve ( M_GIdsKnownMesh.size() );

	std::set<ID>::const_iterator i;

	for (UInt dim = 0; dim < nDimensions; ++dim)
		for ( i = M_GIdsKnownMesh.begin(); i != M_GIdsKnownMesh.end(); ++i )
		{
			dofKnown.push_back ( *i + dim * M_knownField->size()/nDimensions );
		}

	int* pointerToDofs (0);
	if (dofKnown.size() > 0)
	{
		pointerToDofs = &dofKnown[0];
	}

	M_knownInterfaceMap.reset ( new MapEpetra ( -1, static_cast<int> (dofKnown.size() ), pointerToDofs, M_knownField->mapPtr()->commPtr() ) );

	// std::cout << " M_knownInterfaceMap size " << M_knownInterfaceMap->map(Unique)->NumGlobalElements() << std::endl;
	// std::cout << " M_knownInterfaceMap local size " << M_knownInterfaceMap->map(Unique)->NumMyElements() << std::endl;
}

void
Interpolation::buildUnknownInterfaceMap()
{
	std::vector<int> dofunKnown;
	dofunKnown.reserve ( M_GIdsUnknownMesh.size() );

	std::set<ID>::const_iterator i;

	for (UInt dim = 0; dim < nDimensions; ++dim)
		for ( i = M_GIdsUnknownMesh.begin(); i != M_GIdsUnknownMesh.end(); ++i )
		{
			dofunKnown.push_back ( *i + dim * M_unknownField->size()/nDimensions );
		}

	int* pointerToDofs (0);
	if (dofunKnown.size() > 0)
	{
		pointerToDofs = &dofunKnown[0];
	}

	M_unknownInterfaceMap.reset ( new MapEpetra ( -1, static_cast<int> (dofunKnown.size() ), pointerToDofs, M_unknownField->mapPtr()->commPtr() ) );

//	std::cout << " M_unknownInterfaceMap size " << M_unknownInterfaceMap->map(Unique)->NumGlobalElements() << std::endl;
}

void
Interpolation::buildInterpolationOperatorMap()
{
	std::set<ID>::const_iterator ITrow;

	Int numtasks = M_knownField->mapPtr()->commPtr()->NumProc(); // Numero di processi
	int* numInterfaceDof (new int[numtasks]); // vettore lungo tanti quanti sono i processi
	int pid = M_knownField->mapPtr()->commPtr()->MyPID(); // ID processo
    M_pid = pid;
    int numMyElements;

	numMyElements = M_knownInterfaceMap->map (Unique)->NumMyElements(); // numero di elementi sul processo
	numInterfaceDof[pid] = numMyElements; // Ogni processore mette nella propria posizione il numero di elementi di interfaccia che ha

	mapPtr_Type subMap;
	subMap.reset ( new map_Type ( *M_knownInterfaceMap->map (Unique), (UInt) 0, M_knownField->size()/nDimensions) );

	M_numerationInterfaceKnown.reset (new VectorEpetra (*subMap, Unique) );

	for (int j = 0; j < numtasks; ++j)
	{
		M_knownField->mapPtr()->commPtr()->Broadcast ( &numInterfaceDof[j], 1, j);
	}

	for (int j = numtasks - 1; j > 0 ; --j)
	{
		numInterfaceDof[j] = numInterfaceDof[j - 1];
	}
	numInterfaceDof[0] = 0;
	for (int j = 1; j < numtasks ; ++j)
	{
		numInterfaceDof[j] += numInterfaceDof[j - 1];
	}

	UInt l = 0;

	Real M_interface = (UInt) M_knownInterfaceMap->map (Unique)->NumGlobalElements() / nDimensions; // Quanti dof ci sono nella mappa scalare di interfaccia
	for (l = 0, ITrow = M_GIdsKnownMesh.begin(); ITrow != M_GIdsKnownMesh.end() ; ++ITrow)
	{
		if (M_knownInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0)
		{
			(*M_numerationInterfaceKnown) [*ITrow ] = l + (int) (numInterfaceDof[pid] / nDimensions);
			if ( (int) (*M_numerationInterfaceKnown) (*ITrow ) != floor (l + numInterfaceDof[pid] / nDimensions + 0.2) )
			{
				std::cout << "ERROR! the numeration of the coupling map is not correct" << std::endl;
			}
			++l;
		}
	}

    M_numerationInterfaceKnownColumns.resize(numtasks);

    for (int j = 0; j < numtasks ; ++j)
        (M_numerationInterfaceKnownColumns[j]).reset(new vector_Type ( *M_numerationInterfaceKnown, j ) );

	std::vector<int> couplingVector;
	couplingVector.reserve ( (int) (M_knownInterfaceMap->map (Unique)->NumMyElements() ) );

	for ( ITrow = M_GIdsKnownMesh.begin(); ITrow != M_GIdsKnownMesh.end() ; ++ITrow)
	{
		if (M_knownInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0)
		{
			couplingVector.push_back ( (*M_numerationInterfaceKnown) (*ITrow ) );
		}
	}

	M_interpolationOperatorMap.reset (new MapEpetra (-1, static_cast< Int> ( couplingVector.size() ), &couplingVector[0], M_knownField->mapPtr()->commPtr() ) );

	M_interpolationOperatorMapVectorial.reset ( new MapEpetra ( *M_interpolationOperatorMap ) );
	*M_interpolationOperatorMapVectorial += *M_interpolationOperatorMap;
	*M_interpolationOperatorMapVectorial += *M_interpolationOperatorMap;

//	std::cout << "\n\n\n" << M_interpolationOperatorMap->map(Unique)->NumGlobalElements() << "\n\n\n";

	delete [] numInterfaceDof;
}

void
Interpolation::buildProjectionOperatorMap()
{
    std::set<ID>::const_iterator ITrow;

    Int numtasks = M_unknownField->mapPtr()->commPtr()->NumProc(); // Numero di processi
    int* numInterfaceDof (new int[numtasks]); // vettore lungo tanti quanti sono i processi
    int pid = M_unknownField->mapPtr()->commPtr()->MyPID(); // ID processo
    int numMyElements;

    numMyElements = M_unknownInterfaceMap->map (Unique)->NumMyElements(); // numero di elementi sul processo
    numInterfaceDof[pid] = numMyElements; // Ogni processore mette nella propria posizione il numero di elementi di interfaccia che ha

    mapPtr_Type subMap;
    subMap.reset ( new map_Type ( *M_unknownInterfaceMap->map (Unique), (UInt) 0, M_unknownField->size()/nDimensions) );

    M_numerationInterfaceUnknown.reset (new VectorEpetra (*subMap, Unique) );

    for (int j = 0; j < numtasks; ++j)
    {
        M_unknownField->mapPtr()->commPtr()->Broadcast ( &numInterfaceDof[j], 1, j);
    }

    for (int j = numtasks - 1; j > 0 ; --j)
    {
        numInterfaceDof[j] = numInterfaceDof[j - 1];
    }
    numInterfaceDof[0] = 0;
    for (int j = 1; j < numtasks ; ++j)
    {
        numInterfaceDof[j] += numInterfaceDof[j - 1];
    }

    UInt l = 0;

    Real M_interface = (UInt) M_unknownInterfaceMap->map (Unique)->NumGlobalElements() / nDimensions; // Quanti dof ci sono nella mappa scalare di interfaccia
    for (l = 0, ITrow = M_GIdsUnknownMesh.begin(); ITrow != M_GIdsUnknownMesh.end() ; ++ITrow)
    {
        if (M_unknownInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0)
        {
            (*M_numerationInterfaceUnknown) [*ITrow ] = l + (int) (numInterfaceDof[pid] / nDimensions);
            if ( (int) (*M_numerationInterfaceUnknown) (*ITrow ) != floor (l + numInterfaceDof[pid] / nDimensions + 0.2) )
            {
                std::cout << "ERROR! the numeration of the coupling map is not correct" << std::endl;
            }
            ++l;
        }
    }

    std::vector<int> couplingVector;
    couplingVector.reserve ( (int) (M_unknownInterfaceMap->map (Unique)->NumMyElements() ) );

    for ( ITrow = M_GIdsUnknownMesh.begin(); ITrow != M_GIdsUnknownMesh.end() ; ++ITrow)
    {
        if (M_unknownInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0)
        {
            couplingVector.push_back ( (*M_numerationInterfaceUnknown) (*ITrow ) );
        }
    }

    M_projectionOperatorMap.reset (new MapEpetra (-1, static_cast< Int> ( couplingVector.size() ), &couplingVector[0], M_unknownField->mapPtr()->commPtr() ) );

    M_projectionOperatorMapVectorial.reset ( new MapEpetra ( *M_projectionOperatorMap ) );
    *M_projectionOperatorMapVectorial += *M_projectionOperatorMap;
    *M_projectionOperatorMapVectorial += *M_projectionOperatorMap;

    delete [] numInterfaceDof;
}

void
Interpolation::buildRhs()
{
    M_RhsOne.reset (new vector_Type (*M_interpolationOperatorMap) );
    *M_RhsOne += 1;

    M_RhsF1.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF2.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF3.reset (new vector_Type (*M_interpolationOperatorMap) );

    UInt offset = M_knownField->size()/3;

    for(int i = 0; i < M_numerationInterfaceKnown->epetraVector().MyLength(); ++i)
    {
        (*M_RhsF1)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]]
        = (*M_knownField)(M_numerationInterfaceKnown->blockMap().GID(i));

        (*M_RhsF2)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]]
        = (*M_knownField)(M_numerationInterfaceKnown->blockMap().GID(i) + offset);

        (*M_RhsF3)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]]
        = (*M_knownField)(M_numerationInterfaceKnown->blockMap().GID(i) + 2*offset);

//        std::cout << "GID = " << M_numerationInterfaceKnown->blockMap().GID(i)
//        << ", Value =" << (*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] << std::endl;
    }
}

void
Interpolation::interpolateCostantField()
{
    vectorPtr_Type gamma_one;
    gamma_one.reset (new vector_Type (*M_interpolationOperatorMap) );

    // Preconditioner
    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
    precPtr.reset ( precRawPtr );

    LinearSolver solverOne;
    solverOne.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverOne.setParameters ( *M_belosList );
    solverOne.setPreconditioner ( precPtr );

    solverOne.setOperator (M_interpolationOperator);
    solverOne.setRightHandSide ( M_RhsOne );
    solverOne.solve ( gamma_one );

    M_rbf_one.reset (new vector_Type (*M_projectionOperatorMap) );
    M_projectionOperator->multiply (false, *gamma_one, *M_rbf_one);
}

void
Interpolation::expandGammaToOmega_Known(const vectorPtr_Type& vectorOnGamma, vectorPtr_Type& vectorOnOmega)
{
	vectorOnOmega.reset ( new vector_Type ( M_knownField->map(), Unique ) );
	vectorOnOmega->zero();

	UInt offset = vectorOnOmega->map().map(Unique)->NumGlobalElements()/3;
	UInt offsetGamma = vectorOnGamma->map().map(Unique)->NumGlobalElements()/3;

	for(int i = 0; i < M_numerationInterfaceKnown->epetraVector().MyLength(); ++i)
	{
		(*vectorOnOmega)[M_numerationInterfaceKnown->blockMap().GID(i)]
		                  = (*vectorOnGamma)((*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]);

		(*vectorOnOmega)[M_numerationInterfaceKnown->blockMap().GID(i) + offset]
		                  = (*vectorOnGamma)((*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] + offsetGamma );

		(*vectorOnOmega)[M_numerationInterfaceKnown->blockMap().GID(i) + 2*offset]
		                  = (*vectorOnGamma)((*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] + 2*offsetGamma);

		//        std::cout << "GID = " << M_numerationInterfaceKnown->blockMap().GID(i)
		//        << ", Value =" << (*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] << std::endl;
	}
}

void
Interpolation::restrictOmegaToGamma_Known(const vectorPtr_Type& vectorOnOmega, vectorPtr_Type& vectorOnGamma)
{
	vectorOnGamma.reset ( new vector_Type ( *M_interpolationOperatorMapVectorial, Unique ) );
	vectorOnGamma->zero();

	UInt offset = vectorOnOmega->map().map(Unique)->NumGlobalElements()/3;
	UInt offsetGamma = vectorOnGamma->map().map(Unique)->NumGlobalElements()/3;

	for ( int i = 0; i < M_numerationInterfaceKnown->epetraVector().MyLength(); ++i )
	{
		(*vectorOnGamma)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]]
		           = (*vectorOnOmega)(M_numerationInterfaceKnown->blockMap().GID(i));

		(*vectorOnGamma)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] + offsetGamma ]
		           = (*vectorOnOmega)(M_numerationInterfaceKnown->blockMap().GID(i) + offset);

		(*vectorOnGamma)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] + 2*offsetGamma ]
		           = (*vectorOnOmega)(M_numerationInterfaceKnown->blockMap().GID(i) + 2*offset);

		//        std::cout << "GID = " << M_numerationInterfaceKnown->blockMap().GID(i)
		//        << ", Value =" << (*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] << std::endl;
	}
}

void
Interpolation::interpolate()
{
    vectorPtr_Type gamma_f1;
    gamma_f1.reset (new vector_Type (*M_interpolationOperatorMap) );

    vectorPtr_Type gamma_f2;
    gamma_f2.reset (new vector_Type (*M_interpolationOperatorMap) );

    vectorPtr_Type gamma_f3;
    gamma_f3.reset (new vector_Type (*M_interpolationOperatorMap) );

    // Preconditioner
    if ( !M_precBuilt )
    {
    	prec_Type* precRawPtr;
		precRawPtr = new prec_Type;
		precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
		M_precPtr.reset ( precRawPtr );
		M_precBuilt = true;
    }

    // Solving for component x
    LinearSolver solverRBF1;
    solverRBF1.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverRBF1.setParameters ( *M_belosList );
    solverRBF1.setPreconditioner ( M_precPtr );

    solverRBF1.setOperator (M_interpolationOperator);
    solverRBF1.setRightHandSide (M_RhsF1);
    solverRBF1.solve (gamma_f1);

    // Solving for component y
    LinearSolver solverRBF2;
    solverRBF2.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverRBF2.setParameters ( *M_belosList );
    solverRBF2.setPreconditioner ( M_precPtr );

    solverRBF2.setOperator (M_interpolationOperator);
    solverRBF2.setRightHandSide (M_RhsF2);
    solverRBF2.solve (gamma_f2);

    // Solving for component z
    LinearSolver solverRBF3;
    solverRBF3.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverRBF3.setParameters ( *M_belosList );
    solverRBF3.setPreconditioner ( M_precPtr );

    solverRBF3.setOperator (M_interpolationOperator);
    solverRBF3.setRightHandSide (M_RhsF3);
    solverRBF3.solve (gamma_f3);

    vectorPtr_Type rbf_f1;
    rbf_f1.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution1;
    solution1.reset (new vector_Type (*M_projectionOperatorMap) );

    M_projectionOperator->multiply (false, *gamma_f1, *rbf_f1);

    // Rescaling component x
    *solution1 = *rbf_f1;
    *solution1 /= *M_rbf_one;

    vectorPtr_Type rbf_f2;
    rbf_f2.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution2;
    solution2.reset (new vector_Type (*M_projectionOperatorMap) );

    M_projectionOperator->multiply (false, *gamma_f2, *rbf_f2);

    // Rescaling component y
    *solution2 = *rbf_f2;
    *solution2 /= *M_rbf_one;

    vectorPtr_Type rbf_f3;
    rbf_f3.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution3;
    solution3.reset (new vector_Type (*M_projectionOperatorMap) );

    M_projectionOperator->multiply (false, *gamma_f3, *rbf_f3);

    // Rescaling component z
    *solution3 = *rbf_f3;
    *solution3 /= *M_rbf_one;

    // Solution of the interpolation problem

    UInt offset = M_unknownField->size()/3;

    for(int i = 0; i < M_numerationInterfaceUnknown->epetraVector().MyLength(); ++i)
    {
        (*M_unknownField)[M_numerationInterfaceUnknown->blockMap().GID(i)]
        = (*solution1)((*M_numerationInterfaceUnknown)[M_numerationInterfaceUnknown->blockMap().GID(i)]);

        (*M_unknownField)[M_numerationInterfaceUnknown->blockMap().GID(i) + offset]
        = (*solution2)((*M_numerationInterfaceUnknown)[M_numerationInterfaceUnknown->blockMap().GID(i)]);

        (*M_unknownField)[M_numerationInterfaceUnknown->blockMap().GID(i) + 2*offset]
        = (*solution3)((*M_numerationInterfaceUnknown)[M_numerationInterfaceUnknown->blockMap().GID(i)]);

        //        std::cout << "GID = " << M_numerationInterfaceKnown->blockMap().GID(i)
        //        << ", Value =" << (*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] << std::endl;
    }

    // Storing also the solution on gamma

    M_solutionOnGamma.reset ( new vector_Type ( *M_projectionOperatorMapVectorial ) );
    M_solutionOnGamma->zero();

    UInt offsetScalar = M_solutionOnGamma->size()/3;

    M_solutionOnGamma->subset(*solution1, *M_projectionOperatorMap, 0, 0);
    M_solutionOnGamma->subset(*solution2, *M_projectionOperatorMap, 0, offsetScalar );
    M_solutionOnGamma->subset(*solution3, *M_projectionOperatorMap, 0, 2*offsetScalar);
}

void
Interpolation::solution (vectorPtr_Type& Solution)
{
    Solution->zero();
    *Solution += *M_unknownField;
}

void
Interpolation::updateRhs(const vectorPtr_Type& newRhs)
{
    M_RhsF1->zero();
    M_RhsF2->zero();
    M_RhsF3->zero();

    UInt offset = M_knownField->size()/3;

    for(int i = 0; i < M_numerationInterfaceKnown->epetraVector().MyLength(); ++i)
    {
        (*M_RhsF1)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]]
        = (*newRhs)(M_numerationInterfaceKnown->blockMap().GID(i));

        (*M_RhsF2)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]]
        = (*newRhs)(M_numerationInterfaceKnown->blockMap().GID(i) + offset);

        (*M_RhsF3)[(*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)]]
        = (*newRhs)(M_numerationInterfaceKnown->blockMap().GID(i) + 2*offset);

        //        std::cout << "GID = " << M_numerationInterfaceKnown->blockMap().GID(i)
        //        << ", Value =" << (*M_numerationInterfaceKnown)[M_numerationInterfaceKnown->blockMap().GID(i)] << std::endl;
    }
}

} // Namespace LifeV
