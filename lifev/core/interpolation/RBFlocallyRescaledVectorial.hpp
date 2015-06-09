//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief A short description of the file content

    @author Davide Forti <forti@mathicsepc48.epfl.ch>
    @date 14 Mar 2013

    A more detailed description of the file (if necessary)
 */

#ifndef RBFLOCALLYRESCALEDVECTORIAL_H
#define RBFLOCALLYRESCALEDVECTORIAL_H 1

#include <lifev/core/interpolation/RBFInterpolation.hpp>

namespace LifeV {

template <typename mesh_Type>
class RBFlocallyRescaledVectorial: public RBFInterpolation<mesh_Type>
{
public:

    typedef boost::shared_ptr<mesh_Type>                                          meshPtr_Type;

    typedef VectorEpetra                                                          vector_Type;
    typedef boost::shared_ptr<vector_Type >                                       vectorPtr_Type;

    typedef MatrixEpetra<double>                                                  matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                                        matrixPtr_Type;

    typedef std::vector<int>                                                      flagContainer_Type;

    typedef boost::unordered_set<ID>                                              idContainer_Type;

    typedef MapEpetra                                                             map_Type;
    typedef boost::shared_ptr<MapEpetra>                                          mapPtr_Type;

    typedef GhostHandler<mesh_Type>                                               neighboring_Type;
    typedef boost::shared_ptr<neighboring_Type>                                   neighboringPtr_Type;

    typedef LifeV::Preconditioner                                                 basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>                                      basePrecPtr_Type;

    typedef LifeV::PreconditionerIfpack                                           prec_Type;
    typedef boost::shared_ptr<prec_Type>                                          precPtr_Type;

    typedef Teuchos::RCP< Teuchos::ParameterList >                                parameterList_Type;

    RBFlocallyRescaledVectorial();

    virtual ~RBFlocallyRescaledVectorial();

    void setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags );

    void setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList);

    void buildOperators();

    void interpolationOperator();

    void projectionOperator();

    void buildRhs();

    void interpolateCostantField();

    void identifyNodes (meshPtr_Type LocalMesh, std::set<ID>& GID_nodes, vectorPtr_Type CheckVector);

    bool isInside (ID pointMarker, flagContainer_Type Flags);

    double computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID);

    double rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius);

    void interpolate();

    void solution (vectorPtr_Type& Solution);

    void solutionrbf (vectorPtr_Type& Solution_rbf);

    void updateRhs(const vectorPtr_Type& newRhs);

    void getinterpolationOperatorMap(mapPtr_Type& map)
    {
        map.reset(new map_Type(*M_interpolationOperatorMapVectorial) );
        //map.reset(new map_Type(*M_interpolationOperatorMap));
        //*map += *M_interpolationOperatorMap;
        //*map += *M_interpolationOperatorMap;
    }

    void getSolutionOnGamma (vectorPtr_Type& GammaSolution) { GammaSolution.reset(new vector_Type ( *M_solOnGamma ) ); }

    void getprojectionOperatorMap (mapPtr_Type& map) { map.reset(new map_Type(*M_projectionOperatorMap)); }

    void buildUnknownVectorialInterfaceMap();

    // Methods added after changing the maps

    void buildKnownInterfaceMap();

    void buildInterpolationOperatorMap();

private:

    meshPtr_Type        M_fullMeshKnown;
    meshPtr_Type        M_localMeshKnown;
    meshPtr_Type        M_fullMeshUnknown;
    meshPtr_Type        M_localMeshUnknown;
    flagContainer_Type  M_flags;
    vectorPtr_Type      M_knownField;
    vectorPtr_Type      M_unknownField;
    GetPot              M_datafile;
    parameterList_Type  M_belosList;
    std::set<ID>        M_GIdsKnownMesh;
    std::set<ID>        M_GIdsUnknownMesh;
    matrixPtr_Type      M_interpolationOperator;
    matrixPtr_Type      M_projectionOperator;
    vectorPtr_Type      M_RhsF1;
    vectorPtr_Type      M_RhsF2;
    vectorPtr_Type      M_RhsF3;
    vectorPtr_Type      M_RhsOne;
    vectorPtr_Type      M_rbf_one;
    mapPtr_Type         M_interpolationOperatorMap;
    mapPtr_Type         M_interpolationOperatorMapVectorial;
    mapPtr_Type         M_projectionOperatorMap;
    neighboringPtr_Type M_neighbors;
    vectorPtr_Type      M_unknownField_rbf;

    matrixPtr_Type      M_interpolationOperatorVecMap;
    matrixPtr_Type      M_projectionOperatorVecMap;

    vectorPtr_Type      M_solOnGamma;
    mapPtr_Type         M_gammaMapUnknownVectorial;
    UInt                M_links;
    Real                M_radius;

    // Added for maps

    mapPtr_Type         M_knownInterfaceMap;
    vectorPtr_Type      M_numerationInterfaceKnown;
};

template <typename mesh_Type>
RBFlocallyRescaledVectorial<mesh_Type>::RBFlocallyRescaledVectorial()
{}

template <typename mesh_Type>
RBFlocallyRescaledVectorial<mesh_Type>::~RBFlocallyRescaledVectorial()
{}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags )
{
    M_fullMeshKnown = fullMeshKnown;
    M_localMeshKnown = localMeshKnown;
    M_fullMeshUnknown = fullMeshUnknown;
    M_localMeshUnknown = localMeshUnknown;
    M_flags = flags;
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList)
{
    M_knownField.reset( new vector_Type ( *KnownField ) );
    M_unknownField.reset( new vector_Type ( *UnknownField ) );
    M_datafile     = datafile;
    M_belosList    = belosList;
    M_links = M_datafile("interpolation/n_links",1);
}

template <typename Mesh>
void RBFlocallyRescaledVectorial<Mesh>::buildOperators()
{
    this->interpolationOperator();
//    this->projectionOperator();
//    this->buildRhs();
//    this->interpolateCostantField();
}

template <typename Mesh>
void RBFlocallyRescaledVectorial<Mesh>::buildKnownInterfaceMap()
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
}

template <typename Mesh>
void RBFlocallyRescaledVectorial<Mesh>::buildInterpolationOperatorMap()
{
	std::set<ID>::const_iterator ITrow;

	Int numtasks = M_knownField->mapPtr()->commPtr()->NumProc(); // Numero di processi
	int* numInterfaceDof (new int[numtasks]); // vettore lungo tanti quanti sono i processi
	int pid = M_knownField->mapPtr()->commPtr()->MyPID(); // ID processo
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

	delete [] numInterfaceDof;
}

template <typename Mesh>
void RBFlocallyRescaledVectorial<Mesh>::interpolationOperator()
{
    // Identifying dofs to be taken into account
	this->identifyNodes (M_localMeshKnown, M_GIdsKnownMesh, M_knownField);

	// This map will be used to select dofs from the whole vector as right hand side
	buildKnownInterfaceMap();

	// Building map that will be used to close the matrix, this map will be ordered from
	// 1 to the number of DOFs at the interface
	buildInterpolationOperatorMap();

	// Object needed to find neighbors by mesh connectivity
	M_neighbors.reset ( new neighboring_Type ( M_fullMeshKnown, M_localMeshKnown, M_knownField->mapPtr(), M_knownField->mapPtr()->commPtr() ) );
    if (M_flags[0] == -1)
    {
        M_neighbors->setUpNeighbors ();
    }
    else
    {
    	// List of neighbors, for each dof finding its neighbors
    	M_neighbors->createPointPointNeighborsList (M_flags);
    }

    int LocalNodesNumber = M_GIdsKnownMesh.size();

    std::vector<double>   RBF_radius (LocalNodesNumber);
    std::vector<boost::unordered_set<ID> > MatrixGraph (LocalNodesNumber);
    int* GlobalID = new int[LocalNodesNumber];

    int k = 0;
    int Max_entries = 0;

    for (std::set<ID>::iterator it = M_GIdsKnownMesh.begin(); it != M_GIdsKnownMesh.end(); ++it)
    {
    	// Selecting each dof taken into account
    	GlobalID[k] = *it;
        
    	// For each of them, identify the neighbors within a certain number of circles M_links
        neighbors_Type Neighbors;
        Neighbors = M_neighbors->circleNeighbors ( *it, M_links );
        Neighbors.insert ( *it );
        MatrixGraph[k] = Neighbors;
        
        RBF_radius[k] = computeRBFradius ( M_fullMeshKnown, M_fullMeshKnown, Neighbors, GlobalID[k]);
        neighbors_Type NeighborsR;
        NeighborsR = M_neighbors->neighborsWithinRadius ( GlobalID[k], RBF_radius[k] );
        NeighborsR.insert ( GlobalID[k] );
        MatrixGraph[k] = NeighborsR;

        ++k;
    }

    // Matrix for interpolation operator
    M_interpolationOperator.reset (new matrix_Type (*M_interpolationOperatorMap, 100) );
    Real Values;
    std::cout << M_numerationInterfaceKnown->epetraVector();

    vector_Type test(*M_numerationInterfaceKnown,0);

    std::cout << test.epetraVector();

    M_interpolationOperator->zero();

    // Loop over all the dofs taken into account by the interpolation
    for ( int i = 0 ; i < LocalNodesNumber; ++i )
    {
        k = 0;

        if (GlobalID[i] == 9 )
        {
        	std::cout << "\n\n Linea dof 9, che locale e' " << (*M_numerationInterfaceKnown)( GlobalID[i] ) << "\n\n";
        }

        // For each i-th dof, evaluate the rbf between the dof and its neighbors
        for ( boost::unordered_set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
        {
//        	if ( M_numerationInterfaceKnown->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (*it) ) >= 0 )
//        	{
        		Values = rbf (  M_fullMeshKnown->point (GlobalID[i]).x(),
        						M_fullMeshKnown->point (GlobalID[i]).y(),
        						M_fullMeshKnown->point (GlobalID[i]).z(),
        						M_fullMeshKnown->point (*it).x(),
        						M_fullMeshKnown->point (*it).y(),
        						M_fullMeshKnown->point (*it).z(),
        						RBF_radius[i] );

            	if (GlobalID[i] == 9 )
            	{
            		std::cout << "(" << *it << ", " << (*M_numerationInterfaceKnown)(*it) << ", " << Values << ")  ";
            	}

            	//M_interpolationOperator->matrixPtr()->InsertGlobalValues ()

//				M_interpolationOperator->addToCoefficient( (*M_numerationInterfaceKnown)( GlobalID[i] ),
//														   *it,
//														   Values );

//        	}
        }
    }
    
    M_interpolationOperator->globalAssemble();
    M_interpolationOperator->spy("M_interpolationOperator2procs");
    // M_interpolationOperator->exportToHDF5("InputField", "M_interpolationOperator", false);

    delete[] GlobalID;
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::projectionOperator()
{
	/*
    // Identifying dofs to be taken into account
    this->identifyNodes (M_localMeshUnknown, M_GIdsUnknownMesh, M_unknownField);

    // Total number of dofs taken into account
    int LocalNodesNumber = M_GIdsUnknownMesh.size();

    std::vector<double>  RBF_radius (LocalNodesNumber);
    std::vector<boost::unordered_set<ID> > MatrixGraph (LocalNodesNumber);
    int* ElementsPerRow = new int[LocalNodesNumber];
    int* GlobalID = new int[LocalNodesNumber];
    int k = 0;
    int Max_entries = 0;

    // I need to find the closest point in the "known mesh" to use its radius
    double d;
    double d_min;
    int nearestPoint;

    for (boost::unordered_set<ID>::iterator it = M_GIdsUnknownMesh.begin(); it != M_GIdsUnknownMesh.end(); ++it)
    {
        GlobalID[k] = *it;
        d_min = 100;
        for (int j = 0; j <  M_fullMeshKnown->numVertices(); ++j)
        {
            if ( M_flags[0] == -1 || this->isInside (M_fullMeshKnown->point (j).markerID(), M_flags) )
            {
                d = std::sqrt ( ( M_fullMeshKnown->point (j).x() - M_fullMeshUnknown->point (GlobalID[k]).x() ) *
                		        ( M_fullMeshKnown->point (j).x() - M_fullMeshUnknown->point (GlobalID[k]).x() ) +
                                ( M_fullMeshKnown->point (j).y() - M_fullMeshUnknown->point (GlobalID[k]).y() ) *
                                ( M_fullMeshKnown->point (j).y() - M_fullMeshUnknown->point (GlobalID[k]).y() ) +
                                ( M_fullMeshKnown->point (j).z() - M_fullMeshUnknown->point (GlobalID[k]).z() ) *
                                ( M_fullMeshKnown->point (j).z() - M_fullMeshUnknown->point (GlobalID[k]).z() ) );
                if (d < d_min)
                {
                    d_min = d;
                    nearestPoint = M_fullMeshKnown->point (j).id();
                }
            }
        }

        // For each of them, identify the neighbors on the other mesh within a certain number of circles M_links
        neighbors_Type Neighbors;
        Neighbors = M_neighbors->circleNeighbors ( nearestPoint, M_links );
        Neighbors.insert ( nearestPoint );
        
        RBF_radius[k] = computeRBFradius ( M_fullMeshKnown, M_fullMeshUnknown, Neighbors, GlobalID[k]);
        neighbors_Type NeighborsR;
        NeighborsR = M_neighbors->neighborsWithinRadius ( nearestPoint, RBF_radius[k] );
        NeighborsR.insert ( nearestPoint );
        MatrixGraph[k] = NeighborsR;
        
        ElementsPerRow[k] = MatrixGraph[k].size();
        if (ElementsPerRow[k] > Max_entries)
        {
            Max_entries = ElementsPerRow[k];
        }
        ++k;
    }

    M_projectionOperatorMap.reset (new map_Type (-1, LocalNodesNumber, GlobalID, M_unknownField->mapPtr()->commPtr() ) );
    M_projectionOperator.reset (new matrix_Type (*M_projectionOperatorMap, ElementsPerRow) );

    int* Indices = new int[Max_entries];
    double* Values = new double[Max_entries];

    for ( int i = 0 ; i < LocalNodesNumber; ++i )
    {
        k = 0;
        for ( boost::unordered_set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
        {
            Indices[k] = *it;
            Values[k]  = rbf ( M_fullMeshUnknown->point (GlobalID[i]).x(),
                               M_fullMeshUnknown->point (GlobalID[i]).y(),
                               M_fullMeshUnknown->point (GlobalID[i]).z(),
                               M_fullMeshKnown->point (*it).x(),
                               M_fullMeshKnown->point (*it).y(),
                               M_fullMeshKnown->point (*it).z(),
                               RBF_radius[i]);
            ++k;
        }
        M_projectionOperator->matrixPtr()->InsertGlobalValues (GlobalID[i], k, Values, Indices);
    }
    M_projectionOperator->globalAssemble (M_interpolationOperatorMap, M_projectionOperatorMap);
    delete[] Indices;
    delete[] Values;
    delete[] ElementsPerRow;
    delete[] GlobalID;
    */
}

template <typename mesh_Type>
inline double RBFlocallyRescaledVectorial<mesh_Type>::computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID)
{
    double r = 0;
    double r_max = 0;
    for (idContainer_Type::iterator it = Neighbors.begin(); it != Neighbors.end(); ++it)
    {
        r = std::sqrt (  ( MeshGID->point ( GlobalID ).x() - MeshNeighbors->point ( *it ).x() ) * 
                         ( MeshGID->point ( GlobalID ).x() - MeshNeighbors->point ( *it ).x() ) + 
                         ( MeshGID->point ( GlobalID ).y() - MeshNeighbors->point ( *it ).y() ) *
                         ( MeshGID->point ( GlobalID ).y() - MeshNeighbors->point ( *it ).y() ) +
                         ( MeshGID->point ( GlobalID ).z() - MeshNeighbors->point ( *it ).z() ) *
			             ( MeshGID->point ( GlobalID ).z() - MeshNeighbors->point ( *it ).z() ) );
        r_max = ( r > r_max ) ? r : r_max;
    }
    return r_max;
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::buildRhs()
{
	/*
    M_RhsF1.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF2.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF3.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsOne.reset (new vector_Type (*M_interpolationOperatorMap) );

    M_RhsF1->subset (*M_knownField, *M_interpolationOperatorMap, 0, 0);
    M_RhsF2->subset (*M_knownField, *M_interpolationOperatorMap, M_knownField->size()/3, 0);
    M_RhsF3->subset (*M_knownField, *M_interpolationOperatorMap, M_knownField->size()/3*2, 0);
    *M_RhsOne += 1;
    */
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::interpolateCostantField()
{
	/*
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
    */
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::interpolate()
{
	/*
    vectorPtr_Type gamma_f1;
    gamma_f1.reset (new vector_Type (*M_interpolationOperatorMap) );

    vectorPtr_Type gamma_f2;
    gamma_f2.reset (new vector_Type (*M_interpolationOperatorMap) );

    vectorPtr_Type gamma_f3;
    gamma_f3.reset (new vector_Type (*M_interpolationOperatorMap) );


    // Preconditioner
    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
    precPtr.reset ( precRawPtr );

    // Solving for component x
    LinearSolver solverRBF1;
    solverRBF1.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverRBF1.setParameters ( *M_belosList );
    solverRBF1.setPreconditioner ( precPtr );

    solverRBF1.setOperator (M_interpolationOperator);
    solverRBF1.setRightHandSide (M_RhsF1);
    solverRBF1.solve (gamma_f1);

    // Solving for component y
    LinearSolver solverRBF2;
    solverRBF2.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverRBF2.setParameters ( *M_belosList );
    solverRBF2.setPreconditioner ( precPtr );

    solverRBF2.setOperator (M_interpolationOperator);
    solverRBF2.setRightHandSide (M_RhsF2);
    solverRBF2.solve (gamma_f2);

    // Solving for component z
    LinearSolver solverRBF3;
    solverRBF3.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverRBF3.setParameters ( *M_belosList );
    solverRBF3.setPreconditioner ( precPtr );

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

    M_unknownField_rbf.reset (new vector_Type (M_unknownField->map() ) );
    M_unknownField_rbf->subset (*rbf_f1, *M_projectionOperatorMap, 0, 0);
    M_unknownField_rbf->subset (*rbf_f2, *M_projectionOperatorMap, 0, M_unknownField->size()/3);
    M_unknownField_rbf->subset (*rbf_f3, *M_projectionOperatorMap, 0, M_unknownField->size()/3*2);

    // Solution of the interpolation problem
    M_unknownField->subset (*solution1, *M_projectionOperatorMap, 0, 0);
    M_unknownField->subset (*solution2, *M_projectionOperatorMap, 0, M_unknownField->size()/3);
    M_unknownField->subset (*solution3, *M_projectionOperatorMap, 0, M_unknownField->size()/3*2);
	*/

//    /*
//    mapPtr_Type solOnGammaMap;
//
//    solOnGammaMap.reset( new map_Type( *M_projectionOperatorMap ) );
//    *solOnGammaMap += *M_projectionOperatorMap;
//    *solOnGammaMap += *M_projectionOperatorMap;
//    */
//
//    /*
//    M_solOnGamma.reset(new vector_Type( *M_gammaMapUnknownVectorial ) );
//
//    M_solOnGamma->subset (*solution1, *M_projectionOperatorMap, 0, 0);
//    M_solOnGamma->subset (*solution2, *M_projectionOperatorMap, 0, M_unknownField->size()/3);
//    M_solOnGamma->subset (*solution3, *M_projectionOperatorMap, 0, 2*M_unknownField->size()/3);
//    */
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::buildUnknownVectorialInterfaceMap()
{
	/*
    boost::unordered_set<ID> GID_vectorial;
    for ( UInt i = 0; i < M_localMeshUnknown->numVertices(); ++i )
        if ( this->isInside (M_localMeshUnknown->point (i).markerID(), M_flags) )
            if (M_unknownField->blockMap().LID ( static_cast<EpetraInt_Type> (M_localMeshUnknown->point (i).id()) ) != -1)
                for(int nDim = 0; nDim < 3; ++nDim)
                    GID_vectorial.insert (M_localMeshUnknown->point (i).id() + nDim*M_fullMeshUnknown->numVertices());

    int LocalNodesNumber = GID_vectorial.size();
    //std::cout << LocalNodesNumber << std::endl;
    int* GlobalID = new int[LocalNodesNumber];
    int k = 0;

    for (boost::unordered_set<ID>::iterator it = GID_vectorial.begin(); it != GID_vectorial.end(); ++it)
    {
        GlobalID[k] = *it;
        ++k;
    }

    M_gammaMapUnknownVectorial.reset (new map_Type (-1, LocalNodesNumber, GlobalID, M_unknownField->mapPtr()->commPtr() ) );
    delete GlobalID;
    */
}


template <typename mesh_Type>
inline void RBFlocallyRescaledVectorial<mesh_Type>::identifyNodes (meshPtr_Type LocalMesh, std::set<ID>& GID_nodes, vectorPtr_Type CheckVector)
{
    if (M_flags[0] == -1)
    {
        for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
            if (CheckVector->blockMap().LID ( static_cast<EpetraInt_Type> ( LocalMesh->point (i).id() ) ) != -1)
            {
                GID_nodes.insert (LocalMesh->point (i).id() );
            }
    }
    else
    {
	// Checking which points have to be considered. 
	// TODO: here the code could be otimized but currently it is 
	// not since we want to keep the code general, i.e., such
	// that it works for mesh to mesh and also for interface to interface.
	// In the case interface to interface we do not have to loop all over
	// the vertices, but on the boundary faces. 
	// TODO: works just for P1 now since we ask for the markerID of the 
	// points.
        for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
            if ( this->isInside (LocalMesh->point (i).markerID(), M_flags) )
                if (CheckVector->blockMap().LID ( static_cast<EpetraInt_Type> (LocalMesh->point (i).id()) ) != -1)
                {
                    GID_nodes.insert (LocalMesh->point (i).id() );
                }
    }
}

template <typename mesh_Type>
inline bool RBFlocallyRescaledVectorial<mesh_Type>::isInside (ID pointMarker, flagContainer_Type flags)
{
    // Checking if the point with id number ID has to be considered
    for (UInt i = 0; i < flags.size(); ++i)
        if (pointMarker == flags[i])
        {
	    return true;
        }
    return false;
}


template <typename mesh_Type>
inline double RBFlocallyRescaledVectorial<mesh_Type>::rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius)
{
    double distance = std::sqrt ( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2) );
    return (1 - distance / radius) * (1 - distance / radius) * (1 - distance / radius) * (1 - distance / radius) * (4 * distance / radius + 1);
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::updateRhs(const vectorPtr_Type& newRhs)
{
	/*
    M_RhsF1->zero();
    M_RhsF1->subset (*newRhs, *M_interpolationOperatorMap, 0, 0);
    M_RhsF2->zero();
    M_RhsF2->subset (*newRhs, *M_interpolationOperatorMap, newRhs->size()/3, 0);
    M_RhsF3->zero();
    M_RhsF3->subset (*newRhs, *M_interpolationOperatorMap, newRhs->size()/3*2, 0);
    */
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::solution (vectorPtr_Type& Solution)
{
	/*
    Solution->zero();
    *Solution += *M_unknownField;
    */
}

template <typename mesh_Type>
void RBFlocallyRescaledVectorial<mesh_Type>::solutionrbf (vectorPtr_Type& Solution_rbf)
{
	/*
    Solution_rbf = M_unknownField_rbf;
	*/
}

//! Factory create function
template <typename mesh_Type>
inline RBFInterpolation<mesh_Type> * createRBFlocallyRescaledVectorial()
{
    return new RBFlocallyRescaledVectorial< mesh_Type > ();
}
namespace
{
static bool S_registerTriLRV = RBFInterpolation<LifeV::RegionMesh<LinearTriangle > >::InterpolationFactory::instance().registerProduct ( "RBFlocallyRescaledVectorial", &createRBFlocallyRescaledVectorial<LifeV::RegionMesh<LinearTriangle > > );
static bool S_registerTetLRV = RBFInterpolation<LifeV::RegionMesh<LinearTetra > >::InterpolationFactory::instance().registerProduct ( "RBFlocallyRescaledVectorial", &createRBFlocallyRescaledVectorial<LifeV::RegionMesh<LinearTetra > > );
}

} // Namespace LifeV

#endif /* RBFLOCALLYRESCALEDVECTORIAL_H */
