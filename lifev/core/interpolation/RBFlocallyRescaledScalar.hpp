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
    @date 13 Mar 2013

    A more detailed description of the file (if necessary)
 */

#ifndef RBFLOCALLYRESCALEDSCALAR_H
#define RBFLOCALLYRESCALEDSCALAR_H 1

#include <lifev/core/interpolation/RBFInterpolation.hpp>

namespace LifeV
{

template <typename mesh_Type>
class RBFlocallyRescaledScalar: public RBFInterpolation<mesh_Type>
{
public:

    typedef boost::shared_ptr<mesh_Type>                                          meshPtr_Type;

    typedef VectorEpetra                                                          vector_Type;
    typedef boost::shared_ptr<vector_Type >                                       vectorPtr_Type;

    typedef MatrixEpetra<double>                                                  matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                                        matrixPtr_Type;

    typedef std::vector<int>                                                      flagContainer_Type;

    typedef boost::unordered_set<ID>                                                          idContainer_Type;

    typedef MapEpetra                                                             map_Type;
    typedef boost::shared_ptr<MapEpetra>                                          mapPtr_Type;

    typedef GhostHandler<mesh_Type>                                               neighbors_Type;
    typedef boost::shared_ptr<neighbors_Type>                                     neighborsPtr_Type;

    typedef LifeV::Preconditioner                                                 basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>                                      basePrecPtr_Type;

    typedef LifeV::PreconditionerIfpack                                           prec_Type;
    typedef boost::shared_ptr<prec_Type>                                          precPtr_Type;

    typedef Teuchos::RCP< Teuchos::ParameterList >                                parameterList_Type;

    RBFlocallyRescaledScalar();

    virtual ~RBFlocallyRescaledScalar();

    void setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags );

    void setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList);

    void buildOperators();

    void interpolationOperator();

    void projectionOperator();

    void buildRhs();

    void interpolateCostantField();

    void identifyNodes (meshPtr_Type LocalMesh, boost::unordered_set<ID>& GID_nodes, vectorPtr_Type CheckVector);

    bool isInside (ID pointMarker, flagContainer_Type Flags);

    double computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID);

    double rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius);

    void interpolate();

    void solution (vectorPtr_Type& Solution);

    void solutionrbf (vectorPtr_Type& Solution_rbf);

    void updateRhs(vectorPtr_Type newRhs);

    void setRadius ( double radius );

    void getinterpolationOperatorMap(mapPtr_Type& map){ map.reset(new map_Type(*M_interpolationOperatorMap)); }

    void getprojectionOperatorMap(mapPtr_Type& map){ map.reset(new map_Type(*M_projectionOperatorMap)); }

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
    idContainer_Type    M_GIdsKnownMesh;
    idContainer_Type    M_GIdsUnknownMesh;
    matrixPtr_Type      M_interpolationOperator;
    matrixPtr_Type      M_projectionOperator;
    vectorPtr_Type      M_RhsF;
    vectorPtr_Type      M_RhsOne;
    vectorPtr_Type      M_rbf_one;
    mapPtr_Type         M_interpolationOperatorMap;
    mapPtr_Type         M_projectionOperatorMap;
    neighborsPtr_Type   M_neighbors;
    vectorPtr_Type      M_unknownField_rbf;
    double              M_radius;

    std::vector<double> M_kx;
    std::vector<double> M_ky;
    std::vector<double> M_kz;

    std::vector<double> M_ukx;
    std::vector<double> M_uky;
    std::vector<double> M_ukz;
};

template <typename mesh_Type>
RBFlocallyRescaledScalar<mesh_Type>::RBFlocallyRescaledScalar()
{}

template <typename mesh_Type>
RBFlocallyRescaledScalar<mesh_Type>::~RBFlocallyRescaledScalar()
{}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags )
{
    M_fullMeshKnown = fullMeshKnown;
    M_localMeshKnown = localMeshKnown;
    M_fullMeshUnknown = fullMeshUnknown;
    M_localMeshUnknown = localMeshUnknown;
    M_flags = flags;

    M_kx.resize(M_fullMeshKnown->numVertices());
    M_ky.resize(M_fullMeshKnown->numVertices());
    M_kz.resize(M_fullMeshKnown->numVertices());

    for (int j = 0; j <  M_fullMeshKnown->numVertices(); ++j)
    {
    	M_kx[j] =  M_fullMeshKnown->point(j).x();
    	M_ky[j] =  M_fullMeshKnown->point(j).y();
    	M_kz[j] =  M_fullMeshKnown->point(j).z();
    }

    M_ukx.resize(M_fullMeshUnknown->numVertices());
    M_uky.resize(M_fullMeshUnknown->numVertices());
    M_ukz.resize(M_fullMeshUnknown->numVertices());

    for (int j = 0; j <  M_fullMeshUnknown->numVertices(); ++j)
    {
        M_ukx[j] =  M_fullMeshUnknown->point(j).x();
        M_uky[j] =  M_fullMeshUnknown->point(j).y();
        M_ukz[j] =  M_fullMeshUnknown->point(j).z();
    }

}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList)
{
    M_knownField   = KnownField;
    M_unknownField = UnknownField;
    M_datafile     = datafile;
    M_belosList    = belosList;
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::setRadius ( double radius )
{
    M_radius = radius;
}

template <typename Mesh>
void RBFlocallyRescaledScalar<Mesh>::buildOperators()
{
    LifeChrono TimeBuilding;
    TimeBuilding.start();
    this->interpolationOperator();
    this->projectionOperator();
    TimeBuilding.stop();
    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "Time to assembly operators = " << TimeBuilding.diff() << std::endl;

    this->buildRhs();
    this->interpolateCostantField();
}

template <typename Mesh>
void RBFlocallyRescaledScalar<Mesh>::interpolationOperator()
{
    this->identifyNodes (M_localMeshKnown, M_GIdsKnownMesh, M_knownField);
    M_neighbors.reset ( new neighbors_Type ( M_fullMeshKnown, M_localMeshKnown, M_knownField->mapPtr(), M_knownField->mapPtr()->commPtr() ) );
    if (M_flags[0] == -1)
    {
        M_neighbors->setUpNeighbors();
    }
    else
    {
        M_neighbors->createPointPointNeighborsList (M_flags);
    }

    int LocalNodesNumber = M_GIdsKnownMesh.size();

    std::vector<double>   RBF_radius (LocalNodesNumber);
    std::vector<boost::unordered_set<ID> > MatrixGraph (LocalNodesNumber);
    int* ElementsPerRow = new int[LocalNodesNumber];
    int* GlobalID = new int[LocalNodesNumber];
    int k = 0;
    int Max_entries = 0;

    for (boost::unordered_set<ID>::iterator it = M_GIdsKnownMesh.begin(); it != M_GIdsKnownMesh.end(); ++it)
    {
        GlobalID[k] = *it;
        MatrixGraph[k] = M_neighbors->pointPointNeighborsList() [GlobalID[k]];
        MatrixGraph[k].insert (GlobalID[k]);
        RBF_radius[k] = computeRBFradius ( M_fullMeshKnown, M_fullMeshKnown, MatrixGraph[k], GlobalID[k]);
        ElementsPerRow[k] = MatrixGraph[k].size();
        if (ElementsPerRow[k] > Max_entries)
        {
            Max_entries = ElementsPerRow[k];
        }
        ++k;
    }

    M_interpolationOperatorMap.reset (new map_Type (-1, LocalNodesNumber, GlobalID, M_knownField->mapPtr()->commPtr() ) );
    M_interpolationOperator.reset (new matrix_Type (*M_interpolationOperatorMap, ElementsPerRow) );

    int* Indices = new int[Max_entries];
    double* Values = new double[Max_entries];

    for ( int i = 0 ; i < LocalNodesNumber; ++i )
    {
        k = 0;
        for ( boost::unordered_set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
        {
            Indices[k] = *it;
            Values[k]  = rbf ( M_kx[GlobalID[i]], M_ky[GlobalID[i]], M_kz[GlobalID[i]],
                                               M_kx[*it], M_ky[*it], M_kz[*it],
                                               RBF_radius[i]);
            ++k;
        }
        M_interpolationOperator->matrixPtr()->InsertGlobalValues (GlobalID[i], k, Values, Indices);
    }
    M_interpolationOperator->globalAssemble();
    M_interpolationOperator->spy("TEST");
    delete Indices;
    delete Values;
    delete ElementsPerRow;
    delete GlobalID;
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::projectionOperator()
{

    this->identifyNodes (M_localMeshUnknown, M_GIdsUnknownMesh, M_unknownField);

    int LocalNodesNumber = M_GIdsUnknownMesh.size();

    std::vector<double>        RBF_radius (LocalNodesNumber);
    std::vector<boost::unordered_set<ID> > MatrixGraph (LocalNodesNumber);
    int* ElementsPerRow = new int[LocalNodesNumber];
    int* GlobalID = new int[LocalNodesNumber];
    int k = 0;
    int Max_entries = 0;
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
                d = std::sqrt ( pow (M_kx[j] - M_ukx[GlobalID[k]], 2)
                                + pow (M_ky[j] - M_uky[GlobalID[k]], 2)
                                + pow (M_kz[j] - M_ukz[GlobalID[k]], 2) );
                if (d < d_min)
                {
                    d_min = d;
                    nearestPoint = M_fullMeshKnown->point (j).id();
                }
            }
        }
        MatrixGraph[k] = M_neighbors->pointPointNeighborsList() [nearestPoint];
        MatrixGraph[k].insert (nearestPoint);
        RBF_radius[k] = computeRBFradius ( M_fullMeshKnown, M_fullMeshUnknown, MatrixGraph[k], GlobalID[k]);
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
            Values[k]  = rbf ( M_ukx[GlobalID[i]], M_uky[GlobalID[i]], M_ukz[GlobalID[i]],
            				   M_kx[*it], M_ky[*it], M_kz[*it], RBF_radius[i]);
            ++k;
        }
        M_projectionOperator->matrixPtr()->InsertGlobalValues (GlobalID[i], k, Values, Indices);
    }
    M_projectionOperator->globalAssemble (M_interpolationOperatorMap, M_projectionOperatorMap);
    delete Indices;
    delete Values;
    delete ElementsPerRow;
    delete GlobalID;
}

template <typename mesh_Type>
double RBFlocallyRescaledScalar<mesh_Type>::computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID)
{
    double r = 0;
    double r_max = 0;
    for (idContainer_Type::iterator it = Neighbors.begin(); it != Neighbors.end(); ++it)
    {
        r = std::sqrt ( pow ( MeshGID->point ( GlobalID ).x() - MeshNeighbors->point ( *it ).x(), 2 )
                        + pow ( MeshGID->point ( GlobalID ).y() - MeshNeighbors->point ( *it ).y(), 2 )
                        + pow ( MeshGID->point ( GlobalID ).z() - MeshNeighbors->point ( *it ).z(), 2 ) );
        r_max = ( r > r_max ) ? r : r_max;
    }
    return r_max;
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::buildRhs()
{
    M_RhsF.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsOne.reset (new vector_Type (*M_interpolationOperatorMap) );

    M_RhsF->subset (*M_knownField, *M_interpolationOperatorMap, 0, 0);
    *M_RhsOne += 1;
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::interpolateCostantField()
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

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::interpolate()
{
    vectorPtr_Type gamma_f;
    gamma_f.reset (new vector_Type (*M_interpolationOperatorMap) );

    // Preconditioner
    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
    precPtr.reset ( precRawPtr );

    LinearSolver solverRBF;
    solverRBF.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solverRBF.setParameters ( *M_belosList );
    solverRBF.setPreconditioner ( precPtr );

    solverRBF.setOperator (M_interpolationOperator);
    solverRBF.setRightHandSide (M_RhsF);
    solverRBF.solve (gamma_f);

    vectorPtr_Type rbf_f;
    rbf_f.reset (new vector_Type (*M_projectionOperatorMap) );

    vectorPtr_Type solution;
    solution.reset (new vector_Type (*M_projectionOperatorMap) );

    M_projectionOperator->multiply (false, *gamma_f, *rbf_f);

    *solution = *rbf_f;
    *solution /= *M_rbf_one;

    M_unknownField_rbf.reset (new vector_Type (M_unknownField->map() ) );
    M_unknownField_rbf->subset (*rbf_f, *M_projectionOperatorMap, 0, 0);

    M_unknownField->subset (*solution, *M_projectionOperatorMap, 0, 0);
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::identifyNodes (meshPtr_Type LocalMesh, boost::unordered_set<ID>& GID_nodes, vectorPtr_Type CheckVector)
{
    if (M_flags[0] == -1)
    {
        for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
            if (CheckVector->blockMap().LID (LocalMesh->point (i).id() ) != -1)
            {
                GID_nodes.insert (LocalMesh->point (i).id() );
            }
    }
    else
    {
        for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
            if ( this->isInside (LocalMesh->point (i).markerID(), M_flags) )
                if (CheckVector->blockMap().LID (LocalMesh->point (i).id() ) != -1)
                {
                    GID_nodes.insert (LocalMesh->point (i).id() );
                }
    }
}

template <typename mesh_Type>
bool RBFlocallyRescaledScalar<mesh_Type>::isInside (ID pointMarker, flagContainer_Type flags)
{
    int check = 0;
    for (UInt i = 0; i < flags.size(); ++i)
        if (pointMarker == flags[i])
        {
            ++check;
        }
    return (check > 0) ? true : false;
}


template <typename mesh_Type>
double RBFlocallyRescaledScalar<mesh_Type>::rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius)
{
    double distance = sqrt ( pow (x1 - x2, 2) + pow (y1 - y2, 2) + pow (z1 - z2, 2) );
    return pow (1 - distance / radius, 4) * (4 * distance / radius + 1);
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::updateRhs(vectorPtr_Type newRhs)
{
    *M_RhsF *= 0;
    *M_RhsF = *newRhs;
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::solution (vectorPtr_Type& Solution)
{
    Solution = M_unknownField;
}

template <typename mesh_Type>
void RBFlocallyRescaledScalar<mesh_Type>::solutionrbf (vectorPtr_Type& Solution_rbf)
{
    Solution_rbf = M_unknownField_rbf;
}

//! Factory create function
template <typename mesh_Type>
inline RBFInterpolation<mesh_Type> * createRBFlocallyRescaledScalar()
{
    return new RBFlocallyRescaledScalar< mesh_Type > ();
}
namespace
{
static bool S_registerTriLRS = RBFInterpolation<LifeV::RegionMesh<LinearTriangle > >::InterpolationFactory::instance().registerProduct ( "RBFlocallyRescaledScalar", &createRBFlocallyRescaledScalar<LifeV::RegionMesh<LinearTriangle > > );
static bool S_registerTetLRS = RBFInterpolation<LifeV::RegionMesh<LinearTetra > >::InterpolationFactory::instance().registerProduct ( "RBFlocallyRescaledScalar", &createRBFlocallyRescaledScalar<LifeV::RegionMesh<LinearTetra > > );
}

} // Namespace LifeV

#endif /* RBFLOCALLYRESCALEDSCALAR_H */
