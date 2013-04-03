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

#ifndef RBFSCALAR_H
#define RBFSCALAR_H 1

#include <lifev/core/interpolation/RBFInterpolation.hpp>

namespace LifeV {

template <typename mesh_Type>
class RBFscalar: public RBFInterpolation<mesh_Type>
{
public:

    typedef boost::shared_ptr<mesh_Type>                                          meshPtr_Type;

    typedef VectorEpetra                                                          vector_Type;
    typedef boost::shared_ptr<vector_Type >                                       vectorPtr_Type;

    typedef MatrixEpetra<double>                                                  matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                                        matrixPtr_Type;

    typedef std::vector<int>                                                      flagContainer_Type;

    typedef std::set<ID>                                                          idContainer_Type;

    typedef MapEpetra                                                             map_Type;
    typedef boost::shared_ptr<MapEpetra>                                          mapPtr_Type;

    typedef GhostHandler<mesh_Type>                                               neighbors_Type;
    typedef boost::shared_ptr<neighbors_Type>                                     neighborsPtr_Type;

    typedef LifeV::Preconditioner                                                 basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>                                      basePrecPtr_Type;

    typedef LifeV::PreconditionerIfpack                                           prec_Type;
    typedef boost::shared_ptr<prec_Type>                                          precPtr_Type;

    typedef Teuchos::RCP< Teuchos::ParameterList >                                parameterList_Type;

    RBFscalar();

    virtual ~RBFscalar();

    void setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags );

    void setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList);

    void buildOperators();

    void interpolationOperator();

    void projectionOperator();

    void buildRhs();

    void identifyNodes (meshPtr_Type LocalMesh, std::set<ID>& GID_nodes, vectorPtr_Type CheckVector);

    bool isInside (ID pointMarker, flagContainer_Type Flags);

    double rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius);

    void interpolate();

    void solution (vectorPtr_Type& Solution);

    void updateRhs(vectorPtr_Type newRhs);

    void setRadius ( double radius );

    void setBasis (const std::string & basis);

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
    std::string         M_basis;
    int                 M_globalNodesNumber;
};

template <typename mesh_Type>
RBFscalar<mesh_Type>::RBFscalar():
    M_radius(0)
{}

template <typename mesh_Type>
RBFscalar<mesh_Type>::~RBFscalar()
{}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::setup ( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags )
{
    M_fullMeshKnown = fullMeshKnown;
    M_localMeshKnown = localMeshKnown;
    M_fullMeshUnknown = fullMeshUnknown;
    M_localMeshUnknown = localMeshUnknown;
    M_flags = flags;
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList)
{
    M_knownField   = KnownField;
    M_unknownField = UnknownField;
    M_datafile     = datafile;
    M_belosList    = belosList;
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::setRadius ( double radius )
{
    M_radius = radius;
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::setBasis ( const std::string & basis )
{
    M_basis = basis;
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::buildOperators()
{
    LifeChrono TimeBuilding;
    TimeBuilding.start();
    this->interpolationOperator();
    TimeBuilding.stop();
    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "Time to assembly operators = " << TimeBuilding.diff() << std::endl;

    LifeChrono TimeBuildingp;
    TimeBuildingp.start();
    this->projectionOperator();
    TimeBuildingp.stop();
    if(M_knownField->mapPtr()->commPtr()->MyPID()==0)
        std::cout << "Time to assembly operators = " << TimeBuildingp.diff() << std::endl;

    this->buildRhs();

}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::interpolationOperator()
{
    ASSERT(M_radius!=0, "Please set the basis radius using RBFscalar<mesh_Type>::setRadius(double radius)");

    if(M_basis=="TPS"||M_basis=="IMQ"){
        this->identifyNodes (M_localMeshKnown, M_GIdsKnownMesh, M_knownField);
        int LocalNodesNumber = M_GIdsKnownMesh.size();
        M_globalNodesNumber = 0;
        MPI_Allreduce(&LocalNodesNumber,&M_globalNodesNumber,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        int* ElementsPerRow = new int[LocalNodesNumber];
        int* GlobalID = new int[LocalNodesNumber];
        int k = 0;

        for (std::set<ID>::iterator it = M_GIdsKnownMesh.begin(); it != M_GIdsKnownMesh.end(); ++it)
        {
            GlobalID[k] = *it;
            ElementsPerRow[k] = M_globalNodesNumber;
            ++k;
        }

        M_interpolationOperatorMap.reset (new map_Type (-1, LocalNodesNumber, GlobalID, M_knownField->mapPtr()->commPtr() ) );
        M_interpolationOperator.reset (new matrix_Type (*M_interpolationOperatorMap, ElementsPerRow) );

        int* Indices = new int[M_globalNodesNumber];
        double* Values = new double[M_globalNodesNumber];

        for ( int i = 0 ; i < LocalNodesNumber; ++i )
        {
            k = 0;
            for ( int j = 0; j < M_fullMeshKnown->numVertices(); ++j)
            {
                if(isInside(M_fullMeshKnown->point (j).markerID(),M_flags))
                {
                    Indices[k] = j;
                    Values[k]  = rbf ( M_fullMeshKnown->point (GlobalID[i]).x(),
                                       M_fullMeshKnown->point (GlobalID[i]).y(),
                                       M_fullMeshKnown->point (GlobalID[i]).z(),
                                       M_fullMeshKnown->point (j).x(),
                                       M_fullMeshKnown->point (j).y(),
                                       M_fullMeshKnown->point (j).z(),
                                       M_radius);
                    ++k;
                }
            }
            M_interpolationOperator->matrixPtr()->InsertGlobalValues (GlobalID[i], k, Values, Indices);
        }
        M_interpolationOperator->globalAssemble();
        delete Indices;
        delete Values;
        delete ElementsPerRow;
        delete GlobalID;
    }
    if(M_basis=="BW"){
        this->identifyNodes (M_localMeshKnown, M_GIdsKnownMesh, M_knownField);
        M_neighbors.reset ( new neighbors_Type ( M_fullMeshKnown, M_localMeshKnown, M_knownField->mapPtr(), M_knownField->mapPtr()->commPtr() ) );
        if (M_flags[0] == -1)
        {
            M_neighbors->setUp();
        }
        else
        {
            M_neighbors->setUp (M_flags);
        }

        int LocalNodesNumber = M_GIdsKnownMesh.size();

        std::vector<std::set<ID> > MatrixGraph (LocalNodesNumber);
        int* ElementsPerRow = new int[LocalNodesNumber];
        int* GlobalID = new int[LocalNodesNumber];
        int k = 0;
        int Max_entries = 0;

        for (std::set<ID>::iterator it = M_GIdsKnownMesh.begin(); it != M_GIdsKnownMesh.end(); ++it)
        {
            GlobalID[k] = *it;
            MatrixGraph[k] = M_neighbors->createNodeNodeNeighborsMapWithinRadius (M_radius, GlobalID[k]);
            MatrixGraph[k].insert (GlobalID[k]);
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
            for ( std::set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
            {
                Indices[k] = *it;
                Values[k]  = rbf ( M_fullMeshKnown->point (GlobalID[i]).x(),
                                   M_fullMeshKnown->point (GlobalID[i]).y(),
                                   M_fullMeshKnown->point (GlobalID[i]).z(),
                                   M_fullMeshKnown->point (*it).x(),
                                   M_fullMeshKnown->point (*it).y(),
                                   M_fullMeshKnown->point (*it).z(),
                                   M_radius);
                ++k;
            }
            M_interpolationOperator->matrixPtr()->InsertGlobalValues (GlobalID[i], k, Values, Indices);
        }
        M_interpolationOperator->globalAssemble();
        delete Indices;
        delete Values;
        delete ElementsPerRow;
        delete GlobalID;
    }
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::projectionOperator()
{
    this->identifyNodes (M_localMeshUnknown, M_GIdsUnknownMesh, M_unknownField);
    int LocalNodesNumber = M_GIdsUnknownMesh.size();

    if(M_basis=="TPS"||M_basis=="IMQ"){

        int* ElementsPerRow = new int[LocalNodesNumber];
        int* GlobalID = new int[LocalNodesNumber];
        int k = 0;

        for (std::set<ID>::iterator it = M_GIdsUnknownMesh.begin(); it != M_GIdsUnknownMesh.end(); ++it)
        {
            GlobalID[k] = *it;
            ElementsPerRow[k] = M_globalNodesNumber;
            ++k;
        }

        M_projectionOperatorMap.reset (new map_Type (-1, LocalNodesNumber, GlobalID, M_unknownField->mapPtr()->commPtr() ) );
        M_projectionOperator.reset (new matrix_Type (*M_projectionOperatorMap, ElementsPerRow) );

        int* Indices = new int[M_globalNodesNumber];
        double* Values = new double[M_globalNodesNumber];

        for ( int i = 0 ; i < LocalNodesNumber; ++i )
        {
            k = 0;
            for ( int j = 0; j < M_fullMeshKnown->numVertices(); ++j)
            {
                if(isInside(M_fullMeshKnown->point(j).markerID(),M_flags))
                {
                    Indices[k] = j;
                    Values[k]  = rbf ( M_fullMeshUnknown->point (GlobalID[i]).x(),
                                       M_fullMeshUnknown->point (GlobalID[i]).y(),
                                       M_fullMeshUnknown->point (GlobalID[i]).z(),
                                       M_fullMeshKnown->point (j).x(),
                                       M_fullMeshKnown->point (j).y(),
                                       M_fullMeshKnown->point (j).z(),
                                       M_radius);
                    ++k;
                }
            }
            M_projectionOperator->matrixPtr()->InsertGlobalValues (GlobalID[i], k, Values, Indices);
        }
        M_projectionOperator->globalAssemble (M_interpolationOperatorMap, M_projectionOperatorMap);
        delete Indices;
        delete Values;
        delete ElementsPerRow;
        delete GlobalID;
    }

    if(M_basis=="BW"){

        std::vector<std::set<ID> > MatrixGraph (LocalNodesNumber);
        int* ElementsPerRow = new int[LocalNodesNumber];
        int* GlobalID = new int[LocalNodesNumber];
        int k = 0;
        int Max_entries = 0;
        double d;
        double d_min;
        int nearestPoint;

        for (std::set<ID>::iterator it = M_GIdsUnknownMesh.begin(); it != M_GIdsUnknownMesh.end(); ++it)
        {
            GlobalID[k] = *it;
            d_min = 100;
            for (int j = 0; j <  M_fullMeshKnown->numVertices(); ++j)
            {
                if ( M_flags[0] == -1 || this->isInside (M_fullMeshKnown->point (j).markerID(), M_flags) )
                {
                    d = std::sqrt ( pow (M_fullMeshKnown->point (j).x() - M_fullMeshUnknown->point (GlobalID[k]).x(), 2)
                                    + pow (M_fullMeshKnown->point (j).y() - M_fullMeshUnknown->point (GlobalID[k]).y(), 2)
                                    + pow (M_fullMeshKnown->point (j).z() - M_fullMeshUnknown->point (GlobalID[k]).z(), 2) );
                    if (d < d_min)
                    {
                        d_min = d;
                        nearestPoint = M_fullMeshKnown->point (j).id();
                    }
                }
            }

            MatrixGraph[k] = M_neighbors->createNodeNodeNeighborsMapWithinRadius (M_radius, nearestPoint);
            MatrixGraph[k].insert (nearestPoint);
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
            for ( std::set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
            {
                Indices[k] = *it;
                Values[k]  = rbf ( M_fullMeshUnknown->point (GlobalID[i]).x(),
                                   M_fullMeshUnknown->point (GlobalID[i]).y(),
                                   M_fullMeshUnknown->point (GlobalID[i]).z(),
                                   M_fullMeshKnown->point (*it).x(),
                                   M_fullMeshKnown->point (*it).y(),
                                   M_fullMeshKnown->point (*it).z(),
                                   M_radius);
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

}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::buildRhs()
{
    M_RhsF.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF->subset (*M_knownField, *M_interpolationOperatorMap, 0, 0);
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::interpolate()
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

    /*
    vectorPtr_Type solution;
    solution.reset (new vector_Type (*M_projectionOperatorMap) );

    M_projectionOperator->multiply (false, *gamma_f, *solution);
    M_unknownField->subset (*solution, *M_projectionOperatorMap, 0, 0);
    */
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::identifyNodes (meshPtr_Type LocalMesh, std::set<ID>& GID_nodes, vectorPtr_Type CheckVector)
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
bool RBFscalar<mesh_Type>::isInside (ID pointMarker, flagContainer_Type flags)
{
    int check = 0;
    if(flags[0]==-1)
        return true;
    else
    {
        for (UInt i = 0; i < flags.size(); ++i)
            if (pointMarker == flags[i])
            {
                ++check;
            }
        return (check > 0) ? true : false;
    }
}

template <typename mesh_Type>
double RBFscalar<mesh_Type>::rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius)
{
    double distance = sqrt ( pow (x1 - x2, 2) + pow (y1 - y2, 2) + pow (z1 - z2, 2) );

    if(M_basis=="BW")
        return pow (1 - distance / radius, 4) * (4 * distance / radius + 1);
    else if(M_basis=="TPS")
        if (distance == 0)
            return 0;
        else
            return abs(distance/radius)*abs(distance/radius)*log(distance/radius);
    else if(M_basis=="IMQ")
        return 1/sqrt( abs(distance)*abs(distance) + radius*radius);

}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::updateRhs(vectorPtr_Type newRhs)
{
    *M_RhsF *= 0;
    *M_RhsF = *newRhs;
}

template <typename mesh_Type>
void RBFscalar<mesh_Type>::solution (vectorPtr_Type& Solution)
{
    Solution = M_unknownField;
}

//! Factory create function
template <typename mesh_Type>
inline RBFInterpolation<mesh_Type> * createRBFscalar()
{
    return new RBFscalar< mesh_Type > ();
}
namespace
{
static bool S_registerTriS = RBFInterpolation<LifeV::RegionMesh<LinearTriangle > >::InterpolationFactory::instance().registerProduct ( "RBFscalar", &createRBFscalar<LifeV::RegionMesh<LinearTriangle > > );
static bool S_registerTetS = RBFInterpolation<LifeV::RegionMesh<LinearTetra > >::InterpolationFactory::instance().registerProduct ( "RBFscalar", &createRBFscalar<LifeV::RegionMesh<LinearTetra > > );
}

} // Namespace LifeV

#endif /* RBFSCALAR_H */
