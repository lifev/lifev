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
    @date 01 Mar 2013

    A more detailed description of the file (if necessary)
 */

#ifndef RBFINTERPOLATIONSTANDARD_H
#define RBFINTERPOLATIONSTANDARD_H 1

#include <lifev/core/LifeV.hpp>
#include <Epetra_Vector.h>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/GhostHandler.hpp>
#include <Epetra_FECrsMatrix.h>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

namespace LifeV
{

//! RBFInterpolationStandard - Short description of the class
/*!
    @author Davide Forti
    @see Reference to papers (if available)

    Here write a long and detailed description of the class.

    For this purpose you can use a lot of standard HTML code.
    Here there is a list with some useful examples:

    For bold text use: <b>BOLD TEXT</b>

    For empatyze a word type @e followed by the word

    For verbatim a word type @c followed by the word

    For vertical space (empty lines) use: <br>

    For creating list type:
    <ol>
        <li> First element of the enumerated list
        <ul>
             <li> First element of the dotted sublist.
             <li> Second element of the dotted sublist
        </ul>
        <li> Second element of the enumerated list
        <li> Third element of the enumerated list
    </ol>

    For writing a warning type: @warning followed by the description
    of the warning

    It is possible to use a lot of other standard HTML commands.
    Visit http://www.stack.nl/~dimitri/doxygen/htmlcmds.html for
    a detailed list.

    For any other kind of information visit www.doxygen.org.
 */

template <typename Mesh>
class RBFInterpolationStandard
{

public:

    typedef Mesh                                              mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                      meshPtr_Type;

    typedef VectorEpetra                                      vector_Type;
    typedef boost::shared_ptr<vector_Type >                   vectorPtr_Type;

    typedef Epetra_FECrsMatrix                                matrixEpetra_Type;
    typedef boost::shared_ptr<matrixEpetra_Type>              matrixEpetraPtr_Type;

    typedef MatrixEpetra<double>                              matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                    matrixPtr_Type;

    typedef std::vector<int>                                  flagContainer_Type;

    typedef std::set<ID>                                      idContainer_Type;

    typedef MapEpetra                                         map_Type;
    typedef boost::shared_ptr<MapEpetra>                      mapPtr_Type;

    typedef Epetra_Map                                        mapEpetra_Type;
    typedef boost::shared_ptr<mapEpetra_Type>                 mapEpetraPtr_Type;

    typedef GhostHandler<mesh_Type>                           neighbors_Type;
    typedef boost::shared_ptr<neighbors_Type>                 neighborsPtr_Type;

    typedef LifeV::Preconditioner                             basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>                  basePrecPtr_Type;

    typedef LifeV::PreconditionerIfpack                       prec_Type;
    typedef boost::shared_ptr<prec_Type>                      precPtr_Type;

    typedef Teuchos::RCP< Teuchos::ParameterList >            parameterList_Type;

    //! Constructor
    RBFInterpolationStandard ( meshPtr_Type fullMeshKnown,
                               meshPtr_Type localMeshKnown,
                               meshPtr_Type fullMeshUnknown,
                               meshPtr_Type localMeshUnknown,
                               flagContainer_Type flags,
                               double radius);

    //! Destructor
    ~RBFInterpolationStandard() {}


    //! Setup the RBF data
    void setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, Teuchos::RCP< Teuchos::ParameterList > belosList);


    //! Build the RBF Operators, namely the interpolant and the projection ones.
    void buildOperators();

    //! Build the RBF interpolant
    void InterpolationOperator();

    //! Identifies nodes with an assigned markerID
    void identifyNodes (meshPtr_Type LocalMesh, std::set<ID>& GID_nodes, vectorPtr_Type CheckVector);

    //! Check if the point with markerID pointMarker has to be selected
    bool isInside (ID pointMarker, flagContainer_Type Flags);

    //! Evaluate the RBF radius
    double computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID);

    //! Evaluation of the RBF
    double rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius);

    //! Build the projection operator
    void ProjectionOperator();

    //! Prepare Rhs
    void buildRhs();

    //! Manages the solution of the interpolation problem
    void interpolate();

    //! Getter for the solution
    void solution (vectorPtr_Type& Solution);

    //! Getter for the solution obtained by a pure RBF approach
    void solutionrbf (vectorPtr_Type& Solution_rbf);

private:

    meshPtr_Type                M_fullMeshKnown;
    meshPtr_Type                M_localMeshKnown;
    meshPtr_Type                M_fullMeshUnknown;
    meshPtr_Type                M_localMeshUnknown;
    matrixPtr_Type              M_interpolationOperator;
    matrixPtr_Type              M_projectionOperator;
    flagContainer_Type          M_flags;
    vectorPtr_Type              M_knownField;
    vectorPtr_Type              M_unknownField;
    idContainer_Type            M_GIdsKnownMesh;
    idContainer_Type            M_GIdsUnknownMesh;
    vectorPtr_Type              M_RhsF;
    vectorPtr_Type              M_RhsOne;
    mapEpetraPtr_Type           M_interpolationOperatorEpetraMap;
    mapEpetraPtr_Type           M_projectionOperatorEpetraMap;
    mapPtr_Type                 M_interpolationOperatorMap;
    mapPtr_Type                 M_projectionOperatorMap;
    neighborsPtr_Type           M_neighbors;
    vectorPtr_Type              M_unknownField_rbf;
    GetPot                      M_datafile;
    parameterList_Type          M_belosList;

    double                      M_radius;
};

template <typename Mesh>
RBFInterpolationStandard<Mesh>::RBFInterpolationStandard ( meshPtr_Type fullMeshKnown,
                                                           meshPtr_Type localMeshKnown,
                                                           meshPtr_Type fullMeshUnknown,
                                                           meshPtr_Type localMeshUnknown,
                                                           flagContainer_Type flags,
                                                           double radius) :
    M_fullMeshKnown ( fullMeshKnown ),
    M_localMeshKnown ( localMeshKnown ),
    M_fullMeshUnknown ( fullMeshUnknown ),
    M_localMeshUnknown ( localMeshUnknown ),
    M_flags ( flags ),
    M_radius ( radius )
{
}


template <typename Mesh>
void RBFInterpolationStandard<Mesh>::setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList)
{
    M_knownField   = KnownField;
    M_unknownField = UnknownField;
    M_datafile     = datafile;
    M_belosList    = belosList;
}


template <typename Mesh>
void RBFInterpolationStandard<Mesh>::buildOperators()
{
    this->InterpolationOperator();
    this->ProjectionOperator();
    this->buildRhs();
}

template <typename Mesh>
void RBFInterpolationStandard<Mesh>::InterpolationOperator()
{
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
    int TotalNodesNumber = 0;

    MPI_Allreduce (&LocalNodesNumber, &TotalNodesNumber, 1,  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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

    M_interpolationOperatorMap.reset (new map_Type (TotalNodesNumber, LocalNodesNumber, GlobalID, M_knownField->mapPtr()->commPtr() ) );
    M_interpolationOperatorEpetraMap.reset (new mapEpetra_Type (TotalNodesNumber, LocalNodesNumber, GlobalID, 0, * (M_knownField->mapPtr()->commPtr() ) ) );

    matrixEpetraPtr_Type InterpolationOperator;
    InterpolationOperator.reset (new matrixEpetra_Type (Copy, *M_interpolationOperatorEpetraMap, ElementsPerRow) );

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
            //RBF_radius[i]);
            ++k;
        }
        InterpolationOperator->InsertGlobalValues (GlobalID[i], k, Values, Indices);
    }
    InterpolationOperator->FillComplete();

    delete Indices;
    delete Values;
    delete ElementsPerRow;
    delete GlobalID;

    M_interpolationOperator.reset (new matrix_Type (*M_interpolationOperatorMap, InterpolationOperator) );
    M_interpolationOperator->spy ("M_interpolationOperator");

}

template <typename Mesh>
void RBFInterpolationStandard<Mesh>::ProjectionOperator()
{

    this->identifyNodes (M_localMeshUnknown, M_GIdsUnknownMesh, M_unknownField);

    int LocalNodesNumber = M_GIdsUnknownMesh.size();
    int TotalNodesNumber = 0;

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

    MPI_Allreduce (&LocalNodesNumber, &TotalNodesNumber, 1,  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    M_projectionOperatorMap.reset (new map_Type (TotalNodesNumber, LocalNodesNumber, GlobalID, M_unknownField->mapPtr()->commPtr() ) );
    M_projectionOperatorEpetraMap.reset (new mapEpetra_Type (TotalNodesNumber, LocalNodesNumber, GlobalID, 0, * (M_unknownField->mapPtr()->commPtr() ) ) );

    matrixEpetraPtr_Type ProjectionOperator;
    ProjectionOperator.reset (new matrixEpetra_Type (Copy, *M_projectionOperatorEpetraMap, ElementsPerRow) );

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
            //RBF_radius[i]);
            ++k;
        }
        ProjectionOperator->InsertGlobalValues (GlobalID[i], k, Values, Indices);
    }
    ProjectionOperator->FillComplete (*M_interpolationOperatorEpetraMap, *M_projectionOperatorEpetraMap);

    delete Indices;
    delete Values;
    delete ElementsPerRow;
    delete GlobalID;

    M_projectionOperator.reset (new matrix_Type (*M_projectionOperatorMap, ProjectionOperator) );
    M_projectionOperator->spy ("M_projectionOperator.m");
}

template <typename Mesh>
double RBFInterpolationStandard<Mesh>::computeRBFradius (meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID)
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


template <typename Mesh>
void RBFInterpolationStandard<Mesh>::buildRhs()
{
    M_RhsF.reset (new vector_Type (*M_interpolationOperatorMap) );
    M_RhsF->subset (*M_knownField, *M_interpolationOperatorMap, 0, 0);
}


template <typename Mesh>
void RBFInterpolationStandard<Mesh>::interpolate()
{
    vectorPtr_Type gamma_f;
    gamma_f.reset (new vector_Type (*M_interpolationOperatorMap) );

    // Preconditioner
    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
    precPtr.reset ( precRawPtr );

    LinearSolver solver;
    solver.setCommunicator ( M_knownField->mapPtr()->commPtr() );
    solver.setParameters ( *M_belosList );
    solver.setPreconditioner ( precPtr );

    solver.setOperator (M_interpolationOperator);
    solver.setRightHandSide ( M_RhsF );
    solver.solve ( gamma_f );

    vectorPtr_Type solution;
    solution.reset (new vector_Type (*M_projectionOperatorMap) );

    M_projectionOperator->multiply (false, *gamma_f, *solution);

    M_unknownField->subset (*solution, *M_projectionOperatorMap, 0, 0);

}


template <typename Mesh>
void RBFInterpolationStandard<Mesh>::identifyNodes (meshPtr_Type LocalMesh, std::set<ID>& GID_nodes, vectorPtr_Type CheckVector)
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

template <typename Mesh>
bool RBFInterpolationStandard<Mesh>::isInside (ID pointMarker, flagContainer_Type flags)
{
    int check = 0;
    for (UInt i = 0; i < flags.size(); ++i)
        if (pointMarker == flags[i])
        {
            ++check;
        }
    return (check > 0) ? true : false;
}


template <typename Mesh>
double RBFInterpolationStandard<Mesh>::rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius)
{
    double distance = sqrt ( pow (x1 - x2, 2) + pow (y1 - y2, 2) + pow (z1 - z2, 2) );
    return pow (1 - distance / radius, 4) * (4 * distance / radius + 1);
    //return exp ( -pow(distance,2) / 2*pow(radius,2));
}

template <typename Mesh>
void RBFInterpolationStandard<Mesh>::solution (vectorPtr_Type& Solution)
{
    Solution = M_unknownField;
}


template <typename Mesh>
void RBFInterpolationStandard<Mesh>::solutionrbf (vectorPtr_Type& Solution_rbf)
{
    Solution_rbf = M_unknownField_rbf;
}

} // namespace LifeV

#endif // RBF_INTERPOLATIONSTANDARD_HPP
