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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H 1

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <lifev/core/fem/CurrentFEManifold.hpp>

namespace LifeV
{

class Interpolation
{
public:

	typedef RegionMesh<LinearTetra> mesh_Type;
    typedef std::shared_ptr<mesh_Type> meshPtr_Type;
    typedef VectorEpetra vector_Type;
    typedef std::shared_ptr<vector_Type> vectorPtr_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef std::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef MapEpetra map_Type;
    typedef std::shared_ptr<MapEpetra> mapPtr_Type;
    typedef LifeV::Preconditioner basePrec_Type;
    typedef std::shared_ptr<basePrec_Type> basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack prec_Type;
    typedef std::shared_ptr<prec_Type> precPtr_Type;
    typedef Teuchos::RCP< Teuchos::ParameterList > parameterList_Type;
    typedef FESpace<mesh_Type, MapEpetra> FESpace_Type;
    typedef std::shared_ptr<FESpace_Type> FESpacePtr_Type;

    Interpolation();

    ~Interpolation();

    void setup( const GetPot& datafile, parameterList_Type belosList );

    void setFlag( const UInt& flag ) { M_flag = flag; }

    void buildTableDofs_known ( const FESpacePtr_Type& fespace );

    void buildTableDofs_unknown ( const FESpacePtr_Type& fespace );

    void identifyNodes_known ( );

    void identifyNodes_unknown ( );

    void setVectors (const vectorPtr_Type& KnownField, const vectorPtr_Type& UnknownField );

    void buildKnownInterfaceMap();

    void buildUnknownInterfaceMap();

    void buildOperators();

    void interpolationOperator();

    void projectionOperator();

    void expandGammaToOmega_Known(const vectorPtr_Type& vectorOnGamma, vectorPtr_Type& vectorOnOmega);

    void restrictOmegaToGamma_Known(const vectorPtr_Type& vectorOnOmega, vectorPtr_Type& vectorOnGamma);

    void interpolate();

    void solution (vectorPtr_Type& Solution);

    void getSolutionOnGamma (vectorPtr_Type& GammaSolution) { GammaSolution.reset(new vector_Type ( *M_solutionOnGamma ) ); };

    void updateRhs(const vectorPtr_Type& newRhs);

    void getInterpolationOperatorMap ( mapPtr_Type& map ) { map.reset(new map_Type(*M_interpolationOperatorMap) ); };

    void getprojectionOperatorMap (mapPtr_Type& map) { map.reset(new map_Type(*M_projectionOperatorMap)); };

    void getKnownInterfaceMap(mapPtr_Type& map){ map.reset ( new map_Type(*M_knownInterfaceMap) ); };

    void getNumerationInterfaceKnown(vectorPtr_Type& vector){ vector.reset ( new vector_Type ( *M_numerationInterfaceKnown ) ); };

    void getVectorialInterpolationMap ( mapPtr_Type& map ) { map.reset ( new map_Type(*M_interpolationOperatorMapVectorial) ); };

    void free_space ();

    void setMeshSize(const Real& mesh_size ) { M_links = (1.5*mesh_size); };

private:

    void buildInterpolationOperatorMap();

    void buildProjectionOperatorMap();

    void buildRhs();

    void interpolateCostantField();

    double computeRBFradius_known ( const ID& index, std::vector<int> indexes );

    double computeRBFradius_unknown ( const ID& index, std::vector<int> indexes );

    double rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius);

    std::vector<Real>   M_xcoord_known;
    std::vector<Real>   M_ycoord_known;
    std::vector<Real>   M_zcoord_known;
    std::vector<Real>   M_xcoord_unknown;
    std::vector<Real>   M_ycoord_unknown;
    std::vector<Real>   M_zcoord_unknown;

    std::vector<UInt>   M_marker_known;
    std::vector<UInt>   M_marker_unknown;

    UInt 				M_flag;

    std::set<ID>        M_GIdsKnownMesh_common;
    std::set<ID>        M_GIdsUnknownMesh_common;

    std::set<ID>        M_GIdsKnownMesh;
    std::set<ID>        M_GIdsUnknownMesh;

    vectorPtr_Type      M_knownField;
    vectorPtr_Type      M_unknownField;

    mapPtr_Type         M_knownInterfaceMap;
    mapPtr_Type         M_unknownInterfaceMap;

    vectorPtr_Type      M_numerationInterfaceKnown;
    vectorPtr_Type      M_numerationInterfaceUnknown;
    std::vector<vectorPtr_Type> M_numerationInterfaceKnownColumns;
    int M_pid;

    mapPtr_Type         M_interpolationOperatorMap;
    mapPtr_Type         M_interpolationOperatorMapVectorial;

    std::vector< std::vector<int> > M_dof_connectivity_known; // per ogni dof i suoi vicini
    std::vector< std::vector<int> > M_dof_connectivity_unknown;

    matrixPtr_Type      M_interpolationOperator;

    mapPtr_Type         M_projectionOperatorMap;
    mapPtr_Type         M_projectionOperatorMapVectorial;

    matrixPtr_Type      M_projectionOperator;

    vectorPtr_Type      M_RhsF1;
    vectorPtr_Type      M_RhsF2;
    vectorPtr_Type      M_RhsF3;
    vectorPtr_Type      M_RhsOne;
    vectorPtr_Type      M_rbf_one;

    GetPot              M_datafile;
    parameterList_Type  M_belosList;

    vectorPtr_Type      M_solutionOnGamma;
    
    Real                M_links;

    basePrecPtr_Type    M_precPtr;

    bool                M_precBuilt;
};

} // Namespace LifeV

#endif /* INTERPOLATION_H */
