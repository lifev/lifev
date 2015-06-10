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
    @brief FSIData - File containing the implementation of Radial Basis Functions suited for interpolation
                     between non-matching grids

    @author Davide Forti <davide.forti@epfl.ch>
    @date 01-31-2010

    @maintainer Davide Forti <davide.Forti@epfl.ch>
 */

#ifndef RBF_INTERPOLATION_HPP
#define RBF_INTERPOLATION_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/GhostHandler.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

namespace LifeV
{
template <typename mesh_Type>
class RBFInterpolation
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

    typedef FactorySingleton<Factory<RBFInterpolation<mesh_Type>, std::string> >  InterpolationFactory;

    RBFInterpolation();

    virtual ~RBFInterpolation() {}

    virtual void setup( meshPtr_Type fullMeshKnown, meshPtr_Type localMeshKnown, meshPtr_Type fullMeshUnknown, meshPtr_Type localMeshUnknown, flagContainer_Type flags ) {};

    virtual void setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField, GetPot datafile, parameterList_Type belosList) {};

    virtual void setupRBFData (vectorPtr_Type KnownField, vectorPtr_Type UnknownField) {};

    virtual void buildOperators() {};

    virtual void interpolationOperator() {};

    virtual void projectionOperator() {};

    virtual void buildRhs() {};

    virtual void interpolateCostantField() {};

    virtual void identifyNodes (meshPtr_Type LocalMesh, std::set<ID>& GID_nodes, vectorPtr_Type CheckVector) {};

    virtual bool isInside (ID pointMarker, flagContainer_Type Flags) {};

    virtual double computeRBFradius (meshPtr_Type , meshPtr_Type , idContainer_Type , ID ) {};

    virtual void setBasis (const std::string &) {};

    virtual double rbf (double x1, double y1, double z1, double x2, double y2, double z2, double radius) {};

    virtual void interpolate() {};

    virtual void solution (vectorPtr_Type& Solution) {};

    virtual void solutionrbf (vectorPtr_Type & ) {};

    virtual void updateRhs(const vectorPtr_Type& ) {};

    virtual void setRadius ( double ){};

    virtual void approximateInverse ( ) {};

    virtual void getInterpolationOperatorMap(mapPtr_Type&){ };

    virtual void getprojectionOperatorMap(mapPtr_Type& ){ };

    virtual void buildUnknownVectorialInterfaceMap(){};

    // Methods added after changing the maps

    virtual void buildKnownInterfaceMap(){};

    virtual void buildUnknownInterfaceMap(){};
    
    virtual void buildInterpolationOperatorMap(){};
    
    virtual void buildProjectionOperatorMap(){};

    virtual void getKnownInterfaceMap(mapPtr_Type& map){};

    virtual void getNumerationInterfaceKnown(vectorPtr_Type& vector){};

    virtual void getSolutionOnGamma(vectorPtr_Type& ) { };

    virtual void expandGammaToOmega_Known(const vectorPtr_Type&, vectorPtr_Type& ) { };

    virtual void restrictOmegaToGamma_Known(const vectorPtr_Type&, vectorPtr_Type& ) { };

    virtual void getVectorialInterpolationMap ( mapPtr_Type& ) { };

private:

};

template <typename mesh_Type>
RBFInterpolation<mesh_Type>::RBFInterpolation()
{}


} // namespace LifeV

#endif // RBF_INTERPOLATION_HPP
