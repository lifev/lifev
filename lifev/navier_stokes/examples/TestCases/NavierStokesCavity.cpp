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
    @file NavierStokesCavity.cpp
    @brief Navier-Stokes cavity problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 19-08-2011
 */

#include <vector>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/fem/BCFunction.hpp>
#include <lifev/navier_stokes/examples/TestCases/NavierStokesCavity.hpp>
#include <lifev/navier_stokes/examples/TestCases/MeshUtility.hpp>

namespace LifeV
{

// Walls
const Int NavierStokesCavity::LEFTWALL   = 4;
const Int NavierStokesCavity::RIGHTWALL  = 2;
const Int NavierStokesCavity::FRONTWALL  = 1;
const Int NavierStokesCavity::BACKWALL   = 3;
const Int NavierStokesCavity::TOPWALL    = 6;
const Int NavierStokesCavity::BOTTOMWALL = 5;
// Edges
const Int NavierStokesCavity::BOTTOMEDGE1 =  7;
const Int NavierStokesCavity::BOTTOMEDGE2 =  8;
const Int NavierStokesCavity::BOTTOMEDGE3 =  9;
const Int NavierStokesCavity::BOTTOMEDGE4 = 10;
const Int NavierStokesCavity::SIDEEDGE1   = 11;
const Int NavierStokesCavity::SIDEEDGE2   = 12;
const Int NavierStokesCavity::SIDEEDGE3   = 13;
const Int NavierStokesCavity::SIDEEDGE4   = 14;
const Int NavierStokesCavity::TOPEDGE1    = 15;
const Int NavierStokesCavity::TOPEDGE2    = 16;
const Int NavierStokesCavity::TOPEDGE3    = 17;
const Int NavierStokesCavity::TOPEDGE4    = 18;
// Corners
const Int NavierStokesCavity::BOTTOMCORNER1 = 19;
const Int NavierStokesCavity::BOTTOMCORNER2 = 20;
const Int NavierStokesCavity::BOTTOMCORNER3 = 21;
const Int NavierStokesCavity::BOTTOMCORNER4 = 22;
const Int NavierStokesCavity::TOPCORNER1    = 23;
const Int NavierStokesCavity::TOPCORNER2    = 24;
const Int NavierStokesCavity::TOPCORNER3    = 25;
const Int NavierStokesCavity::TOPCORNER4    = 26;

NavierStokesCavity::NavierStokesCavity()
    : NavierStokesProblem< RegionMesh< LinearTetra > >()
{

}

NavierStokesCavity::~NavierStokesCavity()
{

}

void
NavierStokesCavity::showMe ( std::ostream& output ) const
{
    output << "To be implemented";
}

Real
NavierStokesCavity::lidBC ( const LifeV::Real& /*t*/,
                            const LifeV::Real& /*x*/,
                            const LifeV::Real& /*y*/,
                            const LifeV::Real& /*z*/,
                            const LifeV::ID& i)
{
    switch (i)
    {
        case 1:
            return 1.0;
        default:
            return 0.0;
    }
}

Real
NavierStokesCavity::fZero ( const LifeV::Real& /* t */,
                            const LifeV::Real& /* x */,
                            const LifeV::Real& /* y */,
                            const LifeV::Real& /* z */,
                            const LifeV::ID& /* i */ )
{
    return 0.0;
}

void
NavierStokesCavity::mesh ( boost::shared_ptr< RegionMesh<LinearTetra> >& mesh ) const
{
    if ( M_refinement == 0 )
    {
        std::cout << "ERROR: The mesh refinement requested is not valid." << std::endl;
        exit ( 0 );
    }
    std::vector<Real> width (3, -1.0);
    std::vector<Real> shift (3, -1.0);
    std::vector<UInt> numMeshElem (3, M_refinement);
    MeshUtility::fillWithStructuredMesh ( mesh, 1, numMeshElem, false, width, shift );

}

void
NavierStokesCavity::boundaryConditions ( boost::shared_ptr<BCHandler> bcHandler ) const
{
    LifeV::BCFunctionBase uZero (fZero);
    LifeV::BCFunctionBase uLid (lidBC);

    std::vector<LifeV::ID> yComp (1);
    yComp[0] = 1;

    bcHandler->addBC ( "Top"   , TOPWALL      , LifeV::Essential        , LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPEDGE1     , LifeV::EssentialEdges   , LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPEDGE2     , LifeV::EssentialEdges   , LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPEDGE3     , LifeV::EssentialEdges   , LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPEDGE4     , LifeV::EssentialEdges   , LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPCORNER1   , LifeV::EssentialVertices, LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPCORNER2   , LifeV::EssentialVertices, LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPCORNER3   , LifeV::EssentialVertices, LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Top"   , TOPCORNER4   , LifeV::EssentialVertices, LifeV::Full     , uLid , 3     );
    bcHandler->addBC ( "Left"  , LEFTWALL     , LifeV::Essential        , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Left"  , SIDEEDGE1    , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Left"  , SIDEEDGE4    , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Front" , FRONTWALL    , LifeV::Essential        , LifeV::Component, uZero, yComp );
    bcHandler->addBC ( "Right" , RIGHTWALL    , LifeV::Essential        , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Right" , SIDEEDGE2    , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Right" , SIDEEDGE3    , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Back"  , BACKWALL     , LifeV::Essential        , LifeV::Component, uZero, yComp );
    bcHandler->addBC ( "Bottom", BOTTOMWALL   , LifeV::Essential        , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE1  , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE2  , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE3  , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE4  , LifeV::EssentialEdges   , LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER1, LifeV::EssentialVertices, LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER2, LifeV::EssentialVertices, LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER3, LifeV::EssentialVertices, LifeV::Full     , uZero, 3     );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER4, LifeV::EssentialVertices, LifeV::Full     , uZero, 3     );
}

std::string
NavierStokesCavity::name() const
{
    return "Lid-driven cavity";
}

} // namespace LifeV
