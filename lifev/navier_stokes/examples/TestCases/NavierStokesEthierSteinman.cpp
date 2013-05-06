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
    @file NavierStokesEthierSteinman.cpp
    @brief Navier-Stokes Ethier-Steinman problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 19-08-2011
 */

#include <vector>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/fem/BCFunction.hpp>
#include <lifev/navier_stokes/function/RossEthierSteinmanDec.hpp>
#include <lifev/navier_stokes/function/RossEthierSteinmanInc.hpp>
#include <lifev/navier_stokes/examples/TestCases/NavierStokesEthierSteinman.hpp>
#include <lifev/navier_stokes/examples/TestCases/MeshUtility.hpp>

namespace LifeV
{

// Walls
const Int NavierStokesEthierSteinman::LEFTWALL   = 4;
const Int NavierStokesEthierSteinman::RIGHTWALL  = 2;
const Int NavierStokesEthierSteinman::FRONTWALL  = 1;
const Int NavierStokesEthierSteinman::BACKWALL   = 3;
const Int NavierStokesEthierSteinman::TOPWALL    = 6;
const Int NavierStokesEthierSteinman::BOTTOMWALL = 5;
// Edges
const Int NavierStokesEthierSteinman::BOTTOMEDGE1 =  7;
const Int NavierStokesEthierSteinman::BOTTOMEDGE2 =  8;
const Int NavierStokesEthierSteinman::BOTTOMEDGE3 =  9;
const Int NavierStokesEthierSteinman::BOTTOMEDGE4 = 10;
const Int NavierStokesEthierSteinman::SIDEEDGE1   = 11;
const Int NavierStokesEthierSteinman::SIDEEDGE2   = 12;
const Int NavierStokesEthierSteinman::SIDEEDGE3   = 13;
const Int NavierStokesEthierSteinman::SIDEEDGE4   = 14;
const Int NavierStokesEthierSteinman::TOPEDGE1    = 15;
const Int NavierStokesEthierSteinman::TOPEDGE2    = 16;
const Int NavierStokesEthierSteinman::TOPEDGE3    = 17;
const Int NavierStokesEthierSteinman::TOPEDGE4    = 18;
// Corners
const Int NavierStokesEthierSteinman::BOTTOMCORNER1 = 19;
const Int NavierStokesEthierSteinman::BOTTOMCORNER2 = 20;
const Int NavierStokesEthierSteinman::BOTTOMCORNER3 = 21;
const Int NavierStokesEthierSteinman::BOTTOMCORNER4 = 22;
const Int NavierStokesEthierSteinman::TOPCORNER1    = 23;
const Int NavierStokesEthierSteinman::TOPCORNER2    = 24;
const Int NavierStokesEthierSteinman::TOPCORNER3    = 25;
const Int NavierStokesEthierSteinman::TOPCORNER4    = 26;

NavierStokesEthierSteinman::NavierStokesEthierSteinman()
    : NavierStokesProblem()
{
    RossEthierSteinmanUnsteadyDec::setA ( 1.0 );
    RossEthierSteinmanUnsteadyDec::setD ( 1.0 );
    RossEthierSteinmanUnsteadyDec::setViscosity ( M_viscosity );
    RossEthierSteinmanUnsteadyDec::setDensity ( M_density );
}

NavierStokesEthierSteinman::~NavierStokesEthierSteinman()
{

}

bool
NavierStokesEthierSteinman::hasExactSolution() const
{
    return true;
}

NavierStokesProblem::function_Type
NavierStokesEthierSteinman::xexact()
{
    return RossEthierSteinmanUnsteadyDec::xexact;
}

NavierStokesProblem::function_Type
NavierStokesEthierSteinman::uexact()
{
    return RossEthierSteinmanUnsteadyDec::uexact;
}

NavierStokesProblem::function_Type
NavierStokesEthierSteinman::uderexact()
{
    return RossEthierSteinmanUnsteadyDec::uderexact;
}

NavierStokesProblem::function_Type
NavierStokesEthierSteinman::pexact()
{
    return RossEthierSteinmanUnsteadyDec::pexact;
}

void
NavierStokesEthierSteinman::showMe ( std::ostream& output ) const
{
    output << "To be implemented";
}

void
NavierStokesEthierSteinman::setViscosity ( const Real& viscosity )
{
    NavierStokesProblem::setViscosity ( viscosity );
    RossEthierSteinmanUnsteadyDec::setViscosity ( viscosity );
}

void
NavierStokesEthierSteinman::setDensity ( const Real& density )
{
    NavierStokesProblem::setDensity ( density );
    RossEthierSteinmanUnsteadyDec::setDensity ( density );
}

void
NavierStokesEthierSteinman::mesh ( boost::shared_ptr< RegionMesh<LinearTetra> >& mesh ) const
{
    UInt numMeshElem ( M_refinement );
    if ( M_refinement == 0 )
    {
        std::cout << "ERROR: The mesh refinement requested is not valid." << std::endl;
        exit ( 0 );
    }
    MeshUtility::fillWithStructuredMesh ( mesh,
                                          1,
                                          numMeshElem, numMeshElem, numMeshElem,
                                          false,
                                          2.0,   2.0,   2.0,
                                          -1.0,  -1.0,  -1.0 );
}

void
NavierStokesEthierSteinman::boundaryConditions ( boost::shared_ptr<BCHandler> bcHandler ) const
{
    BCFunctionBase uDirichlet ( RossEthierSteinmanUnsteadyDec::uexact );
    BCFunctionBase uNeumann  ( RossEthierSteinmanUnsteadyDec::fNeumann );

    bcHandler->addBC ( "Flux", FRONTWALL, Natural, Full, uNeumann, 3 );

    bcHandler->addBC ( "Top"   , TOPWALL      , LifeV::Essential        , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPEDGE1     , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPEDGE2     , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPEDGE3     , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPEDGE4     , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPCORNER1   , LifeV::EssentialVertices, Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPCORNER2   , LifeV::EssentialVertices, Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPCORNER3   , LifeV::EssentialVertices, Full, uDirichlet, 3 );
    bcHandler->addBC ( "Top"   , TOPCORNER4   , LifeV::EssentialVertices, Full, uDirichlet, 3 );
    bcHandler->addBC ( "Left"  , LEFTWALL     , LifeV::Essential        , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Left"  , SIDEEDGE1    , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Left"  , SIDEEDGE4    , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Right" , RIGHTWALL    , LifeV::Essential        , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Right" , SIDEEDGE2    , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Right" , SIDEEDGE3    , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Back"  , BACKWALL     , LifeV::Essential        , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMWALL   , LifeV::Essential        , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE1  , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE2  , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE3  , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMEDGE4  , LifeV::EssentialEdges   , Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER1, LifeV::EssentialVertices, Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER2, LifeV::EssentialVertices, Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER3, LifeV::EssentialVertices, Full, uDirichlet, 3 );
    bcHandler->addBC ( "Bottom", BOTTOMCORNER4, LifeV::EssentialVertices, Full, uDirichlet, 3 );
}

std::string
NavierStokesEthierSteinman::name() const
{
    return "Ethier-Steinman problem";
}

} // namespace LifeV
