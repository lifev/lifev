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
    @brief Quadrature Rule test

	@author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @author Umberto Villa <uvilla@emory.edu>
    @contributor
    @maintainer Umberto Villa <uvilla@emory.edu>

    @date 02-03-2010

The program tests the degree of exactness (DoE) and the convergence rate (CR)
of all the quadrature rule in 3D (Tetrahedra) or in 2D (Triangles).

The code produce the following output:

quadRuleTetra.txt ==> Show the Degree of Exactness of all the quadrature rules on Tetrahedral elements
quadRuleTetra.plt ==> Show the Convergence Rate of all the quadrature rules on Tetrahedral elements
                               using gnuplot
 */

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <life/lifecore/GetPot.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/quadRule.hpp>
#include <life/lifefem/CurrentFE.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifemesh/structuredMesh3D.hpp>
#include <life/lifefem/DOF.hpp>

#include "test_quadrule.hpp"

int main(int argc, char** argv)
{
    using namespace LifeV;

#ifdef HAVE_MPI
    std::cout << "MPI Initialization" << std::endl;
    MPI_Init( &argc, &argv );
#endif

    bool checkGlobal(true), check(true);

    // All the quadrature rules for tetrahedra
    container_Type allQuadRuleTetra;
    allQuadRuleTetra.reserve(5);
    allQuadRuleTetra.push_back(&quadRuleTetra1pt);
    allQuadRuleTetra.push_back(&quadRuleTetra4pt);
    allQuadRuleTetra.push_back(&quadRuleTetra5pt);
    allQuadRuleTetra.push_back(&quadRuleTetra15pt);
    allQuadRuleTetra.push_back(&quadRuleTetra64pt);

    // All the quadrature rules for triangles
    container_Type allQuadRuleTria;
    allQuadRuleTria.reserve(5);
    allQuadRuleTria.push_back(&quadRuleTria1pt);
    allQuadRuleTria.push_back(&quadRuleTria3pt);
    allQuadRuleTria.push_back(&quadRuleTria4pt);
    allQuadRuleTria.push_back(&quadRuleTria6pt);
    allQuadRuleTria.push_back(&quadRuleTria7pt);

    // All the quadrature rules for segments
    container_Type allQuadRuleSegments;
    allQuadRuleSegments.reserve(3);
    allQuadRuleSegments.push_back(&quadRuleSeg1pt);
    allQuadRuleSegments.push_back(&quadRuleSeg2pt);
    allQuadRuleSegments.push_back(&quadRuleSeg3pt);

    // Check Quadrature Rule on Tetrahedra
    check = quad_check_doe< RegionMesh3D<LinearTetra> >(feTetraP1, geoLinearTetra, allQuadRuleTetra, "quadRuleTetra");
    checkGlobal = checkGlobal & check;
    check = quad_check_cr< RegionMesh3D<LinearTetra> >(feTetraP1, geoLinearTetra, allQuadRuleTetra, "quadRuleTetra");
    checkGlobal = checkGlobal & check;

    // Check the quadrature rules for triangles
    for (constIterator_Type it(allQuadRuleTria.begin()); it != allQuadRuleTria.end(); ++it)
    {
        check = (*it)->degreeOfExactness() == (*it)->checkExactness();
        checkGlobal = checkGlobal & check;
    }

    // Check the quadrature rules for triangles
    for (constIterator_Type it(allQuadRuleSegments.begin()); it != allQuadRuleSegments.end(); ++it)
    {
        check = (*it)->degreeOfExactness() == (*it)->checkExactness();
        checkGlobal = checkGlobal & check;
    }

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif


    if (check)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}//end main


