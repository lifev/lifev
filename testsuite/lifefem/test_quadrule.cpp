/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/* ========================================================

The program tests the degree of exactness (DoE) and the convergence rate (CR)
of all the quadrature rule in 3D (Tetrahedra and Hexahedra) or in 2D (Triangles, and soon
Quadrilaterals).

NOTE: If you want to add a new quadrature rule, it will be automatically tested. No modification of the test
is needed

The code produce the following output:

quadRule_ELEMENT_SHAPE.txt ==> Show the Degree of Exactness of all the quadrature rules on _ELEMENT_SHAPE.
                               _ELEMENT_SHAPE = Tria (for 2d),  Tetra Hexa (for 3d).
quadRule_ELEMENT_SHAPE.plt ==> Show the Convergence Rate of all the quadrature rules on _ELEMENT_SHAPE
                               using gnuplot

\author Umberto Villa <uvilla@emory.edu>
\date 02/03/2010
*/

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <life/lifecore/GetPot.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/quadRule.hpp>
#include <life/lifefem/currentFE.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifemesh/structuredMesh3D.hpp>
#include <life/lifefem/dof.hpp>

#include "test_quadrule.hpp"

int main(int argc, char** argv)
{
    using namespace LifeV;

    #ifdef HAVE_MPI
    std::cout << "MPI Initialization" << std::endl;
    MPI_Init( &argc, &argv );
	#endif

    bool check(true);
    GetPot command_line(argc,argv);

    std::vector<QuadRule*> allQuadRuleTetra;
    allQuadRuleTetra.reserve(5);
    allQuadRuleTetra.push_back(&quadRuleTetra1pt);
    allQuadRuleTetra.push_back(&quadRuleTetra4pt);
    allQuadRuleTetra.push_back(&quadRuleTetra5pt);
    allQuadRuleTetra.push_back(&quadRuleTetra15pt);
    allQuadRuleTetra.push_back(&quadRuleTetra64pt);

    // Check Quadrature Rule on Tetrahedra
    std::string data_file_name = command_line.follow("data", 2, "-f","--file");

    check = quad_check_doe< RegionMesh3D<LinearTetra> >(feTetraP2, geoLinearTetra, allQuadRuleTetra, "quadRuleTetra");

    quad_check_cr< RegionMesh3D<LinearTetra> >(feTetraP2, geoLinearTetra, allQuadRuleTetra, "quadRuleTetra");

    #ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
	#endif


    if(check)
    	return EXIT_SUCCESS;
    else
    	return EXIT_FAILURE;
}//end main


