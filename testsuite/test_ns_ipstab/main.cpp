/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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

/*!

  \author M.A. Fernandez
  \date 01/01/2004

  \brief Main program for solving he Ossen's equation with reaction using
         P1/P1 or P2/P2 FEM with Interior Penalty (IP) stabilization
         (See Burman-Fernandez-Hansbo 2004).

*/

#include "lifeV.hpp"
#include "NavierStokesSolverIP.hpp"
#include "ud_functions.hpp"
#include "picard.hpp"
#include "vectorNorms.hpp"

int main(int argc, char** argv)
{
    using namespace LifeV;

    // Reading from data file
    //
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);

    // Number of boundary conditions for the fluid velocity,
    // solid displacement, and fluid mesh motion
    //
    BC_Handler BCh(1);

    NavierStokesSolverIP< RegionMesh3D<LinearTetra> > fluid(data_file, feTetraP1,
                                                            quadRuleTetra4pt,
                                                            quadRuleTria3pt, BCh);
    fluid.showMe();

    // Boundary conditions for the fluid velocity
    //vector<ID> icomp(1);
    //icomp[0]=1;
    BCFunction_Base u_wall(uexact);
    // BCFunction_Base bcf(fZero);
    BCh.addBC("Wall",   1,  Essential, Full, u_wall,  3);
    //BCh.addBC("Wall",   2,  Essential, Component, bcf,  icomp);
    //BCFunction_Base bcf(afZero);
    //BCFunction_Base in_flow(u2);
    //BCh.addBC("Wall",   1,  Essential, Full, bcf,  3);
    //BCh.addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    //BCh.addBC("Edges",  20, Essential, Full, bcf,     3);
    //========================================================================================
    //  TEMPORAL LOOP
    //========================================================================================
    //
    BCh.showMe();



    // Picard-Aitken iterations: steady version
    //

    Real abstol=1.e6;
    Real reltol=0.0;
    int maxiter=300;
    Real omega=0;

    UInt dim_u = fluid.uDof().numTotalDof();

    Vector fx1(3*dim_u), fx0(3*dim_u), gx1(3*dim_u), gx0(3*dim_u);
    Vector x1(3*dim_u), x0(3*dim_u);

    fx1 = 0.0;
    fx0 = 0.0;
    gx1 = 0.0;
    gx0 = 0.0;
    x1  = 0.0;
    x0  = 0.0;

    // Compute right hand side
    fluid.initialize(u0);
    fluid.timeAdvance(f,0.0);

    int  status = picard(&fluid,maxnorm,fx1,fx0,gx1,gx0,
                         x1,x0,abstol,reltol,maxiter,1,omega);

    if(status == 1) {
        cout << "Inners iterations failed\n";
        exit(1);
    }  else {
        // cout << "End of time "<< time << endl;
        cout << "Number of inner iterations       : " << maxiter << endl;
        fluid.postProcess();
    }

    return 0;

}
