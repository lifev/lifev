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

#include <lifeV.hpp>

#include <boost/timer.hpp>

#include <FSISolver.hpp>

#include "ud_functions.hpp"

typedef boost::shared_ptr<LifeV::FSISolver> fsi_solver_ptr;

/*!
  This routine sets up the problem:

  -# create the standard boundary conditions for the fluid and
     structure problems.

  -# initialize and setup the FSIsolver
 */
fsi_solver_ptr
fsi_setup( GetPot const& data_file )
{
    using namespace LifeV;

    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    BCFunctionBase bcf(fZero);
    BCHandler BCh_mesh;
    BCh_mesh.addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_mesh.addBC("Edges",    20, Essential, Full, bcf,   3);

    // Boundary conditions for the fluid velocity
    BCFunctionBase in_flow(u2);
    BCHandler BCh_u;
    BCh_u.addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_u.addBC("Edges",  20, Essential, Full, bcf,     3);

    // Boundary conditions for the solid displacement
    //    BCh_d.addBC("Interface", 1, Natural, Full, g_wall, 3);
    BCHandler BCh_d;
    BCh_d.addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_d.addBC("Base",      2, Essential, Full, bcf,  3);

    fsi_solver_ptr fsi(  new FSISolver( data_file, BCh_u, BCh_d, BCh_mesh ) );
    fsi->showMe();
    fsi->setSourceTerms( fZero, fZero );
    fsi->initialize( u0, d0, w0 );

    return fsi;
}
/*!
  This routine runs the temporal loop
*/
void
fsi_run( fsi_solver_ptr __fsi )
{
    double dt     = __fsi->timeStep();
    double T      = __fsi->timeEnd();

    boost::timer __overall_timer;

    int __i = 1;
    for (double time=dt; time <= T; time+=dt, ++__i)
    {
        boost::timer __timer;

        __fsi->iterate( time );

        std::cout << "[fsi_run] Iteration " << __i << " was done in : "
                  << __timer.elapsed() << "\n";
    }
    std::cout << "Total computation time = "
              << __overall_timer.elapsed() << "s" << "\n";
}

/*
  This routine checks that the results given by the different
  FSIOperator are actually the same.
*/
void
fsi_check( GetPot const& data_file)
{
    fsi_solver_ptr __fsi = fsi_setup( data_file );

    __fsi->setFSIOperator( "steklovPoincare" );
    __fsi->FSIOperator()->setPreconditioner( LifeV::DIRICHLET_NEUMANN );

    fsi_run( __fsi );

    __fsi->setFSIOperator( "exactJacobian" );

    fsi_run( __fsi );

}
int main(int argc, char** argv)
{
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);

    fsi_solver_ptr __fsi = fsi_setup( data_file );
    fsi_run( __fsi );

    return 0;
}

