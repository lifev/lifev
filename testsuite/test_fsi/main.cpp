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
#include <cassert>

#include <lifeV.hpp>

#include <boost/timer.hpp>

#include <FSISolver.hpp>

#include "ud_functions.hpp"

class Problem
{
public:
    typedef boost::shared_ptr<LifeV::FSISolver> fsi_solver_ptr;

    /*!
      This routine sets up the problem:

      -# create the standard boundary conditions for the fluid and
      structure problems.

      -# initialize and setup the FSIsolver
    */
    Problem( GetPot const& data_file, std::string __oper = "" )
        {
            using namespace LifeV;

            // Boundary conditions for the harmonic extension of the
            // interface solid displacement
            Debug( 10000 ) << "Boundary condition for the harmonic extension\n";
            BCFunctionBase bcf(fZero);
            FSISolver::bchandler_type BCh_mesh( new BCHandler );
            BCh_mesh->addBC("Top",       3, Essential, Full, bcf,   3);
            BCh_mesh->addBC("Base",      2, Essential, Full, bcf,   3);
            BCh_mesh->addBC("Edges",    20, Essential, Full, bcf,   3);

            // Boundary conditions for the fluid velocity
            Debug( 10000 ) << "Boundary condition for the fluid\n";
            BCFunctionBase in_flow(u2);
            FSISolver::bchandler_type BCh_u( new BCHandler );
            BCh_u->addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
            BCh_u->addBC("Edges",  20, Essential, Full, bcf,     3);

            // Boundary conditions for the solid displacement
            Debug( 10000 ) << "Boundary condition for the solid\n";
            FSISolver::bchandler_type BCh_d( new BCHandler );
            BCh_d->addBC("Top",       3, Essential, Full, bcf,  3);
            BCh_d->addBC("Base",      2, Essential, Full, bcf,  3);

            Debug( 10000 ) << "creating FSISolver with operator :  " << __oper << "\n";
            _M_fsi = fsi_solver_ptr(  new FSISolver( data_file, BCh_u, BCh_d, BCh_mesh, __oper ) );
            _M_fsi->showMe();
            _M_fsi->setSourceTerms( fZero, fZero );
            _M_fsi->initialize( u0, d0, w0 );
        }

    fsi_solver_ptr fsiSolver() { return _M_fsi; }

    /*!
      This routine runs the temporal loop
    */
    void
    run( double dt, double T)
        {
            boost::timer __overall_timer;

            int __i = 1;
            for (double time=dt; time <= T; time+=dt, ++__i)
            {
                boost::timer __timer;

                _M_fsi->iterate( time );

                std::cout << "[fsi_run] Iteration " << __i << " was done in : "
                          << __timer.elapsed() << "\n";
            }
            std::cout << "Total computation time = "
                      << __overall_timer.elapsed() << "s" << "\n";

        }

private:

    fsi_solver_ptr _M_fsi;
};

LifeV::Vector
check( GetPot const& data_file,  std::string __oper, LifeV::OperFSPreconditioner __prec = LifeV::NO_PRECONDITIONER )
{
    try
    {
        Problem fsip( data_file, __oper );
        fsip.fsiSolver()->FSIOperator()->setPreconditioner( __prec );

        fsip.run( fsip.fsiSolver()->timeStep(), fsip.fsiSolver()->timeStep() ); // only one iteration

        return fsip.fsiSolver()->displacement();
    }
    catch ( std::exception const& __ex )
    {
        std::cout << "caught exception :  " << __ex.what() << "\n";
    }
}
int main(int argc, char** argv)
{
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);

    LifeV::Vector __ej_disp = check( data_file,  "fixedPoint" );
    LifeV::Vector __sp_disp = check( data_file,  "steklovPoincare",  LifeV::DIRICHLET_NEUMANN );

    std::cout << "norm_2(displacement error) = " << LifeV::norm_2( __sp_disp - __ej_disp ) << "\n";
    return 0;
}

