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

/*
  \author Daniele A. Di Pietro <dipietro@unibg.it>
  \date 2-6-2005
  \brief Two-fluid Navier-Stokes sample problem
*/

#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/elemMat.hpp>

#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bdf.hpp>

#include <life/lifealg/SolverAztec.hpp>

#include <life/lifesolver/LevelSetSolver.hpp>
#include <life/lifesolver/NSSolver2FluidsMixed.hpp>

#include <life/lifefilters/openDX_wrtrs.hpp>

#include <main.hpp>
#include <ud_functions.hpp>

/* #include <ethierSteinman.hpp> */

#define ONE_QUARTER 0
#define P1BUBBLE 0

int main() {
    using namespace LifeV;

    typedef RegionMesh3D<LinearTetra> mesh_type;

    Chrono chrono;

    // Retriving data from data file
    GetPot datafile("data");
    Real t0 = datafile("navier-stokes/time-discretization/inittime", 0.);
    Real T = datafile("navier-stokes/time-discretization/endtime", 1.);
    Real delta_t = datafile("navier-stokes/bdf/delta_t", 0.04);

    std::cout << "==================================================" << std::endl;
    std::cout << "t0         " << t0 << std::endl;
    std::cout << "T          " << T << std::endl;
    std::cout << "==================================================" << std::endl;

    // Boundary conditions

    //EthierSteinmanUnsteady::setParamsFromGetPot( datafile );
    //BCFunctionBase gv(EthierSteinmanUnsteady::uexact);
    BCFunctionBase gv(g);

    BCHandler BCh( 6, BCHandler::HINT_BC_ONLY_ESSENTIAL );
    BCh.addBC("Wall", 1, Essential, Full, gv, 3);
    BCh.addBC("Wall", 2, Essential, Full, gv, 3);
    BCh.addBC("Wall", 3, Essential, Full, gv, 3);
    BCh.addBC("Wall", 4, Essential, Full, gv, 3);
    BCh.addBC("Wall", 5, Essential, Full, gv, 3);
    BCh.addBC("Wall", 6, Essential, Full, gv, 3);

    //BCh.addBC("Wall", 10, Essential, Full, gv, 3);

    // Source term

    srcfct f;

    // Finite element stuff

#if P1BUBBLE
    const RefFE& refFE_u = feTetraP1bubble;
#else
    const RefFE& refFE_u = feTetraP2;
#endif

    const RefFE& refFE_p = feTetraP1;
    const RefFE& refFE_lss = feTetraP2;

    const QuadRule& qr = quadRuleTetra15pt;
    const QuadRule& qr_bd = quadRuleTria4pt;

    // Post-processing

    std::string ofile_root_ls = "./results/2f";
    std::string ofile_root_velocity = "./results/v";
    std::string ofile_root_pressure = "./results/p";

    // Solving the problem

    std::cout << "** NS2F test ** Initializing 2 fluid solver" << std::endl;

    NSSolver2FluidsMixed<mesh_type> NSS("data",
                                        refFE_u,
                                        refFE_p,
                                        BCh,
                                        refFE_lss,
                                        BCh,
                                        qr, qr_bd);

    NSS.initialize(zero, zero, sphere, t0, T);
    NSS.setVerbose();

    // Export initial conditions to OpenDX format

    wr_opendx_header(ofile_root_ls + "0000.dx", NSS.mesh(), NSS.lsDof(), NSS.fe_ls(), "P2" );
    wr_opendx_scalar(ofile_root_ls + "0000.dx", "levelset_ipstab", NSS.lsfunction());
#if P1BUBBLE
    wr_opendx_header(ofile_root_velocity + "0000.dx", NSS.mesh(), NSS.uDof(), NSS.fe_u(), "P1bubble");
#else
    wr_opendx_header(ofile_root_velocity + "0000.dx", NSS.mesh(), NSS.uDof(), NSS.fe_u(), "P2");
#endif
    wr_opendx_vector(ofile_root_velocity + "0000.dx", "u", NSS.velocity(), 3);

    std::cout << "** NS2F test ** Advancing in time" << std::endl;

    UInt current_step = 1;
    UInt steps_after_last_save = 1;
    UInt steps_after_last_reini = 1;

    UInt save_every = datafile("navier-stokes/miscellaneous/save_every", 1);
    UInt reini_every = datafile("navier-stokes/miscellaneous/reini_every", 5);

    for(Real t = t0; t < T; t += delta_t) {
        std::cout << "** NS2F test ** Time: " << t << std::endl;
        NSS.timeAdvance(gravity, 0.);

        // Export solution
        if(steps_after_last_save == save_every) {
            std::cout << "** NS2F test ** Exporting solution to OpenDX format" << std::endl;
            std::ostringstream number;
            number.width(4);
            number.fill('0');
            number << current_step;

            wr_opendx_header(ofile_root_ls + number.str() + ".dx", NSS.mesh(), NSS.lsDof(), NSS.fe_ls(), "P2" );
            wr_opendx_scalar(ofile_root_ls + number.str() + ".dx", "ls", NSS.lsfunction());

            wr_opendx_header(ofile_root_pressure + number.str() + ".dx", NSS.mesh(), NSS.pDof());
            wr_opendx_scalar(ofile_root_pressure + number.str() + ".dx", "p", NSS.pressure());
            
#if P1BUBBLE
            wr_opendx_header(ofile_root_velocity + number.str() + ".dx", NSS.mesh(), NSS.uDof(), NSS.fe_u(), "P1bubble");
#else
            wr_opendx_header(ofile_root_velocity + number.str() + ".dx", NSS.mesh(), NSS.uDof(), NSS.fe_u(), "P2");
#endif
            wr_opendx_vector(ofile_root_velocity + number.str() + ".dx", "u", NSS.velocity(), 3);

            steps_after_last_save = 1;
        } else
            steps_after_last_save++;

        // Re-initialize solution
        if(steps_after_last_reini == reini_every) {
            std::cout << "** NS2F test ** Reinitializing signed distance function" << std::endl;
            NSS.reinitialize("direct");
            steps_after_last_reini = 1;
        } else
            steps_after_last_reini++;

        current_step++;
    }

    std::cout << "** NS2F test ** Done" << std::endl;;

    return 0;
}
