/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
Date: 2004-10-26

Copyright (C) 2004 Politecnico di Milano

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file main_twotubes.cpp
   \author Vincent Martin <vincent.martin@mate.polimi.it>
   \date 2004-10-26
*/

#include <life/lifecore/life.hpp>
#include <life/lifecore/chrono.hpp>
#include "dataOneDModel.hpp"
#include "oneDModelSolver.hpp"
#include "interface2Vessels.hpp"
#include "ud_functions.hpp"
#include <life/lifecore/GetPot.hpp>
#include <sstream>


int main(int argc, char** argv)
{
    using namespace LifeV;

    //! ********** Reading from data file ******************************************
    GetPot command_line(argc,argv);
    /*const char* data_file_name = command_line.follow("data", 2, "-f","--file");
      GetPot data_file(data_file_name);
    */


    GetPot data_file( "datanl" );

    OneDNonLinModelParam onedparamNL( data_file );
    onedparamNL.initParam(4e6);

    GetPot data_file_t2( "datanl2" );
    OneDNonLinModelParam onedparamNL_t2( data_file_t2 );
    onedparamNL_t2.initParam(8e6);


    LinearSimpleParam onedparamLin(data_file);
    LinearSimpleParam onedparamLin_t2(data_file_t2);

    std::cout << "======\n\tNon Linear model tube 1" << std::endl;
    onedparamNL.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;
    std::cout << "======\n\tNon Linear model tube 2" << std::endl;
    onedparamNL_t2.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;
    /*
      std::cout << "======\n\tLinear model tube 1" << std::endl;
      onedparamLin.showMeData(std::cout);
      std::cout << "-----------------------------" << std::endl;
      std::cout << "======\n\tLinear model tube 2" << std::endl;
      onedparamLin_t2.showMeData(std::cout);
      std::cout << "-----------------------------" << std::endl;
    */

    //  OneDModelSolver onedm(data_file, onedparamLin);
    OneDModelSolver onedm   (data_file, onedparamNL);
    OneDModelSolver onedm_t2(data_file_t2, onedparamNL_t2);

    std::cout << "======\n\ttube 1" << std::endl;
    onedm.showMeData();
    std::cout << "======\n\ttube 2" << std::endl;
    onedm_t2.showMeData();
    //  onedm.showMeHandler(cout, 6);

    Interface2Vessels interf_t1_t2( onedm, onedm_t2 );

    // Initialization
    //
    Real dt     = onedm.timestep();
    Real startT = onedm.inittime();
    Real T      = onedm.endtime();

    Real dt_t2     = onedm_t2.timestep();
    Real startT_t2 = onedm_t2.inittime();
    Real T_t2      = onedm_t2.endtime();

    ASSERT_PRE( dt == dt_t2, "Same time step required!");
    ASSERT_PRE( startT == startT_t2, "Same initial time required!");
    ASSERT_PRE( T == T_t2, "Same final time required!");

    ASSERT( onedm.xRight() == onedm_t2.xLeft(),
            "Same interface point required!");

    /*
      Real u1_0 = 0.; //! constant initial condition
      Real u2_0 = 0.; //! constant initial condition
    */
    //! tube 1
    Real u1_0 = onedparamNL.Area0(0); //! constant initial condition
    Real u2_0 = 0.; //! constant initial condition

    std::cout << "initialize tube 1 with constant (u1_0, u2_0)" << std::endl;
    onedm.initialize(u1_0, u2_0);

    //! tube 2
    u1_0 = onedparamNL_t2.Area0(0); //! constant initial condition
    u2_0 = 0.; //! constant initial condition

    std::cout << "initialize tube 2 with constant (u1_0, u2_0)" << std::endl;
    onedm_t2.initialize(u1_0, u2_0);

    std::cout << "startT T dt " << startT << " " <<  T << " " << dt << std::endl;

    char ch;
    std::cout << "Hit return to continue" << std::endl;
    std::cin.get(ch);

    std::cout << "++++++++++++++++++++++++++++++\n\tTemporal loop starting ... " << std::endl;
    // Temporal loop
    //
    Chrono chrono;
    int count = 0;
    for ( Real time = startT + dt ; time <= T ; time += dt ) {
        count++;
        std::cout << "Iteration " <<  count  << ", t = " << time
                  << "s... \n";
        chrono.start();

        //! compute the interface values
        interf_t1_t2.updateInterface2Vessels( onedm, onedm_t2 );

        int cvg_newton = interf_t1_t2.computeInterface2TubesValues();
        Vector bcDir_t1 = interf_t1_t2.BcDir_alpha();
        Vector bcDir_t2 = interf_t1_t2.BcDir_beta();

        std::cout << "bcDir_t1 " << bcDir_t1   << "\nbcDir_t2 " << bcDir_t2  << std::endl;

        //! set the interface values
        onedm.setBCValuesRight( bcDir_t1[0], bcDir_t1[1] );
        onedm_t2.setBCValuesLeft( bcDir_t2[0], bcDir_t2[1] );

        ASSERT_PRE( !cvg_newton,"Newton iteration for interface values computation not achieved.");

        //! tube 1
        onedm.timeAdvance( time );
        onedm.iterate( time , count );

        //! tube 2
        onedm_t2.timeAdvance( time );
        onedm_t2.iterate( time , count );

        chrono.stop();
        std::cout << "Iteration " <<  count  << " computed in " << chrono.diff() << " s.\n" << std::endl;

        if ( data_file( "miscellaneous/show_graceplot", 0 ) )
            onedm.gplot();

    }

    std::cout << "Hit return to close" << std::endl;
    std::cin.get(ch);

    return 0;

}



