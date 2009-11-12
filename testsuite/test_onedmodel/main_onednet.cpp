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
#include <lifemc/lifesolver/dataOneDModel.hpp>
#include <lifemc/lifesolver/oneDModelSolver.hpp>
#include "interface2Vessels.hpp"
#include <life/lifecore/GetPot.hpp>
#include <sstream>

using namespace LifeV;

int main(int argc, char** argv)
{


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::flush;

    MPI_Init(&argc,&argv);

    Epetra_Comm* comm = new Epetra_MpiComm( MPI_COMM_WORLD );

    //    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    std::cout << "ok" << std::endl;
#else
    Epetra_Comm* comm = new Epetra_SerialComm();
#endif;



    //! ********** Reading from data file ******************************************
    GetPot command_line(argc,argv);
    /*string data_file_name = command_line.follow("data", 2, "-f","--file");
      GetPot data_file(data_file_name);
    */
    typedef NonLinearFluxFun1D                            Flux1D;
    typedef NonLinearSourceFun1D                          Source1D;
    typedef OneDNonLinModelParam                          Params1D;
    typedef OneDModelSolver<Params1D, Flux1D, Source1D>   onedsolver_type;
    typedef Interface2Vessels<Params1D, Flux1D, Source1D> interface_type;

    GetPot data_file_t1                   ( "datanl" );
//     OneDNonLinModelParam onedparamNL   ( data_file );
//     onedparamNL.initParam              ( data_file );

    GetPot data_file_t2                ( "datanl2" );
//     OneDNonLinModelParam onedparamNL_t2( data_file_t2 );
//     onedparamNL_t2.initParam           ( data_file );

    const RefFE*    refFE = &feSegP1;
    const QuadRule* qR    = &quadRuleSeg3pt;
    const QuadRule* bdQr  = &quadRuleSeg1pt;

    std::cout << "======\n\tNon Linear model tube 1" << std::endl;

    DataOneDModel data_t1(data_file_t1);
    Params1D params_t1(data_file_t1);
    params_t1.showMeData(std::cout);

    std::cout << "    Building FE Space ... " << std::flush;
    FESpace<RegionMesh, EpetraMap> odFESpace_t1(data_t1.mesh(), *refFE, *qR, *bdQr, 1, *comm);
    std::cout << "ok." << std::endl;

    onedsolver_type onedm_t1(data_t1, params_t1, odFESpace_t1, *comm );

    std::cout << "-----------------------------" << std::endl;

    std::cout << "======\n\tNon Linear model tube 2" << std::endl;

    DataOneDModel data_t2(data_file_t2);
    Params1D params_t2(data_file_t2);
    params_t2.showMeData(std::cout);

    std::cout << "    Building FE Space ... " << std::flush;
    FESpace<RegionMesh, EpetraMap> odFESpace_t2(data_t2.mesh(), *refFE, *qR, *bdQr, 1, *comm);
    std::cout << "ok." << std::endl;

    onedsolver_type onedm_t2(data_t2, params_t2, odFESpace_t2, *comm );

    std::cout << "-----------------------------" << std::endl;

    

    interface_type interf_t1_t2( onedm_t1, onedm_t2 );

    // Initialization
    //
    Real dt_t1     = data_t1.timestep();
    Real startT_t1 = data_t1.inittime();
    Real T_t1      = data_t1.endtime();

    Real dt_t2     = data_t2.timestep();
    Real startT_t2 = data_t2.inittime();
    Real T_t2      = data_t2.endtime();

    ASSERT_PRE( dt_t1     == dt_t2,     "Same time step required!");
    ASSERT_PRE( startT_t1 == startT_t2, "Same initial time required!");
    ASSERT_PRE( T_t1      == T_t2,      "Same final time required!");

    ASSERT( onedm.xRight() == onedm_t2.xLeft(),
            "Same interface point required!");


    //  Real u1_0 = 0.; //! constant initial condition
    //  Real u2_0 = 0.; //! constant initial condition

    //! tube 1
    Real u1_0 = params_t1.Area0(0); //! constant initial condition
    Real u2_0 = 0.;                 //! constant initial condition

    std::cout << "initialize tube 1 with constant (u1_0, u2_0)" << std::endl;
    onedm_t1.initialize(u1_0, u2_0);

    //! tube 2
    u1_0 = params_t2.Area0(0); //! constant initial condition
    u2_0 = 0.; //! constant initial condition

    std::cout << "initialize tube 2 with constant (u1_0, u2_0)" << std::endl;
    onedm_t2.initialize(u1_0, u2_0);


    std::cout << "startT, T,  dt = " << startT_t1 << ", " <<  T_t1 << ", " << dt_t1 << std::endl;

    std::cout << "++++++++++++++++++++++++++++++\n\tTemporal loop starting ... " << std::endl;
    // Temporal loop
    //
    Chrono chrono;
    int count = 0;
    for ( Real time = startT_t1 + dt_t1 ; time <= T_t1 ; time += dt_t1 ) {

        count++;
        std::cout << "Iteration " <<  count  << ", t = " << time
                  << "s... \n";
        chrono.start();

        //! compute the interface values
        interf_t1_t2.updateInterface2Vessels( onedm_t1, onedm_t2 );

	std::cout << "OK" << std::endl;
        int cvg_newton  = interf_t1_t2.computeInterface2TubesValues();
        Vector bcDir_t1 = interf_t1_t2.BcDir_alpha();
        Vector bcDir_t2 = interf_t1_t2.BcDir_beta();
	std::cout << "OK" << std::endl;
        std::cout << "bcDir_t1 " << bcDir_t1   << "\nbcDir_t2 " << bcDir_t2  << std::endl;

        //! set the interface values
        onedm_t1.setBCValuesRight  ( bcDir_t1[0], bcDir_t1[1] );
        onedm_t2.setBCValuesLeft( bcDir_t2[0], bcDir_t2[1] );

        ASSERT_PRE( !cvg_newton,"Newton iteration for interface values computation not achieved.");

        //! tube 1
        onedm_t1.timeAdvance( time );
        onedm_t1.iterate( time , count );

        //! tube 2
        onedm_t2.timeAdvance( time );
        onedm_t2.iterate( time , count );

        chrono.stop();
        std::cout << "Iteration " <<  count  << " computed in " << chrono.diff() << " s.\n" << std::endl;

        if ( data_file_t1( "miscellaneous/show_graceplot", 0 ) )
        {
            onedm_t1.postProcess( time );
            onedm_t2.postProcess( time );
        }

    }
    

#ifdef EPETRA_MPI
    MPI_Finalize();
#endif

    return 0;





}



