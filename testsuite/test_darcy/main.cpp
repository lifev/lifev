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
#include <iostream>
#include "GetPot.hpp"
#include "darcySolver.hpp"
#include "chrono.hpp"

/*

Darcy soler using Mixed Hybrid finite element


usage: darcy       : read the data file "data" and run the solver
darcy -f otherdata : read the data file "otherdata" and run the solver
darcy -h           : read the data file, print help and exit
darcy -i           : read the data file, print the read values and exit

*/
int main(int argc, char** argv)
{
    using namespace LifeV;
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);
    if( command_line.search(2, "-i","--info") ) {
        data_file.print();
        exit(0);
    }

    if( command_line.search( 2 , "-h" , "--help" ) ) {
        std::cout << std::endl << std::endl;
        std::cout <<"usage: darcy              : read the data file 'data' \n";
        std::cout <<"       darcy -f otherdata : read the data file 'otherdata' \n";
        std::cout <<"	    darcy -h           : help and exit\n";
        std::cout <<"       darcy -i           : : read the data file, print the read values and exit\n";
        std::cout << std::endl;
        exit(0);
    }

    enum MeshElementShape{ TetraElt, HexaElt };

    MeshElementShape elemshape = (MeshElementShape)
        data_file( "darcy/discretization/element_shape", TetraElt );

    std::cout << "*** shape " << elemshape << std::endl;


    Chrono chrono;
    //
    std::cout << "*** Initialisation --->" << std::endl;


    //-----------------------------------------------------
    if ( elemshape == HexaElt ) {
        //! case for hexahedric mesh
        chrono.start();

        DarcySolver< RegionMesh3D<LinearHexa> >
            darcyslv(data_file, feHexaRT0, feHexaQ0, feHexaRT0Hyb,
                     feHexaRT0VdotNHyb, feHexaQ1, quadRuleHexa8pt, quadRuleQuad4pt);

        chrono.stop();
        std::cout << "<--- Initialisation done in " << chrono.diff() << "s." << std::endl;

        if(darcyslv.verbose)
            std::cout << "*** Compute the matrix --->" << std::endl;
        chrono.start();
        darcyslv.computeHybridMatrixAndSourceRHS();
        chrono.stop();
        if(darcyslv.verbose)
            std::cout << "<--- matrix computation done in "<< chrono.diff() << "s." << std::endl;
        darcyslv.applyBC();
        //
        if(darcyslv.verbose) std::cout << "*** Resolution of the hybrid system --->\n";
        chrono.start();
        darcyslv.solveDarcy();
        chrono.stop();
        if(darcyslv.verbose)std::cout << "<--- Linear system solved in " << chrono.diff()
                                      << "s." << std::endl << std::endl;

        if(darcyslv.verbose) std::cout << "*** Compute pressure and flux --->" << std::endl;
        chrono.start();
        darcyslv.computePresFlux();
        chrono.stop();
        if(darcyslv.verbose)
            std::cout << "<---  done in " << chrono.diff() << "s.\n" << std::endl;

        if(darcyslv.verbose) std::cout << "*** Postproc --->" << std::endl;
        chrono.start();
        darcyslv.postProcessTraceOfPressureRT0();
        darcyslv.postProcessVelocityRT0();
        darcyslv.postProcessPressureQ0();
        darcyslv.postProcessPressureQ1();
        darcyslv.postProcessVelocityQ1();
        chrono.stop();
        if(darcyslv.verbose)
            std::cout << "<---  done in " << chrono.diff() << "s.\n" << std::endl;
        //-----------------------------------------------------

    }
    else if ( elemshape == TetraElt ) {
        //-----------------------------------------------------
        //! case for tetrahedric mesh
        DarcySolver< RegionMesh3D<LinearTetra> >
            darcyslv( data_file, feTetraRT0, feTetraP0, feTetraRT0Hyb, feTetraRT0VdotNHyb, 
                      feTetraP1, quadRuleTetra15pt, quadRuleTria4pt );
        if(darcyslv.verbose)
            std::cout << "*** Compute the matrix --->" << std::endl;
        chrono.start();
        darcyslv.computeHybridMatrixAndSourceRHS();
        chrono.stop();
        if(darcyslv.verbose)
            std::cout << "<--- matrix computation done in "<< chrono.diff() << "s." << std::endl;
        darcyslv.applyBC();
        //
        if(darcyslv.verbose) std::cout << "*** Resolution of the hybrid system --->\n";
        chrono.start();
        darcyslv.solveDarcy();
        chrono.stop();
        if(darcyslv.verbose)std::cout << "<--- Linear system solved in " << chrono.diff()
                                      << "s." << std::endl << std::endl;

        if(darcyslv.verbose) std::cout << "*** Compute pressure and flux --->" << std::endl;
        chrono.start();
        darcyslv.computePresFlux();
        chrono.stop();
        if(darcyslv.verbose)
            std::cout << "<---  done in " << chrono.diff() << "s.\n" << std::endl;

        if(darcyslv.verbose) std::cout << "*** Postproc --->" << std::endl;
        chrono.start();
        darcyslv.postProcessTraceOfPressureRT0();
        darcyslv.postProcessVelocityRT0();
        darcyslv.postProcessPressureQ0();
        darcyslv.postProcessPressureQ1();
        darcyslv.postProcessVelocityQ1();
        chrono.stop();
        if(darcyslv.verbose)
            std::cout << "<---  done in " << chrono.diff() << "s.\n" << std::endl;

    }
    //-----------------------------------------------------

    return 0;
}
