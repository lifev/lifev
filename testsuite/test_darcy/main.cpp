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

#include <GetPot.hpp>
#include <darcySolverBase.hpp>
#include <chrono.hpp>

#include "user_diffusion.hpp"
#include "user_fct.hpp"

struct printError
{
    void
    operator()( std::string const& __norm,
                LifeV::darcy_unknown_type __t,
                double norm_l2_2,
                double norm_l2_sol_2,
                double norm_l2_diff_2 )
        {
            double normL2     = sqrt( norm_l2_2 );
            double normL2sol  = sqrt( norm_l2_sol_2 );
            double normL2diff = sqrt( norm_l2_diff_2 );

            switch( __t )
            {
                case LifeV::DARCY_PRESSURE_GLOBAL:
                    std::cout << "PRESSURE ERROR (Q0)" << "\n";
                    std::cout << "|| p        ||_{L^2}                   = " << normL2 << "\n";
                    std::cout << "|| p_ex     ||_{L^2}                   = " << normL2sol << "\n";
                    std::cout << "|| p - p_ex ||_{L^2}                   = " << normL2diff<< "\n";
                    std::cout << "|| P - p_ex ||_{L^2} / || p_ex ||_{L^2} = "
                              << normL2diff / normL2sol << "\n" << "\n";

                    std::cout << "SQUARE of PRESSURE ERROR (Q0)" << "\n";
                    std::cout << "|| p       ||^2_{L^2}                   = " << norm_l2_2 << "\n";
                    std::cout << "|| p_ex     ||^2_{L^2}                   = " << norm_l2_sol_2 << "\n";
                    std::cout << "|| p - p_ex ||^2_{L^2}                   = " << norm_l2_diff_2 << "\n";
                    std::cout << "|| P - p_ex ||^2_{L^2} / || p_ex ||^2_{L^2} = "
                              << norm_l2_diff_2 / norm_l2_sol_2 << "\n";
                    break;
                case LifeV::DARCY_PRESSURE:
                    std::cout << "PRESSURE ERROR (Q1)" << "\n";

                    std::cout << "|| P         ||_{" << __norm << "}                   = " << normL2 << "\n";
                    std::cout << "|| exact     ||_{" << __norm << "}                   = " << normL2sol << "\n";
                    std::cout << "|| P - exact ||_{" << __norm << "}                   = " << normL2diff<< "\n";
                    std::cout << "|| P - exact ||_{" << __norm << "}"
                              << "/"
                              << "|| exact ||_{" << __norm << "} = " << normL2diff/normL2sol
                              << "\n";
                    break;
                case LifeV::DARCY_VELOCITY:
                    std::cout << "VELOCITY ERROR (Q1)" << "\n";

                    std::cout << "|| U         ||_{" << __norm << "}                   = " << normL2 << "\n";
                    std::cout << "|| exact     ||_{" << __norm << "}                   = " << normL2sol << "\n";
                    std::cout << "|| U - exact ||_{" << __norm << "}                   = " << normL2diff<< "\n";
                    std::cout << "|| U - exact ||_{" << __norm << "}"
                              << "/"
                              << "|| exact ||_{" << __norm << "} = " << normL2diff/normL2sol
                              << "\n";
                    break;
            }

        }

};

void
checkAnalytical()
{
    using namespace LifeV;

    //! boundary conditions functions
    //! on a cube one might use them as follows:
    BCFunctionBase bcFct1;  //!< low  X
    BCFunctionBase bcFct2;  //!< high X
    BCFunctionBase bcFct3;  //!< low  Y
    BCFunctionBase bcFct4;  //!< high Y
    BCFunctionBase bcFct5;  //!< low  Z
    BCFunctionBase bcFct6;  //!< high Z

    // number of boundary conditions
    int nb_bc;
    // boundary conditions handler
    BCHandler bc;

    /*
      Analytical solution defined in user_fct

      example of mesh: hexahexa10x10x10.mesh

    */
    bcFct1.setFunction(zero); //!< low   X
    bcFct2.setFunction(zero); //!< high  X
    bcFct3.setFunction(zero); //!< low   Y
    bcFct4.setFunction(zero); //!< high  Y
    bcFct5.setFunction(zero); //!< low   Z
    bcFct6.setFunction(zero); //!< high  Z

    nb_bc = 6;
    bc.setNumber(nb_bc);

    bc.addBC("Analytical, real BC",    1, Essential,   Scalar, bcFct1);
    bc.addBC("Analytical, real BC",    2, Essential,   Scalar, bcFct2);
    bc.addBC("Analytical, real BC",    3, Essential,   Scalar, bcFct3);
    bc.addBC("Analytical, real BC",    4, Essential,   Scalar, bcFct4);
    bc.addBC("Analytical, real BC",    5, Essential,   Scalar, bcFct5);
    bc.addBC("Analytical, real BC",    6, Essential,   Scalar, bcFct6);

    boost::shared_ptr<DarcySolverBase> __darcy( FactoryDarcy::instance().createObject( "darcy_hexa" ) );

    __darcy->setBC( bc );
    __darcy->setSourceTerm( SourceAnalyticalFct() );
    __darcy->setup();
    __darcy->solve();
    __darcy->doOnErrorComputation( printError() );
    __darcy->errorL2( LifeV::DARCY_PRESSURE_GLOBAL, AnalyticalSolPres() );
    __darcy->errorL2( LifeV::DARCY_PRESSURE, AnalyticalSolPres() );
    __darcy->errorL2( AnalyticalSolFlux() );

}
/*

Darcy soler using Mixed Hybrid finite element


usage: darcy       : read the data file "data" and run the solver
darcy -f otherdata : read the data file "otherdata" and run the solver
darcy -h           : read the data file, print help and exit
darcy -i           : read the data file, print the read values and exit

*/
int main(int argc, char** argv)
{
    checkAnalytical();
}


//
//
//
#if 0
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

        if(darcyslv.verbose) std::cout << "*** Error wrt analytical solution --->" << std::endl;
        darcyslv.doOnErrorComputation( printError() );
        darcyslv.errorL2( LifeV::DARCY_PRESSURE_GLOBAL, AnalyticalSolPres() );

        if(darcyslv.verbose) std::cout << "*** Postproc --->" << std::endl;
        chrono.start();
        darcyslv.postProcessTraceOfPressureRT0();
        darcyslv.postProcessVelocityRT0();
        darcyslv.postProcessPressureQ0();
        darcyslv.postProcessPressureQ1();
        darcyslv.postProcessVelocityQ1();
        darcyslv.postProcessEnsight();
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

        if(darcyslv.verbose) std::cout << "*** Error wrt analytical solution --->" << std::endl;
        darcyslv.doOnErrorComputation( printError() );
        darcyslv.errorL2( LifeV::DARCY_PRESSURE_GLOBAL, AnalyticalSolPres() );
        darcyslv.errorL2( LifeV::DARCY_PRESSURE, AnalyticalSolPres() );
        darcyslv.errorL2( AnalyticalSolFlux() );

        if(darcyslv.verbose) std::cout << "*** Postproc --->" << std::endl;
        chrono.start();
        darcyslv.postProcessTraceOfPressureRT0();
        darcyslv.postProcessVelocityRT0();
        darcyslv.postProcessPressureQ0();
        darcyslv.postProcessPressureQ1();
        darcyslv.postProcessVelocityQ1();
        darcyslv.postProcessEnsight();
        chrono.stop();
        if(darcyslv.verbose)
            std::cout << "<---  done in " << chrono.diff() << "s.\n" << std::endl;

    }
    //-----------------------------------------------------
    return 0;
}
#endif

//! boundary conditions functions

#if 0
    //! on a cube one might use them as follows:
    BCFunctionBase bcFct1;  //!< low  X
    BCFunctionBase bcFct2;  //!< high X
    BCFunctionBase bcFct3;  //!< low  Y
    BCFunctionBase bcFct4;  //!< high Y
    BCFunctionBase bcFct5;  //!< low  Z
    BCFunctionBase bcFct6;  //!< high Z

    BCFunctionMixte bcFct_rob;  //!< a mixte (or Robin) bc function

    int nb_bc;                //!< number of boundary conditions
    BCHandler bc;            //!< boundary conditions handler

    // define the boundary conditions
    switch( data_file( "darcy/test_case" ,33 ) )
    {
        case 1:
            /*
              Neumann condition (on p) at inlet (ref 1) and outlet (ref 3)
              example of mesh: cylhexa.mesh
            */
            bcFct1.setFunction(g1);
            bcFct2.setFunction(g3);
            nb_bc = 2;
            bc.setNumber(nb_bc);
            bc.addBC("Inlet",      1, Natural,   Scalar, bcFct1);
            bc.addBC("Outlet",  3, Natural, Scalar, bcFct2);
            break;

        case 2:
            /*
              Dirichlet condition (on p) at inlet (ref 1) and outlet (ref 3)
              example of mesh: cylhexa.mesh
            */
            bcFct1.setFunction(g1);
            bcFct2.setFunction(g3);
            nb_bc = 2;
            bc.setNumber(nb_bc);
            bc.addBC("Inlet",      1, Essential,   Scalar, bcFct1);
            bc.addBC("Outlet",  3, Essential, Scalar, bcFct2);
            break;

        case 3:
            /*
              Robin condition (dp/dq + alpha p = 1, alpha=1) at inlet (ref 1)
              and Dirichlet (p=-1) at outlet (ref 3)
              example of mesh: cylhexa.mesh
            */
            bcFct_rob.setFunctions_Mixte(g1, mixte_coeff); //! Robin coeff = 1.
            bcFct2.setFunction(g3);
            nb_bc = 2;
            bc.setNumber(nb_bc);
            bc.addBC("Inlet",   1,     Mixte, Scalar, bcFct_rob);
            bc.addBC("Outlet",  3, Essential, Scalar, bcFct2);
            break;

        case 33:
            /*
              Analytical solution defined in user_fct

              example of mesh: hexahexa10x10x10.mesh

            */
            bcFct1.setFunction(zero); //!< low   X
            bcFct2.setFunction(zero); //!< high  X
            bcFct3.setFunction(zero); //!< low   Y
            bcFct4.setFunction(zero); //!< high  Y
            bcFct5.setFunction(zero); //!< low   Z
            bcFct6.setFunction(zero); //!< high  Z

            nb_bc = 6;
            bc.setNumber(nb_bc);

            bc.addBC("Analytical, real BC",    1, Essential,   Scalar, bcFct1);
            bc.addBC("Analytical, real BC",    2, Essential,   Scalar, bcFct2);
            bc.addBC("Analytical, real BC",    3, Essential,   Scalar, bcFct3);
            bc.addBC("Analytical, real BC",    4, Essential,   Scalar, bcFct4);
            bc.addBC("Analytical, real BC",    5, Essential,   Scalar, bcFct5);
            bc.addBC("Analytical, real BC",    6, Essential,   Scalar, bcFct6);

            break;

        case 4:
            /*
              Robin condition (dp/dq = alpha + beta*p, alpha=2.535e-6, beta=-7.596) at endothel (ref 3)
              and Dirichlet (p=7.4226e-8) at adventitia (ref 1)
              example of mesh: piece-of-tube-wall-tetra.mesh
              diffusion coefficient: 1.0
              This is the scaled problem of plasma filtration through the artery wall with phys.
              realistic data: diff.coeff (=darcy_perm/mu_plasma = 2e-14/0.0072 = 2.78e-12
              RobinBC resulting from the flux balance at the endothel (L_p_end=2.11e-11)
              pressure at the adventitia = 20 mmHg
              (scaling: p* = diff.coeff * p otherwise wrong result - bad condition of Matrices)
            */
            bcFct_rob.setFunctions_Mixte(alpha, beta);
            bcFct2.setFunction(p_adv);
            nb_bc = 2;
            bc.setNumber(nb_bc);
            bc.addBC("Endothel",   3,     Mixte, Scalar, bcFct_rob);
            bc.addBC("Adventitia",  1, Essential, Scalar, bcFct2);
            break;

        default:
            ERROR_MSG("Unknown test case");
    }

#endif

//
// diffusion
//
#if 0
    double xg,yg,zg;

    switch(diffusion_type){
    case 0:
        //-------------------------
        // *** scalar diffusion ***
        //-------------------------
        switch(diffusion_function){
        case 0: // constant diffusion given in the data file
            mass_Hdiv(1./diffusion_scalar,elmatMix,vfe,0,0);
            break;
        case 9: // porous medium with a function
            pfe.barycenter(xg,yg,zg); // coordinate of the barycenter of the current element
            diffusion_scalar = permeability_sd009(xg, yg, zg);
            mass_Hdiv(1./diffusion_scalar,elmatMix,vfe,0,0);
            break;
        case 10: // porous medium with a function
            pfe.barycenter(xg,yg,zg); // coordinate of the barycenter of the current element
            diffusion_scalar = permeability_sd010(xg, yg, zg);
            mass_Hdiv(1./diffusion_scalar,elmatMix,vfe,0,0);
            break;
        default:
            std::cerr << "Unknown function for scalar diffusion. Change physics/diffusion_function in the data file\n"
                      << "\n";
            exit(1);
        }
        break;
    case 1:
        //-------------------------
        // *** tensor diffusion ***
        //-------------------------
        {
            KNM<double> permlower(3,3),invpermea(3,3);
            switch(diffusion_function){
            case 0: // constant diffusion tensor given in the data file
                permlower = diffusion_tensor;
                permlower *= diffusion_scalar;
                /* // the diffusion matrix is divided by the viscosity which sould be 1/diffusion_scalar!
                   Remark: not very optimal since in this case
                   this constant matrix is inverted on each elements.
                   (can be easily improved if needed)
                */
                break;
            case 1: // fibrous medium
                pfe.barycenter(xg,yg,zg); // coordinate of the barycenter of the current element
                permlower = fibrous_permea(1./diffusion_scalar,diffusion_tensor, xg); //! added "1./"
                break;
            default:
                std::cerr << "Unknown function for tensor diffusion. Change physics/diffusion_function in the data file\n"
                          << "\n";
                exit(1);
            }
            //
            // we compute the inverse of permlower
            //
            // PERM <- L and Lt where L Lt is the Cholesky factorization of PERM
            int NBT[1] = {3}; //  dim of tensor permeabilite
            int INFO[1] = {0};
            dpotrf_("L", NBT, permlower, NBT, INFO);
            ASSERT_PRE(!INFO[0],"Lapack factorization of PERM is not achieved.");
            dpotri_("L", NBT, permlower, NBT, INFO);
            ASSERT_PRE(!INFO[0],"Lapack solution of PERM is not achieved.");
            permlower(0,1) = permlower(1,0);
            permlower(0,2) = permlower(2,0);
            permlower(1,2) = permlower(2,1);
            invpermea = permlower;
            //
            mass_Hdiv(invpermea, elmatMix, vfe,0,0); //  modify the (0,0) block of the matrix
            break;
        }
    default:
        std::cerr << "diffusion_type=" << diffusion_type << " ??? \n"
                  << "\n";
        exit(1);
    }
#endif
