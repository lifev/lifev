/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-11-18

  Copyright (C) 2004 EPFL

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
   \file FSISolver.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-18
 */
#include <lifeV.hpp>

#include <FSISolver.hpp>

namespace LifeV
{
FSISolver::FSISolver( GetPot const& data_file,
                      bchandler_type& __bcu, bchandler_type& __bcd, bchandler_type& __bchext,
                      std::string __oper ):
    M_BCh_u( __bcu ),
    M_BCh_d( __bcd ),
    M_BCh_mesh( __bchext ),
    M_fluid( new FSIOperator::fluid_type::value_type (data_file, feTetraP1bubble, feTetraP1,quadRuleTetra64pt,
                                                 quadRuleTria3pt, quadRuleTetra64pt, quadRuleTria3pt,
                                                 *M_BCh_u,*M_BCh_mesh) ),
    M_solid( new FSIOperator::solid_type::value_type (data_file, feTetraP1, quadRuleTetra4pt,
                                                 quadRuleTria3pt, *M_BCh_d) ),
    M_disp(3*M_solid->dDof().numTotalDof()),
    M_velo(3*M_solid->dDof().numTotalDof()),
    M_method( data_file("problem/method"    , "steklovPoincare") ),
    M_maxpf( data_file("problem/maxSubIter", 300) ),
    M_defomega( data_file("problem/defOmega"  , 0.01) ),
    M_abstol( data_file("problem/abstol"  , 1.e-07) ),
    M_reltol( data_file("problem/reltol"  , 1.e-04) ),
    M_etamax( data_file("problem/etamax"  , 1.e-03) ),
    M_linesearch( data_file("problem/linesearch"  , 0) ),
    out_iter("iter"),
    out_res ("res")

{
    Debug( 6220 ) << "FSISolver::FSISolver starts\n";

    M_disp   = ZeroVector( M_disp.size() );
    M_velo   = ZeroVector( M_velo.size() );

    Debug( 6220 ) << "FSISolver::M_disp: " << M_disp.size() << "\n";
    Debug( 6220 ) << "FSISolver::M_velo: " << M_velo.size() << "\n";

    Preconditioner precond  = ( Preconditioner ) data_file("problem/precond"   , DIRICHLET_NEUMANN );

    Debug( 6220 ) << "FSISolver::preconditioner: " << precond << "\n";

    if ( !__oper.empty() )
    {
        M_method = __oper;
    }

    this->setFSIOperator( M_method );

    M_oper->setDataFromGetPot( data_file );

    M_oper->setPreconditioner( precond );

    M_oper->setUpBC();

    Debug( 6220 ) << "FSISolver::FSISolver ends\n";
}
void
FSISolver::setFSIOperator( std::string const& __op )
{
    Debug( 6220 ) << "FSISolver::setFSIOperator with operator " << __op << "\n";

    M_oper = oper_fsi_ptr( FSIFactory::instance().createObject( __op ) );
    M_oper->setFluid( M_fluid );
    M_oper->setSolid( M_solid );

    UInt dim_solid = M_oper->solid().dDof().numTotalDof();
    UInt dim_fluid = M_oper->fluid().uDof().numTotalDof();


    // Passing data from structure to the fluid: solid velocity at the interface velocity
    //
    BCVectorInterface u_wall( M_oper->fluid().wInterpolated(),
                              dim_fluid,
                              M_oper->dofMeshToFluid() );

    // Boundary conditions for the fluid velocity
    M_BCh_u->addBC("Wall",   1,  Essential, Full, u_wall,  3);

    M_oper->setBC( M_BCh_u, M_BCh_d, M_BCh_mesh );
    M_oper->setup();
    Debug( 6220 ) << "FSISolver::setFSIOperator done\n";
}
void
FSISolver::iterate( Real time )
{
    Debug( 6220 ) << "============================================================\n";
    Debug( 6220 ) << "Solving FSI at time " << time << " with FSIOperator: " << M_method  << "\n";

    UInt dim_solid = M_oper->solid().dDof().numTotalDof();
    UInt dim_fluid = M_oper->fluid().uDof().numTotalDof();


    M_oper->fluid().timeAdvance( M_oper->fluid().sourceTerm(), time);
    M_oper->solid().timeAdvance( M_oper->solid().sourceTerm(), time);
    M_oper->setTime(time);

    // displacement prediction

    M_disp   = M_oper->solid().d() + timeStep()*(1.5*M_oper->solid().w() - 0.5*M_velo);

    M_velo = M_oper->solid().w();

    std::cout << "norm( disp   ) init = " << norm_inf(M_disp)   << std::endl;
    std::cout << "norm( velo ) init = " << norm_inf(M_velo) << std::endl;

    int maxiter = M_maxpf;

    // the newton solver
    uint status = 1;
    if ( M_method == "exactJacobian" )
    {
        status = newton(M_disp, *M_oper, norm_inf_adaptor(),
                        M_abstol, M_reltol, maxiter, M_etamax,
                        M_linesearch, out_res, time);
    }
    else
    {
        status = nonLinRichardson(M_disp, *M_oper, norm_inf_adaptor(),
                                  M_abstol, M_reltol, maxiter, M_etamax,
                                  M_linesearch, out_res, time, M_defomega);
    }

    if(status == 1)
    {
        std::ostringstream __ex;
        __ex << "FSISolver::iterate ( " << time << " ) Inners iterations failed to converge\n";
        throw std::logic_error( __ex.str() );
    }
    else
    {
        std::cout << "End of time "<< time << std::endl;
        std::cout << "Number of inner iterations       : "
                  << maxiter << std::endl;
        out_iter << time << " " << maxiter << " "
                 << M_oper->nbEval() << std::endl;

        M_oper->fluid().postProcess();
        M_oper->solid().postProcess();
    }

    Debug( 6220 ) << "FSISolver iteration at time " << time << " done\n";
    Debug( 6220 ) << "============================================================\n";
}
}
