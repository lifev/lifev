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

#include "quasiNewton.hpp"

namespace LifeV
{

quasiNewton::quasiNewton(fluid_type       _fluid,
                         operFS* const    _op):
    M_operFS(_op),
    M_fluid( _fluid ),
    M_BCh_dp( new BCHandler ),
    M_refFE( _fluid->refFEp() ),
    M_dacc( _op->solid().dDof().numTotalDof() ),
    M_Qr( quadRuleTetra4pt ),
    M_bdQr( quadRuleTria3pt ),
    M_dof(_fluid->mesh(), M_refFE),
    M_dim(M_dof.numTotalDof()),
    M_pattC(M_dof),
    M_C(M_pattC),
    M_CAux(M_pattC),
    M_fe(M_refFE, getGeoMap(_fluid->mesh()), M_Qr ),
    M_feBd(M_refFE.boundaryFE(), getGeoMap(_fluid->mesh()  ).boundaryMap(),M_bdQr),
    M_elmatC(M_fe.nbNode, 1, 1),
    M_dp(M_dim),
    M_residual_dp(3*_fluid->uDof().numTotalDof()),
    M_f(M_dim),
    M_computedC(false)
{
}

void quasiNewton::setUpBC(bchandler_type _BCh_dp)
{
    M_BCh_dp = _BCh_dp;
}


void quasiNewton::setLinearSolver(GetPot const &_data_file)
{
    M_linearSolver.setOptionsFromGetPot( _data_file, "reducedfluid/aztec" );
    M_linearSolver.setMatrix( M_C );
}

void quasiNewton::solveReducedLinearFluid()
{
    if (!M_computedC) { // computing laplacien fe matrix (begining of the time step)

        // Initializing matrix
        M_CAux.zeros();

        // Loop on elements
        for ( UInt i = 1; i <= M_fluid->mesh().numVolumes(); i++ )
        {
            M_fe.updateFirstDerivQuadPt( M_fluid->mesh().volumeList( i ) );

            M_elmatC.zero();

            stiff(1.0, M_elmatC, M_fe );
            assemb_mat(M_CAux, M_elmatC, M_fe, M_dof, 0, 0);
        }

        M_BCh_dp->bdUpdate(M_fluid->mesh(), M_feBd, M_dof);
        M_computedC = true;
    }

    M_C = M_CAux;
    M_f = ZeroVector( M_f.size() );

    bcManage(M_C, M_f, M_fluid->mesh(), M_dof, *M_BCh_dp, M_feBd, 1.0, M_operFS->time());

    M_linearSolver.setRecursionLevel( 1 );
    M_dp = ZeroVector( M_dp.size() );

    M_linearSolver.solve( M_dp, M_f, SolverAztec::SAME_PRECONDITIONER );

    evalResidual();
    M_minusdp = -1.0*M_dp;
}


void quasiNewton::evalResidual()
{

    M_residual_dp = ZeroVector(M_residual_dp.size());

    int iBC = M_BCh_dp->getBCbyName("Wall");

    UInt nDofF    = M_feBd.nbNode;
    UInt totalDof = M_dof.numTotalDof();

    const IdentifierNatural* pId;
    ID ibF, icDof, gDof;

    // Number of components involved in this boundary condition
    UInt nComp = (*M_BCh_dp)[iBC].numberOfComponents();

    for ( ID i = 1; i <= (*M_BCh_dp)[iBC].list_size(); ++i )
    {

        // Pointer to the i-th itdentifier in the list
        pId = static_cast< const IdentifierNatural* >( (*M_BCh_dp)[iBC]( i ) );

        // Number of the current boundary face
        ibF = pId->id();

        // Updating face stuff
        M_feBd.updateMeasNormalQuadPt( M_fluid->mesh().boundaryFace( ibF ) );

        // Loop on total Dof per Face
        for ( ID l = 1; l <= nDofF; ++l )
        {
            gDof = pId->bdLocalToGlobal( l );

            // Loop on components involved in this boundary condition
            for ( UInt ic = 0; ic < nComp; ++ic )
            {
                icDof = gDof + ic * totalDof;

                // Loop on quadrature points
                for ( int iq = 0; iq < M_feBd.nbQuadPt; ++iq )
                {
                    // Adding right hand side contribution
                    M_residual_dp[ icDof - 1 ] += M_dp[ gDof ] * M_feBd.phi( int( l - 1 ), iq )
                        * M_feBd.normal( int( ic ), iq )
                        * M_feBd.weightMeas( iq );
                }
            }
        }
    }

    std::cout << "dp residual = " << norm_2(M_residual_dp) << std::endl;
    std::cout << "dp          = " << norm_2(M_dp) << std::endl;

}

}

