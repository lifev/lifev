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

#include "reducedLinFluid.hpp"

namespace LifeV
{
reducedLinFluid::reducedLinFluid(FSIOperator* const _op,
                         fluid_type _fluid,
                         solid_type _solid):
    M_FSIOperator(_op),
    M_fluid( _fluid ),
    M_solid( _solid ),
    M_BCh_dp( new BCHandler ),
    M_dacc( M_solid->dDof().numTotalDof() ),
    M_refFE( M_fluid->refFEp() ),
    M_Qr( quadRuleTetra4pt ),
    M_bdQr( quadRuleTria3pt ),
    M_dof(M_fluid->mesh(), M_refFE),
    M_dim(M_dof.numTotalDof()),
    M_pattC(M_dof),
    M_C(M_pattC),
    M_CAux(M_pattC),
    M_fe(M_refFE, getGeoMap(M_fluid->mesh()), M_Qr ),
    M_feBd(M_refFE.boundaryFE(), getGeoMap(M_fluid->mesh()  ).boundaryMap(), M_bdQr),
    M_elmatC(M_fe.nbNode, 1, 1),
    M_dp(M_dim),
    M_residual_dp(3*M_dim),
    M_f(M_dim),
    M_computedC(false),
    M_computedResidual(false)

{
}

void reducedLinFluid::setUpBC(bchandler_type _BCh_dp)
{
    M_BCh_dp = _BCh_dp;
}


void reducedLinFluid::setLinearSolver(GetPot const &_data_file)
{
    M_linearSolver.setOptionsFromGetPot( _data_file, "reducedfluid/aztec" );
    M_linearSolver.setMatrix( M_C );
}

const Vector& reducedLinFluid::residual()
{
    if (!M_computedResidual)
        evalResidual();

    return M_residual_dp;
}

void reducedLinFluid::solveReducedLinearFluid()
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

    bcManage(M_C, M_f, M_fluid->mesh(), M_dof, *M_BCh_dp, M_feBd, 1.0, M_FSIOperator->time());

    M_linearSolver.setRecursionLevel( 1 );
    M_dp = ZeroVector( M_dp.size() );

    M_linearSolver.solve( M_dp, M_f, SolverAztec::SAME_PRECONDITIONER );

    M_computedResidual = false;

    M_minusdp = -1.0*M_dp;
}


void reducedLinFluid::evalResidual()

{
    M_residual_dp = ZeroVector(M_residual_dp.size());

    bchandler_type BCb = M_FSIOperator->BCh_mesh();

    if (!BCb->bdUpdateDone())
        BCb->bdUpdate(M_fluid->mesh(), M_fluid->feBd_u(), M_fluid->uDof() );

    UInt iBC = BCb->getBCbyName("Interface");
    // Number of local Dof (i.e. nodes) in this face
    UInt nDofF = this->M_feBd.nbNode;

    // Number of total scalar Dof
    UInt totalDof = this->M_dof.numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = (*BCb)[iBC].numberOfComponents();

    const IdentifierNatural* pId;
    ID ibF, idDof, icDof, gDof;

    std::cout << "iBC = " << iBC << std::endl;
    std::cout << "NoC = " << nComp << std::endl;
    std::cout << "dof = " << totalDof << std::endl;

    for ( ID i = 1; i <= (*BCb)[iBC].list_size(); ++i )
    {
        // Pointer to the i-th itdentifier in the list
        pId = static_cast< const IdentifierNatural* >( (*BCb)[iBC]( i ) );

        // Number of the current boundary face
        ibF = pId->id();

        std::cout << "iBF = " << ibF << std::flush << std::endl;
        // Updating face stuff
        M_feBd.updateMeasNormalQuadPt( M_fluid->mesh().boundaryFace( ibF ) );

        std::cout << "normal updated ... " << std::flush << std::endl;
        // Loop on total Dof per Face
        for ( ID l = 1; l <= nDofF; ++l )
        {
            gDof = pId->bdLocalToGlobal( l );

            std::cout << "gDof = " << gDof << std::flush << std::endl;
            // Loop on components involved in this boundary condition
            for ( UInt ic = 0; ic < nComp; ++ic )
            {
                icDof = gDof + ic * totalDof;

                // Loop on quadrature points
                for ( int iq = 0; iq < M_feBd.nbQuadPt; ++iq )
                {
                    std::cout << icDof << " " << gDof << std::flush << std::endl;
                    // Adding right hand side contribution
                    M_residual_dp[ icDof ] += M_dp[gDof]*
                        M_feBd.phi( int( l - 1 ), iq )*
                        M_feBd.normal( int( ic ), iq )*
                        M_feBd.weightMeas( iq );
                }
            }
        }
    }



//     for ( UInt gDof = 0; gDof < M_dp.size(); ++gDof )
//     {
//         M_feBd.updateMeasNormalQuadPt( M_fluid->mesh().boundaryFace( gDof ) );

//         for ( ID l = 1; l <= nDofF; ++l )
//             for ( UInt ic = 0; ic < nComp; ++ic )
//             {
//                 icDof = gDof + ic * totalDof;

//                 // Loop on quadrature points
//                 for ( int iq = 0; iq < this->M_feBd.nbQuadPt; ++iq )
//                 {
//                     M_residual_dp[ icDof ] += M_dp[gDof]*
//                         this->M_feBd.phi( int( l - 1 ), iq )*
//                         this->M_feBd.normal( int( ic ), iq )*
//                         this->M_feBd.weightMeas( iq );
//                 }
//             }
//         std::cout << icDof << "->" << gDof << "->" << M_residual_dp[ icDof ] << std::endl;
//    }

    M_computedResidual = true;

    std::cout << "dp residual = " << norm_2(M_residual_dp) << std::endl;
    std::cout << "dp          = " << norm_2(M_dp) << std::endl;

}

}

