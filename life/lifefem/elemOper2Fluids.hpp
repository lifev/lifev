/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Daniele Antonio Di Pietro <dipietro@unibg.it>
       Date: 2-1-2005

  Copyright (C) 2005 Università degli Studi di Bergamo

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
   \file elemOper2Fluids.hpp
   \author Daniele A. Di Pietro <dipietro@unibg.it>
   \date 2-1-2005
 */

#ifndef _ELEMOPER2FLUIDS_H_
#define _ELEMOPER2FLUIDS_H_

#include "lifeV.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "currentFE.hpp"
#include "currentBdFE.hpp"
#include "dof.hpp"

namespace LifeV {

    /**
       \Two-fluid mass contribution
    */

    void mass_2f(Real coef1, Real coef2, 
                 const ElemVec& ls_fun_loc, const CurrentFE& fe_ls,
                 ElemMat& elmat, const CurrentFE& fe,
                 int iblock, int jblock);

    /**
       \Two-fluid stiffness-strain contribution
    */

    void stiff_strain_2f(Real coef1, Real coef2,
                         const ElemVec& ls_fun_loc, const CurrentFE& fe_ls,
                         ElemMat& elmat, const CurrentFE& fe);


    /**
       \Two-fluid passive advection term
    */

    void advection_2f(const int icoor, Real coef1, Real coef2,
                      const ElemVec& ls_fun_loc, const CurrentFE& fe_ls,
                      ElemMat& elmat, const CurrentFE& fe,
                      int iblock, int jblock);

    // Implementations

    void mass_2f(Real coef1, Real coef2, 
                 const ElemVec& ls_fun_loc, const CurrentFE& fe_ls,
                 ElemMat& elmat, const CurrentFE& fe,
                 int iblock = 0, int jblock = 0) {
        ASSERT_PRE( fe.hasJac(), "Mass matrix needs at least the Jacobian" );
        ElemMat::matrix_view mat = elmat.block(iblock, jblock);
        Real s, coef_s;
        Real ls_fun_on_qn;

        for (int i = 0; i < fe.nbNode; i++) {
            for(int j = 0; j < fe.nbNode; j++) {
                s = 0;
                for (int iq = 0; iq < fe.nbQuadPt; iq++) {
                    // Evaluate level set function on current quadrature node

                    ls_fun_on_qn = 0.;
                    for(int k = 0; k < fe_ls.nbNode; k++)
                        ls_fun_on_qn += ls_fun_loc(k) * fe_ls.phi(k, iq);

                    // Evaluate coefficient on current quadrature node

                    coef_s = coef1 * (ls_fun_on_qn > 0) + coef2 * (ls_fun_on_qn < 0);

                    // Compute local contribution

                    s += coef_s * fe.phi(i, iq) * fe.phi(j, iq) * fe.weightDet(iq);
                }

                mat(i, j) += s;
            }
        }
    }

    void stiff_strain_2f(Real coef1, Real coef2,
                         const ElemVec& ls_fun_loc, const CurrentFE& fe_ls,
                         ElemMat& elmat, const CurrentFE& fe) {
    
        ASSERT_PRE( fe.hasFirstDeriv(),
                    "Stiffness Strain matrix needs at least the first derivatives" );
        Real s, coef_s;
        Real ls_fun_on_qn;

        ElemMat::matrix_type mat_tmp(fe.nbNode, fe.nbNode);

        for ( int i = 0; i < fe.nbNode; ++i ) {
            for ( int j = 0; j < fe.nbNode; ++j ) {
                s = 0;
                for (int iq = 0; iq < fe.nbQuadPt; ++iq) {
                    // Evaluate level set function on current quadrature node

                    ls_fun_on_qn = 0.;
                    for(int k = 0; k < fe_ls.nbNode; k++)
                        ls_fun_on_qn += ls_fun_loc(k) * fe_ls.phi(k, iq);

                    // Evaluate coefficient on current quadrature node

                    coef_s = 0.5 * ( coef1 * (ls_fun_on_qn > 0) + coef2 * (ls_fun_on_qn < 0) );

                    // Compute local contribution
                    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
                        s += coef_s * fe.phiDer( i, icoor, iq ) * fe.phiDer( j, icoor, iq ) * fe.weightDet( iq );
                }
                mat_tmp( i, j ) = s;
            }
        }

        for ( int icoor = 0; icoor < fe.nbCoor; ++icoor ) {
            ElemMat::matrix_view mat = elmat.block( icoor, icoor );
            mat += mat_tmp;
        }

        for ( int icoor = 0; icoor < fe.nbCoor; ++icoor ) {
            for ( int jcoor = 0; jcoor < fe.nbCoor; ++jcoor ) {
                ElemMat::matrix_view mat = elmat.block( icoor, jcoor );
                for ( int i = 0; i < fe.nbNode; ++i ) {
                    for ( int j = 0; j < fe.nbNode; ++j ) {
                        s = 0;
                        for (int iq = 0; iq < fe.nbQuadPt; ++iq) {
                            // Evaluate level set function on current quadrature node

                            ls_fun_on_qn = 0.;
                            for(int k = 0; k < fe_ls.nbNode; k++)
                                ls_fun_on_qn += ls_fun_loc(k) * fe_ls.phi(k, iq);

                            // Evaluate coefficient on current quadrature node

                            coef_s = 0.5 * ( coef1 * (ls_fun_on_qn > 0) + coef2 * (ls_fun_on_qn < 0) );

                            // Compute local contribution
                            s += coef_s * fe.phiDer( i, jcoor, iq ) * fe.phiDer( j, icoor, iq ) * fe.weightDet( iq );
                        }
                        mat( i, j ) += s;
                    }
                }
            }
        }
    }

    void advection_2f(Real coef1, Real coef2,
                      const ElemVec& ls_fun_loc, const CurrentFE& fe_ls,
                      const ElemVec& vel_loc,
                      ElemMat& elmat, const CurrentFE& fe,
                      int iblock = 0, int jblock = 0) {
        ASSERT_PRE( fe.hasJac(), "Mass matrix needs at least the Jacobian" );
        ElemMat::matrix_view mat = elmat.block(iblock, jblock);
        Real s, coef_s, coef_v;
        Real ls_fun_on_qn;

        for (int i = 0; i < fe.nbNode; i++) {
            for(int j = 0; j < fe.nbNode; j++) {
                s = 0;
                for (int iq = 0;iq < fe.nbQuadPt;iq++) {
                    // Evaluate level set function on current quadrature node

                    ls_fun_on_qn = 0.;
                    for(int k = 0; k < fe_ls.nbNode; k++)
                        ls_fun_on_qn += ls_fun_loc(k) * fe_ls.phi(k, iq);

                    // Evaluate coefficient on current quadrature node

                    coef_s = coef1 * (ls_fun_on_qn > 0) + coef2 * (ls_fun_on_qn < 0);

                    // Compute local contribution

                    for(int icoor = 0; icoor < NDIM; icoor++) {
                        coef_v = 0.;

                        // Compute k-th component of velocity vector

                        for(int k = 0; k < fe.nbNode; k++)
                            coef_v += vel_loc.vec() [k + icoor * fe.nbNode] * fe.phi(k, iq);
                        
                        // Add local contribution

                        s += coef_s * coef_v * fe.phiDer(j, icoor, iq) * fe.phi(i, iq) * fe.weightDet(iq);
                    }
                }

                mat(i, j) += s;
            }
        }
    }
}

#endif
