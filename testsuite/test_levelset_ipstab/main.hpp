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
#define N_PHY_UNK 1

#define PI 3.141592653589793283

#include <utility>
#include <algorithm>
#include <functional>
#include <vector>
#include <list>
#include <set>

#include <boost/numeric/ublas/vector.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/elemOper_ext.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/dofDG.hpp>
#include <life/lifemesh/markers.hpp>
#include <life/lifemesh/basisElSh.hpp>

#include <life/lifefem/values.hpp>
#include <life/lifefem/assemb.hpp>

#include <life/lifealg/dataAztec.hpp>

namespace LifeV
{

    /**
       \Analytical velocity field
    */

    class vortex {
    public:
        inline Real operator()(Real& x, Real& y, Real& z, int i) const {
            Real s = 0.;
            switch(i){
            case 1: {
                s = sin(PI * x) * sin(PI * x) * (sin(2 * PI * z) - sin(2 * PI * y));
                break;
            }
            case 2: {
                s = sin(PI * y) * sin(PI * y) * (sin(2 * PI * x) - sin(2 * PI * z));
                break;
            }
            case 3: {
                s = sin(PI * z) * sin(PI * z) * (sin(2 * PI * y) - sin(2 * PI * x));
                break;
            }
            }
            return s;
        }
    };

    class uniform
    {
    public:
        uniform() {};
        static void set( Real ux, Real uy, Real uz )
            {
                M_u[0] = ux;
                M_u[1] = uy;
                M_u[2] = uz;
            }
        inline Real operator() ( Real /* x */, Real /* y */,
                                 Real /* z */, int i ) const
            {
                return M_u[i+1];
            }
        static Real u( Real /* t */, Real /* x */, Real /* y */, Real /* z */,
                       int i )
            {
                return M_u[i+1];
            }
    private:
        static Real M_u[3];
    };

Real uniform::M_u[3];

    /**
       \Analytical expression for signed distance function
    */

    class sphere {
    public:
        static void init( Real x0, Real y0, Real z0, Real r,
                          Real ux, Real uy, Real uz )
            {
                M_x0 = x0;
                M_y0 = y0;
                M_z0 = z0;
                M_r = r;
                M_ux = ux;
                M_uy = uy;
                M_uz = uz;
            }
        static Real phi(const Real& t,
                        const Real& x,
                        const Real& y,
                        const Real& z,
                        const ID& /* i */) {
            Real dx = M_x0 + M_ux * t - x;
            Real dy = M_y0 + M_uy * t - y;
            Real dz = M_z0 + M_uz * t - z;
            return sqrt( dx*dx + dy*dy + dz*dz ) - M_r;
        }
    private:
        static Real M_x0;
        static Real M_y0;
        static Real M_z0;
        static Real M_r;
        static Real M_ux;
        static Real M_uy;
        static Real M_uz;
    };

Real sphere::M_x0 = 0.5;
Real sphere::M_y0 = 0.5;
Real sphere::M_z0 = 0.75;
Real sphere::M_r = 0.2;
Real sphere::M_ux = 1.;
Real sphere::M_uy = 1.;
Real sphere::M_uz = 1.;

    /**
       \Project the analytical velocity field onto finite element space
    */

    template <class DOF, class MESH, class VECTOR, class  FUNCTION, class MATRIX, class SOLVER>
    void projectVelocityField(const MESH& mesh, CurrentFE& fe, const DOF& dof, VECTOR& U, const FUNCTION& u_ex, UInt nc,
                              MATRIX& M, SOLVER& solver) {

        VECTOR b(U.size());
        b = ZeroVector(U.size());

        ElemVec elvec(fe.nbNode, nc);

        for (UInt ic = 0; ic < nc; ic++)
            for (UInt i = 1; i <= mesh.numVolumes(); i++ ){
                fe.updateJacQuadPt(mesh.volumeList(i));

                elvec.zero();
                compute_vec(u_ex, elvec, fe, (int)ic);
                assemb_vec(b, elvec, fe, dof, (int)ic);
            }

        solver.setMatrix(M);
        solver.solve(U, b);
    }

    /**
       \A passive advection element operator to be used with analytical velocity
       \field
    */

    template<typename BETA>
    void passiveAdvection(BETA beta, ElemMat& elmat, CurrentFE& fe, int iblock = 0, int jblock = 0, int nb = 0) {
        Real x, y, z, s;
        Real b_dot_grad;

        ElemMat::matrix_view mat = elmat.block(iblock, jblock);

        for (int i = 0; i < fe.nbNode; i++) {
            for (int j = 0; j < fe.nbNode; j++) {
                s = 0.;
                for (int iq = 0; iq < fe.nbQuadPt; iq++) {
                    // Co-ordinates of current quadrature point
                    x = fe.quadPt(iq, 0);
                    y = fe.quadPt(iq, 1);
                    z = fe.quadPt(iq, 2);

                    b_dot_grad = 0.;
                    for(int icoor = 0; icoor < fe.nbCoor; icoor++)
                        b_dot_grad += beta(icoor, x, y, z) * fe.phiDer(j,icoor,iq) * fe.phi(i,iq);

                    s += fe.weightDet(iq) * b_dot_grad;
                } // iq
                mat(i, j) += s;
            } // j
        } // i
    }
}
