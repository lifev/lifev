/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-
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
  This file implements Jean Frederic Gerbeau's and Marina Vidrascu's
  quasi-newton algorithm. See INRIA RR 4691 for further details.

  Miguel Fernandez miguel.fernandez@inria.fr
  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file reducedLinFluid.hpp
*/

//! previously called quasiNewton in the CVS repository

#ifndef _QNEWTON
#define _QNEWTON


#include <life/lifecore/life.hpp>
#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifesolver/NavierStokesAleSolverPC.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
#include <life/lifefem/regionMesh3D_ALE.hpp>


namespace LifeV
{


class reducedLinFluid{

    typedef FSIOperator::fluid_type     fluid_type;
    typedef FSIOperator::solid_type     solid_type;
    typedef FSIOperator::bchandler_type bchandler_type;

public:


    //! constructor
    reducedLinFluid(FSIOperator* const _op,
                fluid_type _fluid,
                solid_type _solid);

    //! bondary conditions setup

    void setUpBC   (bchandler_type _BCh_dp);
    void setUpInvBC(bchandler_type _BCh_dp_inv);

    //! jacobian computation

    void solveReducedLinearFluid();
    void solveInvReducedLinearFluid();


    const Vector& dacc()    {return M_dacc;}
    void setDacc(Vector const &_vec){M_dacc = _vec;}
    void setComputedMatrix(bool pred){M_computedC = pred;}
    Vector& minusdp() {return M_minusdp;}

    const Vector& residual();
    void evalResidual();

    void setLinearSolver (GetPot const &_data_file);

private:

    //! pointer to FSIOperator
    FSIOperator*     M_FSIOperator;

    //! pointer to fluid
    fluid_type       M_fluid;

    //!pointer to solid
    solid_type       M_solid;

    //! BC for pressure
    bchandler_type   M_BCh_dp;
    bchandler_type   M_BCh_dp_inv;

    //! linearized acceleration
    Vector           M_dacc;


    //
    // laplacian stuff
    //

    //! reference FE
    const RefFE&     M_refFE;

    //! volumic quadrature rule
    const QuadRule&  M_Qr;

    //! boundary quadrature ru;e
    const QuadRule&  M_bdQr;

    //! dof
    Dof              M_dof;

    //! number of total dof
    UInt             M_dim;

    //! pattern of C
    MSRPatt          M_pattC;

    //! global matrix
    MSRMatr<double>  M_C;
    MSRMatr<double>  M_CAux;

    //! Current FE
    CurrentFE        M_fe;
    CurrentBdFE      M_feBd;

    //! elementary matrix and vector
    ElemMat          M_elmatC;

    //! unknown: presure variation and residual
    Vector           M_dp;
    Vector           M_minusdp;
    Vector           M_residual_dp;

    //! right handside
    Vector           M_f;

    //! linear solver
    SolverAztec      M_linearSolver;

    //! computed matrix C?
    bool             M_computedC;

    //! computed residual ?
    bool             M_computedResidual;
};
}

#endif
