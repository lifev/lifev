//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
   @file
   @brief Ionic model based on Fox model.
   @date 04-2013
   @author Marie Dupraz <dupraz.marie@gmail.com>

   @contributors
   @mantainer Marie Dupraz <dupraz.marie@gmail.com>
   @last update 04-2013
  */

#include <lifev/electrophysiology/solver/IonicModels/IonicFox.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <cmath>


namespace LifeV
{
// ===================================================
//! Constructors
// ===================================================
IonicFox::IonicFox()    :
    super       ( 13, 7 ),
    M_ACap      ( 1.534e-4 ),
    M_VMyo      ( 25.84e-6 ),
    M_Vup       ( 0.1 ),
    M_Vsr       ( 2e-6 ),
    M_NaO       ( 138.0 ),
    M_CaO       ( 2000.0 ),
    M_K0        ( 4.0 ),
    M_NaIn      ( 10.0 ),
    M_KIn       ( 149.4 ),
    M_CmdnTot   ( 10.0 ),
    M_CsqnTot   ( 10000.0 ),
    M_KmCmdn    ( 2.0 ),
    M_KmCsqn    ( 600.0 ),
    M_Cm        ( 1.0 ),
    M_F         ( 96.5 ),
    M_T         ( 310.0 ),
    M_R         ( 8.314 ),
    M_GNa       ( 12.8 ),
    M_GKp       ( 0.002216 ),
    M_KmNa      ( 87.5 ),
    M_KmCa      ( 1380.0 ),
    M_kSat      ( 0.2 ),
    M_kNaCa     ( 1500.0 ),
    M_eta       ( 0.35 ),
    M_INaK      ( 0.693 ),
    M_KmNai     ( 10.0 ),
    M_KmK0      ( 1.5 ),
    M_IpCa      ( 0.05 ),
    M_KmPCa     ( 0.05 ),
    M_GCab      ( 0.0003842 ),
    M_GNab      ( 0.0031 ),
    M_KmUp      ( 0.32 ),
    M_PCa       ( 0.0000226 ),
    M_ICaHalf   ( -0.265 ),
    M_GK1       ( 2.8 ),
    M_GKr       ( 0.0136 ),
    M_GKs       ( 0.0245 ),
    M_Gt0       ( 0.23815 ),
    M_PCaK      ( 5.79e-7 ),
    M_Prel      ( 6.0 ),
    M_Pleak     ( 0.000001 ),
    M_KmfCa     ( 0.18 ),
    M_KmK1      ( 13.0 )
{
    M_restingConditions.at (0) = - 94.7;
    M_restingConditions.at (1) = 2.4676e-4;
    M_restingConditions.at (2) = 0.99869;
    M_restingConditions.at (3) = 0.99887;
    M_restingConditions.at (4) = 0.229;
    M_restingConditions.at (5) = 0.0001;
    M_restingConditions.at (6) = 3.742e-5;
    M_restingConditions.at (7) = 1.0;
    M_restingConditions.at (8) = 0.983;
    M_restingConditions.at (9) = 0.0001;
    M_restingConditions.at (10) = 0.942;
    M_restingConditions.at (11) = 0.0472;
    M_restingConditions.at (12) = 320.0;
}

IonicFox::IonicFox ( Teuchos::ParameterList& parameterList ) :
    super       ( 13, 7 )
{
    M_ACap       = parameterList.get ( "areaCap", 1.534e-4 );
    M_VMyo       = parameterList.get ( "volMyo", 25.84e-6 );
    M_Vup        = parameterList.get ( "volUp", 0.1 );
    M_Vsr        = parameterList.get ( "volSR", 2e-6 );
    M_NaO        = parameterList.get ( "concNa0", 138.0 );
    M_CaO        = parameterList.get ( "concCa0", 2000.0 );
    M_K0         = parameterList.get ( "concK0", 4.0 );
    M_NaIn       = parameterList.get ( "concNaIn", 10.0 );
    M_KIn        = parameterList.get ( "concKIn", 149.4 );
    M_CmdnTot    = parameterList.get ( "cmdnTot", 10.0 );
    M_CsqnTot    = parameterList.get ( "csqnTot", 10000.0 );
    M_KmCmdn     = parameterList.get ( "constmCmdn", 2.0 );
    M_KmCsqn     = parameterList.get ( "constmCsqn", 600.0 );
    M_Cm         = parameterList.get ( "capMem", 1.0 );
    M_F          = parameterList.get ( "farad", 96.5 );
    M_T          = parameterList.get ( "temp", 310.0 );
    M_R          = parameterList.get ( "gasConst", 8.314 );
    M_GNa        = parameterList.get ( "maxCondNa", 12.8 );
    M_GKp        = parameterList.get ( "maxCondKp", 0.002216 );
    M_KmNa       = parameterList.get ( "constmNa", 87.5 );
    M_KmCa       = parameterList.get ( "constmCa", 1380.0 );
    M_kSat       = parameterList.get ( "kSat", 0.2 );
    M_kNaCa      = parameterList.get ( "kNaCa", 1500.0 );
    M_eta        = parameterList.get ( "eta", 0.35 );
    M_INaK       = parameterList.get ( "courNaK", 0.693 );
    M_KmNai      = parameterList.get ( "constmNai", 10.0 );
    M_KmK0       = parameterList.get ( "constmK0", 1.5 );
    M_IpCa       = parameterList.get ( "courpCa", 0.05 );
    M_KmPCa      = parameterList.get ( "constmpCa", 0.05 );
    M_GCab       = parameterList.get ( "maxCondCab", 0.0003842 );
    M_GNab       = parameterList.get ( "maxCondNab", 0.0031 );
    M_KmUp       = parameterList.get ( "constmUp", 0.32 );
    M_PCa        = parameterList.get ( "permCa", 0.0000226 );
    M_ICaHalf    = parameterList.get ( "courCaHalf", -0.265 );
    M_GK1        = parameterList.get ( "maxCondK1", 2.8 );
    M_GKr        = parameterList.get ( "maxCondKr", 0.0136 );
    M_GKs        = parameterList.get ( "maxCondKs", 0.0245 );
    M_Gt0        = parameterList.get ( "maxCondt0", 0.23815 );
    M_PCaK       = parameterList.get ( "permCaK", 5.79e-7 );
    M_Prel       = parameterList.get ( "permrel", 6.0 );
    M_Pleak      = parameterList.get ( "permleak", 0.000001 );
    M_KmfCa      = parameterList.get ( "constmfCa", 0.18 );
    M_KmK1       = parameterList.get ( "constmK1", 13.0 );
}

IonicFox::IonicFox ( const IonicFox& model )
{
    M_ACap      = model.M_ACap;
    M_VMyo      = model.M_VMyo;
    M_Vup       = model.M_Vup;
    M_Vsr       = model.M_Vsr;
    M_NaO       = model.M_NaO;
    M_CaO       = model.M_CaO;
    M_K0        = model.M_K0;
    M_NaIn      = model.M_NaIn;
    M_KIn       = model.M_KIn;
    M_CmdnTot   = model.M_CmdnTot;
    M_CsqnTot   = model.M_CsqnTot;
    M_KmCmdn    = model.M_KmCmdn;
    M_KmCsqn    = model.M_KmCsqn;
    M_Cm        = model.M_Cm;
    M_F         = model.M_F;
    M_T         = model.M_T;
    M_R         = model.M_R;
    M_GNa       = model.M_GNa;
    M_GKp       = model.M_GKp;
    M_KmNa      = model.M_KmNa;
    M_KmCa      = model.M_KmCa;
    M_kSat      = model.M_kSat;
    M_kNaCa     = model.M_kNaCa;
    M_eta       = model.M_eta;
    M_INaK      = model.M_INaK;
    M_KmNai     = model.M_KmNai;
    M_KmK0      = model.M_KmK0;
    M_IpCa      = model.M_IpCa;
    M_KmPCa     = model.M_KmPCa;
    M_GCab      = model.M_GCab;
    M_GNab      = model.M_GNab;
    M_KmUp      = model.M_KmUp;
    M_PCa       = model.M_PCa;
    M_ICaHalf   = model.M_ICaHalf;
    M_GK1       = model.M_GK1;
    M_GKr       = model.M_GKr;
    M_GKs       = model.M_GKs;
    M_Gt0       = model.M_Gt0;
    M_PCaK      = model.M_PCaK;
    M_Prel      = model.M_Prel;
    M_Pleak     = model.M_Pleak;
    M_KmfCa     = model. M_KmfCa;
    M_KmK1      = model. M_KmK1;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;
    M_numberOfGatingVariables = model.M_numberOfGatingVariables;
}

// ===================================================
//! Operator
// ===================================================
IonicFox& IonicFox::operator= ( const IonicFox& model )
{
    M_ACap      = model.M_ACap;
    M_VMyo      = model.M_VMyo;
    M_Vup       = model.M_Vup;
    M_Vsr       = model.M_Vsr;
    M_NaO       = model.M_NaO;
    M_CaO       = model.M_CaO;
    M_K0        = model.M_K0;
    M_NaIn      = model.M_NaIn;
    M_KIn       = model.M_KIn;
    M_CmdnTot   = model.M_CmdnTot;
    M_CsqnTot   = model.M_CsqnTot;
    M_KmCmdn    = model.M_KmCmdn;
    M_KmCsqn    = model.M_KmCsqn;
    M_Cm        = model.M_Cm;
    M_F         = model.M_F;
    M_T         = model.M_T;
    M_R         = model.M_R;
    M_GNa       = model.M_GNa;
    M_GKp       = model.M_GKp;
    M_KmNa      = model.M_KmNa;
    M_KmCa      = model.M_KmCa;
    M_kSat      = model.M_kSat;
    M_kNaCa     = model.M_kNaCa;
    M_eta       = model.M_eta;
    M_INaK      = model.M_INaK;
    M_KmNai     = model.M_KmNai;
    M_KmK0      = model.M_KmK0;
    M_IpCa      = model.M_IpCa;
    M_KmPCa     = model.M_KmPCa;
    M_GCab      = model.M_GCab;
    M_GNab      = model.M_GNab;
    M_KmUp      = model.M_KmUp;
    M_PCa       = model.M_PCa;
    M_ICaHalf   = model.M_ICaHalf;
    M_GK1       = model.M_GK1;
    M_GKr       = model.M_GKr;
    M_GKs       = model.M_GKs;
    M_Gt0       = model.M_Gt0;
    M_PCaK      = model.M_PCaK;
    M_Prel      = model.M_Prel;
    M_Pleak     = model.M_Pleak;
    M_KmfCa     = model. M_KmfCa;
    M_KmK1      = model. M_KmK1;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;
    M_numberOfGatingVariables = model.M_numberOfGatingVariables;

    return *this;
}


// ===================================================
//! Methods
// ===================================================
//Only gating variables
void IonicFox::computeGatingRhs ( const std::vector<Real>&  v, std::vector<Real>& rhs )
{
    std::vector<Real> gatingRhs     ( computeLocalGatingRhs (v) );
    std::vector<Real> subSysCaRhs   ( computeLocalSubSysCaRhs (v) );

    std::copy ( gatingRhs.begin(), gatingRhs.end(), rhs.begin() );
    std::copy ( subSysCaRhs.begin(), subSysCaRhs.end() - 1, rhs.begin() + gatingRhs.size() );
}

//Potential and gating variables
void IonicFox::computeRhs (const   std::vector<Real>&  v, std::vector<Real>& rhs )
{
    std::vector<Real> gatingRhs     ( computeLocalGatingRhs (v) );
    std::vector<Real> subSysCaRhs   ( computeLocalSubSysCaRhs (v) );

    rhs[0] = computeLocalPotentialRhs (v );

    if (rhs.size() > M_numberOfEquations)
    {
        std::copy ( gatingRhs.begin(), gatingRhs.end(), rhs.begin() + 1 );
        std::copy ( subSysCaRhs.begin(), subSysCaRhs.end(), rhs.begin() + 1 + gatingRhs.size() );
    }
    else
    {
        std::copy ( gatingRhs.begin(), gatingRhs.end(), rhs.begin() + 1 );
        std::copy ( subSysCaRhs.begin(), subSysCaRhs.end() - 1, rhs.begin() + 1 + gatingRhs.size() );
    }

}

void IonicFox::computeNonGatingRhs ( const std::vector<Real>&  v, std::vector<Real>& rhs )
{
    //    std::vector<Real> subSysCaRhs   ( computeLocalSubSysCaRhs (v) );
    //    std::copy ( subSysCaRhs.begin(), subSysCaRhs.end() - 1, rhs.begin() + M_numberOfGatingVariables );
    std::vector<Real> gatingRhs     ( computeLocalGatingRhs (v) );
    std::vector<Real> subSysCaRhs   ( computeLocalSubSysCaRhs (v) );

    std::copy ( gatingRhs.begin(), gatingRhs.end(), rhs.begin() );
    std::copy ( subSysCaRhs.begin(), subSysCaRhs.end() - 1, rhs.begin() + gatingRhs.size() );
}

Real IonicFox::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    std::vector<Real> fastNa (fastINa (v) );
    std::vector<Real> rapidK (rapidIK (v) );
    std::vector<Real> slowK  (slowIK  (v) );
    std::vector<Real> transOutK (transOutIK (v) );
    std::vector<Real> subSysCaRHS (computeLocalSubSysCaRhs (v) );
    Real INa = fastNa[0];
    Real IKr = rapidK[0];
    Real IKs = slowK[0];
    Real It0 = transOutK[0];
    Real ICa = subSysCaRHS[5];

    return   - ( INa + timeIIK1 (v) + IKr + IKs + It0 + plaIKp (v) + pumpINaK (v)
                 + exINaCa (v) + backINab (v) + backICab (v) + pumpIpCa (v) + ICa ) ;
}

std::vector<Real> IonicFox::computeLocalGatingRhs ( const std::vector<Real>& v )
{
    Real m ( v[1] );
    Real h ( v[2] );
    Real j ( v[3] );

    Real xr ( v[4] );
    Real xks ( v[5] );
    Real xt0 ( v[6] );
    Real yt0 ( v[7] );

    std::vector<Real> gatingRhs (7);
    std::vector<Real> paramNa   ( fastINa (v) );
    std::vector<Real> paramKr   ( rapidIK (v) );
    std::vector<Real> paramKs   ( slowIK (v) );
    std::vector<Real> paramKt0   ( transOutIK (v) );

    gatingRhs[0] = paramNa[1] * ( 1 - m ) - paramNa[2] * m;
    gatingRhs[1] = paramNa[3] * ( 1 - h ) - paramNa[4] * h;
    gatingRhs[2] = paramNa[5] * ( 1 - j ) - paramNa[6] * j;

    gatingRhs[3] = ( paramKr[1] - xr ) / paramKr[2];
    gatingRhs[4] = ( paramKs[1] - xks ) / paramKs[2];
    gatingRhs[5] = paramKt0[1] * ( 1 - xt0 ) - paramKt0[2] * xt0;
    gatingRhs[6] = paramKt0[3] * ( 1 - yt0 ) - paramKt0[4] * yt0;

    return gatingRhs;
}

void IonicFox::computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt )
{
    Real m ( v[1] );
    Real h ( v[2] );
    Real j ( v[3] );

    Real xr ( v[4] );
    Real xks ( v[5] );
    Real xt0 ( v[6] );
    Real yt0 ( v[7] );

    std::vector<Real> paramNa   ( fastINa (v) );
    std::vector<Real> paramKr   ( rapidIK (v) );
    std::vector<Real> paramKs   ( slowIK (v) );
    std::vector<Real> paramKt0   ( transOutIK (v) );

    Real TAU_M = 1.0 / ( paramNa[1] + paramNa[2] );
    Real M_INF = paramNa[1] * TAU_M;
    Real TAU_H = 1.0 / ( paramNa[3] + paramNa[4] );
    Real H_INF = paramNa[3] * TAU_M;
    Real TAU_J = 1.0 / ( paramNa[5] + paramNa[6] );
    Real J_INF = paramNa[5] * TAU_M;

    Real TAU_Xr = paramKr[2];
    Real Xr_INF = paramKr[1];
    Real TAU_Xks = paramKs[2];
    Real Xks_INF = paramKs[1];

    Real TAU_xto = 1.0 / ( paramKt0[1] + paramKt0[2] );
    Real xto_INF = paramKt0[1] * TAU_M;
    Real TAU_yto = 1.0 / ( paramKt0[3] + paramKt0[4] );
    Real yto_INF = paramKt0[3] * TAU_M;

    v[1] = M_INF    - ( M_INF    - m   ) * std::exp (- dt / TAU_M   );
    v[2] = H_INF    - ( H_INF    - h   ) * std::exp (- dt / TAU_H   );
    v[3] = J_INF    - ( J_INF    - j   ) * std::exp (- dt / TAU_J   );
    v[4] = Xr_INF   - ( Xr_INF   - xr  ) * std::exp (- dt / TAU_Xr  );
    v[5] = Xks_INF  - ( Xks_INF  - xks ) * std::exp (- dt / TAU_Xks );
    v[6] = xto_INF  - ( xto_INF  - xt0 ) * std::exp (- dt / TAU_xto );
    v[7] = yto_INF  - ( yto_INF  - yt0 ) * std::exp (- dt / TAU_yto );
}


//! Ca2+ Subsystem

std::vector<Real> IonicFox::computeLocalSubSysCaRhs ( const std::vector<Real>& v )
{
    std::vector<Real> subSysCaRHS (6);

    Real V       ( v[0] );
    Real cCaIn   ( v[11] );
    Real cCaSR   ( v[12] );

    Real gating_d ( v[9] );
    Real gating_f ( v[8] );
    Real gating_f_Ca ( v[10] );

    // Internal Parameters

    Real gamma = 1.0 / ( 1.0 + pow ( ( 2000.0 / cCaSR ), 3) );

    Real jRel  = M_Prel * gating_d * gating_f * gating_f_Ca * ( gamma * cCaSR - cCaIn ) / ( 1.0 + 1.65 * exp (V / 20.0) );
    Real jLeak = M_Pleak * ( cCaSR - cCaIn );
    Real jUp   = M_Vup / ( 1.0 + pow ( ( M_KmUp / cCaIn ) , 2.0 ) );

    Real betaIn   = 1 / ( 1 + M_CmdnTot * M_KmCmdn / ( ( M_KmCmdn + cCaIn ) * ( M_KmCmdn + cCaIn ) ) );
    Real betaSR = 1 / ( 1 + M_CsqnTot * M_KmCsqn / ( ( M_KmCsqn + cCaSR ) * ( M_KmCsqn + cCaSR ) ) );

    Real sd = 1.0 / ( 1.0 + exp ( - ( V + 10.0 ) / 6.24 ) );
    Real td = 1.0 / ( (0.25 * exp (-0.01 * V) / (1.0 + exp (-0.07 * V) ) ) + (0.07 * exp (-0.05 * (V + 40.0) ) / (1.0 + exp (0.05 * (V + 40.0) ) ) ) );

    Real sfca = 1.0 / (1.0 + pow ( (cCaIn / M_KmfCa), 3.0) );
    Real tfca = 30.0;

    Real sf = 1.0 / (1 + exp (V + 12.5) / 5.0);
    Real tf = 30.0 + (200.0 / (1 + exp ( (V + 20.0) / 9.5) ) );

    Real I_Ca_full  = (M_PCa / M_Cm) * (4.0 * ( V * M_F * M_F ) / ( M_R * M_T ) ) *
                      ( (cCaIn * exp (2.0 * ( V * M_F ) / ( M_R * M_T ) ) - 0.341 * M_CaO) /
                        (exp (2.0 * ( V * M_F ) / ( M_R * M_T ) ) - 1.0) );
    Real I_Cal = gating_d * gating_f * gating_f_Ca * I_Ca_full;
    Real I_Ca_K = (M_PCaK / M_Cm) * (gating_d * gating_f * gating_f_Ca / (1.0 + (I_Ca_full / M_ICaHalf) ) ) *
                  (1000.0 * ( V * M_F * M_F ) / ( M_R * M_T ) ) * ( (M_KIn * exp ( ( V * M_F )
                                                                                   / ( M_R * M_T ) ) - M_K0) / (exp ( ( V * M_F ) / ( M_R * M_T ) ) - 1.0) );

    Real I_Ca = I_Cal + I_Ca_K;

    // RHS of the Ca2+ subsystem

    subSysCaRHS[4]  = betaSR * ( jUp - jLeak - jRel) * ( M_VMyo / M_Vsr ); // dCaIn
    subSysCaRHS[3]  = betaIn * ( jRel + jLeak - jUp - (M_ACap * M_Cm / (2 * M_F * M_VMyo) ) * (I_Cal + backICab (v) + M_IpCa - 2 * exINaCa (v) ) ); //dCaSR

    subSysCaRHS[0]  = ( sf - gating_f ) / tf;
    subSysCaRHS[1]  = ( sd - gating_d ) / td;
    subSysCaRHS[2]  = ( sfca - gating_f_Ca ) / tfca;

    subSysCaRHS[5]  = I_Ca;

    return subSysCaRHS;
}

//! Ionic Currents (Luo and Rudy)

// Fast Na+ Current INa
std::vector<Real> IonicFox::fastINa ( const std::vector<Real>& v )
{
    std::vector<Real> fastNa (7);

    Real V   ( v[0] );
    Real m   ( v[1] );
    Real h   ( v[2] );
    Real j   ( v[3] );

    Real potNa = ( M_R * M_T / M_F ) * log ( M_NaO / M_NaIn );

    fastNa[0] = M_GNa * m * m * m * h * j * ( V - potNa );

    Real alpha_m = 0.32 * ( V + 47.13 ) / ( 1.0 - exp ( - 0.1 * ( V + 47.13 ) ) );
    Real beta_m  = 0.08 * exp (- V / 11.0 );

    Real alpha_h = 0.135 * exp ( - ( V + 80.0 ) / 6.8 );
    Real beta_h = 7.5 / ( 1.0 + exp ( - 0.1 * ( V + 11.0 ) ) );
    Real alpha_j = ( 0.175 * exp ( - ( V + 100.0 ) / 23.0 ) ) / ( 1.0 + exp ( 0.15 * ( V + 79.0 ) ) );
    Real beta_j = 0.3 / ( 1.0 + exp ( -0.1 * ( V + 32.0 ) ) );

    fastNa[1] = alpha_m;
    fastNa[2] = beta_m;
    fastNa[3] = alpha_h;
    fastNa[4] = beta_h;
    fastNa[5] = alpha_j;
    fastNa[6] = beta_j;

    return fastNa;
}

// rapid K+ current
std::vector<Real> IonicFox::rapidIK ( const std::vector<Real>& v )
{
    std::vector<Real> rapidK (3);

    Real V   ( v[0] );
    Real xr  ( v[4] );

    Real potK1  = ( M_R * M_T / M_F ) * std::log ( M_K0 / M_KIn );
    Real RV     = 1.0 / ( 1.0 + 2.5 * std::exp ( 0.1 * ( V + 28.0 ) ) );

    Real sxr    = 1.0 / ( 1.0 + std::exp ( -2.182 - 0.1819 * V ) );
    Real txr    = 43.0 + ( 1.0 / ( std::exp ( -5.495 + 0.1691 * V ) + std::exp ( -7.677 - 0.0128 * V ) ) );

    rapidK[0] = M_GKr * xr * RV * std::sqrt ( M_K0 / 4.0 ) * ( V - potK1 );

    rapidK[1] = sxr;
    rapidK[2] = txr;

    return rapidK;
}

// slow K+ current
std::vector<Real> IonicFox::slowIK ( const std::vector<Real>& v )
{


    Real V   ( v[0] );
    Real xs  ( v[5] );

    Real potKs  = ( M_R * M_T / M_F ) * std::log ( ( M_K0 + 0.01833 * M_NaO ) / ( M_KIn + 0.01833 * M_NaIn ) );

    Real sxs    = 1.0 / ( 1.0 + std::exp ( - ( V - 16.0 ) / 13.6 ) );
    Real txs    = 1.0 / ( ( 0.0000719 * ( V - 10.0 ) / ( 1.0 - std::exp ( -0.148 * ( V - 10.0 ) ) ) )
                          + ( 0.000131 * ( V - 10.0 ) / ( std::exp ( 0.0687 * (V - 10.0) ) - 1.0) ) );

    std::vector<Real> slow (3, 0.0);
    slow[0] = M_GKs * xs * xs * ( V - potKs );
    slow[1] = sxs;
    slow[2] = txs;

    return slow;
}

// transient outward K+ current
std::vector<Real> IonicFox::transOutIK ( const std::vector<Real>& v )
{
    std::vector<Real> transOutK (5);

    Real V    ( v[0] );
    Real xto  ( v[6] );
    Real yto  ( v[7] );

    Real potK1    = ( M_R * M_T / M_F ) * log ( M_K0 / M_KIn );

    Real axto     = 0.04516 * exp ( 0.03577 * V );
    Real bxto     = 0.0989 * exp ( -0.06237 * V );
    Real ayto     = ( 0.005415 * exp ( - ( V + 33.5 ) / 5.0 ) ) / ( 1.0 + 0.051335 * exp ( - ( V + 33.5 ) / 5.0 ) );
    Real byto     = ( 0.005415 * exp ( ( V + 33.5 ) / 5.0 ) ) / ( 1.0 + 0.051335 * exp ( ( V + 33.5 ) / 5.0 ) );

    transOutK[0]  = M_Gt0 * xto * yto * ( V - potK1 );

    transOutK[1]  = axto;
    transOutK[2]  = bxto;
    transOutK[3]  = ayto;
    transOutK[4]  = byto;

    return transOutK;
}

// Time-independent K+ Current IK1
Real IonicFox::timeIIK1 ( const std::vector<Real>& v )
{
    Real V   ( v[0] );

    Real potK1      = ( M_R * M_T / M_F ) * log ( M_K0 / M_KIn );
    Real maxCondK1  = M_GK1;
    Real sK1        = 1.0 / ( 2.0 + exp ( ( 1.62 / ( M_R * M_T / M_F ) ) * ( V - potK1 ) ) );
    Real constK1inf = M_K0 / ( M_K0 + M_KmK1 );

    return  maxCondK1 * sK1 * constK1inf * ( V - potK1 );
}

// Plateau K+ current IKp
Real IonicFox::plaIKp ( const std::vector<Real>& v )
{
    Real V   ( v[0] );

    Real potKp   = ( M_R * M_T / M_F ) * log ( M_K0 / M_KIn );
    Real constKp = 1 / ( 1 + exp ( ( 7.488 - V ) / 5.98 ) );

    return M_GKp * constKp * ( V - potKp );
}

// Na+/Ca2+ exchanger current INaCa
Real IonicFox::exINaCa ( const std::vector<Real>& v )
{

    Real V         ( v[0] );
    Real cCaIn     ( v[11] );

    return M_kNaCa * ( 1.0 / ( pow (M_KmNa, 3) + pow (M_NaO, 3) ) ) * ( 1.0 / ( M_KmCa + M_CaO ) )
           * (1.0 / ( 1.0 + M_kSat * exp ( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) ) )
           * ( exp ( M_eta * ( V * M_F ) / ( M_R * M_T ) ) * pow (M_NaIn, 3) * M_CaO -
               exp ( ( M_eta - 1.0 ) * ( V * M_F ) / ( M_R * M_T ) ) * pow (M_NaO, 3) * cCaIn );
}

// Na+/K+ pump INaK
Real IonicFox::pumpINaK ( const std::vector<Real>& v )
{
    Real V   ( v[0] );

    Real sigma = ( 1.0 / 7.0 ) * ( exp ( M_NaO / 67.3 ) - 1.0 );
    Real fNak  = 1.0 / ( 1.0 + 0.1245 * exp ( -0.1 * ( V * M_F ) / ( M_R * M_T ) ) +
                         0.0365 * sigma * exp ( - ( V * M_F ) / ( M_R * M_T ) ) );

    return M_INaK * fNak * ( 1.0 / ( 1.0 + pow ( M_KmNai / M_NaIn, 1.5) ) ) * ( M_K0 / ( M_K0 + M_KmK0 ) );
}

// Sarcolemmal Ca2+ pump current IpCa
Real IonicFox::pumpIpCa ( const std::vector<Real>& v )
{
    return M_IpCa * v[11] / ( M_KmPCa + v[11] );
}

// Ca2+ background current ICab OK
Real IonicFox::backICab ( const std::vector<Real>& v )
{
    Real V      ( v[0] );
    Real cCaIn  ( v[11] );

    Real potCaN = ( ( M_R * M_T ) / ( 2 * M_F ) ) * log ( M_CaO / cCaIn );

    return M_GCab * ( V - potCaN );
}

// Na+ background current INab
Real IonicFox::backINab ( const std::vector<Real>& v)
{
    Real V   ( v[0] );

    Real potNaN = M_R * M_T / M_F * log ( M_NaO / M_NaIn );

    return M_GNab * ( V - potNaN );
}

void IonicFox::showMe()
{
    std::cout << "\n\n\t\tIonicFox Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    //        std::cout << "\n\t\tList of model parameters:\n\n";
    //        std::cout << "areaCap: " << this->areaCap() << std::endl;
    //        std::cout << "volMyo: " << this->volMyo() << std::endl;
    //        std::cout << "volJSR: " << this->volJSR() << std::endl;
    //        std::cout << "volNSR: " << this->volNSR() << std::endl;
    //        std::cout << "volSS: " << this->volSS() << std::endl;
    //        std::cout << "concNa0: " << this->concNa0() << std::endl;
    //        std::cout << "concCa0: " << this->concCa0() << std::endl;
    //        std::cout << "lTrpnTot: " << this->lTrpnTot() << std::endl;
    //        std::cout << "hTrpnTot: " << this->hTrpnTot() << std::endl;
    //        std::cout << "kpHtrpn: " << this->kpHtrpn() << std::endl;
    //        std::cout << "knHtrpn: " << this->knHtrpn() << std::endl;
    //        std::cout << "kpLtrpn: " << this->kpLtrpn() << std::endl;
    //        std::cout << "knLtrpn: " << this->knLtrpn() << std::endl;
    //        std::cout << "cmdnTot: " << this->cmdnTot() << std::endl;
    //        std::cout << "csqnTot: " << this->csqnTot() << std::endl;
    //        std::cout << "constmCmdn: " << this->constmCmdn() << std::endl;
    //        std::cout << "constmCsqn: " << this->constmCsqn() << std::endl;
    //        std::cout << "capMem: " << this->capMem() << std::endl;
    //        std::cout << "farad: " << this->farad() << std::endl;
    //        std::cout << "temp: " << this->temp() << std::endl;
    //        std::cout << "gasConst: " << this->gasConst() << std::endl;
    //        std::cout << "maxCondNa: " << this->maxCondNa() << std::endl;
    //        std::cout << "maxCondKp: " << this->maxCondKp() << std::endl;
    //        std::cout << "permNaK: " << this->permNaK() << std::endl;
    //        std::cout << "KNaCa: " << this->kNaCa() << std::endl;
    //        std::cout << "constmNa: " << this->constmNa() << std::endl;
    //        std::cout << "constmCa: " << this->constmCa() << std::endl;
    //        std::cout << "kSat: " << this->kSat() << std::endl;
    //        std::cout << "eta: " << this->eta() << std::endl;
    //        std::cout << "courNaK: " << this->courNaK() << std::endl;
    //        std::cout << "constmNai: " << this->constmNai() << std::endl;
    //        std::cout << "constmK0: " << this->constmK0() << std::endl;
    //        std::cout << "permNsCa: " << this->permNsCa() << std::endl;
    //        std::cout << "constmNsCa: " << this->constmNsCa() << std::endl;
    //        std::cout << "courpCa: " << this->courpCa() << std::endl;
    //        std::cout << "constmpCa: " << this->constmCa() << std::endl;
    //        std::cout << "maxCondCab: " << this->maxCondCab() << std::endl;
    //        std::cout << "maxCondNab: " << this->maxCondNab() << std::endl;
    //        std::cout << "maxRyRPerm: " << this->maxRyRPerm() << std::endl;
    //        std::cout << "leakRateConst: " << this->leakRateConst() << std::endl;
    //        std::cout << "pumpRateATPase: " << this->pumpRateATPase() << std::endl;
    //        std::cout << "constmUp: " << this->constmUp() << std::endl;
    //        std::cout << "timeConstNsrJsr: " << this->timeConstNsrJsr() << std::endl;
    //        std::cout << "timeConstSubMyo: " << this->timeConstSubMyo() << std::endl;
    //        std::cout << "kAPlus: " << this->kAPlus() << std::endl;
    //        std::cout << "kANeg: " << this->kANeg() << std::endl;
    //        std::cout << "kBPlus: " << this->kBPlus() << std::endl;
    //        std::cout << "kBNeg: " << this->kBNeg() << std::endl;
    //        std::cout << "kCPlus: " << this->kCPlus() << std::endl;
    //        std::cout << "kCNeg: " << this->kCNeg() << std::endl;
    //        std::cout << "coopParamN: " << this->coopParamN() << std::endl;
    //        std::cout << "coopParamM: " << this->coopParamM() << std::endl;
    //        std::cout << "intoOpenSt: " << this->intoOpenSt() << std::endl;
    //        std::cout << "outOpenSt: " << this->outOpenSt() << std::endl;
    //        std::cout << "intoOpenStca: " << this->intoOpenStCa() << std::endl;
    //        std::cout << "outOpentSt2: " << this->outOpenSt2() << std::endl;
    //        std::cout << "modeTParamA: " << this->modeTParamA() << std::endl;
    //        std::cout << "modeTParamB: " << this->modeTParamB() << std::endl;
    //        std::cout << "modeTParamO: " << this->modeTParamO() << std::endl;
    //        std::cout << "permCa: " << this->permCa() << std::endl;
    //        std::cout << "permK: " << this->permK() << std::endl;
    //        std::cout << "courCaHalf: " << this->courCaHalf() << std::endl;


    std::cout << "\n\t\t End of IonicFox Informations\n\n\n";
}


}

