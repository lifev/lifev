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
  @brief Ionic model of ten Tusscher-Panfilov 2006
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher06.hpp>


namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
IonicTenTusscher06::IonicTenTusscher06()  :
    super       (  19 , 18  ),
    flag        ( MCell )
{
    //Membrane capacitance
    M_membraneCapacitance = 2.0;// In the paper they give it as 2 but in their code I cannot find it

    //External concentrations
    Ko = 5.4;
    Cao = 2.0;
    Nao = 140.0;

    //Intracellular volumes
    Vc = 0.016404;
    Vsr = 0.001094;
    Vss = 0.00005468;

    //Calcium buffering dynamics
    Bufc = 0.2;
    Kbufc = 0.001;
    Bufsr = 10.;
    Kbufsr = 0.3;
    Bufss = 0.4;
    Kbufss = 0.00025;

    //Intracellular calcium flux dynamics
    Vmaxup = 0.006375;
    Kup = 0.00025;
    Vrel = 0.102; //40.8;
    k1_ = 0.15;
    k2_ = 0.045;
    k3 = 0.060;
    k4 = 0.005; //0.000015;
    EC = 1.5;
    maxsr = 2.5;
    minsr = 1.;
    Vleak = 0.00036;
    Vxfer = 0.0038;



    //Constants
    R = 8314.472;
    F = 96485.3415;
    T = 310.0;
    RTONF = (R * T) / F;

    //Cellular capacitance
    M_cellularCapacitance = 0.185;

    //Parameters for currents
    //Parameters for IKr
    Gkr = 0.153;
    //Parameters for Iks
    pKNa = 0.03;

    if (flag == MCell)
    {
        Gks = 0.098;
    }
    else
    {
        Gks = 0.932;
    }
    //Parameters for Ik1
    GK1 = 5.405;
    //Parameters for Ito

    if (flag == Endo)
    {
        Gks = 0.073;
    }
    else
    {
        Gto = 0.294;
    }

    //Parameters for INa
    GNa = 14.838;
    //Parameters for IbNa
    GbNa = 0.00029;
    //Parameters for INaK
    KmK = 1.0;
    KmNa = 40.0;
    knak = 2.724;
    //Parameters for ICaL
    GCaL = 0.00003980;
    //Parameters for IbCa
    GbCa = 0.000592;
    //Parameters for INaCa
    knaca = 1000;
    KmNai = 87.5;
    KmCa = 1.38;
    ksat = 0.1;
    n = 0.35;
    //Parameters for IpCa
    GpCa = 0.1238;
    KpCa = 0.0005;
    //Parameters for IpK;
    GpK = 0.0146;

    inverseVcF2 = 1. / (2.*Vc * F);
    inverseVcF = 1. / (Vc * F);
    inversevssF2 = 1. / (2.*Vss * F);

    //V
    M_restingConditions.at (0) = -86.2;
    //m
    M_restingConditions.at (1) = 0.0;
    //h
    M_restingConditions.at (2) = 0.75;
    //j
    M_restingConditions.at (3) = 0.75;
    //d
    M_restingConditions.at (4) = 0.;
    //f
    M_restingConditions.at (5) = 1.;
    //f2
    M_restingConditions.at (6) = 1;
    //fCass
    M_restingConditions.at (7) = 1.0;
    //r
    M_restingConditions.at (8) = 0.;
    //s
    M_restingConditions.at (9) = 1.0;
    //Xr1
    M_restingConditions.at (10) = 0.0;
    //Xr2
    M_restingConditions.at (11) = 1.0;
    //Xs
    M_restingConditions.at (12) =  0.0;
    //Nai
    M_restingConditions.at (13) = 7.67;
    //Ki
    M_restingConditions.at (14) = 138.3;
    //Cai
    M_restingConditions.at (15) = 0.00007;
    //Cass
    M_restingConditions.at (16) = 0.00007;
    //Casr
    M_restingConditions.at (17) = 1.3;
    //Rprime
    M_restingConditions.at (18) = 1.0;

    //  M_restingConditions.at(0) = -86.2;
    //  //m
    //  M_restingConditions.at(1) = 0.00165;
    //  //h
    //  M_restingConditions.at(2) = 0.749;
    //  //j
    //  M_restingConditions.at(3) = 0.6788;
    //  //d
    //  M_restingConditions.at(4) = 3.288e-5;
    //  //f
    //  M_restingConditions.at(5) = 0.7026;
    //  //f2
    //  M_restingConditions.at(6) = 0.9526;
    //  //fCass
    //  M_restingConditions.at(7) = 1.0;
    //  //r
    //  M_restingConditions.at(8) = 2.347e-8;
    //  //s
    //  M_restingConditions.at(9) = 0.999998;
    //  //Xr1
    //  M_restingConditions.at(10) = 0.0165;
    //  //Xr2
    //  M_restingConditions.at(11) = 0.473;
    //  //Xs
    //  M_restingConditions.at(12) =  0.0174;
    //  //Nai
    //  M_restingConditions.at(13) = 7.67;
    //  //Ki
    //  M_restingConditions.at(14) = 138.3;
    //  //Cai
    //  M_restingConditions.at(15) = 0.00007;
    //  //Cass
    //  M_restingConditions.at(16) = 0.00007;
    //  //Casr
    //  M_restingConditions.at(17) = 1.3;
    //  //Rprime
    //  M_restingConditions.at(18) = 0.8978;


}

IonicTenTusscher06::IonicTenTusscher06 ( Teuchos::ParameterList& parameterList     )   :
    super       ( 19, 12 )
{
    knak = parameterList.get ("knak", 0.0 );
    KmNa = parameterList.get ("KmNa", 0.0 );
    KmK = parameterList.get ("KmK", 0.0 );
    knaca = parameterList.get ("knaca", 0.0 );
    KmNai = parameterList.get ("KmNai", 0.0 );
    KmCa = parameterList.get ("KmCa", 0.0 );
    ksat = parameterList.get ("ksat", 0.0 );
    n = parameterList.get ("n", 0.0 );


    Ko = parameterList.get ("Ko", 0.0 );
    Cao = parameterList.get ("Cao", 0.0 );
    Nao = parameterList.get ("Nao", 0.0 );


    Bufc = parameterList.get ("Bufc", 0.0 );
    Kbufc = parameterList.get ("Kbufc", 0.0 );
    Bufsr = parameterList.get ("Bufsr", 0.0 );
    Kbufsr = parameterList.get ("Kbufsr", 0.0 );
    Bufss = parameterList.get ("Bufss", 0.0 );
    Kbufss = parameterList.get ("Kbufss", 0.0 );

    Vmaxup = parameterList.get ("Vmaxup", 0.0 );
    Kup = parameterList.get ("Kup", 0.0 );
    Vrel = parameterList.get ("Vrel", 0.0 );
    k1_ = parameterList.get ("k1", 0.0 );
    k2_ = parameterList.get ("k2", 0.0 );
    k3 = parameterList.get ("k3", 0.0 );
    k4 = parameterList.get ("k4", 0.0 );
    EC = parameterList.get ("EC", 0.0 );
    maxsr = parameterList.get ("maxsr", 0.0 );
    minsr = parameterList.get ("minsr", 0.0 );
    Vleak = parameterList.get ("Vleak", 0.0 );
    Vxfer = parameterList.get ("Vxfer", 0.0 );


    pKNa = parameterList.get ("pKNa", 0.0 );


    M_cellularCapacitance = parameterList.get ("Cm", 0.0 );
    F = parameterList.get ("F", 0.0 );
    R = parameterList.get ("R", 0.0 );
    T = parameterList.get ("T", 0.0 );


    Gkr = parameterList.get ("Gkr", 0.0 );
    Gks = parameterList.get ("Gks", 0.0 );
    GK1 = parameterList.get ("GK1", 0.0 );
    Gto = parameterList.get ("Gto", 0.0 );
    GNa = parameterList.get ("GNa", 0.0 );
    GbNa = parameterList.get ("GbNa", 0.0 );
    GCaL = parameterList.get ("GCal", 0.0 );
    GbCa = parameterList.get ("GbCa", 0.0 );
    GpCa = parameterList.get ("GpCa", 0.0 );
    KpCa = parameterList.get ("KpCa", 0.0 );
    GpK = parameterList.get ("GpK", 0.0 );


    Vc = parameterList.get ("Vc", 0.0 );
    Vsr = parameterList.get ("Vsr", 0.0 );
    Vss = parameterList.get ("Vss", 0.0 );

    computeRTONF();
    computeinverseVcF2();
    computeinverseVcF();
    computeinversevssF2();

    std::map< std::string, WallFlag > WallFlagMap;
    WallFlagMap["Endo"]     = Endo;
    WallFlagMap["Epi"]      = Epi;
    WallFlagMap["MCell"]    = MCell;

    flag = WallFlagMap[ parameterList.get ( "WallFlag", "MCell" ) ];
    //V
    M_restingConditions.at (0) = parameterList.get ("V", -85.423);
    //m
    M_restingConditions.at (1) = parameterList.get ("m", 0.00165);
    //h
    M_restingConditions.at (2) = parameterList.get ("h", 0.749);
    //j
    M_restingConditions.at (3) = parameterList.get ("j", 0.6788);
    //d
    M_restingConditions.at (4) = parameterList.get ("d", 3.288e-5);
    //f
    M_restingConditions.at (5) = parameterList.get ("f", 0.7026);
    //f2
    M_restingConditions.at (6) = parameterList.get ("f2", 0.9526);
    //fCass
    M_restingConditions.at (7) = parameterList.get ("fCass", 0.9942);
    //r
    M_restingConditions.at (8) = parameterList.get ("r", 2.347e-8);
    //s
    M_restingConditions.at (9) = parameterList.get ("s", 0.999998);
    //Xr1
    M_restingConditions.at (10) = parameterList.get ("Xr1", 0.0165);
    //Xr2
    M_restingConditions.at (11) = parameterList.get ("Xr2", 0.473);
    //Xs
    M_restingConditions.at (12) = parameterList.get ("Xs", 0.0174);
    //Nai
    M_restingConditions.at (13) = parameterList.get ("Nai", 10.132);
    //Ki
    M_restingConditions.at (14) = parameterList.get ("Ki", 138.52);
    //Cai
    M_restingConditions.at (15) = parameterList.get ("Cai", 0.000153);
    //Cass
    M_restingConditions.at (16) = parameterList.get ("Cass", -85.423);
    //Casr
    M_restingConditions.at (17) = parameterList.get ("Casr", 4.272);
    //Rprime
    M_restingConditions.at (18) = parameterList.get ("Rprime", 0.8978);
}

IonicTenTusscher06::IonicTenTusscher06 ( const IonicTenTusscher06& model )
{
    knak = model.knak;
    KmNa = model.KmNa;
    KmK = model.KmK;
    knaca = model.knaca;
    KmNai = model.KmNai;
    KmCa = model.KmCa;
    ksat = model.ksat;
    n = model.n;


    Ko = model.Ko;
    Cao = model.Cao;
    Nao = model.Nao;


    Bufc = model.Bufc;
    Kbufc = model.Kbufc;
    Bufsr = model.Bufsr;
    Kbufsr = model.Kbufsr;
    Bufss = model.Bufss;
    Kbufss = model.Kbufss;

    Vmaxup = model.Vmaxup;
    Kup = model.Kup;
    Vrel = model.Vrel;
    k1_ = model.k1_;
    k2_ = model.k2_;
    k3 = model.k3;
    k4 = model.k4;
    EC = model.EC;
    maxsr = model.maxsr;
    minsr = model.minsr;
    Vleak = model.Vleak;
    Vxfer = model.Vxfer;


    pKNa = model.pKNa;

    F = model.F;
    R = model.R;
    T = model.T;


    Gkr = model.Gkr;
    Gks = model.Gks;
    GK1 = model.GK1;
    Gto = model.Gto;
    GNa = model.GNa;
    GbNa = model.GbNa;
    GCaL = model.GCaL;
    GbCa = model.GbCa;
    GpCa = model.GpCa;
    KpCa = model.KpCa;
    GpK = model.GpK;


    Vc = model.Vc;
    Vsr = model.Vsr;
    Vss = model.Vss;

    computeRTONF();
    computeinverseVcF2();
    computeinverseVcF();
    computeinversevssF2();

    flag = model.flag;

    M_numberOfEquations = model.M_numberOfEquations;
    M_numberOfGatingVariables = model.M_numberOfGatingVariables;
    M_restingConditions = model.M_restingConditions;

}

// ===================================================
//! Operator
// ===================================================
IonicTenTusscher06& IonicTenTusscher06::operator= ( const IonicTenTusscher06& model )
{
    knak = model.knak;
    KmNa = model.KmNa;
    KmK = model.KmK;
    knaca = model.knaca;
    KmNai = model.KmNai;
    KmCa = model.KmCa;
    ksat = model.ksat;
    n = model.n;


    Ko = model.Ko;
    Cao = model.Cao;
    Nao = model.Nao;


    Bufc = model.Bufc;
    Kbufc = model.Kbufc;
    Bufsr = model.Bufsr;
    Kbufsr = model.Kbufsr;
    Bufss = model.Bufss;
    Kbufss = model.Kbufss;

    Vmaxup = model.Vmaxup;
    Kup = model.Kup;
    Vrel = model.Vrel;
    k1_ = model.k1_;
    k2_ = model.k2_;
    k3 = model.k3;
    k4 = model.k4;
    EC = model.EC;
    maxsr = model.maxsr;
    minsr = model.minsr;
    Vleak = model.Vleak;
    Vxfer = model.Vxfer;


    pKNa = model.pKNa;

    F = model.F;
    R = model.R;
    T = model.T;


    Gkr = model.Gkr;
    Gks = model.Gks;
    GK1 = model.GK1;
    Gto = model.Gto;
    GNa = model.GNa;
    GbNa = model.GbNa;
    GCaL = model.GCaL;
    GbCa = model.GbCa;
    GpCa = model.GpCa;
    KpCa = model.KpCa;
    GpK = model.GpK;


    Vc = model.Vc;
    Vsr = model.Vsr;
    Vss = model.Vss;

    computeRTONF();
    computeinverseVcF2();
    computeinverseVcF();
    computeinversevssF2();

    flag = model.flag;

    M_numberOfEquations = model.M_numberOfEquations;
    M_numberOfGatingVariables = model.M_numberOfGatingVariables;
    M_restingConditions = model.M_restingConditions;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
void IonicTenTusscher06::computeGatingRhs ( const   std::vector<Real>&  v,
                                            std::vector<Real>& rhs )
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real f2 = v[6];
    Real fcass = v[7];
    Real r = v[8];
    Real s = v[9];
    Real xr1 = v[10];
    Real xr2 = v[11];
    Real xs = v[12];
    Real Nai = v[13];
    Real Ki = v[14];
    Real Cai = v[15];
    Real CaSS = v[16];
    Real CaSR = v[17];
    Real RR = v[18];

    //m
    rhs[0] = dM (V, m);
    //h
    rhs[1] = dH (V, h);
    //j
    rhs[2] = dJ (V, j);
    //d
    rhs[3] = dD (V, d);
    //f
    rhs[4] = dF (V, f);
    //f2
    rhs[5] = dF2 (V, f2);
    //fCass
    rhs[6] = dFCaSS (V, fcass);
    //r
    rhs[7] = dR (V, r);
    //s
    rhs[8] = dS (V, s);
    //Xr1
    rhs[9] = dXr1 (V, xr1);
    //Xr2
    rhs[10] = dXr2 (V, xr2);
    //Xs
    rhs[11] = dXs (V, xs);
    //Nai
    rhs[12] = dNai (V, m, h, j, Nai, Cai);
    //Ki
    rhs[13] = dKi (V, r, s, xr1, xr2, xs, Ki, Nai);
    //Cai
    rhs[14] = dCai (V, Nai, Cai, CaSR, CaSS);
    //Cass
    rhs[15] = dCaSS (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass);
    //Casr
    rhs[16] =  dCaSR (Cai, CaSR, CaSS, RR);
    //Rprime
    rhs[17] = dRR (CaSR, CaSS, RR);
}

void IonicTenTusscher06::computeNonGatingRhs ( const   std::vector<Real>&  /*v*/,
                                               std::vector<Real>& /*rhs*/ )
{
}

void IonicTenTusscher06::computeRhs ( const   std::vector<Real>&  v,
                                      std::vector<Real>& rhs )
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real f2 = v[6];
    Real fcass = v[7];
    Real r = v[8];
    Real s = v[9];
    Real xr1 = v[10];
    Real xr2 = v[11];
    Real xs = v[12];
    Real Nai = v[13];
    Real Ki = v[14];
    Real Cai = v[15];
    Real CaSS = v[16];
    Real CaSR = v[17];
    Real RR = v[18];

    //V
    rhs[0] = - Itot (V, m, h, j, d, f, f2, fcass, r, s, xr1, xr2, xs, Nai, Ki, Cai, CaSS );
    //m
    rhs[1] = dM (V, m);
    //h
    rhs[2] = dH (V, h);
    //j
    rhs[3] = dJ (V, j);
    //d
    rhs[4] = dD (V, d);
    //f
    rhs[5] = dF (V, f);
    //f2
    rhs[6] = dF2 (V, f2);
    //fCass
    rhs[7] = dFCaSS (V, fcass);
    //r
    rhs[8] = dR (V, r);
    //s
    rhs[9] = dS (V, s);
    //Xr1
    rhs[10] = dXr1 (V, xr1);
    //Xr2
    rhs[11] = dXr2 (V, xr2);
    //Xs
    rhs[12] = dXs (V, xs);
    //Nai
    rhs[13] = dNai (V, m, h, j, Nai, Cai);
    //Ki
    rhs[14] = dKi (V, r, s, xr1, xr2, xs, Ki, Nai);
    //Cai
    rhs[15] = dCai (V, Nai, Cai, CaSR, CaSS);
    //Cass
    rhs[16] = dCaSS (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass);
    //Casr
    rhs[17] =  dCaSR (Cai, CaSR, CaSS, RR);
    //Rprime
    rhs[18] = dRR (CaSR, CaSS, RR);

}


Real IonicTenTusscher06::computeLocalPotentialRhs ( const std::vector<Real>& v)
{
    Real dPotential (0.0);

    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real f2 = v[6];
    Real fCass = v[7];
    Real r = v[8];
    Real s = v[9];
    Real Xr1 = v[10];
    Real Xr2 = v[11];
    Real Xs = v[12];
    Real Nai = v[13];
    Real Ki = v[14];
    Real Cai = v[15];
    Real Cass = v[16];

    dPotential = -Itot (V, m, h, j, d, f, f2, fCass, r, s, Xr1, Xr2, Xs, Nai, Ki, Cai, Cass );

    return dPotential;
}


void IonicTenTusscher06::computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt )
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real f2 = v[6];
    Real fcass = v[7];
    Real r = v[8];
    Real s = v[9];
    Real xr1 = v[10];
    Real xr2 = v[11];
    Real xs = v[12];
    Real Nai = v[13];
    Real Ki = v[14];
    Real Cai = v[15];
    Real CaSS = v[16];
    Real CaSR = v[17];
    Real RR = v[18];

    v[1] = M_INF (V) - ( M_INF (V) - m ) * std::exp (- dt / TAU_M (V) );
    v[2] = H_INF (V) - ( H_INF (V) - h ) * std::exp (- dt / TAU_H (V) );
    v[3] = J_INF (V) - ( J_INF (V) - j ) * std::exp (- dt / TAU_J (V) );
    v[4] = D_INF (V) - ( D_INF (V) - d ) * std::exp (- dt / TAU_D (V) );
    v[5] = F_INF (V) - ( F_INF (V) - f ) * std::exp (- dt / TAU_F (V) );
    v[6] = F2_INF (V) - ( F2_INF (V) - f2 ) * std::exp (- dt / TAU_F2 (V) );
    v[7] = FCaSS_INF (CaSS) - ( FCaSS_INF (CaSS) - fcass ) * std::exp ( -dt / TAU_FCaSS (CaSS) );
    v[8] = R_INF (V) - ( R_INF (V) - r ) * std::exp (- dt / TAU_R (V) );
    v[9] = S_INF (V) - ( S_INF (V) - s ) * std::exp (- dt / TAU_S (V) );
    v[10] = Xr1_INF (V) - ( Xr1_INF (V) - xr1 ) * std::exp (- dt / TAU_Xr1 (V) );
    v[11] = Xr2_INF (V) - ( Xr2_INF (V) - xr2 ) * std::exp (- dt / TAU_Xr2 (V) );
    v[12] = Xs_INF (V) - ( Xs_INF (V) - xs ) * std::exp (- dt / TAU_Xs (V) );
    v[13] = solveNai (V, m, h, j, Nai, Cai, dt);
    v[14] = solveKi (V, r, s, xr1, xr2, xs, Nai, Ki, dt);
    v[15] =  solveCai (V, Nai, Cai, CaSR, CaSS, dt);
    v[16] = solveCaSS (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, dt);
    v[17] = solveCaSR (Cai, CaSR, CaSS, RR, dt);
    v[18] = solveRR (CaSR, CaSS, RR, dt);


}

void IonicTenTusscher06::showMe()
{
    std::cout << "\n\n************************************";
    std::cout << "\n\tHi, I'm the Ten Tusscher model";
    std::cout << "\n\t I've so many parameters that I don't think it's a good idea to display them all\n\n";
    std::cout << "\n\tPlease use a getter, or implement this method otherwise.";
    std::cout << "\n************************************\n\n";
}

void IonicTenTusscher06::solveOneStep (std::vector<Real>& v, Real dt)
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real f2 = v[6];
    Real fcass = v[7];
    Real r = v[8];
    Real s = v[9];
    Real xr1 = v[10];
    Real xr2 = v[11];
    Real xs = v[12];
    Real Nai = v[13];
    Real Ki = v[14];
    Real Cai = v[15];
    Real CaSS = v[16];
    Real CaSR = v[17];
    Real RR = v[18];

    computeGatingVariablesWithRushLarsen (v, dt);

    v[0] = solveV (V, m, h, j, d, f, f2, fcass, r, s, xr1, xr2, xs, Nai, Ki, Cai, CaSS, dt);
    v[13] = solveNai (V, m, h, j, Nai, Cai, dt);
    v[14] = solveKi (V, r, s, xr1, xr2, xs, Nai, Ki, dt);
    v[15] =  solveCai (V, Nai, Cai, CaSR, CaSS, dt);
    v[16] = solveCaSS (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, dt);
    v[17] = solveCaSR (Cai, CaSR, CaSS, RR, dt);
    v[18] = solveRR (CaSR, CaSS, RR, dt);




}


void IonicTenTusscher06::showCurrents (std::vector<Real>& v)
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real f2 = v[6];
    Real fcass = v[7];
    Real r = v[8];
    Real s = v[9];
    Real xr1 = v[10];
    Real xr2 = v[11];
    Real xs = v[12];
    Real Nai = v[13];
    Real Ki = v[14];
    Real Cai = v[15];
    Real CaSS = v[16];

    showCurrents ( V,  m,  h,  j,  d,  f,  f2,  fcass, r,  s,  xr1,  xr2,  xs,  Nai,  Ki, Cai,  CaSS);
}

void IonicTenTusscher06::showCurrents (Real V, Real m, Real h, Real j, Real d, Real f, Real f2, Real fcass,
                                       Real r, Real s, Real xr1, Real xr2, Real xs, Real Nai, Real Ki,
                                       Real Cai, Real CaSS)
{
    std::cout << "\nIKr = "   << IKr (V, xr1, xr2, Ki);
    std::cout << "\nIKs = "   << IKs (V, xs, Ki, Nai);
    std::cout << "\nIK1 = "   << IK1 (V, Ki);
    std::cout << "\nIto = "   << Ito (V, r, s, Ki);
    std::cout << "\nINa = "   << INa (V, m, h, j, Nai);
    std::cout << "\nIbNa = "  << IbNa (V, Nai);
    std::cout << "\nICaL = "  << ICaL (V, d, f, f2, fcass, CaSS);
    std::cout << "\nIbCa = "  << IbCa (V, Cai);
    std::cout << "\nINaK = "  << INaK (V, Nai);
    std::cout << "\nINaCa = " << INaCa (V, Nai, Cai);
    std::cout << "\nIpCa = "  << IpCa (Cai);
    std::cout << "\nIKpK = "  << IpK (V, Ki);
    std::cout << "\n";
}


}

