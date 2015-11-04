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
  @brief Ionic model Luo-Rudy I
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#include <lifev/electrophysiology/solver/IonicModels/IonicLuoRudyI.hpp>


namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
IonicLuoRudyI::IonicLuoRudyI()  :
    super       (  8 , 6  ),
    M_ENa       ( 54.4 ),
    M_gNa       ( 23.0 ),
    M_gsi       ( 0.09 ),
    M_K0        (  5.4 ),
    M_EK        (  -77 ),
    M_EK1       ( -87.26),
    M_EKp       ( M_EK1 ),
    M_gKp       ( 0.0183),
    M_gb        ( 0.03921)
{
    M_gK = computeGK (M_K0);
    M_gK1 = computeGK1 (M_K0);

    //V
    M_restingConditions.at (0) = -84.0;
    //m
    M_restingConditions.at (1) = minf ( M_restingConditions.at (0) );
    //h
    M_restingConditions.at (2) = hinf ( M_restingConditions.at (0) );
    //j
    M_restingConditions.at (3) = jinf ( M_restingConditions.at (0) );
    //d
    M_restingConditions.at (4) = dinf ( M_restingConditions.at (0) );
    //f
    M_restingConditions.at (5) = finf ( M_restingConditions.at (0) );
    //X
    M_restingConditions.at (6) = Xinf ( M_restingConditions.at (0) );
    //Ca
    M_restingConditions.at (7) = 2e-4;
}

IonicLuoRudyI::IonicLuoRudyI ( Teuchos::ParameterList& parameterList     )   :
    super       ( 8, 6 )
{
    M_ENa        =  parameterList.get ("ENa", 54.4  );
    M_gNa        =  parameterList.get ("gNa", 23.3  );
    M_gsi        =  parameterList.get ("gsi", 0.09  );
    M_K0        =  parameterList.get ("K0", 5.4  );
    M_EK        =  parameterList.get ("EK", -77.0  );
    M_EK1        =  parameterList.get ("EK1", -87.26  );
    M_EKp        =  parameterList.get ("EKp", -87.26  );
    M_gKp        =  parameterList.get ("gKp", 0.0183  );
    M_gb        =  parameterList.get ("gb", 0.03921  );

    M_gK = computeGK (M_K0);
    M_gK1 = computeGK1 (M_K0);

    //V
    M_restingConditions.at (0) = parameterList.get ("V0", -84.0  );
    //m
    M_restingConditions.at (1) = minf ( M_restingConditions.at (0) );
    //h
    M_restingConditions.at (2) = hinf ( M_restingConditions.at (0) );
    //j
    M_restingConditions.at (3) = jinf ( M_restingConditions.at (0) );
    //d
    M_restingConditions.at (4) = dinf ( M_restingConditions.at (0) );
    //f
    M_restingConditions.at (5) = finf ( M_restingConditions.at (0) );
    //X
    M_restingConditions.at (6) = Xinf ( M_restingConditions.at (0) );
    //Ca
    M_restingConditions.at (7) = parameterList.get ("Ca0", 2e-4  );
}

IonicLuoRudyI::IonicLuoRudyI ( const IonicLuoRudyI& model )
{
    M_ENa        =  model.M_ENa;
    M_gNa        =  model.M_gNa;
    M_gsi        =  model.M_gsi;
    M_K0        =  model.M_K0;
    M_EK        =  model.M_EK;
    M_EK1        =  model.M_EK1;
    M_EKp        =  model.M_EKp;
    M_gKp        =  model.M_gKp;
    M_gb        =  model.M_gb;

    M_gK = computeGK (M_K0);
    M_gK1 = computeGK1 (M_K0);
    M_numberOfEquations = model.M_numberOfEquations;
    M_numberOfGatingVariables = model.M_numberOfGatingVariables;
    M_restingConditions = model.M_restingConditions;

}

// ===================================================
//! Operator
// ===================================================
IonicLuoRudyI& IonicLuoRudyI::operator= ( const IonicLuoRudyI& model )
{
    M_ENa        =  model.M_ENa;
    M_gNa        =  model.M_gNa;
    M_gsi        =  model.M_gsi;
    M_K0        =  model.M_K0;
    M_EK        =  model.M_EK;
    M_EK1        =  model.M_EK1;
    M_EKp        =  model.M_EKp;
    M_gKp        =  model.M_gKp;
    M_gb        =  model.M_gb;

    M_gK = computeGK (M_K0);
    M_gK1 = computeGK1 (M_K0);
    M_numberOfEquations = model.M_numberOfEquations;
    M_numberOfGatingVariables = model.M_numberOfGatingVariables;
    M_restingConditions = model.M_restingConditions;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
void IonicLuoRudyI::computeGatingRhs ( const   std::vector<Real>&  v,
                                       std::vector<Real>& rhs )
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real X = v[6];
    Real Ca = v[7];

    //m
    rhs[0] = dm (V, m);
    //h
    rhs[1] = dh (V, h);
    //j
    rhs[2] = dj (V, j);
    //d
    rhs[3] = dd (V, d);
    //f
    rhs[4] = df (V, f);
    //X
    rhs[5] = dX (V, X);
    //Ca
    rhs[6] = dCa (V, d, f, Ca);
}

void IonicLuoRudyI::computeNonGatingRhs ( const   std::vector<Real>&  v,
                                          std::vector<Real>& rhs )
{
    Real V = v[0];
    Real d = v[4];
    Real f = v[5];
    Real Ca = v[7];

    //Ca
    rhs[0] = dCa (V, d, f, Ca);
}

void IonicLuoRudyI::computeRhs ( const   std::vector<Real>&  v,
                                 std::vector<Real>& rhs )
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real X = v[6];
    Real Ca = v[7];

    //V
    rhs[0] = - Itot (V, m , h , j , d, f, X, Ca);
    //m
    rhs[1] = dm (V, m);
    //h
    rhs[2] = dh (V, h);
    //j
    rhs[3] = dj (V, j);
    //d
    rhs[4] = dd (V, d);
    //f
    rhs[5] = df (V, f);
    //X
    rhs[6] = dX (V, X);
    //Ca
    rhs[7] = dCa (V, d, f, Ca);
}

void IonicLuoRudyI::computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt )
{
    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real X = v[6];

    v[1] = minf (V) - ( minf (V) - m ) * std::exp (- dt / tm (V) );
    v[2] = hinf (V) - ( hinf (V) - h ) * std::exp (- dt / th (V) );
    v[3] = jinf (V) - ( jinf (V) - j ) * std::exp (- dt / tj (V) );
    v[4] = dinf (V) - ( dinf (V) - d ) * std::exp (- dt / td (V) );
    v[5] = finf (V) - ( finf (V) - f ) * std::exp (- dt / tf (V) );
    v[6] = Xinf (V) - ( Xinf (V) - X ) * std::exp (- dt / tX (V) );

}

Real IonicLuoRudyI::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    Real dPotential (0.0);

    Real V = v[0];
    Real m = v[1];
    Real h = v[2];
    Real j = v[3];
    Real d = v[4];
    Real f = v[5];
    Real X = v[6];
    Real Ca = v[7];

    dPotential = - Itot (V, m , h , j , d, f, X, Ca);

    return dPotential;
}


void IonicLuoRudyI::showMe()
{
    std::cout << "\n\n************************************";
    std::cout << "\n\tHi, I'm the Luo Rudy Phase I model";

    std::cout << "\nENa: " << M_ENa;
    std::cout << "\ngNa: " << M_gNa;
    std::cout << "\ngsi: " << M_gsi;
    std::cout << "\nK0: " << M_K0;
    std::cout << "\nEK: " << M_EK;
    std::cout << "\nEK1: " << M_EK1;
    std::cout << "\nEKp: " << M_EKp;
    std::cout << "\ngKp: " << M_gKp;
    std::cout << "\ngb: " << M_gb;
    std::cout << "\n************************************\n\n";
}


}

