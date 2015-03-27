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
  @brief Ionic model of Noble for Purkinje cells
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#include <lifev/electrophysiology/solver/IonicModels/IonicNoblePurkinje.hpp>


namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================


IonicNoblePurkinje::IonicNoblePurkinje()  :
    super       ( 4, 3 ),
    M_gi (0.14),
    M_vNa (40.0),
    M_vK (-100.0),
    M_Cm (12.0),
    M_Itot (0)
{
    M_restingConditions.at (0) = -80.0;
    M_restingConditions.at (1) = mInf (M_restingConditions[0]);
    M_restingConditions.at (2) = nInf (M_restingConditions[0]);
    M_restingConditions.at (3) = hInf (M_restingConditions[0]);

}

IonicNoblePurkinje::IonicNoblePurkinje ( Teuchos::ParameterList& parameterList     )   :
    super       ( 4, 3 )
{
    M_gi        =  parameterList.get ("gi",      0.14     );
    M_vNa        =  parameterList.get ("vNa",      40.0     );
    M_vK        =  parameterList.get ("vK",      -100.0     );
    M_Cm        =  parameterList.get ("Cm",      12.0     );
    M_Itot        =  0;

}

IonicNoblePurkinje::IonicNoblePurkinje ( const IonicNoblePurkinje& model )
{

    M_gi        =  model.M_gi;
    M_vNa        =  model.M_vNa;
    M_vK        =  model.M_vK;
    M_Cm        =  model.M_Cm;
    M_Itot      =  model.M_Itot;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;

}

// ===================================================
//! Operator
// ===================================================
IonicNoblePurkinje& IonicNoblePurkinje::operator= ( const IonicNoblePurkinje& model )
{
    M_gi        =  model.M_gi;
    M_vNa        =  model.M_vNa;
    M_vK        =  model.M_vK;
    M_Cm        =  model.M_Cm;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;


    return      *this;
}


// ===================================================
//! Methods
// ===================================================
void IonicNoblePurkinje::computeGatingRhs ( const   std::vector<Real>&  v,
                                            std::vector<Real>& rhs )
{

    Real V = v[0];
    Real M = v[1];
    Real N = v[2];
    Real H = v[3];

    Real alpham = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.1, -1, -15, -48);
    Real betam  = GeneralFunctionAlphaAndBeta (V, 0, 1, -0.12, -1, 5, -8);
    Real alphah = GeneralFunctionAlphaAndBeta (V, 0.17, -20, 0, 0, 1, -90);
    Real betah  = GeneralFunctionAlphaAndBeta (V, 1, 0, 0, 1, -10, -42);
    Real alphan = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.0001, -1, -10, -50);
    Real betan  = GeneralFunctionAlphaAndBeta (V, 0.002, -80, 0, 0, 1, -90);

    rhs[0] = alpham * (1 - M) - betam * M;
    rhs[1] = alphan * (1 - N) - betan * N;
    rhs[2] = alphah * (1 - H) - betah * H;
}

void IonicNoblePurkinje::computeRhs ( const   std::vector<Real>&  v,
                                      std::vector<Real>& rhs )
{

    Real V = v[0];
    Real M = v[1];
    Real N = v[2];
    Real H = v[3];

    Real alpham = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.1, -1, -15, -48);
    Real betam  = GeneralFunctionAlphaAndBeta (V, 0, 1, -0.12, -1, 5, -8);
    Real alphah = GeneralFunctionAlphaAndBeta (V, 0.17, -20, 0, 0, 1, -90);
    Real betah  = GeneralFunctionAlphaAndBeta (V, 1, 0, 0, 1, -10, -42);
    Real alphan = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.0001, -1, -10, -50);
    Real betan  = GeneralFunctionAlphaAndBeta (V, 0.002, -80, 0, 0, 1, -90);

    Real gK1 = 1.2 * std::exp (- (V + 90.0) / 50.0) + 0.015 * std::exp ( (V + 90.0) / 60.0);
    Real gK2 = 1.2 * N * N * N * N;
    Real gNa = 400 * M * M * M * H + M_gi;

    M_Itot = 1 / M_Cm * (-gNa * (V - M_vNa) - (gK1 + gK2) * (V - M_vK) );
    rhs[0] = 1 / M_Cm * (-gNa * (V - M_vNa) - (gK1 + gK2) * (V - M_vK) );
    rhs[1] = alpham * (1 - M) - betam * M;
    rhs[2] = alphan * (1 - N) - betan * N;
    rhs[3] = alphah * (1 - H) - betah * H;
}


Real IonicNoblePurkinje::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    Real dPotential (0.0);

    Real V = v[0];
    Real M = v[1];
    Real N = v[2];
    Real H = v[3];

    /*
    Real alpham = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.1, -1, -15, -48);
    Real betam  = GeneralFunctionAlphaAndBeta (V, 0, 1, -0.12, -1, 5, -8);
    Real alphah = GeneralFunctionAlphaAndBeta (V, 0.17, -20, 0, 0, 1, -90);
    Real betah  = GeneralFunctionAlphaAndBeta (V, 1, 0, 0, 1, -10, -42);
    Real alphan = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.0001, -1, -10, -50);
    Real betan  = GeneralFunctionAlphaAndBeta (V, 0.002, -80, 0, 0, 1, -90);
     */

    Real gK1 = 1.2 * std::exp (- (V + 90.0) / 50.0) + 0.015 * std::exp ( (V + 90.0) / 60.0);
    Real gK2 = 1.2 * N * N * N * N;
    Real gNa = 400 * M * M * M * H + M_gi;

    M_Itot = 1 / M_Cm * (-gNa * (V - M_vNa) - (gK1 + gK2) * (V - M_vK) );
    dPotential = 1 / M_Cm * (-gNa * (V - M_vNa) - (gK1 + gK2) * (V - M_vK) );

    return dPotential;
}

void IonicNoblePurkinje::computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt )
{
    Real V = v[0];


    Real alpham = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.1, -1, -15, -48);
    Real betam  = GeneralFunctionAlphaAndBeta (V, 0, 1, -0.12, -1, 5, -8);
    Real alphah = GeneralFunctionAlphaAndBeta (V, 0.17, -20, 0, 0, 1, -90);
    Real betah  = GeneralFunctionAlphaAndBeta (V, 1, 0, 0, 1, -10, -42);
    Real alphan = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.0001, -1, -10, -50);
    Real betan  = GeneralFunctionAlphaAndBeta (V, 0.002, -80, 0, 0, 1, -90);

    Real taum = alpham + betam;
    Real taun = alphan + betan;
    Real tauh = alphah + betah;

    Real mInf = alpham / (taum);
    Real nInf = alphan / (taun);
    Real hInf = alphah / (tauh);

    v[1] = mInf + (v[1] - mInf) * exp (-dt * taum);
    v[2] = nInf + (v[2] - nInf) * exp (-dt * taun);
    v[3] = hInf + (v[3] - hInf) * exp (-dt * tauh);

}

void IonicNoblePurkinje::showMe()
{

    std::cout << "\n\tHi, I'm the Noble model\n\t Good luck\n\n";
}


}

