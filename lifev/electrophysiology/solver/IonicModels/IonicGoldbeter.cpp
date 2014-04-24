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
  @brief Intracellular Calcium model from Goldbeter et al. (1990). By "potential" we
  @brief refer to the cytosolic calcium concentration, whereas the gating variable
  @brief represents the sarcoplasmic calcium concentration

  @date 09-2013

  @author Ricardo Ruiz <ricardo.ruizbaier@unil.ch>

  @contributors
  @mantainer Ricardo Ruiz <ricardo.ruizbaier@epfl.ch>
  @last update 09-2013
 */

#include <lifev/electrophysiology/solver/IonicModels/IonicGoldbeter.hpp>


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
// ===================================================
//! Constructors
// ===================================================
IonicGoldbeter::IonicGoldbeter()    :
    super       ( 2  ),
    M_nu1       ( 1.58 ),
    M_nu2       ( 16.0 ),
    M_nu3       ( 91.0 ),
    M_nu4       ( 2.0  ),
    M_nu5       ( 0.2  ),
    M_k1       ( 0.0 ),
    M_k2       ( 1.0 ),
    M_k3       ( 4.0 ),
    M_k4       ( 0.7481)

{
    M_restingConditions.at (0) = 0.1;
    M_restingConditions.at (1) = 1.6;
}

IonicGoldbeter::IonicGoldbeter ( Teuchos::ParameterList& parameterList   )   :
    super       ( 2 )
{
    M_nu1       =  parameterList.get ("nu1", 1.58);
    M_nu2       =  parameterList.get ("nu2", 16.0);
    M_nu3       =  parameterList.get ("nu3", 1.58);
    M_nu4       =  parameterList.get ("nu4", 16.0);
    M_nu5       =  parameterList.get ("nu5", 1.58);
    M_k1        =  parameterList.get ("k1", 0.0);
    M_k2        =  parameterList.get ("k2", 1.0);
    M_k3        =  parameterList.get ("k3", 4.0);
    M_k4        =  parameterList.get ("k4", 0.7481);

}

IonicGoldbeter::IonicGoldbeter ( const IonicGoldbeter& model )
{

    M_nu1       =  model.M_nu1;
    M_nu2       =  model.M_nu2;
    M_nu3       =  model.M_nu3;
    M_nu4       =  model.M_nu4;
    M_nu5       =  model.M_nu5;
    M_k1         =  model.M_k1;
    M_k2         =  model.M_k2;
    M_k3         =  model.M_k3;
    M_k4         =  model.M_k4;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;
}

// ===================================================
//! Operator
// ===================================================
IonicGoldbeter& IonicGoldbeter::operator= ( const IonicGoldbeter& model )
{
    M_nu1       =  model.M_nu1;
    M_nu2       =  model.M_nu2;
    M_nu3       =  model.M_nu3;
    M_nu4       =  model.M_nu4;
    M_nu5       =  model.M_nu5;
    M_k1         =  model.M_k1;
    M_k2         =  model.M_k2;
    M_k3         =  model.M_k3;
    M_k4         =  model.M_k4;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
//Only sarcoplasmic calcium
void IonicGoldbeter::computeGatingRhs (    const   std::vector<Real>&  v,
                                                          std::vector<Real>& rhs )
{

    Real dr = M_nu2 * std::pow ( v[0], 2.0) / (M_k2 + std::pow (v[0], 2.0) ) - M_nu3 * std::pow (v[0], 4.0) * std::pow (v[1], 2.0) /
              ( (M_k3 + std::pow (v[1], 2.0) ) * (M_k4 + std::pow (v[0], 4.0) ) ) - M_nu5 * v[1];

    rhs[0] = dr;

}

//Both cytosolic (V) and sarcoplasmic calcium (r)
void IonicGoldbeter::computeRhs (    const   std::vector<Real>&  v,
                                                    std::vector<Real>& rhs )
{

    Real dr = M_nu2 * std::pow (v[0], 2.0) / (M_k2 + std::pow (v[0], 2.0) ) - M_nu3 * std::pow (v[0], 4.0) * std::pow (v[1], 2.0) /
              ( (M_k3 + std::pow (v[1], 2.0) ) * (M_k4 + std::pow (v[0], 4.0) ) ) - M_nu5 * v[1];
    Real dV = M_nu1 - M_nu2 * std::pow (v[0], 2.0) / (M_k2 + std::pow (v[0], 2.0) ) + M_nu3 * std::pow (v[0], 4.0) * std::pow (v[1], 2.0) /
              ( (M_k3 + std::pow (v[1], 2.0) ) * (M_k4 + std::pow (v[0], 4.0) ) ) - M_nu4 * v[0];

    rhs[0] = dV;
    rhs[1] = dr;

}


Real IonicGoldbeter::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    return (M_nu1 - M_nu2 * std::pow (v[0], 2.0) / (M_k2 + std::pow (v[0], 2.0) ) + M_nu3 * std::pow (v[0], 4.0) * std::pow (v[1], 2.0) /
            ( (M_k3 + std::pow (v[1], 2.0) ) * (M_k4 + std::pow (v[0], 4.0) ) ) - M_nu4 * v[0] );
}




void IonicGoldbeter::showMe()
{
    std::cout << "\n\n\t\tIntracellularCalciumGoldbeter Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "nu1: " << this->Nu1() << std::endl;
    std::cout << "nu2: " << this->Nu2() << std::endl;
    std::cout << "nu3: " << this->Nu3() << std::endl;
    std::cout << "nu4: " << this->Nu4() << std::endl;
    std::cout << "nu5: " << this->Nu5() << std::endl;


    std::cout << "k1: " << this->K1() << std::endl;
    std::cout << "k2: " << this->K2() << std::endl;
    std::cout << "k3: " << this->K3() << std::endl;
    std::cout << "k4: " << this->K4() << std::endl;

    std::cout << "\n\t\t End of IntracellularCalciumGoldbeter Informations\n\n\n";

}


}

