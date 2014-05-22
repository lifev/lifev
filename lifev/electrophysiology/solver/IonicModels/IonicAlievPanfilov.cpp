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
  @brief Ionic model of Aliev-Panfilov
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */

#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>


//#include <Teuchos_RCP.hpp>
//#include <Teuchos_ParameterList.hpp>
//#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.
// ===================================================
//! Constructors
// ===================================================
IonicAlievPanfilov::IonicAlievPanfilov()    :
    super       ( 2   ),
    M_mu1       ( 0.12 ),
    M_mu2       ( 0.3 ),
    M_k         ( 8.0 ),
    M_a         ( 0.1 ),
    M_epsilon   ( 0.01 )
{
    M_restingConditions.at (0) = 0.0;
    M_restingConditions.at (1) = 0.16;
}

IonicAlievPanfilov::IonicAlievPanfilov ( Teuchos::ParameterList& parameterList   )   :
    super       ( 2 )
{
	setup ( parameterList );

    M_restingConditions.at (0) = 0.0;
    M_restingConditions.at (1) = 0.16;
}

IonicAlievPanfilov::IonicAlievPanfilov ( const IonicAlievPanfilov& model )
{

    M_mu1       =  model.M_mu1;
    M_mu2       =  model.M_mu2;
    M_k         =  model.M_k;
    M_a         =  model.M_a;
    M_epsilon   =  model.M_epsilon;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;
}

// ===================================================
//! Operator
// ===================================================
IonicAlievPanfilov& IonicAlievPanfilov::operator= ( const IonicAlievPanfilov& model )
{
    M_mu1       =  model.M_mu1;
    M_mu2       =  model.M_mu2;
    M_k         =  model.M_k;
    M_a         =  model.M_a;
    M_epsilon   =  model.M_epsilon;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;

    return      *this;
}

void IonicAlievPanfilov::setup ( Teuchos::ParameterList& parameterList )
{
    M_mu1       =  parameterList.get ("mu1", 0.12);
    M_mu2       =  parameterList.get ("mu2", 0.3);
    M_k         =  parameterList.get ("k", 8.0);
    M_a         =  parameterList.get ("a", 0.1);
    M_epsilon   =  parameterList.get ("epsilon", 0.01);
}
// ===================================================
//! Methods
// ===================================================
//Only gating variables
void IonicAlievPanfilov::computeGatingRhs (    const   std::vector<Real>&  v,
                                               std::vector<Real>& rhs )
{

    Real dr = - ( M_epsilon + M_mu1 * v[1] / ( M_mu2 + v[0] ) ) * ( v[1] + M_k * v[0] * ( v[0] - M_a  - 1.0 ) );

    rhs[0] = dr;

}

//Potential and gating variables
void IonicAlievPanfilov::computeRhs (    const   std::vector<Real>&  v,
                                         std::vector<Real>& rhs )
{

    Real dr = - ( M_epsilon + M_mu1 * v[1] / ( M_mu2 + v[0] ) ) * ( v[1] + M_k * v[0] * ( v[0] - M_a  - 1.0 ) );
    Real dV = - M_k * v[0] * ( v[0] - M_a ) * ( v[0] - 1.0) - v[0] * v[1];

    rhs[0] = dV;
    rhs[1] = dr;

}


Real IonicAlievPanfilov::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    return ( - M_k * v[0] * ( v[0] - M_a ) * ( v[0] - 1.0) - v[0] * v[1] );
}


void IonicAlievPanfilov::showMe()
{
    std::cout << "\n\n\t\tIonicAlievPanfilov Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "mu1: " << this->Mu1() << std::endl;
    std::cout << "mu2: " << this->Mu2() << std::endl;
    std::cout << "k: " << this->K() << std::endl;
    std::cout << "a: " << this->A() << std::endl;
    std::cout << "epsilon: " << this->Epsilon() << std::endl;

    std::cout << "\n\t\t End of IonicAlievPanfilov Informations\n\n\n";

}


}

