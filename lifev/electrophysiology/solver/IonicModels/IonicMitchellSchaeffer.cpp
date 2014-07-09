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
  @brief Ionic model based on Mitchell-Schaeffer model.
  @date 03-2013
  @author Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

  @contributors
  @mantainer Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
  @last update 03-2013
 */
#include <lifev/electrophysiology/solver/IonicModels/IonicMitchellSchaeffer.hpp>


namespace LifeV
{
// ===================================================
//! Constructors
// ===================================================
IonicMitchellSchaeffer::IonicMitchellSchaeffer()    :
    super       ( 2 ),
    M_vGate     ( 0.13 ),
    M_tauClose  ( 150.0 ),
    M_tauOpen   ( 120.0 ),
    M_tauIn     ( 0.3 ),
    M_tauOut    ( 6.0 )
{
    M_restingConditions.at (0) = 0.0;
    M_restingConditions.at (1) = 1.0;
}

IonicMitchellSchaeffer::IonicMitchellSchaeffer ( Teuchos::ParameterList& parameterList ) :
    super       ( 2 )
{
    M_vGate       =  parameterList.get ("vGate", 0.13);
    M_tauClose    =  parameterList.get ("tauClose", 150.0);
    M_tauOpen     =  parameterList.get ("tauOpen", 120.0);
    M_tauIn       =  parameterList.get ("tauIn", 0.3);
    M_tauOut      =  parameterList.get ("tauOut", 6.0);

    M_restingConditions.at (0) = 0.0;
    M_restingConditions.at (1) = 1.0;
}

IonicMitchellSchaeffer::IonicMitchellSchaeffer ( const IonicMitchellSchaeffer& model )
{
    M_vGate     =  model.M_vGate;
    M_tauClose  =  model.M_tauClose;
    M_tauOpen   =  model.M_tauOpen;
    M_tauIn     =  model.M_tauIn;
    M_tauOut    =  model.M_tauOut;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;
}

// ===================================================
//! Operator
// ===================================================
IonicMitchellSchaeffer& IonicMitchellSchaeffer::operator= ( const IonicMitchellSchaeffer& model )
{
    M_vGate     =  model.M_vGate;
    M_tauClose  =  model.M_tauClose;
    M_tauOpen   =  model.M_tauOpen;
    M_tauIn     =  model.M_tauIn;
    M_tauOut    =  model.M_tauOut;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;

    return *this;
}


// ===================================================
//! Methods
// ===================================================
//Only gating variables
void IonicMitchellSchaeffer::computeGatingRhs ( const std::vector<Real>&  v,
                                                std::vector<Real>& rhs )
{

    rhs[0] = computeLocalGatingRhs ( v );

}

//Potential and gating variables
void IonicMitchellSchaeffer::computeRhs (const   std::vector<Real>&  v,
                                         std::vector<Real>& rhs )
{
    rhs[0] = computeLocalPotentialRhs ( v );
    rhs[1] = computeLocalGatingRhs ( v );
}



Real IonicMitchellSchaeffer::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    return ( - ( v[1] / M_tauIn ) * v[0] * v[0] * ( v[0] - 1 ) - v[0] / M_tauOut );
}

Real IonicMitchellSchaeffer::computeLocalGatingRhs ( const std::vector<Real>& v )
{
    if (v[0] <= M_vGate)
    {
        return  (  ( 1 - v[1] ) / M_tauOpen );
    }
    else
    {
        return ( - v[1] / M_tauClose );
    }
}

void IonicMitchellSchaeffer::showMe()
{
    std::cout << "\n\n\t\tIonicMitchellSchaeffer Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "vGate: " << this->vGate() << std::endl;
    std::cout << "tauClose: " << this->tauClose() << std::endl;
    std::cout << "tauOpen: " << this->tauOpen() << std::endl;
    std::cout << "tauIn: " << this->tauIn() << std::endl;
    std::cout << "tauOut: " << this->tauOut() << std::endl;

    std::cout << "\n\t\t End of IonicMitchellSchaeffer Informations\n\n\n";

}


}

