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
  @brief Ionic model of FitzHugh-Nagumo
  @date 01-2013
  @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

  @contributors
  @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
  @last update 01-2013
 */

#include <lifev/electrophysiology/solver/IonicModels/IonicFitzHughNagumo.hpp>

namespace LifeV
{
// ===================================================
//! Constructors
// ===================================================
IonicFitzHughNagumo::IonicFitzHughNagumo()    :
    super   ( 2   ),
    M_G     ( 1.5 ),
    M_Vth   ( 13. ),
    M_Vp    ( 100. ),
    M_Eta1  ( 4.4 ),
    M_Eta2  ( 0.012 ),
    M_Eta3  ( 1.),
    M_Eta   (M_Eta2 / M_Vp),
    M_Gamma (M_Eta2* M_Eta3)

{
    M_restingConditions.at (0) = 1e-8;
    M_restingConditions.at (1) = 0.3;
}

IonicFitzHughNagumo::IonicFitzHughNagumo ( Teuchos::ParameterList& parameterList   )   :
    super       ( 2 )
{
	setup ( parameterList );

    M_restingConditions.at (0) = 1e-8;
    M_restingConditions.at (1) = 0.3;
}

IonicFitzHughNagumo::IonicFitzHughNagumo ( const IonicFitzHughNagumo& model )
{
    M_G     =   model.M_G;
    M_Vth   =   model.M_Vth;
    M_Vp    =   model.M_Vp;
    M_Eta1  =   model.M_Eta1;
    M_Eta2  =   model.M_Eta2;
    M_Eta3  =   model.M_Eta3;
    M_Eta   =   model.M_Eta;
    M_Gamma =   model.M_Gamma;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;
}

// ===================================================
//! Operator
// ===================================================
IonicFitzHughNagumo& IonicFitzHughNagumo::operator= ( const IonicFitzHughNagumo& model )
{
    M_G     =   model.M_G;
    M_Vth   =   model.M_Vth;
    M_Vp    =   model.M_Vp;
    M_Eta1  =   model.M_Eta1;
    M_Eta2  =   model.M_Eta2;
    M_Eta3  =   model.M_Eta3;
    M_Eta   =   model.M_Eta;
    M_Gamma =   model.M_Gamma;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;

    return      *this;
}


// ===================================================
//! Methods
// ===================================================
void IonicFitzHughNagumo::setup ( Teuchos::ParameterList& parameterList )
{
    M_G     =   parameterList.get ("G", 1.5);
    M_Vth   =   parameterList.get ("Vth", 13.);
    M_Vp    =   parameterList.get ("Vp", 100.);
    M_Eta1  =   parameterList.get ("Eta1", 4.4);
    M_Eta2  =   parameterList.get ("Eta2", 0.012);
    M_Eta3  =   parameterList.get ("Eta3", 1.);
    M_Eta   =   M_Eta2 / M_Vp;
    M_Gamma =   M_Eta2 * M_Eta3;
}


//Only gating variables
void IonicFitzHughNagumo::computeGatingRhs (    const   std::vector<Real>&  v,
                                                std::vector<Real>& rhs )
{
    Real dr = M_Eta * v[0] - M_Gamma * v[1] ;

    rhs[0] = dr;
}

//Potential and gating variables
void IonicFitzHughNagumo::computeRhs (    const   std::vector<Real>&  v,
                                          std::vector<Real>& rhs )
{
    Real dV = - ( M_G * v[0] * ( 1.0 - v[0] / M_Vth ) * ( 1.0 - v[0] / M_Vp ) + M_Eta1 * v[0] * v[1] );
    Real dr = M_Eta * v[0] - M_Gamma * v[1] ;

    rhs[0] = dV;
    rhs[1] = dr;

}
void IonicFitzHughNagumo::computeRhs ( const VectorSmall<2>& v, VectorSmall<2>& rhs)
{
    Real dV = - ( M_G * v[0] * ( 1.0 - v[0] / M_Vth ) * ( 1.0 - v[0] / M_Vp ) + M_Eta1 * v[0] * v[1] );
    Real dr = M_Eta * v[0] - M_Gamma * v[1] ;

    rhs[0] = dV;
    rhs[1] = dr;

}

Real IonicFitzHughNagumo::computeLocalPotentialRhs ( const std::vector<Real>& v)
{
    return ( - ( M_G * v[0] * ( 1.0 - v[0] / M_Vth ) * ( 1.0 - v[0] / M_Vp ) + M_Eta1 * v[0] * v[1] ) );
}

std::vector< std::vector<Real> > IonicFitzHughNagumo::getJac (const std::vector<Real>& v, Real /*h*/)
{
    std::vector< std::vector<Real> > J (2, std::vector<Real> (2, 0.0) );
    J[0][0] = - ( M_G / ( M_Vth * M_Vp ) ) * ( M_Vth * ( M_Vp - 2.0 * v[0] ) + v[0] * ( 3.0 * v[0] - 2.0 * M_Vp ) ) - M_Eta1 * v[1];
    J[0][1] = -M_Eta1 * v[0];
    J[1][0] = M_Eta;
    J[1][1] = -M_Gamma;

    return J;
}

MatrixSmall<2, 2> IonicFitzHughNagumo::getJac (const VectorSmall<2>& v, Real /*h*/)
{
    MatrixSmall<2, 2> J;
    J (0, 0) = - ( M_G / ( M_Vth * M_Vp ) ) * ( M_Vth * ( M_Vp - 2.0 * v[0] ) + v[0] * ( 3.0 * v[0] - 2.0 * M_Vp ) ) - M_Eta1 * v[1];
    J (0, 1) = -M_Eta1 * v[0];
    J (1, 0) = M_Eta;
    J (1, 1) = -M_Gamma;

    return J;
}

matrix_Type IonicFitzHughNagumo::getJac (const vector_Type& v, Real h)
{
    matrix_Type J (v.map(), M_numberOfEquations, false);
    const Int* k = v.blockMap().MyGlobalElements();

    int* Indices = new int[M_numberOfEquations];
    double* Values =  new double[M_numberOfEquations];
    Indices[0] = 0;
    Indices[1] = 1;

    MatrixSmall<2, 2> df;
    VectorSmall<2> y;
    y (0) = v[k[0]];
    y (1) = v[k[1]];

    df = getJac (y, h);

    Values[0] = df (0, 0);
    Values[1] = df (0, 1);
    J.matrixPtr()->InsertGlobalValues (0, M_numberOfEquations, Values, Indices);
    Values[0] = df (1, 0);
    Values[1] = df (1, 1);
    J.matrixPtr()->InsertGlobalValues (1, M_numberOfEquations, Values, Indices);

    J.globalAssemble();

    delete[] Indices;
    delete[] Values;

    return J;
}

void IonicFitzHughNagumo::showMe()
{
    std::cout << "\n\n\t\tIonicFitzHughNagumo Informations\n\n";
    std::cout << "number of unkowns: "  << this->Size() << std::endl;

    std::cout << "\n\t\tList of model parameters:\n\n";
    std::cout << "G    : " << this->G() << std::endl;
    std::cout << "Vth  : " << this->Vth() << std::endl;
    std::cout << "Vp   : " << this->Vp() << std::endl;
    std::cout << "Eta1 : " << this->Eta1() << std::endl;
    std::cout << "Eta2 : " << this->Eta2() << std::endl;
    std::cout << "Eta3 : " << this->Eta3() << std::endl;
    std::cout << "Eta  : " << this->Eta() << std::endl;
    std::cout << "Gamma: " << this->Gamma() << std::endl;

    std::cout << "\n\t\t End of IonicFitzHughNagumo Informations\n\n\n";

}


}

