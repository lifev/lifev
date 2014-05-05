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
  @brief Ionic model of Bueno-Orovio et al.
  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <boost/typeof/typeof.hpp>

namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
IonicMinimalModel::IonicMinimalModel()  :
    super       ( 4, 2     ),
    M_uo        ( 0.    ),
    M_uu        ( 1.58  ),
    M_tetav     ( 0.3   ),
    M_tetaw     ( 0.015 ),
    M_tetavm    ( 0.015 ),
    M_tetao     ( 0.006 ),
    M_tauv1     ( 60.0  ),
    M_tauv2     ( 1150.0),
    M_tauvp     ( 1.4506),
    M_tauw1     ( 70.0  ),
    M_tauw2     ( 20.0  ),
    M_kw        ( 65.0  ),
    M_uw        ( 0.03  ),
    M_tauwp     ( 280.0 ),
    M_taufi     ( 0.11  ),
    M_tauo1     ( 6.0   ),
    M_tauo2     ( 6.0   ),
    M_tauso1    ( 43.0  ),
    M_tauso2    ( 0.2   ),
    M_kso       ( 2.0   ),
    M_uso       ( 0.65  ),
    M_taus1     ( 2.7342),
    M_taus2     ( 3.0   ),
    M_ks        ( 2.0994),
    M_us        ( 0.9087),
    M_tausi     ( 2.8723),
    M_tauwinf   ( 0.07  ),
    M_winfstar  ( 0.94  )
{
    M_restingConditions.at (0) = 0.0;
    M_restingConditions.at (1) = 1.0;
    M_restingConditions.at (2) = 1.0;
    M_restingConditions.at (3) = 0.021553043080281;

}

IonicMinimalModel::IonicMinimalModel ( Teuchos::ParameterList& parameterList     )   :
    super       ( 4 )
{
    M_uo        =  parameterList.get ("uo",      1000.0     );
    M_uu        =  parameterList.get ("uu",      1.61    );
    M_tetav     =  parameterList.get ("tetav",   0.30    );
    M_tetaw     =  parameterList.get ("tetaw",   0.130   );
    M_tetavm    =  parameterList.get ("tetavm",  0.10    );
    M_tetao     =  parameterList.get ("tetao",   0.005   );
    M_tauv1     =  parameterList.get ("tauv1",   80.0    );
    M_tauv2     =  parameterList.get ("tauv2",   1.4506  );
    M_tauvp     =  parameterList.get ("tauvp",   1.4506  );
    M_tauw1     =  parameterList.get ("tauw1",   70.0    );
    M_tauw2     =  parameterList.get ("tauw2",   8.0     );
    M_kw        =  parameterList.get ("kw",      200.0   );
    M_uw        =  parameterList.get ("uw",      0.016   );
    M_tauwp     =  parameterList.get ("tauwp",   280.0   );
    M_taufi     =  parameterList.get ("taufi",   0.078   );
    M_tauo1     =  parameterList.get ("tauo1",   410.0   );
    M_tauo2     =  parameterList.get ("tauo2",   7.0     );
    M_tauso1    =  parameterList.get ("tauso1",  91.0    );
    M_tauso2    =  parameterList.get ("tauso2",  0.8     );
    M_kso       =  parameterList.get ("kso",     2.1     );
    M_uso       =  parameterList.get ("uso",     0.6     );
    M_taus1     =  parameterList.get ("taus1",   2.7342  );
    M_taus2     =  parameterList.get ("taus2",   4.0     );
    M_ks        =  parameterList.get ("ks",      2.0994  );
    M_us        =  parameterList.get ("us",      0.9087  );
    M_tausi     =  parameterList.get ("tausi",   3.3849  );
    M_tauwinf   =  parameterList.get ("tauwinf", 0.01    );
    M_winfstar  =  parameterList.get ("winfstar", 0.5     );

    M_restingConditions.at (0) = 0.0;
    M_restingConditions.at (1) = 1.0;
    M_restingConditions.at (2) = 1.0;
    M_restingConditions.at (3) = 0.021553043080281;
}

IonicMinimalModel::IonicMinimalModel ( const IonicMinimalModel& model )
{

    M_uo        =  model.M_uo;
    M_uu        =  model.M_uu;
    M_tetav     =  model.M_tetav;
    M_tetaw     =  model.M_tetaw;
    M_tetavm    =  model.M_tetavm;
    M_tetao     =  model.M_tetao;
    M_tauv1     =  model.M_tauv1;
    M_tauv2     =  model.M_tauv2;
    M_tauvp     =  model.M_tauvp;
    M_tauw1     =  model.M_tauw1;
    M_tauw2     =  model.M_tauw2;
    M_kw        =  model.M_kw;
    M_uw        =  model.M_uw;
    M_tauwp     =  model.M_tauwp;
    M_taufi     =  model.M_taufi;
    M_tauo1     =  model.M_tauo1;
    M_tauo2     =  model.M_tauo2;
    M_tauso1    =  model.M_tauso1;
    M_tauso2    =  model.M_tauso2;
    M_kso       =  model.M_kso;
    M_uso       =  model.M_uso;
    M_taus1     =  model.M_taus1;
    M_taus2     =  model.M_taus2;
    M_ks        =  model.M_ks;
    M_us        =  model.M_us;
    M_tausi     =  model.M_tausi;
    M_tauwinf   =  model.M_tauwinf;
    M_winfstar  =  model.M_winfstar;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;

}

// ===================================================
//! Operator
// ===================================================
IonicMinimalModel& IonicMinimalModel::operator= ( const IonicMinimalModel& model )
{
    M_uo        =  model.M_uo;
    M_uu        =  model.M_uu;
    M_tetav     =  model.M_tetav;
    M_tetaw     =  model.M_tetaw;
    M_tetavm    =  model.M_tetavm;
    M_tetao     =  model.M_tetao;
    M_tauv1     =  model.M_tauv1;
    M_tauv2     =  model.M_tauv2;
    M_tauvp     =  model.M_tauvp;
    M_tauw1     =  model.M_tauw1;
    M_tauw2     =  model.M_tauw2;
    M_kw        =  model.M_kw;
    M_uw        =  model.M_uw;
    M_tauwp     =  model.M_tauwp;
    M_taufi     =  model.M_taufi;
    M_tauo1     =  model.M_tauo1;
    M_tauo2     =  model.M_tauo2;
    M_tauso1    =  model.M_tauso1;
    M_tauso2    =  model.M_tauso2;
    M_kso       =  model.M_kso;
    M_uso       =  model.M_uso;
    M_taus1     =  model.M_taus1;
    M_taus2     =  model.M_taus2;
    M_ks        =  model.M_ks;
    M_us        =  model.M_us;
    M_tausi     =  model.M_tausi;
    M_tauwinf   =  model.M_tauwinf;
    M_winfstar  =  model.M_winfstar;

    M_numberOfEquations = model.M_numberOfEquations;
    M_restingConditions = model.M_restingConditions;


    return      *this;
}


// ===================================================
//! Methods
// ===================================================
void IonicMinimalModel::computeGatingRhs ( const   std::vector<Real>&  v,
                                           std::vector<Real>& rhs )
{

    Real U = v[0];
    Real V = v[1];
    Real W = v[2];
    Real S = v[3];

    Real tauvm = ( 1.0 - Heaviside ( U - M_tetavm ) ) * M_tauv1 + Heaviside ( U - M_tetavm ) * M_tauv2;
    Real tauwm = M_tauw1 + ( M_tauw2  - M_tauw1  ) * ( 1.0 + std::tanh ( M_kw  * ( U - M_uw  ) ) ) / 2.0;
    Real taus  = ( 1.0 - Heaviside ( U - M_tetaw ) ) * M_taus1 + Heaviside ( U - M_tetaw ) * M_taus2;

    Real vinf  = Heaviside ( M_tetavm - U );
    Real winf  = ( 1.0 - Heaviside ( U - M_tetao ) ) * ( 1.0 - U / M_tauwinf ) + Heaviside ( U - M_tetao ) * M_winfstar;

    rhs[0] = ( 1.0 - Heaviside ( U - M_tetav ) ) * ( vinf - V ) / tauvm - Heaviside ( U - M_tetav ) * V / M_tauvp;
    rhs[1] = ( 1.0 - Heaviside ( U - M_tetaw ) ) * ( winf - W ) / tauwm - Heaviside ( U - M_tetaw ) * W / M_tauwp;
    rhs[2] = ( ( 1.0 + std::tanh ( M_ks * ( U - M_us ) ) ) / 2.0 - S ) / taus;

}

void IonicMinimalModel::computeRhs ( const   std::vector<Real>&  v,
                                     std::vector<Real>& rhs )
{

    Real U = v[0];
    Real V = v[1];
    Real W = v[2];
    Real S = v[3];

    Real tauvm = ( 1.0 - Heaviside ( U - M_tetavm ) ) * M_tauv1 + Heaviside ( U - M_tetavm ) * M_tauv2;
    Real tauwm = M_tauw1 + ( M_tauw2  - M_tauw1  ) * ( 1.0 + std::tanh ( M_kw  * ( U - M_uw  ) ) ) / 2.0;
    Real tauso = M_tauso1 + ( M_tauso2 - M_tauso1 ) * ( 1.0 + std::tanh ( M_kso * ( U - M_uso ) ) ) / 2.0;
    Real taus  = ( 1.0 - Heaviside ( U - M_tetaw ) ) * M_taus1 + Heaviside ( U - M_tetaw ) * M_taus2;
    Real tauo  = ( 1.0 - Heaviside ( U - M_tetao ) ) * M_tauo1 + Heaviside ( U - M_tetao ) * M_tauo2;

    Real vinf  = Heaviside ( M_tetavm - U );
    Real winf  = ( 1.0 - Heaviside ( U - M_tetao ) ) * ( 1.0 - U / M_tauwinf ) + Heaviside ( U - M_tetao ) * M_winfstar;

    Real Jfi   = - V * Heaviside ( U - M_tetav ) * ( U - M_tetav ) * ( M_uu - U ) / M_taufi;
    Real Jso   = ( U - M_uo ) * ( 1.0 - Heaviside ( U - M_tetaw )  ) / tauo + Heaviside ( U - M_tetaw ) / tauso;
    Real Jsi   = - Heaviside ( U - M_tetaw ) * W * S / M_tausi;

    rhs[0] = - ( Jfi + Jso + Jsi );
    rhs[1] = ( 1.0 - Heaviside ( U - M_tetav ) ) * ( vinf - V ) / tauvm - Heaviside ( U - M_tetav ) * V / M_tauvp;
    rhs[2] = ( 1.0 - Heaviside ( U - M_tetaw ) ) * ( winf - W ) / tauwm - Heaviside ( U - M_tetaw ) * W / M_tauwp;
    rhs[3] = ( ( 1.0 + std::tanh ( M_ks * ( U - M_us ) ) ) / 2.0 - S ) / taus;

}


void IonicMinimalModel::computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt )
{
    std::cout << "\n\nRush Larsen method, for minimal model not implemented!!!\n";
    std::cout << "\n\nI will use Forward Euler!!!\n";
    std::vector<Real> rhs(3,0.0);
    computeGatingRhs ( v, rhs );
    v[1] += dt * rhs[0];
    v[2] += dt * rhs[1];
    v[3] += dt * rhs[2];
}


Real IonicMinimalModel::computeLocalPotentialRhs ( const std::vector<Real>& v )
{
    Real dPotential (0.0);

    Real U = v[0];
    Real V = v[1];
    Real W = v[2];
    Real S = v[3];

    Real tauso = M_tauso1 + ( M_tauso2 - M_tauso1 ) * ( 1.0 + std::tanh ( M_kso * ( U - M_uso ) ) ) / 2.0;
    Real tauo  = ( 1.0 - Heaviside ( U - M_tetao ) ) * M_tauo1 + Heaviside ( U - M_tetao ) * M_tauo2;

    Real Jfi   = - V * Heaviside ( U - M_tetav ) * ( U - M_tetav ) * ( M_uu - U ) / M_taufi;
    Real Jso   = ( U - M_uo ) * ( 1.0 - Heaviside ( U - M_tetaw )  ) / tauo + Heaviside ( U - M_tetaw ) / tauso;
    Real Jsi   = - Heaviside ( U - M_tetaw ) * W * S / M_tausi;

    dPotential = - ( Jfi + Jso  + Jsi );
    //    std::cout.precision(15);
    //    std::cout << "\nValue: " << Jso;
    return dPotential;
}


void IonicMinimalModel::showMe()
{

    std::cout << "\n\tHi, I'm the minimal model\n\t See you soon\n\n";

    std::cout << "\nuo " <<       M_uo      ;
    std::cout << "\nuu " <<       M_uu     ;
    std::cout << "\ntetav " <<    M_tetav     ;
    std::cout << "\ntetaw " <<    M_tetaw    ;
    std::cout << "\ntetavm " <<   M_tetavm     ;
    std::cout << "\ntetao " <<    M_tetao    ;
    std::cout << "\ntauv1 " <<    M_tauv1      ;
    std::cout << "\ntauv2 " <<    M_tauv2   ;
    std::cout << "\ntauvp " <<    M_tauvp   ;
    std::cout << "\ntauw1 " <<    M_tauw1     ;
    std::cout << "\ntauw2 " <<    M_tauw2      ;
    std::cout << "\nkw " <<       M_kw    ;
    std::cout << "\nuw " <<       M_uw    ;
    std::cout << "\ntauwp " <<    M_tauwp    ;
    std::cout << "\ntaufi " <<    M_taufi    ;
    std::cout << "\ntauo1 " <<    M_tauo1   ;
    std::cout << "\ntauo2 " <<    M_tauo2      ;
    std::cout << "\ntauso1 " <<   M_tauso1     ;
    std::cout << "\ntauso2 " <<   M_tauso2      ;
    std::cout << "\nkso " <<     M_kso       ;
    std::cout << "\nuso " <<      M_uso      ;
    std::cout << "\ntaus1 " <<    M_taus1    ;
    std::cout << "\ntaus2 " <<   M_taus2      ;
    std::cout << "\nks " <<       M_ks   ;
    std::cout << "\nus " <<       M_us   ;
    std::cout << "\ntausi " <<    M_tausi   ;
    std::cout << "\ntauwinf " <<  M_tauwinf     ;
    std::cout << "\nwinfstar " <<   M_winfstar    << "\n"  ;
}

void IonicMinimalModel::computePotentialRhsSVI ( const std::vector<vectorPtr_Type>& v,
                                                 std::vector<vectorPtr_Type>&        rhs,
                                                 FESpace<mesh_Type, MapEpetra>&  uFESpace)
{
    typedef ETFESpace<mesh_Type, MapEpetra, 3, 1> ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type> ETFESpacePtr_Type;

    * (rhs[0]) *= 0.0;

    if (uFESpace.mapPtr() -> commPtr() -> MyPID() == 0)
    {
        std::cout << "\nMinimal Model: Assembling SVI using ETA!\n";
    }

    ETFESpacePtr_Type spaceScalar (
        new ETFESpace_Type (uFESpace.mesh(), &feTetraP1, uFESpace.mapPtr() -> commPtr()  ) );

    {
        using namespace ExpressionAssembly;

        boost::shared_ptr<MMTanhFunctor> tanh (new MMTanhFunctor);
        boost::shared_ptr<MMHFunctor> H (new MMHFunctor);
        boost::shared_ptr<MMSV> sv (new MMSV);

        BOOST_AUTO_TPL (U, value (spaceScalar, * (v[0] ) ) );
        BOOST_AUTO_TPL (V, value (spaceScalar, * (v[1] ) ) );
        BOOST_AUTO_TPL (W, value (spaceScalar, * (v[2] ) ) );
        BOOST_AUTO_TPL (S, value (spaceScalar, * (v[3] ) ) );
        BOOST_AUTO_TPL (Iapp, value (spaceScalar, *M_appliedCurrentPtr ) );


        BOOST_AUTO_TPL (tauso, M_tauso1 + ( M_tauso2 - M_tauso1 ) * ( 1.0 + eval (tanh, M_kso * ( U - M_uso ) ) ) / 2.0);
        BOOST_AUTO_TPL (tauo, ( 1.0 - eval (H, U - M_tetao ) ) * M_tauo1 + eval (H, U - M_tetao ) * M_tauo2);
        BOOST_AUTO_TPL (Jfi, value (-1.0) * V * eval (H, U - M_tetav ) * ( U - M_tetav ) * ( M_uu - U ) / M_taufi);
        BOOST_AUTO_TPL (Jso,  ( U - M_uo ) * ( 1.0 - eval (H, U - M_tetaw )  ) / tauo + eval (H, U - M_tetaw ) / tauso);
        BOOST_AUTO_TPL (Jsi, value (-1.0) * eval (H, U - M_tetaw ) * W * S / M_tausi);

        integrate ( elements ( uFESpace.mesh() ),
        		    uFESpace.qr(),
                    spaceScalar,
                    ( Iapp - ( Jfi +  Jso + Jsi )  ) * phi_i )  >> rhs.at (0);
    }

    rhs.at (0) -> globalAssemble();
}

}

