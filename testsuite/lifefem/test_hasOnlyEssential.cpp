/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2004-12-16

  Copyright (C) 2004 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file test_hasOnlyEssential.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-12-16
*/
#include <bcHandler.hpp>

using namespace LifeV;

Real zero(const Real& t, const Real& x, const Real& y, const Real& z,
          const ID& id)
{
    return 0;
}

int main(int argc, char** argv)
{

    BCFunctionBase bcf(zero);

    std::vector<ID> compX;
    compX.push_back(1);

    std::vector<ID> compY;
    compY.push_back(2);

    std::vector<ID> compZ;
    compZ.push_back(3);

    std::vector<ID> compXY;
    compXY.push_back(1);
    compXY.push_back(2);

    std::vector<ID> compXZ;
    compXZ.push_back(1);
    compXZ.push_back(3);

    std::vector<ID> compYZ;
    compYZ.push_back(2);
    compYZ.push_back(3);

    std::vector<ID> compXYZ;
    compXYZ.push_back(1);
    compXYZ.push_back(2);
    compXYZ.push_back(3);

    std::auto_ptr<BCHandler> bcH(0);

    // 1     tests with one marker
    // 1.1   mode Full
    std::cout << "full" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->hasOnlyEssential();

    // 1.2   mode Component
    // 1.2.1 one component
    std::cout << "X" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compX );
    bcH->hasOnlyEssential();

    std::cout << "Y" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compY );
    bcH->hasOnlyEssential();

    std::cout << "Z" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    // 1.2.2 two components
    std::cout << "X,Y" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compX );
    bcH->addBC( "test", 1, Essential, Component, bcf, compY );
    bcH->hasOnlyEssential();

    std::cout << "X,Z" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compX );
    bcH->addBC( "test", 1, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "Y,Z" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compY );
    bcH->addBC( "test", 1, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "XY" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compXY );
    bcH->hasOnlyEssential();

    std::cout << "XZ" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compXZ );
    bcH->hasOnlyEssential();

    std::cout << "YZ" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compYZ );
    bcH->hasOnlyEssential();

    // 1.2.3 three components
    std::cout << "X,Y,Z" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Component, bcf, compX );
    bcH->addBC( "test", 1, Essential, Component, bcf, compY );
    bcH->addBC( "test", 1, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "X,YZ" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Component, bcf, compX );
    bcH->addBC( "test", 1, Essential, Component, bcf, compYZ );
    bcH->hasOnlyEssential();

    std::cout << "XZ,Y" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Component, bcf, compXZ );
    bcH->addBC( "test", 1, Essential, Component, bcf, compY );
    bcH->hasOnlyEssential();

    std::cout << "XY,Z" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Component, bcf, compXY );
    bcH->addBC( "test", 1, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "XYZ" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Component, bcf, compXYZ );
    bcH->hasOnlyEssential();

    // 1.3 modes Normal and Tangential
    std::cout << "N" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Normal, bcf );
    bcH->hasOnlyEssential();

    std::cout << "T" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Tangential, bcf );
    bcH->hasOnlyEssential();

    std::cout << "N,T" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Normal, bcf );
    bcH->addBC( "test", 1, Essential, Tangential, bcf );
    bcH->hasOnlyEssential();

    // 1.4 modes Component, Normal and Tangential
    std::cout << "X,N" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compX );
    bcH->addBC( "test", 1, Essential, Normal, bcf );
    bcH->hasOnlyEssential();

    std::cout << "Y,N" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compY );
    bcH->addBC( "test", 1, Essential, Normal, bcf );
    bcH->hasOnlyEssential();

    std::cout << "Z,N" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compZ );
    bcH->addBC( "test", 1, Essential, Normal, bcf );
    bcH->hasOnlyEssential();

    std::cout << "X,T" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compX );
    bcH->addBC( "test", 1, Essential, Tangential, bcf );
    bcH->hasOnlyEssential();

    std::cout << "Y,T" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compY );
    bcH->addBC( "test", 1, Essential, Tangential, bcf );
    bcH->hasOnlyEssential();

    std::cout << "Z,T" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Component, bcf, compZ );
    bcH->addBC( "test", 1, Essential, Tangential, bcf );
    bcH->hasOnlyEssential();

    // 2     tests with two markers
    // 2.1   mode Full/Full
    std::cout << "full / full" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Full, bcf, nDimensions );
    bcH->hasOnlyEssential();

    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Natural, Full, bcf, nDimensions );
    bcH->hasOnlyEssential();

    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Natural, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Full, bcf, nDimensions );
    bcH->hasOnlyEssential();

    // 2.2   mode Full/Component
    // 2.2.1 full + one component
    std::cout << "full / X" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compX );
    bcH->hasOnlyEssential();

    std::cout << "full / Y" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compY );
    bcH->hasOnlyEssential();

    std::cout << "full / Z" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    // 2.2.2 full + two components
    std::cout << "full / X,Y" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compX );
    bcH->addBC( "test", 2, Essential, Component, bcf, compY );
    bcH->hasOnlyEssential();

    std::cout << "full / X,Z" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compX );
    bcH->addBC( "test", 2, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "full / Y,Z" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compY );
    bcH->addBC( "test", 2, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "full / XY" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compXY );
    bcH->hasOnlyEssential();

    std::cout << "full / XZ" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compXZ );
    bcH->hasOnlyEssential();

    std::cout << "full / YZ" << std::endl;
    bcH.reset( new BCHandler );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compYZ );
    bcH->hasOnlyEssential();

    // 2.2.3 full + three components
    std::cout << "full / X,Y,Z" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compX );
    bcH->addBC( "test", 2, Essential, Component, bcf, compY );
    bcH->addBC( "test", 2, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "full / X,YZ" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compX );
    bcH->addBC( "test", 2, Essential, Component, bcf, compYZ );
    bcH->hasOnlyEssential();

    std::cout << "full / XZ,Y" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compXZ );
    bcH->addBC( "test", 2, Essential, Component, bcf, compY );
    bcH->hasOnlyEssential();

    std::cout << "full / XY,Z" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compXY );
    bcH->addBC( "test", 2, Essential, Component, bcf, compZ );
    bcH->hasOnlyEssential();

    std::cout << "full / XYZ" << std::endl;
    bcH.reset( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL ) );
    bcH->addBC( "test", 1, Essential, Full, bcf, nDimensions );
    bcH->addBC( "test", 2, Essential, Component, bcf, compXYZ );
    bcH->hasOnlyEssential();

    return 0;
}
