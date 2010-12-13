/* -*- mode: c++ -*- */
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

#include <lifeconfig.h>

#include <lifemc/lifesolver/RobinInterface.hpp>

namespace LifeV
{

// ===================================================
//! Public Methods
// ===================================================

void RobinInterface::setRobinData(const GetPot& data, const std::string& section)
{
    M_alphas=data((section + "/alphas").data(), 0.);
    M_alphaf=data((section + "/alphaf").data(), 0.);
}


void RobinInterface::applyRobinCoupling( std::vector<BlockInterface::matrixPtr_Type> blockVector)
{
    M_robinPart.reset(new BlockInterface::matrix_Type(M_robinCoupling->getMap(), 0));
    for( UInt ITBlock = 0; ITBlock < blockVector.size(); ++ITBlock )
        applyRobinCoupling( blockVector[ITBlock] );
}


// ===================================================
//! Protected Methods
// ===================================================

void RobinInterface::applyRobinCoupling( BlockInterface::matrixPtr_Type block)
{
    BlockInterface::matrixPtr_Type tmpMatrix(new BlockInterface::matrix_Type(M_robinCoupling->getMap(), 0));
    int err = EpetraExt::MatrixMatrix::
              Multiply( *M_robinCoupling->getMatrixPtr(),
                        false,
                        *block->getMatrixPtr(),
                        false,
                        *tmpMatrix->getMatrixPtr()
                      );
    tmpMatrix->GlobalAssemble();
    ASSERT(!err, "Error in multiplication");
    *M_robinPart += *tmpMatrix;
}

} // Namespace LifeV
