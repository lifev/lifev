//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 09 Jul 2010

    A more detailed description of the file (if necessary)
 */

#include <lifemc/lifesolver/RobinInterface.hpp>

namespace LifeV
{

void RobinInterface::setRobinData(const GetPot& data, const std::string& section)
{
    M_alphas=data((section + "/alphas").data(), 0.);
    M_alphaf=data((section + "/alphaf").data(), 0.);
}


void RobinInterface::applyRobinCoupling( std::vector<BlockInterface::matrix_ptrtype> blockVector)
{
    M_robinPart.reset(new BlockInterface::matrix_type(M_robinCoupling->getMap(), 0));
    for ( UInt ITBlock = 0; ITBlock < blockVector.size(); ++ITBlock )
        applyRobinCoupling( blockVector[ITBlock] );
}

void RobinInterface::applyRobinCoupling( BlockInterface::matrix_ptrtype block)
{
    BlockInterface::matrix_ptrtype tmpMatrix(new BlockInterface::matrix_type(M_robinCoupling->getMap(), 0));
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
