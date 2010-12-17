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

#include <life/lifecore/life.hpp>

#include <BlockInterface.hpp>

namespace LifeV
{


// ===================================================
//! Public Methods
// ===================================================

void BlockInterface::couplingMatrix(matrixPtr_Type & bigMatrix,
                                    Int flag,
                                    const std::vector<fespacePtr_Type>& problem,
                                    const std::vector<UInt>& offset,
                                    const std::map<ID, ID>& locDofMap,
                                    const vectorPtr_Type& numerationInterface,
                                    const Real& timeStep,
                                    const Real& value) // not working with non-matching grids
{// coupling from 1 to 31, working as chmod
    if ( flag-31 > 0 )//recursive call
    {
        Int subFlag(flag);
        subFlag -= 31;
        std::vector<fespacePtr_Type> newProblem(problem.begin()+3, problem.end());
        std::vector<UInt> newOffset(offset.begin()+3, offset.end());

        couplingMatrix( bigMatrix, subFlag, newProblem, newOffset, locDofMap, numerationInterface, timeStep, value);
    }
    std::map<ID, ID>::const_iterator ITrow;
    UInt interface(numerationInterface->map().map(Unique)->NumGlobalElements());
    UInt totalSize(offset[0]+problem[0]->map().map(Unique)->NumGlobalElements());

    if (flag-16>=0)//coupling the mesh in FSI
    {
        //UInt interface( numerationInterface->getMap().getMap( Unique )->NumGlobalElements() );
        UInt solidFluidInterface( offset[2] );
        std::map<ID, ID>::const_iterator ITrow;
        for ( UInt dim = 0; dim < nDimensions; ++dim )
        {
            for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow )
            {
                if ( numerationInterface->map().map(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    bigMatrix->addToCoefficient(solidFluidInterface + ITrow->first + dim*problem[2]->dof().numTotalDof() - 1, offset[0] + ITrow->second-1 + dim* problem[0]->dof().numTotalDof(), (-value)*timeStep/**1.e-2*//*scaling of the solid matrix*/ );
                }
            }
        }
        flag -= 16;
    }

    Int newFlag(flag);

    for (UInt dim = 0; dim < nDimensions; ++dim)//coupling F-S in FSI
    {
        for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow, newFlag = flag )
        {
            if (numerationInterface->map().map(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                if (newFlag-8>=0)//right low
                {
                    bigMatrix->addToCoefficient( offset[0] + ITrow->second-1 + dim* problem[0]->dof().numTotalDof(),(int)(*numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*interface + totalSize, value );//right low
                    newFlag -= 8;
                }
                if (newFlag-4>=0)// right up
                {
                    bigMatrix->addToCoefficient( offset[1] + ITrow->first-1 + dim* problem[1]->dof().numTotalDof(), (int)(*numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*interface + totalSize, -value );//right up
                    newFlag -= 4;
                }
                if (newFlag-2>=0)//low left
                {
                    bigMatrix->addToCoefficient( (int)(*numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*interface + totalSize, (ITrow->first)-1 + dim* problem[1]->dof().numTotalDof(), value);//low left
                    newFlag -= 2;
                }
                if (newFlag-1>=0)//low right
                    bigMatrix->addToCoefficient( (int)(*numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*interface + totalSize, (offset[0] + ITrow->second)-1 + dim* problem[0]->dof().numTotalDof(), -value);//low right

                bigMatrix->addToCoefficient( (int)(*numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*interface + totalSize , (int)(*numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*interface + totalSize, 0.0);
            }
        }
    }
}

void BlockInterface::applyBoundaryConditions(const Real& time)
{
    ASSERT( M_bch.size() == M_blocks.size(), "The number of BChandlers must correspond to the number of blocks" )
    for (ID i=0; i<M_bch.size(); ++i)
    {
        applyBoundaryConditions(time, i);
    }
}


void BlockInterface::applyBoundaryConditions(const Real& time, const UInt i)
{
    M_blocks[i]->openCrsMatrix();
    if ( !M_bch[i]->bcUpdateDone() )
    {
        M_bch[i]->bcUpdate( *M_FESpace[i]->mesh(), M_FESpace[i]->feBd(), M_FESpace[i]->dof() );
        M_bch[i]->setOffset(M_offset[i]);
    }
    //M_blocks[i]->GlobalAssemble();
    bcManageMatrix( *M_blocks[i] , *M_FESpace[i]->mesh(), M_FESpace[i]->dof(), *M_bch[i], M_FESpace[i]->feBd(), 1., time);
}

void BlockInterface::setConditions( std::vector<bchandlerPtr_Type>& vec )
{
    M_bch = vec;
}

void
BlockInterface::setSpaces(std::vector<fespacePtr_Type>& vec )
{
    M_FESpace = vec;
}

void
BlockInterface::setOffsets(UInt blocks, ...)
{
    using namespace std;
    va_list arguments;
    va_start(arguments, blocks);
    M_offset.resize(0);

    for (ID i=0; i<M_blocks.size(); ++i)
        M_offset.push_back(va_arg(arguments, UInt));
    va_end(arguments);
}


void
BlockInterface::robinCoupling( matrixPtr_Type& matrix,
                               Real&  alphaf,
                               Real&  alphas,
                               UInt  coupling,
                               const BlockInterface::fespacePtr_Type& FESpace1,
                               const UInt& /*offset1*/,
                               const BlockInterface::fespacePtr_Type& FESpace2,
                               const UInt& offset2,
                               const std::map<ID, ID>& locDofMap,
                               const BlockInterface::vectorPtr_Type& numerationInterface ) // not working with non-matching grids
{//coupling: flag from 1 to 4 working as chmod
    UInt interface(numerationInterface->map().map(Unique)->NumGlobalElements());
    std::map<ID, ID>::const_iterator ITrow;
    UInt totalSize(offset2+FESpace2->map().map(Unique)->NumGlobalElements());


    if (((Int)coupling)-4 >= 0)
    {
        for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if (numerationInterface->map().map(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                for (UInt dim = 0; dim < nDimensions; ++dim)
                {
                    matrix->addToCoefficient( (int)(*numerationInterface)[ITrow->second ] - 1 + dim*interface + totalSize, (offset2 + ITrow->second)-1 + dim* FESpace2->dof().numTotalDof(), alphas);//low right
                }
            }
        }
        coupling -= 4;
    }
    if (((Int)coupling-2) >= 0)
    {
        for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if (numerationInterface->map().map(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                for (UInt dim = 0; dim < nDimensions; ++dim)
                {
                    matrix->addToCoefficient( ITrow->first-1 + dim* FESpace1->dof().numTotalDof(), (int)(*numerationInterface)[ITrow->second ] - 1 + dim*interface + totalSize, alphaf );//right up
                }
            }
        }
        coupling -= 2;
    }
    if (((Int)coupling)-1 >= 0)
    {
        for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if (numerationInterface->map().map(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                for (UInt dim = 0; dim < nDimensions; ++dim)
                {
                    matrix->addToCoefficient( (int)(*numerationInterface)[ITrow->second ] - 1 + dim*interface + totalSize, (ITrow->first)-1 + dim* FESpace1->dof().numTotalDof(), alphas);//low left
                }
            }
        }
    }
}

void BlockInterface::addToBlock( const matrixPtr_Type& Mat, UInt position)
{
    *Mat += *M_blocks[position];
    M_blocks[position] = Mat;
}

void BlockInterface::push_back_oper( BlockInterface& Oper)
{
//     M_blocks.reserve(M_blocks.size()+Oper.getBlockVector().size());
//     M_bch.reserve(M_bch.size()+Oper.getBChVector().size());
//     M_FESpace.reserve(M_FESpace.size()+Oper.getFESpaceVector().size());
//     M_offset.reserve(M_offset.size()+Oper.getOffestVector().size());

    M_blocks.insert(M_blocks.end(), Oper.getBlockVector().begin(), Oper.getBlockVector().end());
    M_bch.insert(M_bch.end(), Oper.getBChVector().begin(), Oper.getBChVector().end());
    M_FESpace.insert(M_FESpace.end(), Oper.getFESpaceVector().begin(), Oper.getFESpaceVector().end());
    M_offset.insert(M_offset.end(), Oper.getOffsetVector().begin(), Oper.getOffsetVector().end());
}

} // Namespace LifeV
