/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-11-29

  Copyright (C) 2010 EPFL

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
   \file PreconditionerPCD.cpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-11-29
 */

#include <vector>
#include "PreconditionerPCD.hpp"
#include <lifemc/lifearray/MatrixBlock.hpp>
#include <lifemc/lifearray/MatrixBlockView.hpp>
#include <lifemc/lifearray/MatrixBlockUtils.hpp>

namespace LifeV {

PreconditionerPCD::PreconditionerPCD():
    ComposedPreconditioner(),
    M_precType(""),
    M_velocityBlockSize(-1),
    M_pressureBlockSize(-1),
    M_timestep(1.0),
    M_viscosity(1.0),
    M_density(1.0)
{
    M_uFESpace.reset();
    M_pFESpace.reset();
    M_beta.reset();
}

PreconditionerPCD::~PreconditionerPCD()
{}

void PreconditionerPCD::createList( list_type&         list,
                                    const GetPot&      dataFile,
                                    const std::string& section,
                                    const std::string& subSection )
{
    createPCDList( list, dataFile, section, subSection );
}

void PreconditionerPCD::createPCDList( list_type&         list,
                                       const GetPot&      dataFile,
                                       const std::string& section,
                                       const std::string& subsection )
{
    bool displayList = dataFile((section + "/displayList").data(),     false);

    std::string precType = dataFile((section + "/prectype").data(),"PCD");
    list.set("prectype", precType);

    //std::string value1 = dataFile((section + "/" + subsection + "/subsection1/value1").data(), "string example");
    //int         value2 = dataFile((section + "/" + subsection + "/subsection2/value2").data(), 1);
    //double      value3 = dataFile((section + "/" + subsection + "/subsection3/value3").data(), 1.0);
    //list.set("subsection1: value1", value1);
    //list.set("subsection2: value2", value2);
    //list.set("subsection3: value3", value3);

    int velocityBlockSize = dataFile((section + "/" + subsection + "/blocks/velocity_block_size").data(), -1);
    int pressureBlockSize = dataFile((section + "/" + subsection + "/blocks/pressure_block_size").data(), -1);
    std::cout << "section: " << section << std::endl;

    list.set("blocks: velocity block size", velocityBlockSize);
    list.set("blocks: pressure block size", pressureBlockSize);

    if (displayList) list.print(std::cout);
}

Real PreconditionerPCD::Condest()
{
    return 0.0;
}

void PreconditionerPCD::updateBeta(const vector_type& beta)
{
    *M_beta = beta;
}

int PreconditionerPCD::buildPreconditioner(operator_type& oper)
{
    if((M_uFESpace.get()==NULL)||(M_pFESpace.get()==NULL))
    {
        std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
        exit(0);
    }

    std::vector<UInt> blockNumRows(2,0);
    blockNumRows[0] = M_velocityBlockSize;
    blockNumRows[1] = M_pressureBlockSize;
    std::vector<UInt> blockNumColumns(blockNumRows);

    bool inversed(true);
    bool notInversed(false);
    //bool transposed(true);
    bool notTransposed(false);

    map_type map(oper->getMap());

    // Getting the block structure of A
    // / F Bt \
    // \ B C  /
    std::cout << "Getting the structure of A" << std::endl;
    MatrixBlockView F,Bt,B,C;
    //oper.getMatrixBlockView(0,0,F );
    F.setup(0,0,blockNumRows[0],blockNumColumns[0],*oper);
    //oper.getMatrixBlockView(0,1,Bt);
    Bt.setup(0,blockNumColumns[0],blockNumRows[0],blockNumColumns[1],*oper);
    //oper.getMatrixBlockView(1,0,B );
    B.setup(blockNumRows[0],0,blockNumRows[1],blockNumColumns[0],*oper);
    //oper.getMatrixBlockView(1,1,C );
    C.setup(blockNumRows[0],blockNumColumns[0],blockNumRows[1],blockNumColumns[1],*oper);

    // PCD:
    // / F Bt \   / I  0 \ / I Bt \ / F 0 \
    // \ 0 -S / = \ 0 -S / \ 0 I  / \ 0 I /

    // PCD^-1:
    // / F  Bt \^-1   / F^-1 0 \ / I -Bt \ / I  0    \
    // \ 0 -S  /    = \ 0    I / \ 0  I  / \ 0 -S^-1 /

    // Getting the block structure of B
    MatrixBlockView B11,B12,B21,B22;

    // Building the block
    // / I  0 \   / I  0       \   / I  0  \ / I 0     \ / I 0  \
    // \ 0 -S / = \ 0 -ApFp^-1 / = \ 0 -Ap / \ 0 Fp^-1 / \ 0 Mp /
    std::cout << "P1a" << std::endl;
    boost::shared_ptr<matrix_type> P1a(new matrix_type( map ));
    *P1a *= 0.0;
    P1a->setBlockStructure(blockNumRows,blockNumColumns);
    P1a->getMatrixBlockView(0,0,B11);
    P1a->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::createIdentityBlock(B11);
    M_adrPressureAssembler.addDiffusion(P1a,-1.0,B22.firstRowIndex(),B22.firstColumnIndex());
    P1a->GlobalAssemble();
    P1a->spy("p1a");
    boost::shared_ptr<parent_matrix_type> p1a = P1a;
    //push_back(p1a,notInversed,notTransposed);

    std::cout << "P1b" << std::endl;
    boost::shared_ptr<matrix_type> P1b(new matrix_type( map ));
    *P1b *= 0.0;
    P1b->setBlockStructure(blockNumRows,blockNumColumns);
    P1b->getMatrixBlockView(0,0,B11);
    P1b->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::createIdentityBlock(B11);
    M_adrPressureAssembler.addDiffusion(P1b,-M_viscosity/M_density,B22.firstRowIndex(),B22.firstColumnIndex());
	M_adrPressureAssembler.addAdvection(P1b,*M_beta,B22.firstRowIndex(),B22.firstColumnIndex());
	M_adrPressureAssembler.addMass(P1b,1.0/M_timestep,B22.firstRowIndex(),B22.firstColumnIndex());
    P1b->GlobalAssemble();
    P1b->spy("p1b");
    boost::shared_ptr<parent_matrix_type> p1b = P1b;
    //boost::shared_ptr<src_matrix_type> p1b(P1b->getMatrixPtr());
    //push_back(p1b,inversed,notTransposed);

    std::cout << "P1c" << std::endl;
    boost::shared_ptr<matrix_type> P1c(new matrix_type( map ));
    *P1c *= 0.0;
    P1c->setBlockStructure(blockNumRows,blockNumColumns);
    P1c->getMatrixBlockView(0,0,B11);
    MatrixBlockUtils::createIdentityBlock(B11);
	M_adrPressureAssembler.addMass(P1c,1.0,B22.firstRowIndex(),B22.firstColumnIndex());
	P1c->GlobalAssemble();
    P1c->spy("p1c");
    boost::shared_ptr<parent_matrix_type> p1c = P1c;
    //push_back(p1c,notInversed,notTransposed);

    // Building the block
    // / F 0 \
    // \ 0 I /
    std::cout << "P2" << std::endl;
    boost::shared_ptr<matrix_type> P2(new matrix_type(map));
    P2->setBlockStructure(blockNumRows,blockNumColumns);
    P2->getMatrixBlockView(0,0,B11);
    P2->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::copyBlock(F,B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P2->GlobalAssemble();
    P2->spy("p2");
    boost::shared_ptr<parent_matrix_type> p2 = P2;
    //push_back(p2,notInversed,notTransposed);

    // Building the block
    // / I Bt \
    // \ 0 I  /
    std::cout << "P3" << std::endl;
    boost::shared_ptr<matrix_type> P3(new matrix_type(map));
    P3->setBlockStructure(blockNumRows,blockNumColumns);
    P3->getMatrixBlockView(0,0,B11);
    P3->getMatrixBlockView(0,1,B12);
    P3->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::createIdentityBlock(B11);
    MatrixBlockUtils::copyBlock(Bt,B12);
    MatrixBlockUtils::createIdentityBlock(B22);
    P3->GlobalAssemble();
    P3->spy("p3");
    boost::shared_ptr<parent_matrix_type> p3 = P3;
    //push_back(p3,notInversed,notTransposed);

    std::cout << "All the blocks are built" << std::endl;

    return ( EXIT_SUCCESS );
}

int PreconditionerPCD::numBlocksRows() const
{
    return 2;
}

int PreconditionerPCD::numBlocksCols() const
{
    return 2;
}

void PreconditionerPCD::setDataFromGetPot( const GetPot& dataFile,
                                           const std::string& section )
{
    createPCDList(M_List, dataFile, section, "PCD");

    M_velocityBlockSize = this->M_List.get("blocks: velocity block size", -1);
    M_pressureBlockSize = this->M_List.get("blocks: pressure block size", -1);
    M_precType          = this->M_List.get("prectype", "PCD");
}

void PreconditionerPCD::setFESpace(FESpace_ptr uFESpace,FESpace_ptr pFESpace){
    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_adrPressureAssembler.setup(pFESpace,uFESpace); // p,beta=u
    M_adrVelocityAssembler.setup(uFESpace,uFESpace); // u,beta=u
    M_beta.reset(new vector_type(M_uFESpace->map() + M_pFESpace->map()));
    *M_beta *= 0;
}

void PreconditionerPCD::setTimestep(const Real& timestep)
{
    M_timestep = timestep;
}

void PreconditionerPCD::setViscosity(const Real& viscosity)
{
    M_viscosity = viscosity;
}

void PreconditionerPCD::setDensity(const Real& density)
{
    M_density = density;
}

} // namespace LifeV
