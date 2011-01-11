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
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>
#include <lifemc/lifearray/MatrixBlock.hpp>
#include <lifemc/lifearray/MatrixBlockView.hpp>
#include <lifemc/lifearray/MatrixBlockUtils.hpp>

namespace LifeV {

PreconditionerPCD::PreconditionerPCD(const  boost::shared_ptr<Epetra_Comm>& comm):
    PreconditionerComposition(comm),
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

void PreconditionerPCD::createParametersList( list_Type&         list,
                                              const GetPot&      dataFile,
                                              const std::string& section,
                                              const std::string& subSection )
{
    createPCDList( list, dataFile, section, subSection );
}

void PreconditionerPCD::createPCDList( list_Type&         list,
                                       const GetPot&      dataFile,
                                       const std::string& section,
                                       const std::string& subsection )
{
    bool displayList = dataFile((section + "/displayList").data(),     false);

    std::string precType = dataFile((section + "/prectype").data(),"PCD");
    list.set("prectype", precType);

    int velocityBlockSize = dataFile((section + "/" + subsection + "/blocks/velocity_block_size").data(), -1);
    int pressureBlockSize = dataFile((section + "/" + subsection + "/blocks/pressure_block_size").data(), -1);

    list.set("blocks: velocity block size", velocityBlockSize);
    list.set("blocks: pressure block size", pressureBlockSize);

    if (displayList) list.print(std::cout);
}

Real PreconditionerPCD::condest()
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

    // Make sure that an operator exists
    initializeOperator();

    std::vector<UInt> blockNumRows(2,0);
    blockNumRows[0] = M_velocityBlockSize;
    blockNumRows[1] = M_pressureBlockSize;
    std::vector<UInt> blockNumColumns(blockNumRows);

    bool inversed(true);
    //bool notInversed(false);
    //bool transposed(true);
    bool notTransposed(false);

    map_type map(oper->map());

    // Getting the block structure of A
    // / F Bt \
    // \ B C  /
    std::cout << std::endl << "     >Getting the structure of A" << std::endl;
    MatrixBlockView F,Bt,B,C;
    //oper.getMatrixBlockView(0,0,F );
    F.setup(0,0,blockNumRows[0],blockNumColumns[0],*oper);
    //F.showMe();

    //oper.getMatrixBlockView(0,1,Bt);
    Bt.setup(0,blockNumColumns[0],blockNumRows[0],blockNumColumns[1],*oper);
    //Bt.showMe();

    //oper.getMatrixBlockView(1,0,B );
    B.setup(blockNumRows[0],0,blockNumRows[1],blockNumColumns[0],*oper);
    //B.showMe();

    //oper.getMatrixBlockView(1,1,C );
    C.setup(blockNumRows[0],blockNumColumns[0],blockNumRows[1],blockNumColumns[1],*oper);
    //C.showMe();

    // PCD:
    // / F Bt \   / I  0 \ / I Bt \ / F 0 \
    // \ 0 -S / = \ 0 -S / \ 0 I  / \ 0 I /

    // PCD^-1:
    // / F  Bt \^-1   / F^-1 0 \ / I -Bt \ / I  0    \
    // \ 0 -S  /    = \ 0    I / \ 0  I  / \ 0 -S^-1 /

    // Getting the block structure of B
    MatrixBlockView B11,B12,B21,B22;

    // Building the block
    // / F^-1 0 \
    // \ 0    I /
    std::cout << " P1" << std::endl;
    boost::shared_ptr<matrix_type> P1(new matrix_type(map));
    P1->setBlockStructure(blockNumRows,blockNumColumns);
    P1->getMatrixBlockView(0,0,B11);
    P1->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::copyBlock(F,B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P1->globalAssemble();
    P1->spy("p1");
    boost::shared_ptr<parent_matrix_type> p1 = P1;
    M_precForBlock1.reset(new PreconditionerML());
    M_precForBlock1->setDataFromGetPot(M_dataFile,M_section);
    pushBack(p1,M_precForBlock1,inversed,notTransposed);

    // Building the block
    // / I -Bt \
    // \ 0  I  /
    std::cout << " P2" << std::endl;
    boost::shared_ptr<matrix_type> P2(new matrix_type(map));
    P2->setBlockStructure(blockNumRows,blockNumColumns);
    P2->getMatrixBlockView(0,0,B11);
    P2->getMatrixBlockView(0,1,B12);
    P2->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::copyBlock(Bt,B12);
    (*P2) *= -1;
    MatrixBlockUtils::createIdentityBlock(B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P2->globalAssemble();
    P2->spy("p2");
    boost::shared_ptr<parent_matrix_type> p2 = P2;
    pushBack(p2,inversed,notTransposed);

    // Building the block
    // / I  0    \   / I  0            \   / I  0     \ / I  0  \ / I 0     \
    // \ 0 -S^-1 / = \ 0 -Mp^-1FpAp^-1 / = \ 0  Mp^-1 / \ 0 -Fp / \ 0 Ap^-1 /
    std::cout << " P3a" << std::endl;
    boost::shared_ptr<matrix_type> P3a(new matrix_type( map ));
    *P3a *= 0.0;
    P3a->setBlockStructure(blockNumRows,blockNumColumns);
    P3a->getMatrixBlockView(0,0,B11);
    P3a->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::createIdentityBlock(B11);
    M_adrPressureAssembler.addMass(P3a,1.0,B22.firstRowIndex(),B22.firstColumnIndex());
    P3a->globalAssemble();
    P3a->spy("p3a");
    boost::shared_ptr<parent_matrix_type> p3a = P3a;
    M_precForBlock2.reset(new PreconditionerIfpack());
    //prec->createParametersList(prec->parametersList(), M_dataFile, M_section, "Ifpack");
    M_precForBlock2->setDataFromGetPot(M_dataFile,M_section);
    pushBack(p3a,M_precForBlock2,inversed,notTransposed);

    std::cout << " P3b" << std::endl;
    boost::shared_ptr<matrix_type> P3b(new matrix_type( map ));
    *P3b *= 0.0;
    P3b->setBlockStructure(blockNumRows,blockNumColumns);
    P3b->getMatrixBlockView(0,0,B11);
    P3b->getMatrixBlockView(1,1,B22);
    M_adrPressureAssembler.addDiffusion(P3b,-M_viscosity/M_density,B22.firstRowIndex(),B22.firstColumnIndex());
	M_adrPressureAssembler.addAdvection(P3b,*M_beta,B22.firstRowIndex(),B22.firstColumnIndex());
	M_adrPressureAssembler.addMass(P3b,1.0/M_timestep,B22.firstRowIndex(),B22.firstColumnIndex());
	(*P3b) *= -1;
    MatrixBlockUtils::createIdentityBlock(B11);
    P3b->globalAssemble();
    P3b->spy("p3b");
    boost::shared_ptr<parent_matrix_type> p3b = P3b;
    pushBack(p3b,inversed,notTransposed);

    std::cout << " P3c" << std::endl;
    boost::shared_ptr<matrix_type> P3c(new matrix_type( map ));
    *P3c *= 0.0;
    P3c->setBlockStructure(blockNumRows,blockNumColumns);
    P3c->getMatrixBlockView(0,0,B11);
    P3c->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::createIdentityBlock(B11);
    M_adrPressureAssembler.addDiffusion(P3c,-1.0,B22.firstRowIndex(),B22.firstColumnIndex());
    P3c->globalAssemble();
    P3c->spy("p3c");
    boost::shared_ptr<parent_matrix_type> p3c = P3c;
    M_precForBlock3.reset(new PreconditionerIfpack());
    M_precForBlock3->setDataFromGetPot(M_dataFile,M_section);
    pushBack(p3c,M_precForBlock3,inversed,notTransposed);

    // Only for debug purposes (Jacobi preconditioner)
    /*
    std::cout << "P1a" << std::endl;
    boost::shared_ptr<matrix_type> P1a(new matrix_type( map ));
    *P1a *= 0.0;
    P1a->setBlockStructure(blockNumRows,blockNumColumns);
    P1a->getMatrixBlockView(0,0,B11);
    P1a->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::createInvDiagBlock(F,B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P1a->globalAssemble();
    boost::shared_ptr<parent_matrix_type> p1a = P1a;
    pushBack(p1a,inversed,notTransposed);
    p1a->spy("diagA");
    oper->spy("A");
    */


    this->M_preconditionerCreated = true;

    std::cout << "[DEBUG] number of operators: " << numOperators() << std::endl;

    std::cout << "     >All the blocks are built" << std::endl;
    std::cout << "     >";
    return ( EXIT_SUCCESS );
}

void PreconditionerPCD::resetPreconditioner()
{
    M_precForBlock1.reset();
    M_precForBlock2.reset();
    M_precForBlock3.reset();
    PreconditionerComposition::resetPreconditioner();
}

int PreconditionerPCD::numBlocksRows() const
{
    return 2;
}

int PreconditionerPCD::numBlocksColumns() const
{
    return 2;
}

void PreconditionerPCD::setDataFromGetPot( const GetPot& dataFile,
                                           const std::string& section )
{
    M_dataFile   = dataFile;
    M_section    = section;
    createPCDList(M_list, dataFile, section, "PCD");

    M_velocityBlockSize = this->M_list.get("blocks: velocity block size", -1);
    M_pressureBlockSize = this->M_list.get("blocks: pressure block size", -1);
    M_precType          = this->M_list.get("prectype", "PCD");
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
