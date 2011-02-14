/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2011-01-24

  Copyright (C) 2011 EPFL

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
   \file PreconditionerSIMPLE.cpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2011-01-24
 */

#include <vector>
#include "PreconditionerSIMPLE.hpp"
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <lifemc/lifearray/MatrixBlock.hpp>
#include <lifemc/lifearray/MatrixBlockView.hpp>
#include <lifemc/lifearray/MatrixBlockUtils.hpp>

namespace LifeV {

PreconditionerSIMPLE::PreconditionerSIMPLE(const  boost::shared_ptr<Epetra_Comm>& comm):
    PreconditionerComposition(comm),
    M_velocityBlockSize(-1),
    M_pressureBlockSize(-1),
    M_timestep(1.0),
    M_viscosity(1.0),
    M_density(1.0),
    M_dampingFactor(0.5),
    M_SIMPLEType("SIMPLE")
{
    M_uFESpace.reset();
    M_pFESpace.reset();
    M_beta.reset();
}

PreconditionerSIMPLE::~PreconditionerSIMPLE()
{}

void PreconditionerSIMPLE::createParametersList( list_Type&         list,
                                                 const GetPot&      dataFile,
                                                 const std::string& section,
                                                 const std::string& subSection )
{
    createSIMPLEList( list, dataFile, section, subSection );
}

void PreconditionerSIMPLE::createSIMPLEList( list_Type&         list,
                                             const GetPot&      dataFile,
                                             const std::string& section,
                                             const std::string& subsection )
{
    bool displayList = dataFile((section + "/displayList").data(),     false);

    std::string precType = dataFile((section + "/prectype").data(),"SIMPLE");
    list.set("prectype", precType);

    int velocityBlockSize = dataFile((section + "/" + subsection + "/blocks/velocity_block_size").data(), -1);
    int pressureBlockSize = dataFile((section + "/" + subsection + "/blocks/pressure_block_size").data(), -1);
    string SIMPLEType = dataFile((section + "/" + subsection + "/SIMPLE_type").data(), "SIMPLE");

    list.set("blocks: velocity block size", velocityBlockSize);
    list.set("blocks: pressure block size", pressureBlockSize);
    list.set("SIMPLE Type", SIMPLEType);

    if (displayList) list.print(std::cout);
}

Real PreconditionerSIMPLE::condest()
{
    return 0.0;
}

void PreconditionerSIMPLE::updateBeta(const vector_type& beta)
{
    *M_beta = beta;
}

int PreconditionerSIMPLE::buildPreconditioner(operator_type& oper)
{
    if((M_uFESpace.get()==NULL)||(M_pFESpace.get()==NULL))
    {
        std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
        exit(0);
    }

    bool verbose(false);
    if(M_comm->MyPID() == 0) verbose = true;

    // Make sure that an operator exists
    initializeOperator();

    std::vector<UInt> blockNumRows(2,0);
    blockNumRows[0] = M_velocityBlockSize;
    blockNumRows[1] = M_pressureBlockSize;
    std::vector<UInt> blockNumColumns(blockNumRows);

    bool inversed(true);
    bool notInversed(false);
    //bool transposed(true);
    bool notTransposed(false);

    map_type map(oper->map());
    //oper->spy("A");

    LifeChrono timer;

    // Getting the block structure of A
    // / F Bt \
    // \ B C  /
    if(verbose) std::cout << std::endl << "      >Getting the structure of A... ";
    timer.start();
    MatrixBlockView F,Bt,B,C;
    //oper.getMatrixBlockView(0,0,F );
    F.setup(0,0,blockNumRows[0],blockNumColumns[0],*oper);

    //oper.getMatrixBlockView(0,1,Bt);
    Bt.setup(0,blockNumColumns[0],blockNumRows[0],blockNumColumns[1],*oper);

    //oper.getMatrixBlockView(1,0,B );
    B.setup(blockNumRows[0],0,blockNumRows[1],blockNumColumns[0],*oper);

    //oper.getMatrixBlockView(1,1,C );
    C.setup(blockNumRows[0],blockNumColumns[0],blockNumRows[1],blockNumColumns[1],*oper);

    if(verbose) std::cout << "done in " << timer.diff() << " s." << std::endl;

    // SIMPLE:
    // P = / F  0 \ / I D^-1Bt \
    //     \ B -S / \ 0 alphaI /

    // Getting the block structure of B
    MatrixBlockView B11,B12,B21,B22,B22base;

    // Building the First block
    // / F  0 \ = / F 0 \/ I 0 \/ I  0 \
    // \ B -S /   \ 0 I /\ B I /\ 0 -S /

    // / F 0 \
    // \ 0 I /
    if(verbose) std::cout << " P1a" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P1a(new matrix_type(map));
    P1a->setBlockStructure(blockNumRows,blockNumColumns);
    P1a->getMatrixBlockView(0,0,B11);
    P1a->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::copyBlock(F,B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P1a->globalAssemble();
    P1a->spy("p1a");
    boost::shared_ptr<parent_matrix_type> p1a = P1a;
    super_PtrType precForBlock1;
    precForBlock1.reset(new PreconditionerIfpack());
    precForBlock1->setDataFromGetPot(M_dataFile,M_section);
    this->pushBack(p1a,precForBlock1,notInversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    // / I 0 \
    // \ B I /
    if(verbose) std::cout << " P1b" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P1b(new matrix_type(map));
    P1b->setBlockStructure(blockNumRows,blockNumColumns);
    P1b->getMatrixBlockView(0,0,B11);
    P1b->getMatrixBlockView(1,0,B21);
    P1b->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::copyBlock(B,B21);
    (*P1b) *= -1;
    MatrixBlockUtils::createIdentityBlock(B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P1b->globalAssemble();
    P1b->spy("p1b");
    boost::shared_ptr<parent_matrix_type> p1b = P1b;
    this->pushBack(p1b,inversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    // / I  0 \
    // \ 0 -S /
    if(verbose) std::cout << " P1c" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P1c(new matrix_type(map));

    boost::shared_ptr<matrix_type> BBlockMat(new matrix_type(map));
    BBlockMat->setBlockStructure(blockNumRows,blockNumColumns);
    BBlockMat->getMatrixBlockView(1,0,B21);
    MatrixBlockUtils::copyBlock(B,B21);
    BBlockMat->globalAssemble();
    boost::shared_ptr<matrix_type> invDBlockMat(new matrix_type(map));
    invDBlockMat->setBlockStructure(blockNumRows,blockNumColumns);
    invDBlockMat->getMatrixBlockView(0,0,B11);
    if(M_SIMPLEType == "SIMPLE")
    {
        MatrixBlockUtils::createInvDiagBlock(F,B11);
        //MatrixBlockUtils::copyBlock(F,B11);
    }
    else if(M_SIMPLEType == "SIMPLEC")
    {
        MatrixBlockUtils::createInvLumpedBlock(F,B11);
    }
    *invDBlockMat *= -1.0;
    invDBlockMat->globalAssemble();
    boost::shared_ptr<matrix_type> tmpResultMat(new matrix_type(map));
    BBlockMat->multiply(false,
                        *invDBlockMat,false,
                        *tmpResultMat,true);
    BBlockMat.reset();
    invDBlockMat.reset();
    boost::shared_ptr<matrix_type> BtBlockMat(new matrix_type(map));
    BtBlockMat->setBlockStructure(blockNumRows,blockNumColumns);
    BtBlockMat->getMatrixBlockView(0,1,B12);
    MatrixBlockUtils::copyBlock(Bt,B12);
    BtBlockMat->globalAssemble();
    tmpResultMat->multiply(false,
                           *BtBlockMat,false,
                           *P1c,false);
    BtBlockMat.reset();
    tmpResultMat.reset();

    P1c->setBlockStructure(blockNumRows,blockNumColumns);
    P1c->getMatrixBlockView(0,0,B11);
    MatrixBlockUtils::createIdentityBlock(B11);
    P1c->globalAssemble();
    P1c->spy("p1c");
    boost::shared_ptr<parent_matrix_type> p1c = P1c;
    super_PtrType precForBlock2;
    precForBlock2.reset(new PreconditionerIfpack());
    precForBlock2->setDataFromGetPot(M_dataFile,M_section);
    this->pushBack(p1c,precForBlock2,notInversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    // Building the First block
    // / I  D^-1Bt \ = / D^-1 0      \/ I Bt \/ D 0 \
    // \ 0  alphaI /   \ 0    alphaI /\ 0 I  /\ 0 I /

    // / D^-1 0     \
    // \ 0    aphaI /
    if(verbose) std::cout << " P2a" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P2a(new matrix_type( map ));
    *P2a *= 0.0;
    P2a->setBlockStructure(blockNumRows,blockNumColumns);
    P2a->getMatrixBlockView(0,0,B11);
    P2a->getMatrixBlockView(1,1,B22);
    if(M_SIMPLEType == "SIMPLE")
    {
        MatrixBlockUtils::createDiagBlock(F,B11);
    }
    else if(M_SIMPLEType == "SIMPLEC")
    {
        MatrixBlockUtils::createLumpedBlock(F,B11);
    }
    MatrixBlockUtils::createScalarBlock(B22,1/M_dampingFactor);
    P2a->globalAssemble();
    P2a->spy("p2a");
    boost::shared_ptr<parent_matrix_type> p2a = P2a;
    this->pushBack(p2a,inversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    // / I Bt \
    // \ 0 I  /
    if(verbose) std::cout << " P2b" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P2b(new matrix_type(map));
    P2b->setBlockStructure(blockNumRows,blockNumColumns);
    P2b->getMatrixBlockView(0,0,B11);
    P2b->getMatrixBlockView(0,1,B12);
    P2b->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::copyBlock(Bt,B12);
    (*P2b) *= -1;
    MatrixBlockUtils::createIdentityBlock(B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P2b->globalAssemble();
    P2b->spy("p2b");
    boost::shared_ptr<parent_matrix_type> p2b = P2b;
    this->pushBack(p2b,inversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    // / D 0 \
    // \ 0 I /
    if(verbose) std::cout << " P2c" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P2c(new matrix_type( map ));
    *P2c *= 0.0;
    P2c->setBlockStructure(blockNumRows,blockNumColumns);
    P2c->getMatrixBlockView(0,0,B11);
    P2c->getMatrixBlockView(1,1,B22);
    if(M_SIMPLEType == "SIMPLE")
    {
        MatrixBlockUtils::createInvDiagBlock(F,B11);
    }
    else if(M_SIMPLEType == "SIMPLEC")
    {
        MatrixBlockUtils::createInvLumpedBlock(F,B11);
    }
    MatrixBlockUtils::createIdentityBlock(B22);
    P2c->globalAssemble();
    P2c->spy("p2c");
    boost::shared_ptr<parent_matrix_type> p2c = P2c;
    this->pushBack(p2c,inversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    this->M_preconditionerCreated = true;

    if(verbose) std::cout << "[DEBUG] number of operators: " << numOperators() << std::endl;

    if(verbose) std::cout << "      >All the blocks are built" << std::endl
                          << "      >";
    return ( EXIT_SUCCESS );
}
int PreconditionerSIMPLE::numBlocksRows() const
{
    return 2;
}

int PreconditionerSIMPLE::numBlocksColumns() const
{
    return 2;
}

void PreconditionerSIMPLE::setDataFromGetPot( const GetPot& dataFile,
                                           const std::string& section )
{
    M_dataFile   = dataFile;
    M_section    = section;
    createSIMPLEList(M_list, dataFile, section, "SIMPLE");

    M_velocityBlockSize = this->M_list.get("blocks: velocity block size", -1);
    M_pressureBlockSize = this->M_list.get("blocks: pressure block size", -1);
    M_precType          = this->M_list.get("prectype", "SIMPLE");
    M_SIMPLEType        = this->M_list.get("SIMPLE Type", "SIMPLE");
}

void PreconditionerSIMPLE::setFESpace(FESpace_ptr uFESpace,FESpace_ptr pFESpace){
    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_adrPressureAssembler.setup(pFESpace,uFESpace); // p,beta=u
    M_adrVelocityAssembler.setup(uFESpace,uFESpace); // u,beta=u
    M_beta.reset(new vector_type(M_uFESpace->map() + M_pFESpace->map()));
    *M_beta *= 0;
}

void PreconditionerSIMPLE::setTimestep(const Real& timestep)
{
    M_timestep = timestep;
}

void PreconditionerSIMPLE::setViscosity(const Real& viscosity)
{
    M_viscosity = viscosity;
}

void PreconditionerSIMPLE::setDensity(const Real& density)
{
    M_density = density;
}

void PreconditionerSIMPLE::setDampingFactor(const Real& dampingFactor)
{
    M_dampingFactor = dampingFactor;
}

} // namespace LifeV
