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
#include <life/lifecore/LifeChrono.hpp>
#include <lifemc/lifearray/MatrixBlock.hpp>
#include <lifemc/lifearray/MatrixBlockView.hpp>
#include <lifemc/lifearray/MatrixBlockUtils.hpp>

#define shuttleworth

#ifdef shuttleworth
#include <EpetraExt_MatrixMatrix.h>
#endif

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
    oper->spy("A");

    LifeChrono timer;

    // Getting the block structure of A
    // / F Bt \
    // \ B C  /
    if(verbose) std::cout << std::endl << "      >Getting the structure of A... ";
    timer.start();
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
    if(verbose) std::cout << "done in " << timer.diff() << " s." << std::endl;

    // PCD:
    // / F Bt \   / I  0 \ / I Bt \ / F 0 \
    // \ 0 -S / = \ 0 -S / \ 0 I  / \ 0 I /

    // PCD^-1:
    // / F  Bt \^-1   / F^-1 0 \ / I -Bt \ / I  0    \
    // \ 0 -S  /    = \ 0    I / \ 0  I  / \ 0 -S^-1 /

    // Getting the block structure of B
    MatrixBlockView B11,B12,B21,B22,B22base;

    // Building the block
    // / I  0 \   / I  0       \   / I  0  \ / I 0     \ / I 0  \
    // \ 0 -S / = \ 0 -ApFp^-1 / = \ 0 -Ap / \ 0 Fp^-1 / \ 0 Mp /
    if(verbose) std::cout << " Building Fp" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> PFp(new matrix_type( map ));
    *PFp *= 0.0;
    PFp->setBlockStructure(blockNumRows,blockNumColumns);
    PFp->getMatrixBlockView(1,1,B22);
    M_adrPressureAssembler.addDiffusion(PFp,-M_viscosity/M_density,B22.firstRowIndex(),B22.firstColumnIndex());
    M_adrPressureAssembler.addAdvection(PFp,*M_beta,B22.firstRowIndex(),B22.firstColumnIndex());
    M_adrPressureAssembler.addMass(PFp,1.0/M_timestep,B22.firstRowIndex(),B22.firstColumnIndex());
    PFp->globalAssemble();
    PFp->spy("pFp");
    boost::shared_ptr<parent_matrix_type> pFp = PFp;
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    if(verbose) std::cout << " Building Ap" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> PAp(new matrix_type( map ));
    *PAp *= 0.0;
    PAp->setBlockStructure(blockNumRows,blockNumColumns);
    PAp->getMatrixBlockView(0,0,B11);
    PAp->getMatrixBlockView(1,1,B22);
#ifdef shuttleworth

    Epetra_CrsMatrix* tmpCrsMatrix(NULL);
    tmpCrsMatrix = PAp->matrixPtr().get();
    EpetraExt::MatrixMatrix::Add(*(pFp->matrixPtr()), false, 0.5*M_density/M_viscosity,
                                 *(pFp->matrixPtr()), true,  0.5*M_density/M_viscosity,
                                 tmpCrsMatrix);

    MatrixBlockUtils::createScalarBlock(B11,1);
#else
    M_adrPressureAssembler.addDiffusion(PAp,-1.0,B22.firstRowIndex(),B22.firstColumnIndex());
    MatrixBlockUtils::createIdentityBlock(B11);
#endif
    PAp->globalAssemble();
    boost::shared_ptr<parent_matrix_type> pAp = PAp;
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    if(verbose) std::cout << " Building Mp" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> PMp(new matrix_type( map ));
    *PMp *= 0.0;
    PMp->setBlockStructure(blockNumRows,blockNumColumns);
    PMp->getMatrixBlockView(0,0,B11);
    PMp->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::createIdentityBlock(B11);
    M_adrPressureAssembler.addMass(PMp,1.0,B22.firstRowIndex(),B22.firstColumnIndex());
    PMp->globalAssemble();
    boost::shared_ptr<parent_matrix_type> pMp = PMp;
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    if(verbose) std::cout << " P1a" << std::endl;
    timer.start();
    //pAp->spy("p1a");
    super_PtrType precForBlock1;
    precForBlock1.reset(new PreconditionerIfpack());
    precForBlock1->setDataFromGetPot(M_dataFile,M_section);
    this->pushBack(pAp,precForBlock1,notInversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    if(verbose) std::cout << " P1b" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P1b(new matrix_type( map ));
    P1b->setBlockStructure(blockNumRows,blockNumColumns);
    *P1b += *PFp;
    P1b->getMatrixBlockView(0,0,B11);
    MatrixBlockUtils::createIdentityBlock(B11);
    P1b->globalAssemble();
    boost::shared_ptr<parent_matrix_type> p1b = P1b;
    //p1b->spy("p1b");
    this->pushBack(p1b,inversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    if(verbose) std::cout << " P1c" << std::endl;
    timer.start();
    //pMp->spy("p1c");
    super_PtrType precForBlock2;
    precForBlock2.reset(new PreconditionerML());
    //prec->createParametersList(prec->parametersList(), M_dataFile, M_section, "Ifpack");
    precForBlock2->setDataFromGetPot(M_dataFile,M_section);
    this->pushBack(pMp,precForBlock2,notInversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    // Building the block (the block is inversed)
    // / I -Bt \
    // \ 0  I  /
    if(verbose) std::cout << " P2" << std::endl;
    timer.start();
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
    //P2->spy("p2");
    boost::shared_ptr<parent_matrix_type> p2 = P2;
    this->pushBack(p2,inversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    // Building the block
    // / F 0 \
    // \ 0 I /
    if(verbose) std::cout << " P3" << std::endl;
    timer.start();
    boost::shared_ptr<matrix_type> P3(new matrix_type(map));
    P3->setBlockStructure(blockNumRows,blockNumColumns);
    P3->getMatrixBlockView(0,0,B11);
    P3->getMatrixBlockView(1,1,B22);
    MatrixBlockUtils::copyBlock(F,B11);
    MatrixBlockUtils::createIdentityBlock(B22);
    P3->globalAssemble();
    //P3->spy("p3");
    boost::shared_ptr<parent_matrix_type> p3 = P3;
    super_PtrType precForBlock3;
    precForBlock3.reset(new PreconditionerIfpack());
    precForBlock3->setDataFromGetPot(M_dataFile,M_section);
    this->pushBack(p3,precForBlock3,notInversed,notTransposed);
    if(verbose) std::cout << " done in " << timer.diff() << " s." << std::endl;

    this->M_preconditionerCreated = true;

    if(verbose) std::cout << "[DEBUG] number of operators: " << numOperators() << std::endl;

    if(verbose) std::cout << "      >All the blocks are built" << std::endl
                          << "      >";
    return ( EXIT_SUCCESS );
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
