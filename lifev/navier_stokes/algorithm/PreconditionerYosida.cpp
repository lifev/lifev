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
    @brief PreconditionerYosida

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 13-10-2011
 */

#include <vector>
#include "PreconditionerYosida.hpp"
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerML2.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

namespace LifeV
{

PreconditionerYosida::PreconditionerYosida ( const  boost::shared_ptr<Epetra_Comm>& comm ) :
    PreconditionerComposition ( comm ),
    M_velocityBlockSize       ( -1 ),
    M_pressureBlockSize       ( -1 ),
    M_timestep                ( 1.0 )
{

    M_uFESpace.reset();
    M_pFESpace.reset();

}

PreconditionerYosida::~PreconditionerYosida()
{

}

void
PreconditionerYosida::createParametersList ( list_Type&         list,
                                             const GetPot&      dataFile,
                                             const std::string& section,
                                             const std::string& subsection )
{
    const bool verbose ( M_comm->MyPID() == 0 );

    bool displayList = dataFile ( ( section + "/displayList" ).data(), false);

    std::string precType = dataFile ( ( section + "/prectype" ).data(), "Yosida" );
    list.set ( "prectype", precType );

    std::string fluidPrec = dataFile ( ( section + "/" + subsection + "/subprecs/fluid_prec" ).data(), "ML" );
    list.set ( "subprecs: fluid prec", fluidPrec );
    std::string fluidPrecDataSection = dataFile ( ( section + "/" + subsection + "/subprecs/fluid_prec_data_section" ).data(), "" );
    list.set ( "subprecs: fluid prec data section", ( fluidPrecDataSection ).data() );

    std::string schurPrec = dataFile ( ( section + "/" + subsection + "/subprecs/schur_prec" ).data(), "ML" );
    list.set ( "subprecs: Schur prec", schurPrec );
    std::string schurPrecDataSection = dataFile ( ( section + "/" + subsection + "/subprecs/schur_prec_data_section" ).data(), "" );
    list.set ( "subprecs: Schur prec data section", ( schurPrecDataSection ).data() );

    if ( displayList && verbose )
    {
        std::cout << "Yosida parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        list.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

Real
PreconditionerYosida::condest()
{
    return 0.0;
}

int
PreconditionerYosida::buildPreconditioner ( matrixPtr_Type& oper )
{
    if ( M_velocityBlockSize < 0 || M_pressureBlockSize < 0 )
    {
        std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
        exit ( 0 );
    }

    const bool verbose ( M_comm->MyPID() == 0 );

    // Make sure that the preconditioner is reset
    this->resetPreconditioner();

    std::vector<UInt> blockNumRows ( 2, 0 );
    blockNumRows[0] = M_velocityBlockSize;
    blockNumRows[1] = M_pressureBlockSize;
    std::vector<UInt> blockNumColumns ( blockNumRows );

    const bool inversed ( true );
    const bool notInversed ( false );
    const bool notTransposed ( false );

    map_Type map ( oper->map() );

    LifeChrono timer;

    /*
     * Getting the block structure of A
     * / F Bt \
     * \ B C  /
     */
    if ( verbose )
    {
        std::cout << "      >Getting the structure of A... ";
    }
    timer.start();
    matrixBlockView_Type F, Bt, B, C, M;
    //oper.blockView( 0, 0, F );
    F.setup ( 0, 0, blockNumRows[0], blockNumColumns[0], oper.get() );

    //oper.blockView( 0, 1, Bt );
    Bt.setup ( 0, blockNumColumns[0], blockNumRows[0], blockNumColumns[1], oper.get() );

    //oper.blockView( 1, 0, B );
    B.setup ( blockNumRows[0], 0, blockNumRows[1], blockNumColumns[0], oper.get() );

    //oper.blockView( 1, 1, C );
    C.setup ( blockNumRows[0], blockNumColumns[0], blockNumRows[1], blockNumColumns[1], oper.get() );

    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * Yosida:
     * P = / F  0 \ / I F^-1Bt \
     *     \ B -S / \ 0   I    /
     *     where S=dt*B*M_lumped*Bt
     */

    // Getting the block structure of B
    matrixBlockView_Type B11, B12, B21, B22, B22base;

    /*
     * Building the First block
     * / F  0 \ = / F 0 \/ I 0 \/ I  0 \
     * \ B -S /   \ 0 I /\ B I /\ 0 -S /
     *
     * / F 0 \
     * \ 0 I /
     */
    if ( verbose )
    {
        std::cout << "       Block 1 ( F )" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P1a ( new matrixBlock_Type ( map ) );
    P1a->setBlockStructure ( blockNumRows, blockNumColumns );
    P1a->blockView ( 0, 0, B11 );
    P1a->blockView ( 1, 1, B22 );
    MatrixEpetraStructuredUtility::copyBlock ( F, B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
    P1a->globalAssemble();
    boost::shared_ptr<matrix_Type> p1a = P1a;
    superPtr_Type precForBlock1 ( PRECFactory::instance().createObject ( M_fluidPrec ) );
    precForBlock1->setDataFromGetPot ( M_dataFile, M_fluidDataSection );
    if ( M_fluidPrec == "ML2" )
    {
        PreconditionerML2* tmpPrecPtr = dynamic_cast<PreconditionerML2*> ( precForBlock1.get() );
        tmpPrecPtr->setFESpace ( M_uFESpace, M_pFESpace );
    }
    this->pushBack ( p1a, precForBlock1, notInversed, notTransposed );
    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * / I 0 \
     * \ B I /
     */
    if ( verbose )
    {
        std::cout << "       Block 1 ( B )" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P1b ( new matrixBlock_Type ( map ) );
    P1b->setBlockStructure ( blockNumRows, blockNumColumns );
    P1b->blockView ( 0, 0, B11 );
    P1b->blockView ( 1, 0, B21 );
    P1b->blockView ( 1, 1, B22 );
    MatrixEpetraStructuredUtility::copyBlock ( B, B21 );
    ( *P1b ) *= -1;
    MatrixEpetraStructuredUtility::createIdentityBlock ( B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
    P1b->globalAssemble();
    boost::shared_ptr<matrix_Type> p1b = P1b;
    this->pushBack ( p1b, inversed, notTransposed );
    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * / I  0 \
     * \ 0 -S /
     */
    if ( verbose )
    {
        std::cout << "       Block 1 ( Schur )" << std::endl;
    }
    timer.start();
    // Computing the mass matrix
    boost::shared_ptr<matrixBlock_Type> massMat ( new matrixBlock_Type ( map ) );
    massMat->setBlockStructure ( blockNumRows, blockNumColumns );
    massMat->blockView ( 0, 0, M );
    M_adrVelocityAssembler.addMass ( massMat, 1.0 / M_timestep, M.firstRowIndex(), M.firstColumnIndex() );
    massMat->globalAssemble();

    boost::shared_ptr<matrixBlock_Type> P1c ( new matrixBlock_Type ( map ) );

    boost::shared_ptr<matrixBlock_Type> BBlockMat ( new matrixBlock_Type ( map ) );
    BBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    BBlockMat->blockView ( 1, 0, B21 );
    MatrixEpetraStructuredUtility::copyBlock ( B, B21 );
    BBlockMat->globalAssemble();
    boost::shared_ptr<matrixBlock_Type> invLumpedMassBlockMat ( new matrixBlock_Type ( map ) );
    invLumpedMassBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    invLumpedMassBlockMat->blockView ( 0, 0, B11 );
    MatrixEpetraStructuredUtility::createInvDiagBlock ( M, B11 );
    massMat.reset();               // Free memory
    *invLumpedMassBlockMat *= -1.0;
    invLumpedMassBlockMat->globalAssemble();
    boost::shared_ptr<matrixBlock_Type> tmpResultMat ( new matrixBlock_Type ( map ) );
    BBlockMat->multiply ( false,
                          *invLumpedMassBlockMat, false,
                          *tmpResultMat, true );
    BBlockMat.reset();             // Free memory
    invLumpedMassBlockMat.reset(); // Free memory
    boost::shared_ptr<matrixBlock_Type> BtBlockMat ( new matrixBlock_Type ( map ) );
    BtBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    BtBlockMat->blockView ( 0, 1, B12 );
    MatrixEpetraStructuredUtility::copyBlock ( Bt, B12 );
    BtBlockMat->globalAssemble();
    tmpResultMat->multiply ( false,
                             *BtBlockMat, false,
                             *P1c, false );
    BtBlockMat.reset();
    tmpResultMat.reset();

    P1c->setBlockStructure ( blockNumRows, blockNumColumns );
    P1c->blockView ( 0, 0, B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B11 );
    P1c->globalAssemble();
    boost::shared_ptr<matrix_Type> p1c = P1c;
    superPtr_Type precForBlock2 ( PRECFactory::instance().createObject ( M_schurPrec ) );
    precForBlock2->setDataFromGetPot ( M_dataFile, M_schurDataSection );
    this->pushBack ( p1c, precForBlock2, notInversed, notTransposed );
    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * Building the Second block
     * / I  F^-1Bt \ = / F^-1 0 \/ I Bt \/ F 0 \
     * \ 0    I    /   \ 0    I /\ 0 I  /\ 0 I /
     */
    if ( verbose )
    {
        std::cout << "       Block 2 ( F^-1, I )" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2a ( new matrixBlock_Type ( map ) );
    *P2a *= 0.0;
    P2a->setBlockStructure ( blockNumRows, blockNumColumns );
    P2a->blockView ( 0, 0, B11 );
    P2a->blockView ( 1, 1, B22 );
    MatrixEpetraStructuredUtility::copyBlock ( F, B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
    P2a->globalAssemble();
    boost::shared_ptr<matrix_Type> p2a = P2a;
    this->pushBack ( p2a, inversed, notTransposed );
    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * / I Bt \
     * \ 0 I  /
     */
    if ( verbose )
    {
        std::cout << "       Block 2 ( Bt )" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2b ( new matrixBlock_Type ( map ) );
    P2b->setBlockStructure ( blockNumRows, blockNumColumns );
    P2b->blockView ( 0, 0, B11 );
    P2b->blockView ( 0, 1, B12 );
    P2b->blockView ( 1, 1, B22 );
    MatrixEpetraStructuredUtility::copyBlock ( Bt, B12 );
    ( *P2b ) *= -1; // We inverse already the block
    MatrixEpetraStructuredUtility::createIdentityBlock ( B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
    P2b->globalAssemble();
    boost::shared_ptr<matrix_Type> p2b = P2b;
    this->pushBack ( p2b, inversed, notTransposed );
    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * / F 0 \
     * \ 0 I /
     */
    if ( verbose )
    {
        std::cout << "       Block2 ( F )" << std::endl;
    }
    timer.start();
    // Here we recycle the first block that we have built (it is the same)
    this->pushBack ( p1a, precForBlock1, notInversed, notTransposed );
    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    this->M_preconditionerCreated = true;

    if ( verbose ) std::cout << "      >All the blocks are built" << std::endl
                                 << "      >";
    return ( EXIT_SUCCESS );
}

int
PreconditionerYosida::numBlocksRows() const
{
    return 2;
}

int
PreconditionerYosida::numBlocksColumns() const
{
    return 2;
}

void
PreconditionerYosida::setDataFromGetPot ( const GetPot& dataFile,
                                          const std::string& section )
{
    M_dataFile   = dataFile;
    this->createParametersList ( M_list, dataFile, section, "Yosida" );
    this->setParameters ( M_list );
}

void
PreconditionerYosida::setParameters ( Teuchos::ParameterList& list )
{
    M_precType         = list.get ( "prectype", "Yosida" );

    M_fluidPrec        = list.get ( "subprecs: fluid prec", "ML" );
    M_fluidDataSection = list.get ( "subprecs: fluid prec data section", "" );

    M_schurPrec        = list.get ( "subprecs: Schur prec", "ML" );
    M_schurDataSection = list.get ( "subprecs: Schur prec data section", "" );
}

void
PreconditionerYosida::setFESpace ( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace )
{
    // We setup the size of the blocks
    M_velocityBlockSize = uFESpace->fieldDim() * uFESpace->dof().numTotalDof();
    M_pressureBlockSize = pFESpace->dof().numTotalDof();

    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_adrVelocityAssembler.setup ( uFESpace, uFESpace ); // u,beta=u
}

void
PreconditionerYosida::setTimestep ( const Real& timestep )
{
    M_timestep = timestep;
}

} // namespace LifeV
