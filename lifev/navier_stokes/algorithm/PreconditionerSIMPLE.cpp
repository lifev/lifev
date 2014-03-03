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
    @brief PreconditionerSIMPLE

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 24-01-2011
 */

#include <vector>
#include "PreconditionerSIMPLE.hpp"
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
// #include <lifev/core/algorithm/PreconditionerML2.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

namespace LifeV
{

PreconditionerSIMPLE::PreconditionerSIMPLE ( boost::shared_ptr<Epetra_Comm> comm ) :
    PreconditionerComposition ( comm ),
    M_velocityBlockSize       ( -1 ),
    M_pressureBlockSize       ( -1 ),
    M_dampingFactor           ( 0.5 ),
    M_SIMPLEType              ( "SIMPLE" )
{

}

PreconditionerSIMPLE::~PreconditionerSIMPLE()
{

}

void
PreconditionerSIMPLE::createParametersList ( list_Type&         list,
                                             const GetPot&      dataFile,
                                             const std::string& section,
                                             const std::string& subsection )
{
    const bool verbose ( M_comm->MyPID() == 0 );

    bool displayList = dataFile ( ( section + "/displayList" ).data(), false);

    std::string precType = dataFile ( ( section + "/prectype" ).data(), "SIMPLE" );
    list.set ( "prectype", precType );
    string SIMPLEType = dataFile ( ( section + "/" + subsection + "/SIMPLE_type" ).data(), "SIMPLE" );

    std::string fluidPrec = dataFile ( ( section + "/" + subsection + "/subprecs/fluid_prec" ).data(), "ML" );
    list.set ( "subprecs: fluid prec", fluidPrec );
    std::string fluidPrecDataSection = dataFile ( ( section + "/" + subsection + "/subprecs/fluid_prec_data_section" ).data(), "" );
    list.set ( "subprecs: fluid prec data section", ( fluidPrecDataSection ).data() );

    std::string schurPrec = dataFile ( ( section + "/" + subsection + "/subprecs/schur_prec" ).data(), "ML" );
    list.set ( "subprecs: Schur prec", schurPrec );
    std::string schurPrecDataSection = dataFile ( ( section + "/" + subsection + "/subprecs/schur_prec_data_section" ).data(), "" );
    list.set ( "subprecs: Schur prec data section", ( schurPrecDataSection ).data() );

    list.set ( "SIMPLE Type", SIMPLEType );

    if ( displayList && verbose )
    {
        std::cout << "SIMPLE parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        list.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

Real
PreconditionerSIMPLE::condest()
{
    return 0.0;
}

int
PreconditionerSIMPLE::buildPreconditioner ( matrixPtr_Type& oper )
{
    if ( M_velocityBlockSize < 0 || M_pressureBlockSize < 0 )
    {
        std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
        exit ( 0 );
    }

    // Make sure that the preconditioner is reset
    this->resetPreconditioner();

    const bool verbose ( M_comm->MyPID() == 0 );

    std::vector<UInt> blockNumRows ( 2, 0 );
    blockNumRows[0] = M_velocityBlockSize;
    blockNumRows[1] = M_pressureBlockSize;
    std::vector<UInt> blockNumColumns ( blockNumRows );

    const bool inversed ( true );
    const bool notInversed ( false );
    const bool notTransposed ( false );

    map_Type map ( oper->map() );
    //oper->spy( "A" );

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
    matrixBlockView_Type F, Bt, B, C;
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
     * SIMPLE:
     * P = / F  0 \ / I D^-1Bt \
     *     \ B -S / \ 0 alphaI /
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
        std::cout << "       Block 1 (F)" << std::endl;
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
//    if ( M_fluidPrec == "ML2" )
//    {
//        PreconditionerML2* tmpPrecPtr = dynamic_cast<PreconditionerML2*> ( precForBlock1.get() );
//        tmpPrecPtr->setFESpace ( M_uFESpace, M_pFESpace );
//    }
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
        std::cout << "       Block 1 (B)" << std::endl;
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
        std::cout << "       Block 1 (Schur)" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P1c ( new matrixBlock_Type ( map ) );

    boost::shared_ptr<matrixBlock_Type> BBlockMat ( new matrixBlock_Type ( map ) );
    BBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    BBlockMat->blockView ( 1, 0, B21 );
    MatrixEpetraStructuredUtility::copyBlock ( B, B21 );
    BBlockMat->globalAssemble();
    boost::shared_ptr<matrixBlock_Type> invDBlockMat ( new matrixBlock_Type ( map ) );
    invDBlockMat->setBlockStructure ( blockNumRows, blockNumColumns );
    invDBlockMat->blockView ( 0, 0, B11 );
    if ( M_SIMPLEType == "SIMPLE" )
    {
        MatrixEpetraStructuredUtility::createInvDiagBlock ( F, B11 );
    }
    else if ( M_SIMPLEType == "SIMPLEC" )
    {
        MatrixEpetraStructuredUtility::createInvLumpedBlock ( F, B11 );
    }
    *invDBlockMat *= -1.0;
    invDBlockMat->globalAssemble();
    boost::shared_ptr<matrixBlock_Type> tmpResultMat ( new matrixBlock_Type ( map ) );
    BBlockMat->multiply ( false,
                          *invDBlockMat, false,
                          *tmpResultMat, true );
    BBlockMat.reset();
    invDBlockMat.reset();
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
     * / I D^-1Bt \ = / D^-1   0    \/ I Bt \/ D 0 \
     * \ 0 alphaI /   \ 0    alphaI /\ 0 I  /\ 0 I /
     */
    if ( verbose )
    {
        std::cout << "       Block 2 (D^-1,alpha I)" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2a ( new matrixBlock_Type ( map ) );
    *P2a *= 0.0;
    P2a->setBlockStructure ( blockNumRows, blockNumColumns );
    P2a->blockView ( 0, 0, B11 );
    P2a->blockView ( 1, 1, B22 );
    if ( M_SIMPLEType == "SIMPLE" )
    {
        MatrixEpetraStructuredUtility::createDiagBlock ( F, B11 );
    }
    else if ( M_SIMPLEType == "SIMPLEC" )
    {
        MatrixEpetraStructuredUtility::createLumpedBlock ( F, B11 );
    }
    MatrixEpetraStructuredUtility::createScalarBlock ( B22, 1 / M_dampingFactor );
    P2a->globalAssemble();
    boost::shared_ptr<matrix_Type> p2a = P2a;
    this->pushBack ( p2a, inversed, notTransposed );
    if ( verbose )
    {
        std::cout << "       done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * / I -Bt \
     * \ 0  I  /
     */
    if ( verbose )
    {
        std::cout << "       Block 2 (Bt)" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2b ( new matrixBlock_Type ( map ) );
    P2b->setBlockStructure ( blockNumRows, blockNumColumns );
    P2b->blockView ( 0, 0, B11 );
    P2b->blockView ( 0, 1, B12 );
    P2b->blockView ( 1, 1, B22 );
    MatrixEpetraStructuredUtility::copyBlock ( Bt, B12 );
    //( *P2b ) *= -1; // We inverse already the block
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
     * / D 0 \
     * \ 0 I /
     */
    if ( verbose )
    {
        std::cout << "       Block2 (D)" << std::endl;
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2c ( new matrixBlock_Type ( map ) );
    *P2c *= 0.0;
    P2c->setBlockStructure ( blockNumRows, blockNumColumns );
    P2c->blockView ( 0, 0, B11 );
    P2c->blockView ( 1, 1, B22 );
    if ( M_SIMPLEType == "SIMPLE" )
    {
        MatrixEpetraStructuredUtility::createInvDiagBlock ( F, B11 );
    }
    else if ( M_SIMPLEType == "SIMPLEC" )
    {
        MatrixEpetraStructuredUtility::createInvLumpedBlock ( F, B11 );
    }
    MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
    P2c->globalAssemble();
    boost::shared_ptr<matrix_Type> p2c = P2c;
    this->pushBack ( p2c, inversed, notTransposed );
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
PreconditionerSIMPLE::numBlocksRows() const
{
    return 2;
}

int
PreconditionerSIMPLE::numBlocksColumns() const
{
    return 2;
}

void
PreconditionerSIMPLE::setDataFromGetPot ( const GetPot& dataFile,
                                          const std::string& section )
{
    M_dataFile   = dataFile;
    this->createParametersList ( M_list, dataFile, section, "SIMPLE" );
    this->setParameters ( M_list );
}

void
PreconditionerSIMPLE::setParameters ( Teuchos::ParameterList& list )
{
    M_precType         = list.get ( "prectype", "SIMPLE" );
    M_SIMPLEType       = list.get ( "SIMPLE Type", "SIMPLE" );

    M_fluidPrec        = list.get ( "subprecs: fluid prec", "ML" );
    M_fluidDataSection = list.get ( "subprecs: fluid prec data section", "" );

    M_schurPrec        = list.get ( "subprecs: Schur prec", "ML" );
    M_schurDataSection = list.get ( "subprecs: Schur prec data section", "" );
}

void
PreconditionerSIMPLE::setFESpace ( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace )
{
    // We setup the size of the blocks
    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_velocityBlockSize = uFESpace->fieldDim() * uFESpace->dof().numTotalDof();
    M_pressureBlockSize = pFESpace->dof().numTotalDof();
}

void
PreconditionerSIMPLE::setDampingFactor ( const Real& dampingFactor )
{
    M_dampingFactor = dampingFactor;
}

} // namespace LifeV
