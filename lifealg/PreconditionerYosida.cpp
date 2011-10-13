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
    @brief PreconditionerPressureCorrection

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 26-05-2011
 */

#include <vector>
#include "PreconditionerPressureCorrection.hpp"
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifefem/BCManage.hpp>
#include <lifemc/lifearray/MatrixBlock.hpp>
#include <lifemc/lifearray/MatrixBlockView.hpp>
#include <lifemc/lifearray/MatrixBlockUtils.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <EpetraExt_MatrixMatrix.h>
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV {

PreconditionerPressureCorrection::PreconditionerPressureCorrection( const  boost::shared_ptr<Epetra_Comm>& comm ):
    PreconditionerComposition ( comm ),
    M_velocityBlockSize       ( -1 ),
    M_pressureBlockSize       ( -1 ),
    M_timestep                ( 1.0 ),
    M_viscosity               ( 1.0 ),
    M_density                 ( 1.0 )
{

}

PreconditionerPressureCorrection::~PreconditionerPressureCorrection()
{

}

void
PreconditionerPressureCorrection::createParametersList( list_Type&         list,
                                                        const GetPot&      dataFile,
                                                        const std::string& section,
                                                        const std::string& subSection )
{
    createPressureCorrectionList( list, dataFile, section, subSection, M_comm->MyPID() == 0 );
}

void
PreconditionerPressureCorrection::createPressureCorrectionList( list_Type&         list,
                                                                const GetPot&      dataFile,
                                                                const std::string& section,
                                                                const std::string& subsection,
                                                                const bool& verbose )
{
    bool displayList = dataFile( ( section + "/displayList" ).data(), false );

    std::string precType = dataFile( ( section + "/prectype" ).data(), "PCD" );
    list.set( "prectype", precType );

    std::string fluidPrec = dataFile( ( section + "/" + subsection + "/subprecs/fluid_prec" ).data(), "ML" );
    list.set( "subprecs: fluid prec", fluidPrec );
    std::string fluidPrecDataSection = dataFile( ( section + "/" + subsection + "/subprecs/fluid_prec_data_section" ).data(), "" );
    list.set( "subprecs: fluid prec data section", ( section + "/" + subsection+"/subprecs/"+fluidPrecDataSection ).data() );

    std::string pressureMassPrec = dataFile( ( section + "/" + subsection + "/subprecs/pressure_mass_prec" ).data(), "ML" );
    list.set( "subprecs: pressure mass prec", pressureMassPrec );
    std::string pressureMassPrecDataSection = dataFile( ( section + "/" + subsection + "/subprecs/pressure_mass_prec_data_section" ).data(), "" );
    list.set( "subprecs: pressure mass prec data section", ( section + "/" + subsection+"/subprecs/"+pressureMassPrecDataSection ).data() );

    if ( displayList && verbose ) list.print( std::cout );
}

Real
PreconditionerPressureCorrection::condest()
{
    return 0.0;
}

int
PreconditionerPressureCorrection::buildPreconditioner( matrixPtr_Type& oper )
{
    if ( ( M_uFESpace.get() == NULL ) || ( M_pFESpace.get() == NULL ) )
    {
        std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
        exit( 0 );
    }

    bool verbose( false );
    if ( M_comm->MyPID() == 0 ) verbose = true;

    // Make sure that an operator exists
    initializeOperator();

    std::vector<UInt> blockNumRows( 2, 0 );
    blockNumRows[0] = M_velocityBlockSize;
    blockNumRows[1] = M_pressureBlockSize;
    std::vector<UInt> blockNumColumns( blockNumRows );

    bool inversed( true );
    bool notInversed( false );
    //bool transposed( true );
    bool notTransposed( false );

    map_Type map( oper->map() );
    //oper->spy( "A" );

    LifeChrono timer;

    /*
     * Getting the block structure of A
     * / F Bt \
     * \ B C  /
     */
    if ( verbose ) std::cout << std::endl << "      >Getting the structure of A... ";
    timer.start();
    MatrixBlockView F, Bt, B, C;
    //oper.getMatrixBlockView( 0, 0, F );
    F.setup( 0, 0, blockNumRows[0], blockNumColumns[0], *oper );

    //oper.getMatrixBlockView( 0, 1, Bt );
    Bt.setup( 0, blockNumColumns[0], blockNumRows[0], blockNumColumns[1], *oper );

    //oper.getMatrixBlockView( 1, 0, B );
    B.setup( blockNumRows[0], 0, blockNumRows[1], blockNumColumns[0], *oper );

    //oper.getMatrixBlockView( 1, 1, C );
    C.setup( blockNumRows[0], blockNumColumns[0], blockNumRows[1], blockNumColumns[1], *oper );
    if ( verbose ) std::cout << "done in " << timer.diff() << " s." << std::endl;

    /*
     * Pressure Correction:
     *     / F 0 \ / I 0 \ / dt*Mu^-1  0 \ / I  Bt \ / Mu/dt 0 \
     * P = \ 0 I / \ B I / \   0      -S / \ 0  I  / \   0   I /
     *
     * Pressure Correction^-1:
     *        / dt*Mu^-1 0 \ / I -Bt \ / Mu/dt  0   \ /  I 0 \ / F^-1 0 \
     * P^-1 = \   0      I / \ 0  I  / \ 0     -S^1 / \ -B I / \ 0    I /
     */

    // Getting the block structure of B
    MatrixBlockView B11, B12, B21, B22, B22base;
    /*
     * Building the block
     * / F 0 \
     * \ 0 I /
     */
     if ( verbose ) std::cout << "       Fluid block... ";
     timer.start();
     boost::shared_ptr<matrixBlock_Type> P3( new matrixBlock_Type( map ) );
     P3->setBlockStructure( blockNumRows, blockNumColumns );
     P3->getMatrixBlockView( 0, 0, B11 );
     P3->getMatrixBlockView( 1, 1, B22 );
     MatrixBlockUtils::copyBlock( F, B11 );
     MatrixBlockUtils::createIdentityBlock( B22 );
     P3->globalAssemble();

     boost::shared_ptr<matrix_Type> p3 = P3;
     superPtr_Type precForBlock3( PRECFactory::instance().createObject( M_fluidPrec ) );
     precForBlock3->setDataFromGetPot( M_dataFile, M_fluidPrecDataSection );
     this->pushBack( p3,precForBlock3, notInversed, notTransposed );
     if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

     /*
      * Building the block (the block is inversed)
      * /  I 0 \
      * \ -B I /
      */
     if ( verbose ) std::cout << "       Divergence block... ";
     timer.start();
     boost::shared_ptr<matrixBlock_Type> P2( new matrixBlock_Type( map ) );
     P2->setBlockStructure( blockNumRows, blockNumColumns );
     P2->getMatrixBlockView( 0, 0, B11 );
     P2->getMatrixBlockView( 1, 0, B21 );
     P2->getMatrixBlockView( 1, 1, B22 );
     MatrixBlockUtils::copyBlock( B, B21 );
     ( *P2 ) *= -1;
     MatrixBlockUtils::createIdentityBlock( B11 );
     MatrixBlockUtils::createIdentityBlock( B22 );
     P2->globalAssemble();

     boost::shared_ptr<matrix_Type> p2 = P2;
     this->pushBack( p2, inversed, notTransposed );
     if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

     /*
      * Building the block
      * / Mu/dt 0 \
      * \   0   I /
      */
      if ( verbose ) std::cout << "       Fluid block... ";
      timer.start();
      boost::shared_ptr<matrixBlock_Type> P3( new matrixBlock_Type( map ) );
      P3->setBlockStructure( blockNumRows, blockNumColumns );
      P3->getMatrixBlockView( 0, 0, B11 );
      P3->getMatrixBlockView( 1, 1, B22 );
      MatrixBlockUtils::copyBlock( F, B11 );
      MatrixBlockUtils::createIdentityBlock( B22 );
      P3->globalAssemble();

      boost::shared_ptr<matrix_Type> p3 = P3;
      superPtr_Type precForBlock3( PRECFactory::instance().createObject( M_fluidPrec ) );
      precForBlock3->setDataFromGetPot( M_dataFile, M_fluidPrecDataSection );
      this->pushBack( p3,precForBlock3, notInversed, notTransposed );
      if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

    /*
     * Building the block
     * / I    0   \
     * \ 0 -Mp/nu /
     */
    if ( verbose ) std::cout << "       Pressure mass block...";
    timer.start();
    boost::shared_ptr<matrixBlock_Type> PMp( new matrixBlock_Type( map ) );
    *PMp *= 0.0;
    PMp->setBlockStructure( blockNumRows, blockNumColumns );
    PMp->getMatrixBlockView( 0, 0, B11 );
    PMp->getMatrixBlockView( 1, 1, B22 );
    MatrixBlockUtils::createIdentityBlock( B11 );
    M_adrPressureAssembler.addMass( PMp, -M_density/M_viscosity, B22.firstRowIndex(), B22.firstColumnIndex() );
    PMp->globalAssemble();
    boost::shared_ptr<matrix_Type> pMp = PMp;
    superPtr_Type precForBlock2( PRECFactory::instance().createObject( M_pressureMassPrec ) );
    precForBlock2->setDataFromGetPot( M_dataFile, M_pressureMassPrecDataSection );
    this->pushBack( pMp,precForBlock2, notInversed, notTransposed );
    if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

    /*
     * Building the block (the block is inversed)
     * / I -Bt \
     * \ 0  I  /
     */
    if ( verbose ) std::cout << "       Gradient block... ";
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2( new matrixBlock_Type( map ) );
    P2->setBlockStructure( blockNumRows,blockNumColumns );
    P2->getMatrixBlockView( 0, 0, B11 );
    P2->getMatrixBlockView( 0, 1, B12 );
    P2->getMatrixBlockView( 1, 1, B22 );
    MatrixBlockUtils::copyBlock( Bt, B12 );
    ( *P2 ) *= -1;
    MatrixBlockUtils::createIdentityBlock( B11 );
    MatrixBlockUtils::createIdentityBlock( B22 );
    P2->globalAssemble();

    boost::shared_ptr<matrix_Type> p2 = P2;
    this->pushBack( p2, inversed, notTransposed );
    if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

    /*
     * Building the block
     * / dt*Mu^-1 0 \
     * \   0      I /
     */
     if ( verbose ) std::cout << "       Fluid block... ";
     timer.start();
     boost::shared_ptr<matrixBlock_Type> P3( new matrixBlock_Type( map ) );
     P3->setBlockStructure( blockNumRows, blockNumColumns );
     P3->getMatrixBlockView( 0, 0, B11 );
     P3->getMatrixBlockView( 1, 1, B22 );
     MatrixBlockUtils::copyBlock( F, B11 );
     MatrixBlockUtils::createIdentityBlock( B22 );
     P3->globalAssemble();

     boost::shared_ptr<matrix_Type> p3 = P3;
     superPtr_Type precForBlock3( PRECFactory::instance().createObject( M_fluidPrec ) );
     precForBlock3->setDataFromGetPot( M_dataFile, M_fluidPrecDataSection );
     this->pushBack( p3,precForBlock3, notInversed, notTransposed );
     if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;


    this->M_preconditionerCreated = true;

    if ( verbose ) std::cout << "      >All the blocks are built" << std::endl
                             << "      >";
    return ( EXIT_SUCCESS );
}

int
PreconditionerPressureCorrection::numBlocksRows() const
{
    return 2;
}

int
PreconditionerPressureCorrection::numBlocksColumns() const
{
    return 2;
}

void
PreconditionerPressureCorrection::setDataFromGetPot( const GetPot& dataFile,
                                                     const std::string& section )
{
    M_dataFile   = dataFile;
    createPressureCorrectionList( M_list, dataFile, section, "PressureCorrection" );

    M_precType                    = this->M_list.get( "prectype", "PressureCorrection" );

    M_fluidPrec                   = this->M_list.get( "subprecs: fluid prec", "ML" );
    M_fluidPrecDataSection        = this->M_list.get( "subprecs: fluid prec data section", "" );

    M_pressureMassPrec            = this->M_list.get( "subprecs: pressure mass prec", "ML" );
    M_pressureMassPrecDataSection = this->M_list.get( "subprecs: pressure mass prec data section", "" );
}

void
PreconditionerPressureCorrection::setFESpace( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace )
{
    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_adrVelocityAssembler.setup( pFESpace, uFESpace ); // p,beta=u

    // We setup the size of the blocks
    M_velocityBlockSize = M_uFESpace->fieldDim() * M_uFESpace->dof().numTotalDof();
    M_pressureBlockSize = M_pFESpace->dof().numTotalDof();
}

void
PreconditionerPressureCorrection::setTimestep( const Real& timestep )
{
    M_timestep = timestep;
}

void
PreconditionerPressureCorrection::setViscosity( const Real& viscosity )
{
    M_viscosity = viscosity;
}

void
PreconditionerPressureCorrection::setDensity( const Real& density )
{
    M_density = density;
}

} // namespace LifeV
