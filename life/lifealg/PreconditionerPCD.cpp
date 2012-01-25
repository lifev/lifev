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
    @brief PreconditionerPCD

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 29-11-2010
 */

#include <vector>
#include "PreconditionerPCD.hpp"
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>
#include <life/lifealg/PreconditionerML2.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifefem/BCManage.hpp>
#include <life/lifearray/MatrixEpetraStructured.hpp>
#include <life/lifearray/MatrixEpetraStructuredView.hpp>
#include <life/lifearray/MatrixEpetraStructuredUtility.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <EpetraExt_MatrixMatrix.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV {

PreconditionerPCD::PreconditionerPCD( const  boost::shared_ptr<Epetra_Comm>& comm ):
    PreconditionerComposition    ( comm ),
    M_velocityBlockSize          ( -1 ),
    M_pressureBlockSize          ( -1 ),
    M_timestep                   ( 1.0 ),
    M_viscosity                  ( 1.0 ),
    M_density                    ( 1.0 ),
    M_pressureBoundaryConditions ( "none" ),
    M_pressureLaplacianOperator  ( "standard" ),
    M_useLumpedPressureMass      ( false ),
    M_setApBoundaryConditions    ( false ),
    M_setFpBoundaryConditions    ( false ),
    M_fullFactorization          ( false )
{
    M_uFESpace.reset();
    M_pFESpace.reset();
    M_beta.reset();
}

PreconditionerPCD::~PreconditionerPCD()
{

}

void
PreconditionerPCD::createParametersList( list_Type&         list,
                                         const GetPot&      dataFile,
                                         const std::string& section,
                                         const std::string& subSection )
{
    createPCDList( list, dataFile, section, subSection, M_comm->MyPID() == 0 );
}

void
PreconditionerPCD::createPCDList( list_Type&         list,
                                  const GetPot&      dataFile,
                                  const std::string& section,
                                  const std::string& subsection,
                                  const bool&        verbose )
{
    bool displayList = dataFile( ( section + "/displayList" ).data(), false );

    std::string precType = dataFile( ( section + "/prectype" ).data(),"PCD" );
    list.set( "prectype", precType );

    std::string fluidPrec = dataFile( ( section + "/" + subsection + "/subprecs/fluid_prec" ).data(), "ML" );
    list.set( "subprecs: fluid prec", fluidPrec );
    std::string fluidPrecDataSection = dataFile( ( section + "/" + subsection + "/subprecs/fluid_prec_data_section" ).data(), "" );
    list.set( "subprecs: fluid prec data section", ( fluidPrecDataSection ).data() );

    std::string pressureLaplacianPrec = dataFile( ( section + "/" + subsection + "/subprecs/pressure_laplacian_prec" ).data(), "ML" );
    list.set( "subprecs: pressure laplacian prec", pressureLaplacianPrec );
    std::string pressureLaplacianPrecDataSection = dataFile( ( section + "/" + subsection + "/subprecs/pressure_laplacian_prec_data_section" ).data(), "" );
    list.set( "subprecs: pressure laplacian prec data section", ( pressureLaplacianPrecDataSection ).data() );

    std::string pressureMassPrec = dataFile( ( section + "/" + subsection + "/subprecs/pressure_mass_prec" ).data(), "ML" );
    list.set( "subprecs: pressure mass prec", pressureMassPrec );
    std::string pressureMassPrecDataSection = dataFile( ( section + "/" + subsection + "/subprecs/pressure_mass_prec_data_section" ).data(), "" );
    list.set( "subprecs: pressure mass prec data section", ( pressureMassPrecDataSection ).data() );

    std::string pressureBoundaryConditions = dataFile( ( section + "/" + subsection + "/pressure_boundary_conditions").data(), "none" );
    list.set( "pressure boundary conditions", pressureBoundaryConditions );

    std::string pressureLaplacianOperator = dataFile( ( section + "/" + subsection + "/pressure_laplacian_operator").data(), "standard" );
    list.set( "pressure laplacian operator", pressureLaplacianOperator );

    bool useLumpedPressureMass = dataFile( ( section + "/" + subsection + "/use_lumped_pressure_mass" ).data(), false );
    list.set( "use lumped pressure mass", useLumpedPressureMass );

    bool setApBoundaryConditions = dataFile( ( section + "/" + subsection + "/set_Ap_boundary_conditions" ).data(), false );
    list.set( "set Ap boundary conditions", setApBoundaryConditions );

    bool setFpBoundaryConditions = dataFile( ( section + "/" + subsection + "/set_Fp_boundary_conditions" ).data(), false );
    list.set( "set Fp boundary conditions", setFpBoundaryConditions );

    bool fullFactorization = dataFile( ( section + "/" + subsection + "/full_factorization" ).data(), false );
    list.set( "full factorization", fullFactorization );

    if ( displayList && verbose ) list.print( std::cout );
}

Real
PreconditionerPCD::condest()
{
    return 0.0;
}

void
PreconditionerPCD::updateBeta( const vector_Type& beta )
{
    *M_beta = beta;
}

int
PreconditionerPCD::buildPreconditioner( matrixPtr_type& oper )
{
    if ( ( M_uFESpace.get() == NULL ) || ( M_pFESpace.get() == NULL ) )
    {
        std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
        exit( -1 );
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
    matrixBlockView_Type F, Bt, B, C;
    //oper.blockView( 0, 0, F );
    F.setup( 0 ,0, blockNumRows[0], blockNumColumns[0], oper.get() );
    //F.showMe();

    //oper.blockView( 0, 1, Bt );
    Bt.setup( 0, blockNumColumns[0], blockNumRows[0], blockNumColumns[1], oper.get() );
    //Bt.showMe();

    //oper.blockView( 1, 0, B );
    B.setup( blockNumRows[0], 0, blockNumRows[1], blockNumColumns[0], oper.get() );
    //B.showMe();

    //oper.blockView( 1, 1, C );
    C.setup( blockNumRows[0], blockNumColumns[0], blockNumRows[1], blockNumColumns[1], oper.get() );
    //C.showMe();
    if ( verbose ) std::cout << "done in " << timer.diff() << " s." << std::endl;

    // Getting the block structure of B
    matrixBlockView_Type B11, B12, B21, B22, B22base;

    /*
     * PCD:
     * / F Bt \   / I  0 \ / I Bt \ / F 0 \
     * \ 0 -S / = \ 0 -S / \ 0 I  / \ 0 I /
     *
     * PCD^-1:
     * / F  Bt \^-1   / F^-1 0 \ / I -Bt \ / I  0    \
     * \ 0 -S  /    = \ 0    I / \ 0  I  / \ 0 -S^-1 /
     */

    boost::shared_ptr<matrix_Type> p3;
    superPtr_Type precForBlock3;
    if( M_fullFactorization == true )
    {
        /*
         * Building the block
         * / I      0 \   / F 0 \ / I 0 \ / F^-1 0 \
         * \ BF^-1  I / = \ 0 I / \ B I / \ 0    I /
         */

        /*
         * Building the block
         * / F 0 \
         * \ 0 I /
         */
        if ( verbose ) std::cout << " Fluid block... ";
        timer.start();
        boost::shared_ptr<matrixBlock_Type> P3( new matrixBlock_Type( map ) );
        P3->setBlockStructure( blockNumRows, blockNumColumns );
        P3->blockView( 0, 0, B11 );
        P3->blockView( 1, 1, B22 );
        MatrixEpetraStructuredUtility::copyBlock( F, B11 );
        MatrixEpetraStructuredUtility::createIdentityBlock( B22 );
        P3->globalAssemble();
        //P3->spy( "p3" );
        p3 = P3;
        precForBlock3.reset( PRECFactory::instance().createObject( M_fluidPrec ) );
        precForBlock3->setDataFromGetPot( M_dataFile, M_fluidPrecDataSection );
        if( M_fluidPrec == "ML2" )
        {
            PreconditionerML2* tmpPrecPtr = dynamic_cast<PreconditionerML2*>( precForBlock3.get() );
            tmpPrecPtr->setFESpace( M_uFESpace, M_pFESpace );
        }
        this->pushBack( p3,precForBlock3, notInversed, notTransposed );
        if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

        /*
         * Building the block (the block is inversed)
         * /  I  0 \
         * \ -B  I /
         */
        if ( verbose ) std::cout << " Divergence block... ";
        timer.start();
        boost::shared_ptr<matrixBlock_Type> P2e( new matrixBlock_Type( map ) );
        P2e->setBlockStructure( blockNumRows, blockNumColumns );
        P2e->blockView( 0, 0, B11 );
        P2e->blockView( 1, 0, B21 );
        P2e->blockView( 1, 1, B22 );
        MatrixEpetraStructuredUtility::copyBlock( B, B21 );
        ( *P2e ) *= -1;
        MatrixEpetraStructuredUtility::createIdentityBlock( B11 );
        MatrixEpetraStructuredUtility::createIdentityBlock( B22 );
        P2e->globalAssemble();
        //P2->spy( "p2" );
        boost::shared_ptr<matrix_Type> p2e = P2e;
        this->pushBack( p2e, inversed, notTransposed );
        if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

        /*
         * Building the block
         * / F^-1 0 \
         * \ 0    I /
         */
        if ( verbose ) std::cout << " Fluid block... ";
        this->pushBack( p3, inversed, notTransposed );
        if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * Building the block
     * / I  0 \   / I     0      \   / I  0  \ / I 0     \ / I 0  \
     * \ 0 -S / = \ 0 -ApFp^-1Mp / = \ 0 -Ap / \ 0 Fp^-1 / \ 0 Mp /
     */
    if ( verbose ) std::cout << " Building Fp... ";
    timer.start();
    boost::shared_ptr<matrixBlock_Type> PFp( new matrixBlock_Type( map ) );
    *PFp *= 0.0;
    PFp->setBlockStructure( blockNumRows, blockNumColumns );
    PFp->blockView( 1, 1, B22 );
    M_adrPressureAssembler.addDiffusion( PFp, -M_viscosity/M_density, B22.firstRowIndex(), B22.firstColumnIndex() );
    M_adrPressureAssembler.addAdvection( PFp, *M_beta, B22.firstRowIndex(), B22.firstColumnIndex() );
    M_adrPressureAssembler.addMass( PFp, 1.0/M_timestep, B22.firstRowIndex(), B22.firstColumnIndex() );
    PFp->globalAssemble();
    boost::shared_ptr<matrix_Type> pFp = PFp;
    if ( verbose ) std::cout << "done in " << timer.diff() << " s." << std::endl;

    if ( verbose ) std::cout << " Building Ap";
    timer.start();
    boost::shared_ptr<matrixBlock_Type> PAp( new matrixBlock_Type( map ) );
    *PAp *= 0.0;
    PAp->setBlockStructure( blockNumRows, blockNumColumns );
    PAp->blockView( 0, 0, B11 );
    PAp->blockView( 1, 1, B22 );

    if ( M_pressureLaplacianOperator == "symmetric" )
    {
        if ( verbose ) std::cout << " (Symm(Fp) version)... ";
        Epetra_CrsMatrix* tmpCrsMatrix( NULL );
        tmpCrsMatrix = PAp->matrixPtr().get();
        EpetraExt::MatrixMatrix::Add( *( pFp->matrixPtr() ), false, 0.5*M_density/M_viscosity,
                                      *( pFp->matrixPtr() ), true,  0.5*M_density/M_viscosity,
                                      tmpCrsMatrix );

        MatrixEpetraStructuredUtility::createScalarBlock( B11, 1.0 );
    }
    else
    {
        if ( verbose ) std::cout << "... ";
        M_adrPressureAssembler.addDiffusion( PAp, 1.0, B22.firstRowIndex(), B22.firstColumnIndex() );
        MatrixEpetraStructuredUtility::createIdentityBlock( B11 );
    }
    PAp->globalAssemble();
    boost::shared_ptr<matrix_Type> pAp = PAp;
    if ( verbose ) std::cout << "done in " << timer.diff() << " s." << std::endl;


    if ( verbose ) std::cout << " Building Mp";
    timer.start();
    boost::shared_ptr<matrixBlock_Type> PMp( new matrixBlock_Type( map ) );
    *PMp *= 0.0;
    PMp->setBlockStructure( blockNumRows, blockNumColumns );
    PMp->blockView( 0, 0, B11 );
    PMp->blockView( 1, 1, B22 );
    MatrixEpetraStructuredUtility::createIdentityBlock( B11 );
    if ( M_useLumpedPressureMass )
    {
        if ( verbose ) std::cout << " (Lumped version)... ";
        boost::shared_ptr<matrixBlock_Type> tmpMass( new matrixBlock_Type( map ) );
        M_adrPressureAssembler.addMass( tmpMass, -1.0, B22.firstRowIndex(), B22.firstColumnIndex() );
        tmpMass->globalAssemble();
        tmpMass->setBlockStructure( blockNumRows, blockNumColumns );
        tmpMass->blockView( 1, 1, B11 );
        MatrixEpetraStructuredUtility::createInvLumpedBlock( B11, B22 );
        tmpMass.reset();

    }
    else
    {
        if ( verbose ) std::cout << "... ";
        M_adrPressureAssembler.addMass( PMp, -1.0, B22.firstRowIndex(), B22.firstColumnIndex() );
    }
    PMp->globalAssemble();
    boost::shared_ptr<matrix_Type> pMp = PMp;
    if ( verbose ) std::cout << "done in " << timer.diff() << " s." << std::endl;

    if ( M_pressureBoundaryConditions == "first_dof_dirichlet" )
    {
        UInt firstIndex = M_pFESpace->map().map( Unique )->MaxMyGID() + B22.firstRowIndex();
        pAp->diagonalize( firstIndex, 1.0 );
        pFp->diagonalize( firstIndex, 1.0 );
        //pMp->diagonalize( firstIndex, 1.0 );
    }
    else if ( M_pressureBoundaryConditions == "dirichlet_to_dirichlet" )
    {
        // Loop on boundary conditions
        for ( ID i = 0; i < M_bcHandlerPtr->size(); ++i )
        {
            if ( M_bcHandlerPtr->operator[]( i ).type() == Essential )
            {
                for ( ID j = 0; j < M_bcHandlerPtr->operator[]( i ).list_size(); ++j )
                {
                    UInt myId = M_bcHandlerPtr->operator[]( i )[j]->id() + B22.firstRowIndex();
                    if ( M_setApBoundaryConditions ) pAp->diagonalize( myId, 1.0);
                    if ( M_setFpBoundaryConditions ) pFp->diagonalize( myId, 1.0);
                    //pMp->diagonalize( myId, 1.0 );
                }
            }
        }
    }
    else if ( M_pressureBoundaryConditions == "neumann_to_dirichlet" )
    {
        // Loop on boundary conditions
        for ( ID i = 0; i < M_bcHandlerPtr->size(); ++i )
        {
            if( M_bcHandlerPtr->operator[]( i ).type() == Natural )
            {
                for ( ID j = 0; j < M_bcHandlerPtr->operator[]( i ).list_size(); ++j )
                {
                    UInt myId = M_bcHandlerPtr->operator[]( i )[j]->id() + B22.firstRowIndex();
                    if ( M_setApBoundaryConditions ) pAp->diagonalize( myId, 1.0 );
                    if ( M_setFpBoundaryConditions ) pFp->diagonalize( myId, 1.0 );
                    //pMp->diagonalize( myId, 1.0 );
                }
            }
        }
    }
    else if ( M_pressureBoundaryConditions == "dirichlet_at_inflow" )
    {
        // Loop on boundary conditions
        for ( ID i = 0; i < M_bcHandlerPtr->size(); ++i )
        {
            if ( M_bcHandlerPtr->operator[]( i ).flag() == 1 )
            {
                for ( ID j = 0; j < M_bcHandlerPtr->operator[]( i ).list_size(); ++j )
                {
                    UInt myId = M_bcHandlerPtr->operator[]( i )[j]->id() + B22.firstRowIndex();
                    if ( M_setApBoundaryConditions ) pAp->diagonalize( myId, 1.0 );
                    if ( M_setFpBoundaryConditions ) pFp->diagonalize( myId, 1.0 );
                }
            }
        }
    }
    else if ( M_pressureBoundaryConditions == "robin_at_inflow" )
    {
        BCHandler bcHandler;

        // Offset to set the BC for the pressure
        bcHandler.setOffset( B22.firstRowIndex() );

        // Loop on boundary conditions
        for ( ID i = 0; i < M_bcHandlerPtr->size(); ++i )
        {
            if ( verbose ) std::cout << "BCHandler name: " << M_bcHandlerPtr->operator[]( i ).name() << std::endl;

            UInt numComponents = M_bcHandlerPtr->operator[]( i ).numberOfComponents();
            numComponents = 1;

            if ( M_bcHandlerPtr->operator[]( i ).flag() == 1 )
            {
                vector_Type convVelocity( *M_beta, Repeated );
                vector_Type robinRHS( M_uFESpace->map(), Repeated );
                /*
                BCVector uRobin( robinRHS, M_uFESpace->dof().numTotalDof(), 0 );
                uRobin.setRobinCoeffVector( convVelocity );
                uRobin.setBetaCoeff( M_viscosity/M_density );

                bcHandler.addBC( M_bcHandlerPtr->operator[]( i ).name(),
                                 M_bcHandlerPtr->operator[]( i ).flag(),
                                 Robin,
                                 Normal,
                                 uRobin,
                                 numComponents );
                */
            }
            else
            {
                bcHandler.addBC( M_bcHandlerPtr->operator[]( i ).name(),
                                 M_bcHandlerPtr->operator[]( i ).flag(),
                                 M_bcHandlerPtr->operator[]( i ).type(),
                                 M_bcHandlerPtr->operator[]( i ).mode(),
                                 const_cast<BCFunctionBase&>( *( M_bcHandlerPtr->operator[]( i ).pointerToFunctor() ) ),
                                 numComponents );
            }
        }
        if ( M_setApBoundaryConditions ) bcManageMatrix( *pAp, *M_uFESpace->mesh(), M_uFESpace->dof(), bcHandler, M_uFESpace->feBd(), 1.0, 0.0 );
        if ( M_setFpBoundaryConditions ) bcManageMatrix( *pFp, *M_uFESpace->mesh(), M_uFESpace->dof(), bcHandler, M_uFESpace->feBd(), 1.0, 0.0 );

    }
    if ( verbose ) std::cout << " Pressure BC type = " << M_pressureBoundaryConditions << std::endl;
    if ( ( verbose )&&( M_setApBoundaryConditions ) ) std::cout << " BC imposed on Ap" << std::endl;
    if ( ( verbose )&&( M_setFpBoundaryConditions ) ) std::cout << " BC imposed on Fp" << std::endl;

    if ( verbose ) std::cout << " Schur block (a)... ";
    timer.start();
    //pAp->spy( "p1a" );
    superPtr_Type precForBlock1( PRECFactory::instance().createObject( M_pressureLaplacianPrec ) );
    precForBlock1->setDataFromGetPot( M_dataFile, M_pressureLaplacianPrecDataSection );
    this->pushBack( pAp, precForBlock1, notInversed, notTransposed );
    if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

    if ( verbose ) std::cout << " Schur block (b)... ";
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P1b( new matrixBlock_Type( map ) );
    P1b->setBlockStructure( blockNumRows, blockNumColumns );
    *P1b += *PFp;
    P1b->blockView( 0, 0, B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock( B11 );
    P1b->globalAssemble();
    boost::shared_ptr<matrix_Type> p1b = P1b;
    //p1b->spy( "p1b" );
    this->pushBack( p1b, inversed,notTransposed );
    if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

    if ( verbose ) std::cout << " Schur block (c)... ";
    timer.start();
    //pMp->spy( "p1c" );
    if ( M_useLumpedPressureMass )
    {
        this->pushBack( pMp, inversed, notTransposed );
    }
    else
    {
        superPtr_Type precForBlock2( PRECFactory::instance().createObject( M_pressureMassPrec ) );
        precForBlock2->setDataFromGetPot( M_dataFile, M_pressureMassPrecDataSection );
        this->pushBack( pMp,precForBlock2, notInversed, notTransposed );
    }
    if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

    /*
     * Building the block (the block is inversed)
     * / I -Bt \
     * \ 0  I  /
     */
    if ( verbose ) std::cout << " Gradient block... ";
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2( new matrixBlock_Type( map ) );
    P2->setBlockStructure( blockNumRows, blockNumColumns );
    P2->blockView( 0, 0, B11 );
    P2->blockView( 0, 1, B12 );
    P2->blockView( 1, 1, B22 );
    MatrixEpetraStructuredUtility::copyBlock( Bt, B12 );
    ( *P2 ) *= -1;
    MatrixEpetraStructuredUtility::createIdentityBlock( B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock( B22 );
    P2->globalAssemble();
    //P2->spy( "p2" );
    boost::shared_ptr<matrix_Type> p2 = P2;
    this->pushBack( p2, inversed, notTransposed );
    if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;

    /*
     * Building the block
     * / F 0 \
     * \ 0 I /
     */
    if( M_fullFactorization == true )
    {
        precForBlock3.reset( PRECFactory::instance().createObject( M_fluidPrec ) );
        precForBlock3->setDataFromGetPot( M_dataFile, M_fluidPrecDataSection );
        this->pushBack( p3,precForBlock3, notInversed, notTransposed );
        if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;
    }
    else
    {
        if ( verbose ) std::cout << " Fluid block... ";
        timer.start();
        boost::shared_ptr<matrixBlock_Type> P3( new matrixBlock_Type( map ) );
        P3->setBlockStructure( blockNumRows, blockNumColumns );
        P3->blockView( 0, 0, B11 );
        P3->blockView( 1, 1, B22 );
        MatrixEpetraStructuredUtility::copyBlock( F, B11 );
        MatrixEpetraStructuredUtility::createIdentityBlock( B22 );
        P3->globalAssemble();
        //P3->spy( "p3" );
        p3 = P3;
        precForBlock3.reset( PRECFactory::instance().createObject( M_fluidPrec ) );
        precForBlock3->setDataFromGetPot( M_dataFile, M_fluidPrecDataSection );
        if( M_fluidPrec == "ML2" )
        {
            PreconditionerML2* tmpPrecPtr = dynamic_cast<PreconditionerML2*>( precForBlock3.get() );
            tmpPrecPtr->setFESpace( M_uFESpace, M_pFESpace );
        }
        this->pushBack( p3,precForBlock3, notInversed, notTransposed );
        if ( verbose ) std::cout << " done in " << timer.diff() << " s." << std::endl;
    }

    this->M_preconditionerCreated = true;

    if ( verbose ) std::cout << "      >All the blocks are built" << std::endl
                          << "      >";
    return ( EXIT_SUCCESS );
}

int
PreconditionerPCD::numBlocksRows() const
{
    return 2;
}

int
PreconditionerPCD::numBlocksColumns() const
{
    return 2;
}

void PreconditionerPCD::setDataFromGetPot( const GetPot& dataFile,
                                           const std::string& section )
{
    M_dataFile   = dataFile;
    createPCDList( M_list, dataFile, section, "PCD" );

    M_precType                         = this->M_list.get( "prectype", "PCD" );

    M_fluidPrec                        = this->M_list.get( "subprecs: fluid prec","ML" );
    M_fluidPrecDataSection             = this->M_list.get( "subprecs: fluid prec data section", "" );

    M_pressureLaplacianPrec            = this->M_list.get( "subprecs: pressure laplacian prec", "ML");
    M_pressureLaplacianPrecDataSection = this->M_list.get( "subprecs: pressure laplacian prec data section", "" );

    M_pressureMassPrec                 = this->M_list.get( "subprecs: pressure mass prec", "ML" );
    M_pressureMassPrecDataSection      = this->M_list.get( "subprecs: pressure mass prec data section", "" );

    M_pressureBoundaryConditions       = this->M_list.get( "pressure boundary conditions", "none" );
    M_pressureLaplacianOperator        = this->M_list.get( "pressure laplacian operator", "standard" );

    M_useLumpedPressureMass            = this->M_list.get( "use lumped pressure mass", false );
    M_setApBoundaryConditions          = this->M_list.get( "set Ap boundary conditions", false );
    M_setFpBoundaryConditions          = this->M_list.get( "set Fp boundary conditions", false );
    M_fullFactorization                = this->M_list.get( "full factorization", false );
}

void
PreconditionerPCD::setFESpace( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace )
{
    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_adrPressureAssembler.setup( pFESpace, uFESpace ); // p,beta=u
    M_adrVelocityAssembler.setup( uFESpace, uFESpace ); // u,beta=u
    M_beta.reset( new vector_Type( M_uFESpace->map() + M_pFESpace->map() ) );
    *M_beta *= 0;

    // We setup the size of the blocks
    M_velocityBlockSize = M_uFESpace->fieldDim() * M_uFESpace->dof().numTotalDof();
    M_pressureBlockSize = M_pFESpace->dof().numTotalDof();
}

void
PreconditionerPCD::setBCHandler( BCHandlerPtr_Type bchPtr )
{
    M_bcHandlerPtr = bchPtr;
}

void
PreconditionerPCD::setTimestep( const Real& timestep )
{
    M_timestep = timestep;
}

void
PreconditionerPCD::setViscosity( const Real& viscosity )
{
    M_viscosity = viscosity;
}

void
PreconditionerPCD::setDensity( const Real& density )
{
    M_density = density;
}

} // namespace LifeV
