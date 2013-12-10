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
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
// #include <lifev/core/algorithm/PreconditionerML2.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/BCBase.hpp>
#include <lifev/core/fem/BCIdentifier.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>
#include <lifev/core/array/VectorBlockStructure.hpp>


// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <EpetraExt_MatrixMatrix.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{

PreconditionerPCD::PreconditionerPCD ( boost::shared_ptr<Epetra_Comm> comm ) :
    PreconditionerComposition    ( comm ),
    M_velocityBlockSize          ( -1 ),
    M_pressureBlockSize          ( -1 ),
    M_timestep                   ( 1.0 ),
    M_viscosity                  ( 1.0 ),
    M_density                    ( 1.0 ),
    M_pressureLaplacianOperator  ( "standard" ),
    M_pressureMassOperator       ( "lumped" ),
    M_setApBoundaryConditions    ( false ),
    M_setFpBoundaryConditions    ( false ),
    M_setMpBoundaryConditions    ( false ),
    M_fullFactorization          ( false ),
    M_useStiffStrain             ( false ),
    M_enableTransient            ( true ),
    M_divergenceCoeff            ( 1.0 ),
    M_recomputeNormalVectors     ( true ),
    M_schurOperatorReverseOrder  ( false ),
    M_inflowBoundaryType         ( "Robin" ),
    M_outflowBoundaryType        ( "Neumann" ),
    M_characteristicBoundaryType ( "Neumann" )
{
    M_uFESpace.reset();
    M_pFESpace.reset();
    M_beta.reset();
}

PreconditionerPCD::~PreconditionerPCD()
{

}

void
PreconditionerPCD::createParametersList ( list_Type&         list,
                                          const GetPot&      dataFile,
                                          const std::string& section,
                                          const std::string& subsection )
{
    const bool verbose ( M_comm->MyPID() == 0 );

    bool displayList = dataFile ( ( section + "/displayList" ).data(), false );

    std::string precType = dataFile ( ( section + "/prectype" ).data(), "PCD" );
    list.set ( "prectype", precType );

    std::string fluidPrec = dataFile ( ( section + "/" + subsection + "/subprecs/fluid_prec" ).data(), "ML" );
    list.set ( "subprecs: fluid prec", fluidPrec );
    std::string fluidPrecDataSection = dataFile ( ( section + "/" + subsection + "/subprecs/fluid_prec_data_section" ).data(), "" );
    list.set ( "subprecs: fluid prec data section", ( fluidPrecDataSection ).data() );

    std::string pressureLaplacianPrec = dataFile ( ( section + "/" + subsection + "/subprecs/pressure_laplacian_prec" ).data(), "ML" );
    list.set ( "subprecs: pressure laplacian prec", pressureLaplacianPrec );
    std::string pressureLaplacianPrecDataSection = dataFile ( ( section + "/" + subsection + "/subprecs/pressure_laplacian_prec_data_section" ).data(), "" );
    list.set ( "subprecs: pressure laplacian prec data section", ( pressureLaplacianPrecDataSection ).data() );

    std::string pressureMassPrec = dataFile ( ( section + "/" + subsection + "/subprecs/pressure_mass_prec" ).data(), "ML" );
    list.set ( "subprecs: pressure mass prec", pressureMassPrec );
    std::string pressureMassPrecDataSection = dataFile ( ( section + "/" + subsection + "/subprecs/pressure_mass_prec_data_section" ).data(), "" );
    list.set ( "subprecs: pressure mass prec data section", ( pressureMassPrecDataSection ).data() );

    std::string pressureBoundaryConditions = dataFile ( ( section + "/" + subsection + "/pressure_boundary_conditions").data(), "none" );
    list.set ( "pressure boundary conditions", pressureBoundaryConditions );

    std::string pressureLaplacianOperator = dataFile ( ( section + "/" + subsection + "/pressure_laplacian_operator").data(), "standard" );
    list.set ( "pressure laplacian operator", pressureLaplacianOperator );

    std::string pressureMassOperator = dataFile ( ( section + "/" + subsection + "/pressure_mass_operator").data(), "lumped" );
    list.set ( "pressure mass operator", pressureMassOperator );

    bool setApBoundaryConditions = dataFile ( ( section + "/" + subsection + "/set_Ap_boundary_conditions" ).data(), false );
    list.set ( "set Ap boundary conditions", setApBoundaryConditions );

    bool setFpBoundaryConditions = dataFile ( ( section + "/" + subsection + "/set_Fp_boundary_conditions" ).data(), false );
    list.set ( "set Fp boundary conditions", setFpBoundaryConditions );

    bool setMpBoundaryConditions = dataFile ( ( section + "/" + subsection + "/set_Mp_boundary_conditions" ).data(), false );
    list.set ( "set Mp boundary conditions", setMpBoundaryConditions );

    bool fullFactorization = dataFile ( ( section + "/" + subsection + "/full_factorization" ).data(), false );
    list.set ( "full factorization", fullFactorization );

    bool useStiffStrain = dataFile ( ( section + "/" + subsection + "/use_StiffStrain" ).data(), false );
    list.set ( "use stiff strain", useStiffStrain );

    bool enableTransient = dataFile ( ( section + "/" + subsection + "/enable_transient" ).data(), true );
    list.set ( "enable transient", enableTransient );

    bool useMinusDivergence = dataFile ( ( section + "/" + subsection + "/use_minus_divergence" ).data(), false );
    list.set ( "use minus divergence", useMinusDivergence );

    bool recomputeNormalVectors = dataFile ( ( section + "/" + subsection + "/recompute_normal_vectors" ).data(), true );
    list.set ( "recompute normal vectors", recomputeNormalVectors );

    bool schurOperatorReverseOrder = dataFile ( ( section + "/" + subsection + "/Schur_operator_reverse_order" ).data(), false );
    list.set ( "Schur operator reverse order", schurOperatorReverseOrder );

    std::string inflowBoundaryType = dataFile ( ( section + "/" + subsection + "/inflow_boundary_type" ).data(), "Robin" );
    list.set ( "inflow boundary type", inflowBoundaryType );

    std::string outflowBoundaryType = dataFile ( ( section + "/" + subsection + "/outflow_boundary_type" ).data(), "Neumann" );
    list.set ( "outflow boundary type", outflowBoundaryType );

    std::string characteristicBoundaryType = dataFile ( ( section + "/" + subsection + "/characteristic_boundary_type" ).data(), "Neumann" );
    list.set ( "characteristic boundary type", characteristicBoundaryType );

    if ( displayList && verbose )
    {
        std::cout << "PCD parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        list.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

Real
PreconditionerPCD::condest()
{
    return 0.0;
}

void
PreconditionerPCD::updateBeta ( const vector_Type& beta )
{
    *M_beta = beta;
}

int
PreconditionerPCD::buildPreconditioner ( matrixPtr_type& oper )
{
    if ( ( M_uFESpace.get() == NULL ) || ( M_pFESpace.get() == NULL ) )
    {
        std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
        exit ( -1 );
    }

    const bool verbose ( M_comm->MyPID() == 0 );

    // Make sure that the preconditioner is reset
    this->resetPreconditioner();

    std::vector<UInt> blockNumRows ( 2, 0 );
    blockNumRows[0] = M_velocityBlockSize;
    blockNumRows[1] = M_pressureBlockSize;
    std::vector<UInt> blockNumColumns ( blockNumRows );
    VectorBlockStructure vectorStructure;
    vectorStructure.setBlockStructure ( blockNumColumns );

    const bool inversed ( true );
    const bool notInversed ( false );
    const bool notTransposed ( false );

    map_Type map ( oper->map() );
    map_Type velocityMap ( M_uFESpace->map() );
    map_Type pressureMap ( M_pFESpace->map() );

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
    F.setup ( 0 , 0, blockNumRows[0], blockNumColumns[0], oper.get() );

    //oper.blockView( 0, 1, Bt );
    Bt.setup ( 0, blockNumColumns[0], blockNumRows[0], blockNumColumns[1], oper.get() );

    //oper.blockView( 1, 0, B );
    B.setup ( blockNumRows[0], 0, blockNumRows[1], blockNumColumns[0], oper.get() );

    //oper.blockView( 1, 1, C );
    C.setup ( blockNumRows[0], blockNumColumns[0], blockNumRows[1], blockNumColumns[1], oper.get() );
    if ( verbose )
    {
        std::cout << "done in " << timer.diff() << " s." << std::endl;
    }

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
    boost::shared_ptr<matrixBlock_Type> P3;
    if ( M_fluidPrec == "LinearSolver" )
    {
        MatrixEpetraStructuredUtility::createMatrixFromBlock ( F, P3, velocityMap, true );
    }
    else
    {
        P3.reset ( new matrixBlock_Type ( map ) );
        P3->setBlockStructure ( blockNumRows, blockNumColumns );
        P3->blockView ( 0, 0, B11 );
        P3->blockView ( 1, 1, B22 );
        MatrixEpetraStructuredUtility::copyBlock ( F, B11 );
        MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
        P3->globalAssemble();
    }
    p3 = P3;

    if ( M_fullFactorization == true )
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
        if ( verbose )
        {
            std::cout << "      >Fluid block... ";
        }
        timer.start();
        precForBlock3.reset ( PRECFactory::instance().createObject ( M_fluidPrec ) );
        precForBlock3->setDataFromGetPot ( M_dataFile, M_fluidPrecDataSection );
//        if ( M_fluidPrec == "ML2" )
//        {
//            PreconditionerML2* tmpPrecPtr = dynamic_cast<PreconditionerML2*> ( precForBlock3.get() );
//            tmpPrecPtr->setFESpace ( M_uFESpace, M_pFESpace );
//        }
        if ( M_fluidPrec == "LinearSolver" )
        {
            this->pushBack ( p3, precForBlock3, vectorStructure, 0, oper->map(), notInversed, notTransposed );
        }
        else
        {
            this->pushBack ( p3, precForBlock3, notInversed, notTransposed );
        }
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }

        /*
         * Building the block (the block is inversed)
         * /  I  0 \
         * \ -B  I /
         */
        if ( verbose )
        {
            std::cout << "      >Divergence block... ";
        }
        timer.start();
        boost::shared_ptr<matrixBlock_Type> P2e ( new matrixBlock_Type ( map ) );
        P2e->setBlockStructure ( blockNumRows, blockNumColumns );
        P2e->blockView ( 0, 0, B11 );
        P2e->blockView ( 1, 0, B21 );
        P2e->blockView ( 1, 1, B22 );
        MatrixEpetraStructuredUtility::copyBlock ( B, B21 );
        ( *P2e ) *= -1;
        MatrixEpetraStructuredUtility::createIdentityBlock ( B11 );
        MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
        P2e->globalAssemble();
        boost::shared_ptr<matrix_Type> p2e = P2e;
        this->pushBack ( p2e, inversed, notTransposed );
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }

        /*
         * Building the block
         * / F^-1 0 \
         * \ 0    I /
         */
        if ( verbose )
        {
            std::cout << "      >Fluid block... ";
        }
        this->pushBack ( p3, inversed, notTransposed );
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }
    }

    /*
     * Building the block
     * / I  0 \   / I     0      \   / I 0  \ / I 0     \ / I  0  \
     * \ 0 -S / = \ 0 -ApFp^-1Mp / = \ 0 Ap / \ 0 Fp^-1 / \ 0 -Mp /
     */
    UInt ApOffset ( 0 );
    UInt FpOffset ( 0 );
    UInt MpOffset ( 0 );

    if ( verbose )
    {
        std::cout << "      >Building Fp... ";
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> PFp;
    PFp.reset ( new matrixBlock_Type ( map ) );
    PFp->setBlockStructure ( blockNumRows, blockNumColumns );
    *PFp *= 0.0;
    PFp->blockView ( 0, 0, B11 );
    PFp->blockView ( 1, 1, B22 );
    MatrixEpetraStructuredUtility::createScalarBlock ( B11, 1.0 );
    if ( M_useStiffStrain )
    {
        M_adrPressureAssembler.addStiffStrain ( PFp, M_viscosity / M_density, B22.firstRowIndex(), B22.firstColumnIndex() );
    }
    else
    {
        M_adrPressureAssembler.addDiffusion ( PFp, M_viscosity / M_density, B22.firstRowIndex(), B22.firstColumnIndex() );
    }
    M_adrPressureAssembler.addAdvection ( PFp, *M_beta, B22.firstRowIndex(), B22.firstColumnIndex() );
    if ( M_enableTransient )
    {
        M_adrPressureAssembler.addMass ( PFp, 1.0 / M_timestep, B22.firstRowIndex(), B22.firstColumnIndex() );
    }
    boost::shared_ptr<matrix_Type> pFp = PFp;
    FpOffset = B22.firstRowIndex();
    if ( verbose )
    {
        std::cout << "done in " << timer.diff() << " s." << std::endl;
    }

    if ( verbose )
    {
        std::cout << "      >Building Ap";
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> PAp;
    if ( M_pressureLaplacianPrec == "LinearSolver" )
    {
        PAp.reset ( new matrixBlock_Type ( pressureMap ) );
        *PAp *= 0.0;
        PAp->blockView ( 0, 0, B22 );
        ApOffset = B22.firstRowIndex();
    }
    else
    {
        PAp.reset ( new matrixBlock_Type ( map ) );
        PAp->setBlockStructure ( blockNumRows, blockNumColumns );
        *PAp *= 0.0;
        PAp->blockView ( 0, 0, B11 );
        PAp->blockView ( 1, 1, B22 );
        MatrixEpetraStructuredUtility::createScalarBlock ( B11, 1.0 );
        ApOffset = B22.firstRowIndex();
    }

    if ( M_pressureLaplacianOperator == "symmetric" )
    {
        if ( verbose )
        {
            std::cout << " (Symm(Fp) version)... ";
        }
        Epetra_CrsMatrix* tmpCrsMatrix ( NULL );
        tmpCrsMatrix = PAp->matrixPtr().get();
        EpetraExt::MatrixMatrix::Add ( * ( pFp->matrixPtr() ), false, 0.5 * M_density / M_viscosity,
                                       * ( pFp->matrixPtr() ), true,  0.5 * M_density / M_viscosity,
                                       tmpCrsMatrix );
        *PAp *= M_divergenceCoeff;
    }
    else if ( M_pressureLaplacianOperator == "BBt" )
    {
        if ( verbose )
        {
            std::cout << " (BBt version)... ";
        }
        boost::shared_ptr<matrixBlock_Type> BMat ( new matrixBlock_Type ( map ) );
        BMat->setBlockStructure ( blockNumRows, blockNumColumns );
        MatrixEpetraStructuredUtility::copyBlock ( B, * ( BMat->block ( 1, 0 ) ) );
        BMat->globalAssemble();
        boost::shared_ptr<matrixBlock_Type> BtMat ( new matrixBlock_Type ( map ) );
        BtMat->setBlockStructure ( blockNumRows, blockNumColumns );
        MatrixEpetraStructuredUtility::copyBlock ( Bt, * ( BtMat->block ( 0, 1 ) ) );
        BtMat->globalAssemble();
        boost::shared_ptr<matrixBlock_Type> BBtMat ( new matrixBlock_Type ( map ) );
        BBtMat->setBlockStructure ( blockNumRows, blockNumColumns );
        BMat->multiply ( false,
                         *BtMat , false,
                         *BBtMat, true );
        MatrixEpetraStructuredUtility::copyBlock ( * ( BBtMat->block ( 1, 1 ) ), B22 );
    }
    else if ( M_pressureLaplacianOperator == "BinvDBt" )
    {
        if ( verbose )
        {
            std::cout << " (BinvDBt version)... ";
        }
        // Create B
        boost::shared_ptr<matrixBlock_Type> BMat ( new matrixBlock_Type ( map ) );
        BMat->setBlockStructure ( blockNumRows, blockNumColumns );
        MatrixEpetraStructuredUtility::copyBlock ( B, * ( BMat->block ( 1, 0 ) ) );
        BMat->globalAssemble();

        // Create the inverse of the diagonal mass matrix D
        boost::shared_ptr<matrixBlock_Type> tmpVelocityMass ( new matrixBlock_Type ( map ) );
        tmpVelocityMass->setBlockStructure ( blockNumRows, blockNumColumns );
        M_adrVelocityAssembler.addMass ( tmpVelocityMass, 1.0 );
        tmpVelocityMass->globalAssemble();
        boost::shared_ptr<matrixBlock_Type> invDMat ( new matrixBlock_Type ( map ) );
        invDMat->setBlockStructure ( blockNumRows, blockNumColumns );
        MatrixEpetraStructuredUtility::createInvDiagBlock ( * ( tmpVelocityMass->block ( 0, 0 ) ), * ( invDMat->block ( 0, 0 ) ) );
        invDMat->globalAssemble();
        tmpVelocityMass.reset(); // Free the memory

        // Compute BD^-1
        boost::shared_ptr<matrixBlock_Type> BinvDMat ( new matrixBlock_Type ( map ) );
        BinvDMat->setBlockStructure ( blockNumRows, blockNumColumns );
        BMat->multiply ( false,
                         *invDMat , false,
                         *BinvDMat, true );
        invDMat.reset(); // Free the memory
        BMat.reset();

        // Compute BD^-1Bt
        boost::shared_ptr<matrixBlock_Type> BtMat ( new matrixBlock_Type ( map ) );
        BtMat->setBlockStructure ( blockNumRows, blockNumColumns );
        MatrixEpetraStructuredUtility::copyBlock ( Bt, * ( BtMat->block ( 0, 1 ) ) );
        BtMat->globalAssemble();
        boost::shared_ptr<matrixBlock_Type> BBtMat ( new matrixBlock_Type ( map ) );
        BBtMat->setBlockStructure ( blockNumRows, blockNumColumns );
        BinvDMat->multiply ( false,
                             *BtMat , false,
                             *BBtMat, true );
        BinvDMat.reset(); // Free the memory
        BtMat.reset();

        // Export the matrix
        MatrixEpetraStructuredUtility::copyBlock ( * ( BBtMat->block ( 1, 1 ) ), B22 );
    }
    else
    {
        if ( verbose )
        {
            std::cout << "... ";
        }
        M_adrPressureAssembler.addDiffusion ( PAp, M_divergenceCoeff, B22.firstRowIndex(), B22.firstColumnIndex() );
    }
    boost::shared_ptr<matrix_Type> pAp = PAp;
    if ( verbose )
    {
        std::cout << "done in " << timer.diff() << " s." << std::endl;
    }


    if ( verbose )
    {
        std::cout << "      >Building Mp";
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> PMp;
    if ( M_pressureMassPrec == "LinearSolver" && M_pressureMassOperator == "standard" )
    {
        PMp.reset ( new matrixBlock_Type ( pressureMap ) );
        *PMp *= 0.0;
        PMp->blockView ( 0, 0, B22 );
        MpOffset = 0;
    }
    else
    {
        PMp.reset ( new matrixBlock_Type ( map ) );
        PMp->setBlockStructure ( blockNumRows, blockNumColumns );
        *PMp *= 0.0;
        PMp->blockView ( 0, 0, B11 );
        PMp->blockView ( 1, 1, B22 );
        MatrixEpetraStructuredUtility::createScalarBlock ( B11, 1.0 );
        MpOffset = B22.firstRowIndex();
    }
    if ( M_pressureMassOperator == "lumped" )
    {
        if ( verbose )
        {
            std::cout << " (Lumped version)... ";
        }
        boost::shared_ptr<matrixBlock_Type> tmpMass ( new matrixBlock_Type ( map ) );
        M_adrPressureAssembler.addMass ( tmpMass, -1.0, B22.firstRowIndex(), B22.firstColumnIndex() );
        tmpMass->globalAssemble();
        tmpMass->setBlockStructure ( blockNumRows, blockNumColumns );
        tmpMass->blockView ( 1, 1, B11 );
        MatrixEpetraStructuredUtility::createInvLumpedBlock ( B11, B22 );
        tmpMass.reset();
    }
    if ( M_pressureMassOperator == "diagonal" )
    {
        if ( verbose )
        {
            std::cout << " (diagonal version)... ";
        }
        boost::shared_ptr<matrixBlock_Type> tmpMass ( new matrixBlock_Type ( map ) );
        M_adrPressureAssembler.addMass ( tmpMass, -1.0, B22.firstRowIndex(), B22.firstColumnIndex() );
        tmpMass->globalAssemble();
        tmpMass->setBlockStructure ( blockNumRows, blockNumColumns );
        tmpMass->blockView ( 1, 1, B11 );
        MatrixEpetraStructuredUtility::createInvDiagBlock ( B11, B22 );
        tmpMass.reset();
    }
    else
    {
        if ( verbose )
        {
            std::cout << "... ";
        }
        M_adrPressureAssembler.addMass ( PMp, -1.0, B22.firstRowIndex(), B22.firstColumnIndex() );
    }
    boost::shared_ptr<matrix_Type> pMp = PMp;
    if ( verbose )
    {
        std::cout << "done in " << timer.diff() << " s." << std::endl;
    }

    // This option is the one used by Elman & al. in the original PCD:
    // * Dirichlet boundary conditions at inflows
    // * Neumann boundary conditions at outflows
    // * Neumann boundary conditions on characteristic boundaries

    // This option is the one used by Elman & al. in the New PCD:
    // * Robin boundary conditions at inflows
    // * Neumann boundary conditions at outflows
    // * Neumann boundary conditions on characteristic boundaries
    this->setBCByBoundaryType ( pAp, ApOffset,
                                pFp, FpOffset,
                                pMp, MpOffset );

    // For enclosed flow where only characteristics BC are imposed,
    // the choice is correct by default, i.e. Neumann BC.
    if ( ( verbose ) && ( M_setApBoundaryConditions ) )
    {
        std::cout << "      >BC imposed on Ap" << std::endl;
    }
    if ( ( verbose ) && ( M_setFpBoundaryConditions ) )
    {
        std::cout << "      >BC imposed on Fp" << std::endl;
    }
    if ( ( verbose ) && ( M_setMpBoundaryConditions ) )
    {
        std::cout << "      >BC imposed on Mp" << std::endl;
    }

    if ( M_schurOperatorReverseOrder )
    {
        if ( verbose )
        {
            std::cout << "      >Use reverse order for the Schur approximation" << std::endl;
        }

        if ( verbose )
        {
            std::cout << "      >Schur block (a)... ";
        }
        timer.start();
        if ( M_pressureMassOperator != "standard" )
        {
            this->pushBack ( pMp, inversed, notTransposed );
        }
        else
        {
            superPtr_Type precForBlock2 ( PRECFactory::instance().createObject ( M_pressureMassPrec ) );
            precForBlock2->setDataFromGetPot ( M_dataFile, M_pressureMassPrecDataSection );
            if ( M_pressureMassPrec == "LinearSolver" )
            {
                this->pushBack ( pMp, precForBlock2, vectorStructure, 1, oper->map(), notInversed, notTransposed );
            }
            else
            {
                this->pushBack ( pMp, precForBlock2, notInversed, notTransposed );
            }
        }
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }

        if ( verbose )
        {
            std::cout << "      >Schur block (b)... ";
        }
        timer.start();
        boost::shared_ptr<matrix_Type> p1b = PFp;
        this->pushBack ( p1b, inversed, notTransposed );
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }

        if ( verbose )
        {
            std::cout << "      >Schur block (c)... ";
        }
        timer.start();
        superPtr_Type precForBlock1 ( PRECFactory::instance().createObject ( M_pressureLaplacianPrec ) );
        precForBlock1->setDataFromGetPot ( M_dataFile, M_pressureLaplacianPrecDataSection );
        if ( M_pressureLaplacianPrec == "LinearSolver" )
        {
            this->pushBack ( pAp, precForBlock1, vectorStructure, 1, oper->map(), notInversed, notTransposed );
        }
        else
        {
            this->pushBack ( pAp, precForBlock1, notInversed, notTransposed );
        }
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }
    }
    else
    {
        if ( verbose )
        {
            std::cout << "      >Schur block (a)... ";
        }
        timer.start();
        superPtr_Type precForBlock1 ( PRECFactory::instance().createObject ( M_pressureLaplacianPrec ) );
        precForBlock1->setDataFromGetPot ( M_dataFile, M_pressureLaplacianPrecDataSection );
        if ( M_pressureLaplacianPrec == "LinearSolver" )
        {
            this->pushBack ( pAp, precForBlock1, vectorStructure, 1, oper->map(), notInversed, notTransposed );
        }
        else
        {
            this->pushBack ( pAp, precForBlock1, notInversed, notTransposed );
        }
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }

        if ( verbose )
        {
            std::cout << "      >Schur block (b)... ";
        }
        timer.start();
        boost::shared_ptr<matrix_Type> p1b = PFp;
        this->pushBack ( p1b, inversed, notTransposed );
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }

        if ( verbose )
        {
            std::cout << "      >Schur block (c)... ";
        }
        timer.start();
        if ( M_pressureMassOperator != "standard" )
        {
            this->pushBack ( pMp, inversed, notTransposed );
        }
        else
        {
            superPtr_Type precForBlock2 ( PRECFactory::instance().createObject ( M_pressureMassPrec ) );
            precForBlock2->setDataFromGetPot ( M_dataFile, M_pressureMassPrecDataSection );
            if ( M_pressureMassPrec == "LinearSolver" )
            {
                this->pushBack ( pMp, precForBlock2, vectorStructure, 1, oper->map(), notInversed, notTransposed );
            }
            else
            {
                this->pushBack ( pMp, precForBlock2, notInversed, notTransposed );
            }
        }
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }
    }

    /*
     * Building the block (the block is inversed)
     * / I -Bt \
     * \ 0  I  /
     */
    if ( verbose )
    {
        std::cout << "      >Gradient block... ";
    }
    timer.start();
    boost::shared_ptr<matrixBlock_Type> P2 ( new matrixBlock_Type ( map ) );
    P2->setBlockStructure ( blockNumRows, blockNumColumns );
    P2->blockView ( 0, 0, B11 );
    P2->blockView ( 0, 1, B12 );
    P2->blockView ( 1, 1, B22 );
    MatrixEpetraStructuredUtility::copyBlock ( Bt, B12 );
    ( *P2 ) *= -1;
    MatrixEpetraStructuredUtility::createIdentityBlock ( B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );
    P2->globalAssemble();
    boost::shared_ptr<matrix_Type> p2 = P2;
    this->pushBack ( p2, inversed, notTransposed );
    if ( verbose )
    {
        std::cout << " done in " << timer.diff() << " s." << std::endl;
    }

    /*
     * Building the block
     * / F 0 \
     * \ 0 I /
     */
    if ( M_fullFactorization == true )
    {
        if ( M_fluidPrec == "LinearSolver" )
        {
            this->pushBack ( p3, precForBlock3, vectorStructure, 0, oper->map(), notInversed, notTransposed );
        }
        else
        {
            this->pushBack ( p3, precForBlock3, notInversed, notTransposed );
        }
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }
    }
    else
    {
        if ( verbose )
        {
            std::cout << "      >Fluid block... ";
        }
        timer.start();
        precForBlock3.reset ( PRECFactory::instance().createObject ( M_fluidPrec ) );
        precForBlock3->setDataFromGetPot ( M_dataFile, M_fluidPrecDataSection );
//        if ( M_fluidPrec == "ML2" )
//        {
//            PreconditionerML2* tmpPrecPtr = dynamic_cast<PreconditionerML2*> ( precForBlock3.get() );
//            tmpPrecPtr->setFESpace ( M_uFESpace, M_pFESpace );
//        }
        if ( M_fluidPrec == "LinearSolver" )
        {
            this->pushBack ( p3, precForBlock3, vectorStructure, 0, oper->map(), notInversed, notTransposed );
        }
        else
        {
            this->pushBack ( p3, precForBlock3, notInversed, notTransposed );
        }
        if ( verbose )
        {
            std::cout << " done in " << timer.diff() << " s." << std::endl;
        }
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

void
PreconditionerPCD::setDataFromGetPot ( const GetPot& dataFile,
                                       const std::string& section )
{
    M_dataFile   = dataFile;
    this->createParametersList ( M_list, dataFile, section, "PCD" );
    this->setParameters ( M_list );
}

void
PreconditionerPCD::setParameters ( Teuchos::ParameterList& list )
{
    M_precType                         = list.get ( "prectype", "PCD" );

    M_fluidPrec                        = list.get ( "subprecs: fluid prec", "ML" );
    M_fluidPrecDataSection             = list.get ( "subprecs: fluid prec data section", "" );

    M_pressureLaplacianPrec            = list.get ( "subprecs: pressure laplacian prec", "ML");
    M_pressureLaplacianPrecDataSection = list.get ( "subprecs: pressure laplacian prec data section", "" );

    M_pressureMassPrec                 = list.get ( "subprecs: pressure mass prec", "ML" );
    M_pressureMassPrecDataSection      = list.get ( "subprecs: pressure mass prec data section", "" );

    M_pressureLaplacianOperator        = list.get ( "pressure laplacian operator", "standard" );

    M_pressureMassOperator             = list.get ( "pressure mass operator", "lumped" );
    M_setApBoundaryConditions          = list.get ( "set Ap boundary conditions", false );
    M_setFpBoundaryConditions          = list.get ( "set Fp boundary conditions", false );
    M_setMpBoundaryConditions          = list.get ( "set Mp boundary conditions", false );
    M_fullFactorization                = list.get ( "full factorization", false );
    M_useStiffStrain                   = list.get ( "use stiff strain", false );
    M_enableTransient                  = list.get ( "enable transient", true );
    bool useMinusDivergence            = list.get ( "use minus divergence", false );
    this->setUseMinusDivergence ( useMinusDivergence );
    M_recomputeNormalVectors           = list.get ( "recompute normal vectors", true );
    M_schurOperatorReverseOrder        = list.get ( "Schur operator reverse order", false );
    M_inflowBoundaryType               = list.get ( "inflow boundary type", "Robin" );
    M_outflowBoundaryType              = list.get ( "outflow boundary type", "Neumann" );
    M_characteristicBoundaryType       = list.get ( "characteristic boundary type", "Neumann" );
}

void
PreconditionerPCD::setFESpace ( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace )
{
    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_adrPressureAssembler.setup ( pFESpace, uFESpace ); // p,beta=u
    M_adrVelocityAssembler.setup ( uFESpace, uFESpace ); // u,beta=u
    M_beta.reset ( new vector_Type ( M_uFESpace->map() + M_pFESpace->map() ) );
    *M_beta *= 0;

    // We setup the size of the blocks
    M_velocityBlockSize = M_uFESpace->fieldDim() * M_uFESpace->dof().numTotalDof();
    M_pressureBlockSize = M_pFESpace->dof().numTotalDof();
}

void
PreconditionerPCD::setBCHandler ( BCHandlerPtr_Type bchPtr )
{
    M_bcHandlerPtr = bchPtr;
}

void
PreconditionerPCD::setTimestep ( const Real& timestep )
{
    M_timestep = timestep;
}

void
PreconditionerPCD::setViscosity ( const Real& viscosity )
{
    M_viscosity = viscosity;
}

void
PreconditionerPCD::setDensity ( const Real& density )
{
    M_density = density;
}

void
PreconditionerPCD::setUseMinusDivergence ( const bool& useMinusDivergence )
{
    if ( useMinusDivergence )
    {
        M_divergenceCoeff = -1.0;
    }
    else
    {
        M_divergenceCoeff = 1.0;
    }
}

void
PreconditionerPCD::computeNormalVectors()
{
    bool verbose ( false );
    if ( M_comm->MyPID() == 0 )
    {
        verbose = true;
    }

    LifeChrono timer;

    if ( verbose )
    {
        std::cout << "      >Computing the normal vectors on the boundary... ";
    }
    timer.start();

    //-----------------------------------------------------
    // STEP 1: Calculating the normals
    //-----------------------------------------------------

    vector_Type repNormals ( M_pFESpace->map() + M_pFESpace->map() + M_pFESpace->map(), Repeated );
    UInt numTotalDofs ( M_pFESpace->dof().numTotalDof() );

    //Loop on the Faces
    UInt numBoundaryFacets ( M_pFESpace->mesh()->numBoundaryFacets() );
    for ( UInt iFace = 0; iFace < numBoundaryFacets; ++iFace )
    {
        //Update the currentBdFE with the face data
        M_pFESpace->feBd().update ( M_pFESpace->mesh()->boundaryFacet ( iFace ), UPDATE_NORMALS | UPDATE_W_ROOT_DET_METRIC );
        UInt nDofF = M_pFESpace->feBd().nbFEDof();

        //For each node on the face
        for ( UInt icheck = 0; icheck < nDofF; ++icheck )
        {
            ID idf = M_pFESpace->dof().localToGlobalMapByBdFacet ( iFace, icheck );

            //If the point is on this processor
            if ( repNormals.isGlobalIDPresent ( idf ) )
            {
                // ID flag = M_flags[idf];

                //if the normal is not already calculated
                //and the marker correspond to the flag of the point
                // if ((flag == M_uFESpace->mesh().boundaryFacet( iFace ).markerID())||(flag == 0))
                // {
                //Warning: the normal is taken in the first Gauss point
                //since the normal is the same over the triangle
                //(not true in the case of quadratic and bilinear maps)
                Real nx ( M_pFESpace->feBd().normal ( 0, 0 ) );
                Real ny ( M_pFESpace->feBd().normal ( 1, 0 ) );
                Real nz ( M_pFESpace->feBd().normal ( 2, 0 ) );

                //We get the area
                Real area ( M_pFESpace->feBd().measure() );

                //We update the normal component of the boundary point
                ( repNormals ) [idf]                  += nx * area;
                ( repNormals ) [idf +   numTotalDofs] += ny * area;
                ( repNormals ) [idf + 2 * numTotalDofs] += nz * area;
                // }
            }
        }
    }

    //-----------------------------------------------------
    // STEP 2: Gathering the data from others processors
    //-----------------------------------------------------

    M_normalVectors.reset ( new vector_Type ( repNormals, Unique ) );

    //-----------------------------------------------------
    // STEP 3: Normalizing the vectors
    //-----------------------------------------------------

    //We obtain the ID of the element
    Int NumMyElements = M_normalVectors->map().map ( Unique )->NumMyElements();
    std::vector<Int> MyGlobalElements ( NumMyElements );
    M_normalVectors->map().map ( Unique )->MyGlobalElements ( & ( MyGlobalElements[0] ) );

    //We normalize the normal
    Real norm;
    UInt id;

    //Need to run only over the first third of MyGlobalElements
    //(the larger values are the y and z components)
    for ( Int i ( 0 ); i < NumMyElements / static_cast<Int> ( nDimensions ); ++i )
    {
        id = MyGlobalElements[i];
        Real nx ( ( *M_normalVectors ) [id]                  );
        Real ny ( ( *M_normalVectors ) [id +   numTotalDofs] );
        Real nz ( ( *M_normalVectors ) [id + 2 * numTotalDofs] );
        norm = sqrt ( nx * nx + ny * ny + nz * nz );

        // If the norm is exactly 0 it means that the point is an interior point
        if ( norm != 0.0 )
        {
            ( *M_normalVectors ) [id]                  /= norm;
            ( *M_normalVectors ) [id +   numTotalDofs] /= norm;
            ( *M_normalVectors ) [id + 2 * numTotalDofs] /= norm;
        }
    }

    if ( verbose )
    {
        std::cout << " done in " << timer.diff() << " s." << std::endl;
    }
}

PreconditionerPCD::vectorPtr_Type
PreconditionerPCD::computeRobinCoefficient()
{
    if ( M_normalVectors.get() == 0 || M_recomputeNormalVectors )
    {
        this->computeNormalVectors();
    }

    vectorPtr_Type robinCoeffVector ( new vector_Type ( M_pFESpace->map(), Unique ) );

    //We obtain the ID of the element
    Int numVelocityTotalDofs ( M_uFESpace->dof().numTotalDof() );
    Int numPressureTotalDofs ( M_pFESpace->dof().numTotalDof() );
    Int NumMyElements = robinCoeffVector->map().map ( Unique )->NumMyElements();
    std::vector<Int> MyGlobalElements ( NumMyElements );
    robinCoeffVector->map().map ( Unique )->MyGlobalElements ( & ( MyGlobalElements[0] ) );

    //We compute (beta*n)/nu
    UInt id;
    Real invNu = -M_density / M_viscosity;

    //Need to run only over the first third of MyGlobalElements
    //(the larger values are the y and z components)
    for ( Int i ( 0 ); i < NumMyElements; ++i )
    {
        id = MyGlobalElements[i];
        Real x ( ( *M_normalVectors ) [id]                          * ( *M_beta ) [id]                          );
        Real y ( ( *M_normalVectors ) [id +   numPressureTotalDofs] * ( *M_beta ) [id +   numVelocityTotalDofs] );
        Real z ( ( *M_normalVectors ) [id + 2 * numPressureTotalDofs] * ( *M_beta ) [id + 2 * numVelocityTotalDofs] );
        ( *robinCoeffVector ) [id] = ( x + y + z ) * invNu;
    }

    return robinCoeffVector;
}

void
PreconditionerPCD::setBCByBoundaryType ( matrixPtr_type Ap, UInt ApOffset,
                                         matrixPtr_type Fp, UInt FpOffset,
                                         matrixPtr_type Mp, UInt MpOffset )
{
    std::string boundaryType ( "none" );
    Real indicator ( 0.0 );

    vector_Type    robinCoeffVector ( *computeRobinCoefficient(), Repeated );
    vector_Type    robinRHS ( M_pFESpace->map(), Repeated );
    BCVector uRobin ( robinRHS, M_pFESpace->dof().numTotalDof(), 0 );
    uRobin.setRobinCoeffVector ( robinCoeffVector );
    uRobin.setBetaCoeff ( 0 );
    BCFunctionBase uZero ( fZero );

    // <-- DEBUG
    /*{
    Int NumMyElements = M_pFESpace->map().map( Unique )->NumMyElements();
    std::vector<Int> MyGlobalElements( NumMyElements );
    M_pFESpace->map().map( Unique )->MyGlobalElements(&MyGlobalElements[0]);
    ID idof(0);

    UInt numPressureDofs = M_pFESpace->dof().numTotalDof();

    VTKWriter writer;
    writer.setLabel( "Debug normal and boundary classification" );
    writer.addScalarField( "BCType" );
    writer.addVectorField( "Normal_vectors" );

    Real x, y, z;
    UInt boundaryTypeId;
    for ( Int i(0); i<NumMyElements; ++i )
    {
        idof = MyGlobalElements[i];

        if( robinCoeffVector[idof] > 0.0 )
        {
            // Inflow
            boundaryTypeId = 1;
        }
        else if( robinCoeffVector[idof] < 0.0 )
        {
            // Outflow
            boundaryTypeId = 2;
        }
        else if( robinCoeffVector[idof] == 0.0 )
        {
            // Characteristic
            boundaryTypeId = 3;
        }
        writer.addToScalarField( "BCType", boundaryTypeId );
        writer.addToScalarField( "Robin_coeff", robinCoeffVector[idof] );
        writer.addToVectorField( "Normal_vectors",
                                 (*M_normalVectors)[idof],
                                 (*M_normalVectors)[idof +   numPressureDofs],
                                 (*M_normalVectors)[idof + 2*numPressureDofs] );

        for ( UInt j( 0 ); j < M_pFESpace->mesh()->pointList.size(); ++j )
        {
            if ( idof == M_pFESpace->mesh()->pointList[j].id() )
            {
                x = M_pFESpace->mesh()->pointList[j].x();
                y = M_pFESpace->mesh()->pointList[j].y();
                z = M_pFESpace->mesh()->pointList[j].z();
                j = M_pFESpace->mesh()->pointList.size();
            }
        }
        writer.addPoint( x, y, z );
    }
    std::string fileName( "debugInfos_" );
    std::ostringstream ossMyPid;
    ossMyPid << M_pFESpace->map().comm().MyPID();
    fileName.append( ossMyPid.str() );
    writer.write( fileName );
    }*/
    // --> DEBUG END

    BCHandler bcHandler;

    // BC for the inflow, outflow, and characteristic boundaries
    BCBase* currentBoundary;
    BCBase* inflowBoundary, *outflowBoundary, *characteristicBoundary;

    // inflow
    if ( M_inflowBoundaryType == "Dirichlet" )
    {
        bcHandler.addBC ( BCBase ( "inflow",
                                   1,
                                   Essential,
                                   Full,
                                   uZero,
                                   1 ) );

    }
    else if ( M_inflowBoundaryType == "Neumann" )
    {
        bcHandler.addBC ( BCBase ( "inflow",
                                   1,
                                   Natural,
                                   Full,
                                   uZero,
                                   1 ) );
    }
    else if ( M_inflowBoundaryType == "Robin" )
    {
        bcHandler.addBC ( BCBase ( "inflow",
                                   1,
                                   Robin,
                                   Full,
                                   uRobin,
                                   1 ) );
    }

    // outflow
    if ( M_outflowBoundaryType == "Dirichlet" )
    {
        bcHandler.addBC ( BCBase ( "outflow",
                                   2,
                                   Essential,
                                   Full,
                                   uZero,
                                   1 ) );
    }
    else if ( M_outflowBoundaryType == "Neumann" )
    {
        bcHandler.addBC ( BCBase ( "outflow",
                                   2,
                                   Natural,
                                   Full,
                                   uZero,
                                   1 ) );
    }
    else if ( M_outflowBoundaryType == "Robin" )
    {
        bcHandler.addBC ( BCBase ( "outflow",
                                   2,
                                   Robin,
                                   Full,
                                   uRobin,
                                   1 ) );
    }

    // characteristic
    if ( M_characteristicBoundaryType == "Dirichlet" )
    {
        bcHandler.addBC ( BCBase ( "characteristic",
                                   3,
                                   Essential,
                                   Full,
                                   uZero,
                                   1 ) );
    }
    else if ( M_characteristicBoundaryType == "Neumann" )
    {
        bcHandler.addBC ( BCBase ( "characteristic",
                                   3,
                                   Natural,
                                   Full,
                                   uZero,
                                   1 ) );
    }
    else if ( M_characteristicBoundaryType == "Robin" )
    {
        bcHandler.addBC ( BCBase ( "characteristic",
                                   3,
                                   Robin,
                                   Full,
                                   uRobin,
                                   1 ) );
    }

    inflowBoundary         = & ( bcHandler.findBCWithFlag ( 1 ) );
    outflowBoundary        = & ( bcHandler.findBCWithFlag ( 2 ) );
    characteristicBoundary = & ( bcHandler.findBCWithFlag ( 3 ) );

    bcFlag_Type elementMarker; //will store the marker of the element

    // Loop on boundary faces
    for ( ID iBoundaryElement = 0 ; iBoundaryElement < M_pFESpace->mesh()->numBoundaryFacets(); ++iBoundaryElement )
    {
        // construction of localToGlobalMapOnBElem (this part should be moved in DOF.hpp)
        M_pFESpace->feBd().update ( M_pFESpace->mesh()->boundaryFacet ( iBoundaryElement ), UPDATE_W_ROOT_DET_METRIC ); // updating finite element information
        elementMarker = M_pFESpace->mesh()->boundaryFacet ( iBoundaryElement ).markerID(); // We keep the element marker


        //vector containing the local to global map on each element
        std::vector<ID> localToGlobalMapOnBElem = M_pFESpace->dof().localToGlobalMapOnBdFacet ( iBoundaryElement );

        // We classify and build bc objects:
        // a) inflow         where -(w*n)/nu > 0
        // b) outflow        where -(w*n)/nu < 0
        // c) characteristic where -(w*n)/nu = 0
        indicator = 0.0;
        for (ID lDof = 0; lDof < localToGlobalMapOnBElem.size(); lDof++)
        {
            ID gDof = localToGlobalMapOnBElem[lDof];
            indicator += robinCoeffVector[gDof];
        }
        if ( indicator > 0.0 )      // inflow
        {
            boundaryType = M_inflowBoundaryType;
            currentBoundary = inflowBoundary;
        }
        else if ( indicator < 0.0 ) // outflow
        {
            boundaryType = M_outflowBoundaryType;
            currentBoundary = outflowBoundary;
        }
        else if ( indicator == 0.0 ) // characteristic
        {
            boundaryType = M_characteristicBoundaryType;
            currentBoundary = characteristicBoundary;
        }

        // Setup the boundary conditions
        if ( boundaryType == "Dirichlet" )
        {
            for (ID lDof = 0; lDof < localToGlobalMapOnBElem.size(); lDof++)
            {
                ID gDof = localToGlobalMapOnBElem[lDof]; // global DOF

                //providing Essential boundary conditions with user defined functions
                currentBoundary->addBCIdentifier ( new BCIdentifierEssential ( gDof, 0.0, 0.0, 0.0 ) );
            }
        }
        else if ( boundaryType == "Robin" )
        {
            //providing Robin boundary conditions with global DOFs on element
            currentBoundary->addBCIdentifier ( new BCIdentifierNatural ( iBoundaryElement, localToGlobalMapOnBElem ) );
        }
        else if ( boundaryType == "Neumann" )
        {
            // For Neumann BC we do not need to do anything
        }
    }

    // Copying the IDs
    inflowBoundary->copyIdSetIntoIdVector();
    outflowBoundary->copyIdSetIntoIdVector();
    characteristicBoundary->copyIdSetIntoIdVector();

    // We do not update the BCHandler!

    if ( M_setApBoundaryConditions )
    {
        bcHandler.setOffset ( ApOffset );
        bcManageMatrix ( *Ap, *M_pFESpace->mesh(), M_pFESpace->dof(), bcHandler, M_pFESpace->feBd(), 1.0, 0.0 );
    }
    if ( M_setFpBoundaryConditions )
    {
        bcHandler.setOffset ( FpOffset );
        bcManageMatrix ( *Fp, *M_pFESpace->mesh(), M_pFESpace->dof(), bcHandler, M_pFESpace->feBd(), 1.0, 0.0 );
    }
    if ( M_setMpBoundaryConditions )
    {
        bcHandler.setOffset ( MpOffset );
        bcManageMatrix ( *Mp, *M_pFESpace->mesh(), M_pFESpace->dof(), bcHandler, M_pFESpace->feBd(), 1.0, 0.0 );
    }
    Ap->globalAssemble();
    Fp->globalAssemble();
    Mp->globalAssemble();
}

} // namespace LifeV
