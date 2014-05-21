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


#include <lifev/core/LifeV.hpp>

#include <lifev/fsi/solver/FSIMonolithicGI.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDND.hpp>
#include <lifev/fsi/solver/MonolithicBlockMatrixRN.hpp>

namespace LifeV
{

// ===================================================
//  Constructors and Descructor
// ===================================================
FSIMonolithicGI::FSIMonolithicGI() :
    super_Type              (),
    M_mapWithoutMesh        (),
    M_uk                    (),
    M_interface             (0),
    M_meshBlock             (),
    M_shapeDerivativesBlock (),
    M_solidDerBlock         ()
{
}

// ===================================================
//  Public Methods
// ===================================================
void FSIMonolithicGI::setup ( const GetPot& dataFile )
{
    super_Type::setup ( dataFile );
}

void FSIMonolithicGI::setupFluidSolid ( UInt const fluxes )
{
    super_Type::setupFluidSolid ( fluxes );
    UInt offset = M_monolithicMap->map ( Unique )->NumGlobalElements();
    M_mapWithoutMesh.reset ( new MapEpetra ( *M_monolithicMap ) );

    *M_monolithicMap += M_mmFESpace->map();

    M_interface = M_monolithicMatrix->interface();

    M_beta.reset ( new vector_Type (M_uFESpace->map() ) );
    M_rhs.reset (new vector_Type (*M_monolithicMap) );
    M_rhsFull.reset (new vector_Type (*M_monolithicMap) );
    if (M_data->dataFluid()->useShapeDerivatives() )
    {
        M_shapeDerivativesBlock.reset (new matrix_Type (*M_monolithicMap) );
    }
    M_uk.reset (new vector_Type (*M_monolithicMap) );

    M_meshMotion.reset (new meshMotion_Type (*M_mmFESpace,
                                             M_epetraComm,
                                             *M_monolithicMap,
                                             offset) );

    M_fluid.reset     (new fluid_Type (M_data->dataFluid(),
                                       *M_uFESpace,
                                       *M_pFESpace,
                                       *M_mmFESpace,
                                       M_epetraComm,
                                       *M_monolithicMap,
                                       fluxes) );
    M_solid.reset (new solid_Type() );

    M_solid->setup (M_data->dataSolid(),
                    M_dFESpace,
                    M_dETFESpace,
                    M_epetraComm,
                    M_dFESpace->mapPtr(),
                    UInt (0)
                   );
}

void
FSIMonolithicGI::buildSystem()
{
    super_Type::buildSystem();
    M_meshMotion->computeMatrix();
}

void
FSIMonolithicGI::evalResidual ( vector_Type&       res,
                                const vector_Type& disp,
                                const UInt          iter )
{
    // disp here is the current solution guess (u,p,ds,df)

    res = 0.;//this is important. Don't remove it!

    M_uk.reset (new vector_Type ( disp ) );

    UInt offset ( M_solidAndFluidDim + nDimensions * M_interface );

    vectorPtr_Type meshDisp ( new vector_Type (M_mmFESpace->map() ) );

    meshDisp->subset (disp, offset); //if the conv. term is to be condidered implicitly

    vectorPtr_Type mmRep ( new vector_Type (*meshDisp, Repeated ) );

    moveMesh ( *mmRep );

    //here should use extrapolationFirstDerivative instead of velocity
    if (iter == 0)
    {
        M_ALETimeAdvance->updateRHSFirstDerivative (M_data->dataFluid()->dataTime()->timeStep() );
    }
    vector_Type meshVelocityRepeated ( this->M_ALETimeAdvance->firstDerivative ( *meshDisp ), Repeated );
    vector_Type interpolatedMeshVelocity (this->M_uFESpace->map() );

    interpolateVelocity ( meshVelocityRepeated, interpolatedMeshVelocity );
    vectorPtr_Type fluid ( new vector_Type ( M_uFESpace->map() ) );
    M_beta->subset ( disp, 0 );

    *M_beta -= interpolatedMeshVelocity; // convective term, u^(n+1) - w^(n+1)

    assembleSolidBlock ( iter, disp );
    assembleFluidBlock ( iter, disp );
    assembleMeshBlock ( iter );

    *M_rhsFull = *M_rhs;

    M_monolithicMatrix->setRobin ( M_robinCoupling, M_rhsFull );
    M_precPtr->setRobin (M_robinCoupling, M_rhsFull);

    if (!M_monolithicMatrix->set() )
    {
        M_BChs.push_back (M_BCh_d);
        M_BChs.push_back (M_BCh_u);

        M_FESpaces.push_back (M_dFESpace);
        M_FESpaces.push_back (M_uFESpace);

        M_BChs.push_back (M_BCh_mesh);
        M_FESpaces.push_back (M_mmFESpace);

        M_monolithicMatrix->push_back_matrix (M_solidBlockPrec, false);
        M_monolithicMatrix->push_back_matrix (M_fluidBlock, true);
        M_monolithicMatrix->push_back_matrix (M_meshBlock, false);
        M_monolithicMatrix->setConditions (M_BChs);
        M_monolithicMatrix->setSpaces (M_FESpaces);
        M_monolithicMatrix->setOffsets (3, M_offset, 0, M_solidAndFluidDim + nDimensions * M_interface);
        M_monolithicMatrix->coupler (M_monolithicMap, M_dofStructureToFluid->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), M_solidTimeAdvance->coefficientFirstDerivative ( 0 ), M_solid->rescaleFactor() );
        M_monolithicMatrix->coupler ( M_monolithicMap, M_dofStructureToFluid->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), M_solidTimeAdvance->coefficientFirstDerivative ( 0 ), M_solid->rescaleFactor(), 2);
    }
    else
    {
        M_monolithicMatrix->replace_matrix (M_solidBlockPrec, 0);
        M_monolithicMatrix->replace_matrix (M_fluidBlock, 1);
        M_monolithicMatrix->replace_matrix (M_meshBlock, 2);
    }

    M_monolithicMatrix->blockAssembling();
    super_Type::checkIfChangedFluxBC ( M_monolithicMatrix );

    // formulation matrix * vector (i.e. linear elastic )
    // todo: pass to boolean for nonlinear structures
    if ( ! (M_data->dataSolid()->lawType().compare ("linear") ) )
    {
        applyBoundaryConditions();
    }

    M_monolithicMatrix->GlobalAssemble();

    super_Type::evalResidual ( disp, M_rhsFull, res, false );

    //case for nonlinear laws which are formulated in the residual form
    if ( ! ( M_data->dataSolid()->lawType().compare ( "nonlinear" ) ) )
    {
        res += *M_meshBlock * disp;

        if ( !M_BCh_u->bcUpdateDone() )
        {
            M_BCh_u->bcUpdate ( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
        }
        M_BCh_d->setOffset ( M_offset );
        if ( !M_BCh_d->bcUpdateDone() )
        {
            M_BCh_d->bcUpdate ( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
        }
        M_BCh_mesh->setOffset ( M_solidAndFluidDim + nDimensions * M_interface );
        if ( !M_BCh_mesh->bcUpdateDone() )
        {
            M_BCh_mesh->bcUpdate ( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );
        }

        M_monolithicMatrix->applyBoundaryConditions ( dataFluid()->dataTime()->time() /*, M_rhsFull*/);
        M_monolithicMatrix->GlobalAssemble();

        bcManageResidual ( res, *M_rhsFull, disp, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_u,
                           M_uFESpace->feBd(), M_data->dataFluid()->dataTime()->time(), 1. );

        // below sol is repeated by BCManageResidual
        bcManageResidual ( res, *M_rhsFull, disp, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d,
                           M_dFESpace->feBd(), M_data->dataSolid()->dataTime()->time(), 1. );
        bcManageResidual ( res, *M_rhsFull, disp, *M_mmFESpace->mesh(), M_mmFESpace->dof(), *M_BCh_mesh,
                           M_mmFESpace->feBd(), M_data->dataFluid()->dataTime()->time(), 1. );
        res -= *M_rhsFull;
    }
}

void FSIMonolithicGI::applyBoundaryConditions()
{
    if ( !M_BCh_u->bcUpdateDone() )
        M_BCh_u->bcUpdate ( *M_uFESpace->mesh(),
                            M_uFESpace->feBd(),
                            M_uFESpace->dof() );
    M_BCh_d->setOffset ( M_offset );
    if ( !M_BCh_d->bcUpdateDone() )
        M_BCh_d->bcUpdate ( *M_dFESpace->mesh(),
                            M_dFESpace->feBd(),
                            M_dFESpace->dof() );
    M_BCh_mesh->setOffset ( M_solidAndFluidDim + nDimensions * M_interface );
    if ( !M_BCh_mesh->bcUpdateDone() )
        M_BCh_mesh->bcUpdate ( *M_mmFESpace->mesh(),
                               M_mmFESpace->feBd(),
                               M_mmFESpace->dof() );

    M_monolithicMatrix->applyBoundaryConditions (dataFluid()->dataTime()->time(), M_rhsFull);
}

void FSIMonolithicGI::setALEVectorInStencil ( const vectorPtr_Type& fluidDisp,
                                              const UInt iter,
                                              const bool lastVector)
{

    //The fluid timeAdvance has size = orderBDF because it is seen as an equation frist order in time.
    //We need to add the solidVector to the fluidVector in the fluid TimeAdvance because we have the
    //extrapolation on it.
    if ( ( iter < M_fluidTimeAdvance->size() - 1 ) && !lastVector )
    {
        vectorPtr_Type harmonicSolutionRestartTime (new vector_Type ( *M_monolithicMap, Unique, Zero ) );

        *harmonicSolutionRestartTime *= 0.0;

        UInt givenOffset ( M_solidAndFluidDim + nDimensions * M_interface );
        harmonicSolutionRestartTime->subset (*fluidDisp, fluidDisp->map(), (UInt) 0, givenOffset );

        //We sum the vector in the first element of fluidtimeAdvance
        * ( M_fluidTimeAdvance->stencil() [ iter + 1 ] ) += *harmonicSolutionRestartTime;
    }

    if ( !lastVector )
    {
        //The shared_pointer for the vectors has to be trasformed into a pointer to VectorEpetra
        //That is the type of pointers that are used in TimeAdvance
        vector_Type* normalPointerToALEVector ( new vector_Type (*fluidDisp) );
        (M_ALETimeAdvance->stencil() ).push_back ( normalPointerToALEVector );
    }
    else
    {
        vectorPtr_Type harmonicSolutionRestartTime (new vector_Type ( *M_monolithicMap, Unique, Zero ) );

        *harmonicSolutionRestartTime *= 0.0;

        UInt givenOffset ( M_solidAndFluidDim + nDimensions * M_interface );
        harmonicSolutionRestartTime->subset (*fluidDisp, fluidDisp->map(), (UInt) 0, givenOffset );

        //We sum the vector in the first element of fluidtimeAdvance
        * ( M_fluidTimeAdvance->stencil() [ 0 ] ) += *harmonicSolutionRestartTime;
    }

}



//============ Protected Methods ===================

void FSIMonolithicGI::setupBlockPrec()
{

    //The following part accounts for a possibly nonlinear structure model, should not be run when linear
    //elasticity is used

    // case of exponential and neohookean model
    // todo: pass to boolean variable for Nonlinear models ( i.e. for vector formulation )
    // if ( M_data->dataSolid()->getUseExactJacobian() && ( M_data->dataSolid()->solidType().compare( "exponential" )
    //                           && M_data->dataSolid()->solidType().compare( "neoHookean" ) ) )
    // {
    //     M_solid->material()->updateJacobianMatrix ( *M_uk * M_solid->rescaleFactor(),
    //                                                 dataSolid(),
    //                                                 M_solid->mapMarkersVolumes(),
    //                                                 M_solid->mapMarkersIndexes(),
    //                                                 M_solid->displayerPtr() ); // computing the derivatives if nonlinear (comment this for inexact Newton);
    //     M_solidBlockPrec.reset ( new matrix_Type ( *M_monolithicMap,
    //                                                1 ) );
    //     *M_solidBlockPrec += *M_solid->massMatrix();
    //     *M_solidBlockPrec += *M_solid->material()->jacobian();
    //     M_solidBlockPrec->globalAssemble();
    //     *M_solidBlockPrec *= M_solid->rescaleFactor();

    //     M_monolithicMatrix->replace_matrix ( M_solidBlockPrec, 0 );
    // }

    M_monolithicMatrix->blockAssembling();

    if ( M_data->dataFluid()->useShapeDerivatives() )
    {
        *M_shapeDerivativesBlock *= 0.;
        M_shapeDerivativesBlock->openCrsMatrix();
        shapeDerivatives ( M_shapeDerivativesBlock );
        M_shapeDerivativesBlock->globalAssemble();
        M_monolithicMatrix->addToGlobalMatrix ( M_shapeDerivativesBlock );
    }

    if ( M_data->dataFluid()->useShapeDerivatives() || M_data->dataSolid()->getUseExactJacobian() )
    {
        if ( !M_BCh_u->bcUpdateDone() )
            M_BCh_u->bcUpdate ( *M_uFESpace->mesh(),
                                M_uFESpace->feBd(),
                                M_uFESpace->dof() );
        M_BCh_d->setOffset ( M_offset );
        if ( !M_BCh_d->bcUpdateDone() )
            M_BCh_d->bcUpdate ( *M_dFESpace->mesh(),
                                M_dFESpace->feBd(),
                                M_dFESpace->dof() );
        M_BCh_mesh->setOffset ( M_solidAndFluidDim + nDimensions * M_interface );
        if ( !M_BCh_mesh->bcUpdateDone() )
            M_BCh_mesh->bcUpdate ( *M_mmFESpace->mesh(),
                                   M_mmFESpace->feBd(),
                                   M_mmFESpace->dof() );

        M_monolithicMatrix->applyBoundaryConditions ( dataFluid()->dataTime()->time() );
        M_monolithicMatrix->GlobalAssemble();
    }

    super_Type::setupBlockPrec();

    if ( M_precPtr->blockVector().size() < 3 )
    {
        M_precPtr->push_back_matrix ( M_meshBlock,
                                      false );
        M_precPtr->setConditions ( M_BChs );
        M_precPtr->setSpaces ( M_FESpaces );
        M_precPtr->setOffsets ( 3, M_offset, 0,  M_solidAndFluidDim + nDimensions * M_interface );
        M_precPtr->coupler ( M_monolithicMap, M_dofStructureToFluid->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep() , M_solidTimeAdvance->coefficientFirstDerivative ( 0 ), M_solid->rescaleFactor(), 2);

        if (M_data->dataFluid()->useShapeDerivatives() )
        {
            boost::shared_ptr<MatrixEpetra<Real> > staticCast = boost::static_pointer_cast<MatrixEpetra<Real> > (M_shapeDerivativesBlock);
            M_precPtr->push_back_coupling ( staticCast );
        }
    }
    else
    {
        //M_precPtr->replace_matrix( M_solidBlockPrec, 0 );
        //M_precPtr->replace_matrix( M_fluidBlock, 1 );
        M_precPtr->replace_matrix ( M_meshBlock, 2 );

        if ( M_data->dataFluid()->useShapeDerivatives() )
        {
            M_precPtr->replace_coupling ( M_shapeDerivativesBlock, 2 );
        }
    }
}

void FSIMonolithicGI::shapeDerivatives ( FSIOperator::fluid_Type::matrixPtr_Type sdMatrix )
{
    Real alpha = M_fluidTimeAdvance->coefficientFirstDerivative ( 0 ) / M_data->dataFluid()->dataTime()->timeStep();
    vectorPtr_Type rhsNew (new vector_Type (*M_monolithicMap) );
    vector_Type un (M_uFESpace->map() );
    vector_Type uk (M_uFESpace->map() + M_pFESpace->map() );

    vectorPtr_Type meshVelRep ( new vector_Type ( M_mmFESpace->map(), Repeated ) );

    *meshVelRep = M_ALETimeAdvance->firstDerivative();

    //When this class is used, the convective term is used implictly
    un.subset ( *M_uk, 0 );

    uk.subset ( *M_uk, 0 );
    vector_Type veloFluidMesh ( M_uFESpace->map(), Repeated );
    this->transferMeshMotionOnFluid ( *meshVelRep, veloFluidMesh );

    //The last two flags are consistent with the currect interface.
    //When this class is used, they should not be changed.
    M_fluid->updateShapeDerivatives ( *sdMatrix, alpha,
                                      un,
                                      uk,
                                      veloFluidMesh, //(xk-xn)/dt (FI), or (xn-xn-1)/dt (CE)//Repeated
                                      M_solidAndFluidDim + M_interface * nDimensions,
                                      *M_mmFESpace,
                                      true /*This flag tells the method to consider the velocity of the domain implicitly*/,
                                      true /*This flag tells the method to consider the convective term implicitly */ );
}

void FSIMonolithicGI::assembleMeshBlock ( UInt /*iter*/ )
{

    M_meshBlock.reset ( new matrix_Type ( *M_monolithicMap ) );
    M_meshMotion->addSystemMatrixTo ( M_meshBlock );
    M_meshBlock->globalAssemble();
    UInt offset ( M_solidAndFluidDim + nDimensions * M_interface );
    std::map< ID, ID >::const_iterator ITrow;
    std::map< ID, ID > locdofmap ( M_dofStructureToFluid->localDofMap() );

    /******************alternative way************************/
    //     BCFunctionBase bcf(fZero);
    //     fluidBchandlerPtr_Type BCh(new fluidBchandler_Type() );
    //     BCh->addBC("Interface", 1, Essential, Full,
    //                bcf, 3);

    //     BCh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);

    //     if ( !BCh->bcUpdateDone() )
    //         BCh->bcUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

    //     bcManage( *M_meshBlock, *M_rhsFull, *M_mmFESpace->mesh(), M_mmFESpace->dof(), *BCh, M_mmFESpace->feBd(), 1., dataFluid()->dataTime()->time());
    /********************************************************/

    matrixPtr_Type diagFilter (new matrix_Type (*M_monolithicMap) );
    diagFilter->insertZeroDiagonal ();
    for ( ID dim = 0; dim < nDimensions; ++dim )
    {
        for ( ITrow = locdofmap.begin(); ITrow != locdofmap.end(); ++ITrow )
        {
            int i = ITrow->first;
            int iRow =  i + offset + dim * M_mmFESpace->dof().numTotalDof();
            if ( M_meshBlock->mapPtr()->map (Unique)->MyGID (iRow) )
            {
                M_meshBlock->diagonalize ( iRow, 1. );
            }
            else
            {
                double myValues[1];
                myValues[0] = -1;
                int myIndices[1];
                myIndices[0] = iRow;
                diagFilter->matrixPtr()->SumIntoGlobalValues ( iRow, 1, myValues, myIndices );
            }
        }
    }

    diagFilter->globalAssemble();
    // Processor informations
    Int  numLocalEntries = diagFilter->matrixPtr()->RowMap().NumMyElements();
    Int* globalEntries   = diagFilter->matrixPtr()->RowMap().MyGlobalElements();
    UInt globalRowIndex (0);

    // Source informations handlers
    double* diagValue;
    Int* indices;
    Int srcRow (0);
    Int controlValue (0); // This value should be one since it is a diagonal matrix

    for (Int i (0); i < numLocalEntries; ++i)
    {
        // Collecting the data from the source
        globalRowIndex = globalEntries[i];

        // Get the data of the row
        srcRow = diagFilter->matrixPtr()->LRID (static_cast<EpetraInt_Type> (globalRowIndex) );
        diagFilter->matrixPtr()->ExtractMyRowView (srcRow, controlValue, diagValue, indices);

        ASSERT ( controlValue == 1, "Error: The matrix should be diagonal.");
        if (diagValue[0] < 0.0)
        {
            M_meshBlock->diagonalize ( globalRowIndex, 1. );
        }
    }
}

// ===================================================
//  Products registration
// ===================================================
bool FSIMonolithicGI::S_register = BlockPrecFactory::instance().registerProduct ( "AdditiveSchwarzGI",
                                   &MonolithicBlockMatrix::createAdditiveSchwarz )
                                   && BlockPrecFactory::instance().registerProduct ( "ComposedDNGI",
                                           &MonolithicBlockComposedDN::createComposedDNGI )
                                   && MonolithicBlockMatrix::Factory_Type::instance().registerProduct ( "AdditiveSchwarzGI",
                                           &MonolithicBlockMatrix::createAdditiveSchwarz )
                                   && MonolithicBlockMatrix::Factory_Type::instance().registerProduct ( "AdditiveSchwarzRNGI",
                                           &MonolithicBlockMatrixRN::createAdditiveSchwarzRN )
                                   && FSIOperator::FSIFactory_Type::instance().registerProduct ( "monolithicGI",
                                           &FSIMonolithicGI::instantiate )
                                   && BlockPrecFactory::instance().registerProduct ( "ComposedDNDGI",
                                           &MonolithicBlockComposedDND::createComposedDNDGI )
                                   && BlockPrecFactory::instance().registerProduct ( "ComposedDND2GI",
                                           &MonolithicBlockComposedDND::createComposedDND2GI )
                                   && BlockPrecFactory::instance().registerProduct ( "ComposedDN2GI",
                                           &MonolithicBlockComposedDN::createComposedDN2GI );

}
