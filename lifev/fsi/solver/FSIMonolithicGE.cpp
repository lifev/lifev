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
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 26 Jul 2010

    A more detailed description of the file (if necessary)
 */

#include <lifev/core/LifeV.hpp>

#include <lifev/fsi/solver/FSIMonolithicGE.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================


void FSIMonolithicGE::setupFluidSolid ( UInt const fluxes )
{
    super_Type::setupFluidSolid ( fluxes );

    M_meshMotion.reset (new FSIOperator::meshMotion_Type (*M_mmFESpace,
                                                          M_epetraComm) );

    M_fluid.reset (new FSIOperator::fluid_Type (M_data->dataFluid(),
                                                *M_uFESpace,
                                                *M_pFESpace,
                                                M_epetraComm,
                                                *M_monolithicMap,
                                                fluxes) );

    M_rhs.reset (new vector_Type (*this->M_monolithicMap) );
    M_rhsFull.reset (new vector_Type (*this->M_monolithicMap) );
    M_beta.reset  (new vector_Type (M_uFESpace->map() ) );

    M_solid.reset (new solid_Type() );

    M_solid->setup (M_data->dataSolid(),
                    M_dFESpace,
                    M_dETFESpace,
                    M_epetraComm,
                    M_dFESpace->mapPtr(),
                    UInt (0)
                   );
}




void FSIMonolithicGE::setupDOF()
{
    M_bcvStructureDispToHarmonicExtension.reset ( new  BCVectorInterface );
    super_Type::setupDOF();
}

void
FSIMonolithicGE::setupSystem( )
{
    super_Type::setupSystem();
    M_meshMotion->setUp ( M_dataFile );
}

void
FSIMonolithicGE::updateSystem()
{
    super_Type::updateSystem();
}


void
FSIMonolithicGE::evalResidual ( vector_Type&       res,
                                const vector_Type& disp,
                                const UInt          iter )
{
    // disp here is the current solution guess (u,p,ds)
    // disp is already "extrapolated", the main is doing it.

    if (iter == 0)
    {

        // Solve HE
        this->iterateMesh (disp);

        // Update displacement

        M_ALETimeAdvance->updateRHSFirstDerivative (M_data->dataFluid()->dataTime()->timeStep() );

        vector_Type meshDispRepeated ( M_meshMotion->disp(), Repeated );
        this->moveMesh (meshDispRepeated);

        //here should use extrapolationFirstDerivative instead of velocity
        vector_Type meshVelocityRepeated ( this->M_ALETimeAdvance->firstDerivative (  M_meshMotion->disp() ), Repeated );
        vector_Type interpolatedMeshVelocity (this->M_uFESpace->map() );

        interpolateVelocity ( meshVelocityRepeated, interpolatedMeshVelocity );
        // maybe we should use disp here too...
        M_fluidTimeAdvance->extrapolation (*M_beta);
        *M_beta -= interpolatedMeshVelocity; // convective term, u^* - w^*

        // in MonolithicGI here it used M_uk, which comes from disp
        assembleSolidBlock (iter, disp);
        assembleFluidBlock (iter, disp);
        *M_rhsFull = *M_rhs;

        applyBoundaryConditions();
    }
    super_Type::evalResidual ( disp,  M_rhsFull, res, M_diagonalScale);
}

void
FSIMonolithicGE::iterateMesh (const vector_Type& disp)
{
    vector_Type lambdaFluid (*M_interfaceMap, Unique);

    monolithicToInterface (lambdaFluid, disp);

    lambdaFluid *= (M_solid->rescaleFactor() );

    this->setLambdaFluid (lambdaFluid); // it must be _disp restricted to the interface

    M_meshMotion->iterate (*M_BCh_mesh);

}

void FSIMonolithicGE::applyBoundaryConditions( )
{

    if ( !M_BCh_u->bcUpdateDone() )
    {
        M_BCh_u->bcUpdate ( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    }
    M_BCh_d->setOffset (M_offset);
    if ( !M_BCh_d->bcUpdateDone() )
    {
        M_BCh_d->bcUpdate ( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    }

    M_monolithicMatrix->setRobin ( M_robinCoupling, M_rhsFull );
    M_precPtr->setRobin (M_robinCoupling, M_rhsFull);

    if (!this->M_monolithicMatrix->set() )
    {
        M_BChs.push_back (M_BCh_d);
        M_BChs.push_back (M_BCh_u);
        M_FESpaces.push_back (M_dFESpace);
        M_FESpaces.push_back (M_uFESpace);

        M_monolithicMatrix->push_back_matrix (M_solidBlockPrec, false);
        M_monolithicMatrix->push_back_matrix (M_fluidBlock, true);
        M_monolithicMatrix->setConditions (M_BChs);
        M_monolithicMatrix->setSpaces (M_FESpaces);
        M_monolithicMatrix->setOffsets (2, M_offset, 0);
        M_monolithicMatrix->coupler (M_monolithicMap, M_dofStructureToFluid->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), M_solidTimeAdvance->coefficientFirstDerivative ( 0 ), M_solid->rescaleFactor() );
    }
    else
    {
        M_monolithicMatrix->replace_matrix (M_fluidBlock, 1);
        M_monolithicMatrix->replace_matrix (M_solidBlockPrec, 0);
    }

    super_Type::checkIfChangedFluxBC ( M_monolithicMatrix );

    M_monolithicMatrix->blockAssembling();
    M_monolithicMatrix->applyBoundaryConditions (dataFluid()->dataTime()->time(), M_rhsFull);

    M_monolithicMatrix->GlobalAssemble();
    //M_monolithicMatrix->matrix()->spy("M");
}

void FSIMonolithicGE::setALEVectorInStencil (const vectorPtr_Type& fluidDisp, const UInt /*iter*/, const bool /*lastVector*/)
{

    //ALE problem
    //The shared_pointer for the vectors has to be transformed into a pointer to VectorEpetra
    //That is the type of pointers that are used in TimeAdvance
    vector_Type* normalPointerToALEVector ( new vector_Type (*fluidDisp) );
    (M_ALETimeAdvance->stencil() ).push_back ( normalPointerToALEVector );

}

void FSIMonolithicGE::updateSolution ( const vector_Type& solution )
{
    super_Type::updateSolution ( solution );

    //This updateRHSFirstDerivative has to be done before the shiftRight
    //In fact it updates the right hand side of the velocity using the
    //previous times. The method velocity() uses it and then, the computation
    //of the velocity is done using the current time and the previous times.
    //M_ALETimeAdvance->updateRHSFirstDerivative ( M_data->dataFluid()->dataTime()->timeStep() );
    M_ALETimeAdvance->shiftRight ( this->M_meshMotion->disp() );
}


// ===================================================
//! Products registration
// ===================================================

bool FSIMonolithicGE::S_register = FSIFactory_Type::instance().registerProduct ( "monolithicGE", &FSIMonolithicGE::instantiate )  &&
                                   BlockPrecFactory::instance().registerProduct ("ComposedDNND"  , &MonolithicBlockComposedDNND::createComposedDNND) &&
                                   BlockPrecFactory::instance().registerProduct ("AdditiveSchwarz"  , &MonolithicBlockMatrix::createAdditiveSchwarz) &&
                                   MonolithicBlockMatrix::Factory_Type::instance().registerProduct ("AdditiveSchwarz"  , &MonolithicBlockMatrix::createAdditiveSchwarz ) &&
                                   BlockPrecFactory::instance().registerProduct ("AdditiveSchwarzRN"  , &MonolithicBlockMatrixRN::createAdditiveSchwarzRN ) &&
                                   MonolithicBlockMatrix::Factory_Type::instance().registerProduct ("AdditiveSchwarzRN"  , &MonolithicBlockMatrixRN::createAdditiveSchwarzRN ) &&
                                   BlockPrecFactory::instance().registerProduct ("ComposedDN"  , &MonolithicBlockComposedDN::createComposedDN ) &&
                                   BlockPrecFactory::instance().registerProduct ("ComposedDN2"  , &MonolithicBlockComposedDN::createComposedDN2 );

} // Namespace LifeV
