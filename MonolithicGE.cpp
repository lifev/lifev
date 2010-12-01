//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 26 Jul 2010

    A more detailed description of the file (if necessary)
 */

#include <MonolithicGE.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================


void MonolithicGE::setupFluidSolid( UInt const fluxes )
{
    super::setupFluidSolid( fluxes );
    M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                            M_epetraComm));

    M_fluid.reset(new FSIOperator::fluid_raw_type(M_data->dataFluid(),
                                                  *M_uFESpace,
                                                  *M_pFESpace,
                                                  M_epetraComm,
                                                  *M_monolithicMap));

    //             if (isLinearFluid())// to be implemented
    //                 M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
    //                                                                    *M_uFESpace,
    //                                                                    *M_pFESpace,
    //                                                                    *M_epetraComm));
    M_un.reset (new vector_type(*this->M_monolithicMap));
    M_rhs.reset(new vector_type(*this->M_monolithicMap));
    M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
    M_beta.reset  (new vector_type(M_uFESpace->map()));

    M_solid.reset(solid_raw_type::StructureSolverFactory::instance().createObject( M_data->dataSolid()->solidType() ));

    M_solid->setup(M_data->dataSolid(),
                   M_dFESpace,
                   M_epetraComm,
                   M_monolithicMap,
                   M_offset
                  );

    //             if (isLinearSolid())// to be implemented with the offset
    //                 M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
    //                                                                    *M_dFESpace,
    //
    //                                                      *M_epetraComm));
}


void
MonolithicGE::evalResidual( vector_type&       res,
                            const vector_type& disp,
                            const UInt          iter )
{

    if ((iter==0)|| !this->M_data->dataFluid()->isSemiImplicit())
    {
        // Solve HE
        iterateMesh(disp);

        // Update displacement
        M_meshMotion->updateDispDiff();

        M_beta.reset(new vector_type(M_uFESpace->map()));
        vector_type meshDispDiff( M_meshMotion->disp(), Repeated );

        this->moveMesh(meshDispDiff);//initialize the mesh position with the total displacement

        meshDispDiff=M_meshMotion->dispDiff();//repeating the mesh dispDiff
        this->interpolateVelocity(meshDispDiff, *this->M_beta);

        *this->M_beta /= -M_data->dataFluid()->dataTime()->getTimeStep(); //mesh velocity w

        vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
        fluid->subset(*M_un, (UInt)0);
        *this->M_beta += *fluid/*M_un*/;//relative velocity beta=un-w

        //M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap));

        assembleFluidBlock(iter, M_un);
        assembleSolidBlock(iter, M_un);

        applyBoundaryConditions();
    }
    super::evalResidual( disp,  M_rhsFull, res, M_diagonalScale);
}


void MonolithicGE::applyBoundaryConditions( )
{

    if ( !M_BCh_u->bdUpdateDone() )
        M_BCh_u->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    M_BCh_d->setOffset(M_offset);
    if ( !M_BCh_d->bdUpdateDone() )
        M_BCh_d->bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );

    M_monolithicMatrix->setRobin( M_robinCoupling, M_rhsFull );
    M_precPtr->setRobin(M_robinCoupling, M_rhsFull);

    if (!this->M_monolithicMatrix->set())
    {
        M_BChs.push_back(M_BCh_d);
        M_BChs.push_back(M_BCh_u);
        M_FESpaces.push_back(M_dFESpace);
        M_FESpaces.push_back(M_uFESpace);

        M_monolithicMatrix->push_back_matrix(M_solidBlockPrec, false);
        M_monolithicMatrix->push_back_matrix(M_fluidBlock, true);
        M_monolithicMatrix->setConditions(M_BChs);
        M_monolithicMatrix->setSpaces(M_FESpaces);
        M_monolithicMatrix->setOffsets(2, M_offset, 0);
        M_monolithicMatrix->coupler(M_monolithicMap, M_dofStructureToHarmonicExtension->locDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->getTimeStep());
    }
    else
    {
        M_monolithicMatrix->replace_matrix(M_fluidBlock, 1);
        M_monolithicMatrix->replace_matrix(M_solidBlockPrec, 0);
    }

    M_monolithicMatrix->blockAssembling();
    M_monolithicMatrix->applyBoundaryConditions(dataFluid()->dataTime()->getTime(), M_rhsFull);

    M_monolithicMatrix->GlobalAssemble();
    M_monolithicMatrix->getMatrix()->spy("M");

}

void MonolithicGE::setupDOF()
{
    M_bcvStructureDispToHarmonicExtension.reset( new  BCVectorInterface );
    super::setupDOF();
}

void
MonolithicGE::setupSystem( )
{
    super::setupSystem();
    M_meshMotion->setUp( M_dataFile );
}

void
MonolithicGE::updateSystem( )
{
    super::updateSystem();

    // Set displacement for solid RHS
    setDispSolid( *M_un );
}


void
MonolithicGE::iterateMesh(const vector_type& disp)
{
    vector_type lambdaFluid(*M_interfaceMap, Unique);

    monolithicToInterface(lambdaFluid, disp);

    lambdaFluid *= (M_data->dataFluid()->dataTime()->getTimeStep()*(M_solid->rescaleFactor()));

    this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

    M_meshMotion->iterate(*M_BCh_mesh);

}

bool MonolithicGE::reg = FSIFactory::instance().registerProduct( "monolithicGE", &MonolithicGE::createM )  &&
                         BlockPrecFactory::instance().registerProduct("ComposedDNND"  , &ComposedDNND::createComposedDNND) &&
                         BlockPrecFactory::instance().registerProduct("AdditiveSchwarz"  , &BlockMatrix::createAdditiveSchwarz) &&
                         BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarz"  , &BlockMatrix::createAdditiveSchwarz ) &&
                         BlockPrecFactory::instance().registerProduct("AdditiveSchwarzRN"  , &BlockMatrixRN::createAdditiveSchwarzRN ) &&
                         BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarzRN"  , &BlockMatrixRN::createAdditiveSchwarzRN ) &&
                         BlockPrecFactory::instance().registerProduct("ComposedDN"  , &ComposedDN::createComposedDN ) &&
                         BlockPrecFactory::instance().registerProduct("ComposedDN2"  , &ComposedDN::createComposedDN2 );

} // Namespace LifeV
