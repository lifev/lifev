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

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================


void MonolithicGE::setupFluidSolid()
{
    super::setupFluidSolid( );
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
MonolithicGE::assembleFluidBlock(UInt iter)
{
    double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();//mesh velocity w
    matrix_ptrtype newMatrix(new matrix_type(*M_monolithicMap));
    M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, newMatrix );
    this->M_fluid->updateStab(*newMatrix);
    newMatrix->GlobalAssemble();

    M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
    *M_fluidBlock += *newMatrix;


        if(iter==0)
        {
            M_nbEval = 0; // new time step
            M_resetPrec=true;
            *this->M_rhs               += M_fluid->matrMass()*M_bdf->time_der( M_data->dataFluid()->dataTime()->getTimeStep() );
            couplingRhs(this->M_rhs, M_un);

            if (!M_restarts)
            {
                this->M_solid->updateVel();
                M_restarts = false;
            }
            updateSolidSystem(this->M_rhs);
        }

        *M_rhsFull = *M_rhs;

}

namespace
{

BlockInterface* createComposedDNND(){
    const Int couplingsDNND[] = { 8, 4, 2, 8, 1, 2 };
    const ComposedBlockOper::Block order[] = { ComposedBlockOper::fluid, ComposedBlockOper::solid};
    const std::vector<Int> couplingVectorDNND(couplingsDNND, couplingsDNND+6);
    const std::vector<ComposedBlockOper::Block> orderVector(order, order+6);
    return new ComposedDNND(couplingVectorDNND, orderVector);
}

BlockInterface* createComposedNN()
{
    const ComposedBlockOper::Block order[] = {  ComposedBlockOper::fluid, ComposedBlockOper::solid};
    const Int couplingsNN[] = { 8, 4, 1, 2};
    const std::vector<Int> couplingVectorNN(couplingsNN, couplingsNN+4);
    const std::vector<ComposedBlockOper::Block> orderVector(order, order+4);
    return new ComposedNN( couplingVectorNN, orderVector );
}

BlockMatrix*    createAdditiveSchwarz()
{
return new BlockMatrix(15);
}

BlockMatrix*    createAdditiveSchwarzRN()
{
return new BlockMatrixRN(15);
}

BlockInterface* createComposedDN()
{
    const ComposedBlockOper::Block order[] = { ComposedBlockOper::solid, ComposedBlockOper::fluid};
    const Int couplingsDN[] = { 0, 7};
    const std::vector<Int> couplingVectorDN(couplingsDN, couplingsDN+2);
    const std::vector<ComposedBlockOper::Block> orderVector(order, order+2);
    return new ComposedDN(couplingVectorDN, orderVector);
}

BlockInterface* createComposedDN2()
{
    const ComposedBlockOper::Block order[] = { ComposedBlockOper::fluid, ComposedBlockOper::solid};
    const Int couplingsDN2[] = { 8, 6};
    const std::vector<Int> couplingVectorDN2(couplingsDN2, couplingsDN2+2);
    const std::vector<ComposedBlockOper::Block> orderVector(order, order+2);
    return new ComposedDN(couplingVectorDN2, orderVector);
}

FSIOperator* createM(){ return new MonolithicGE(); }
}

bool MonolithicGE::reg = FSIFactory::instance().registerProduct( "monolithicGE", &createM ) &&
    BlockPrecFactory::instance().registerProduct("ComposedDNND"  , &createComposedDNND) &&
    BlockPrecFactory::instance().registerProduct("AdditiveSchwarz"  , &createAdditiveSchwarz) &&
    BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarz"  , &createAdditiveSchwarz ) &&
    BlockPrecFactory::instance().registerProduct("AdditiveSchwarzRN"  , &createAdditiveSchwarzRN ) &&
    BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarzRN"  , &createAdditiveSchwarzRN ) &&
    BlockPrecFactory::instance().registerProduct("ComposedDN"  , &createComposedDN ) &&
    BlockPrecFactory::instance().registerProduct("ComposedDN2"  , &createComposedDN2 );

} // Namespace LifeV
