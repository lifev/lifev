/* -*- mode: c++ -*-

   This file is part of the LifeV library

   Author(s): Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   Date: 2008-09-17

   Copyright (C) 2008

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
   \file (monolithicGI)
   \author crosetto <Paolo Crosetto>
   \date (17/09/2008)
*/

#include <life/lifesolver/VenantKirchhofSolver.hpp>
//#include <life/lifesolver/NonLinearVenantKirchhofSolver.hpp>

#include <lifemc/lifesolver/MonolithicGI.hpp>
#include <lifemc/lifesolver/ComposedDN.hpp>
#include <lifemc/lifesolver/ComposedDND.hpp>
#include <lifemc/lifesolver/BlockMatrixRN.hpp>

namespace LifeV
{

MonolithicGI::MonolithicGI():
    super(),
    M_mapWithoutMesh(),
    M_uk(),
    M_domainVelImplicit(true),
    M_convectiveTermDer(true),
    M_interface(0),
    M_meshBlock(),
    M_shapeDerivativesBlock()
{}

void
MonolithicGI::setupFluidSolid( UInt const fluxes )
{
    super::setupFluidSolid( fluxes );
    UInt offset = M_monolithicMap->getMap(Unique)->NumGlobalElements();
    M_mapWithoutMesh.reset(new EpetraMap(*M_monolithicMap));

    *this->M_monolithicMap += this->M_mmFESpace->map();
    /* OBSOLETE
       if(M_data->dataFluid()->useShapeDerivatives())
       {
       M_epetraOper.reset( new Epetra_FullMonolithic(this));
       M_solid->setOperator(*M_epetraOper);
       }*/
    //std::cout<<"map global elements : "<<M_monolithicMap->getMap(Unique)->NumGlobalElements()<<std::endl;
    M_interface=M_monolithicMatrix->getInterface();

    vector_type u0(*this->M_monolithicMap);
    M_bdf.reset(new BdfT<vector_type>(M_data->dataFluid()->dataTime()->getBDF_order()));
    M_bdf->initialize_unk(u0);
    this->M_rhs.reset(new vector_type(*this->M_monolithicMap));
    this->M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
    if(M_data->dataFluid()->useShapeDerivatives())
        M_shapeDerivativesBlock.reset(new matrix_type(*M_monolithicMap));
    M_uk.reset (new vector_type(*this->M_monolithicMap));
    M_un.reset (new vector_type(*this->M_monolithicMap));

    M_meshMotion.reset(new meshmotion_raw_type(*M_mmFESpace,
                                               M_epetraComm,
                                               *M_monolithicMap,
                                               offset));
    M_fluid.reset     (new fluid_raw_type(M_data->dataFluid(),
                                          *M_uFESpace,
                                          *M_pFESpace,
                                          *M_mmFESpace,
                                          M_epetraComm,
                                          *M_monolithicMap));
    M_solid.reset(solid_raw_type::StructureSolverFactory::instance().createObject( M_data->dataSolid()->solidType( ) ));
    M_solid->setup(M_data->dataSolid(),
                    M_dFESpace,
                    M_epetraComm,
                    M_monolithicMap,
                    M_offset
                   );
}

void
MonolithicGI::updateSystem()
{
    //M_meshMotion->dispOld() is at time n-1 !!
    //M_meshMotion->updateSystem();

    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    vector_ptrtype meshDispDiff(new vector_type(M_mmFESpace->map()));
    meshDispDiff->subset(*M_uk, offset); //if the conv. term is to be condidered implicitly
    M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.`
    super::updateSystem();
    M_un.reset(new vector_type(*M_uk));
}

void
MonolithicGI::buildSystem ()
{
    super::buildSystem();
    M_meshMotion->computeMatrix();
}

void
MonolithicGI::evalResidual( vector_type&       res,
                            const vector_type& disp,
                            const UInt          iter )
{
    setDispSolid(disp);
    if(iter > 0)
    {
        this->M_solid->updateVel();
    }
    M_uk.reset(new vector_type( disp ));
    this->M_beta.reset( new vector_type(M_uFESpace->map()) );
    UInt offset( M_solidAndFluidDim + nDimensions*M_interface );

    vector_ptrtype meshDispDiff( new vector_type(M_mmFESpace->map()) );
    vector_ptrtype meshDispOld( new vector_type(M_mmFESpace->map()) );

    meshDispDiff->subset(disp, offset); //if the conv. term is to be condidered implicitly

    meshDispOld->subset(*M_un, offset);

    //meshDispDiff->subset(*M_uk, offset); //if the mesh motion is at the previous nonlinear step (FP) in the convective term
    //meshDispDiff->subset(*M_un, offset); //if we linearize in a semi-implicit way
    M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.
    double alpha = 1/M_data->dataFluid()->dataTime()->getTimeStep();
    vector_type mmRep(*meshDispDiff, Repeated);// just to repeat dispDiff. No way witout copying?
    this->moveMesh(mmRep);// re-initialize the mesh points
    *meshDispDiff -= *meshDispOld;//relative displacement
    if(!M_domainVelImplicit)
    {
        meshDispDiff=meshDispOld;// at time n /*->subset(*M_un, offset)*/; //if the mesh motion is at the previous time step in the convective term
        *meshDispDiff -= M_meshMotion->dispOld();//at time n-1
    }
    *meshDispDiff *= -alpha;// -w, mesh velocity
    mmRep = *meshDispDiff;

    this->interpolateVelocity(mmRep, *this->M_beta);
    //            *this->M_beta *= -alpha; //HE solution scaled!
    vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
    if(!M_convectiveTermDer)
        fluid->subset(*M_un/**M_unOld*/, 0);
    else
        fluid->subset(disp, 0);
    *this->M_beta += *fluid/*M_un or disp, it could be also M_uk in a FP strategy*/;

//      if(iter == 0)
//      {
// //         M_solid->updateSystem();
//      }
//      else
//      {
//          //         M_solid->computeMatrix( disp, 1.);
//      }

    assembleSolidBlock( iter, M_uk );
    assembleFluidBlock( iter, M_uk );
    assembleMeshBlock ( iter );

    applyBoundaryConditions();

    super::evalResidual( disp, M_rhsFull, res, false );
}

void
MonolithicGI::applyBoundaryConditions()
{
    M_monolithicMatrix->setRobin( M_robinCoupling, M_rhsFull );
    M_precPtr->setRobin(M_robinCoupling, M_rhsFull);

    if(!M_monolithicMatrix->set())
    {
        M_BChs.push_back(M_BCh_d);
        M_BChs.push_back(M_BCh_u);
        M_FESpaces.push_back(M_dFESpace);
        M_FESpaces.push_back(M_uFESpace);

        M_BChs.push_back(M_BCh_mesh);
        M_FESpaces.push_back(M_mmFESpace);

        M_monolithicMatrix->push_back_matrix(M_solidBlockPrec, false);
        M_monolithicMatrix->push_back_matrix(M_fluidBlock, true);
        M_monolithicMatrix->push_back_matrix(M_meshBlock, false);
        M_monolithicMatrix->setConditions(M_BChs);
        M_monolithicMatrix->setSpaces(M_FESpaces);
        M_monolithicMatrix->setOffsets(3, M_offset, 0, M_solidAndFluidDim + nDimensions*M_interface);
        M_monolithicMatrix->coupler(M_monolithicMap, M_dofStructureToHarmonicExtension->locDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->getTimeStep());
    }
    else
    {
        M_monolithicMatrix->replace_matrix(M_solidBlockPrec, 0);
        M_monolithicMatrix->replace_matrix(M_fluidBlock, 1);
        M_monolithicMatrix->replace_matrix(M_meshBlock, 2);
    }

    M_monolithicMatrix->blockAssembling();

    if ( !M_BCh_u->bdUpdateDone() )
        M_BCh_u->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    M_BCh_d->setOffset(M_offset);
    if ( !M_BCh_d->bdUpdateDone() )
        M_BCh_d->bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    M_BCh_mesh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);
    if ( !M_BCh_mesh->bdUpdateDone() )
        M_BCh_mesh->bdUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

    M_monolithicMatrix->applyBoundaryConditions(dataFluid()->dataTime()->getTime(), M_rhsFull);
    M_monolithicMatrix->GlobalAssemble();
    //M_monolithicMatrix->getMatrix()->spy("FM");
}


int MonolithicGI::setupBlockPrec( )
{
    super::setupBlockPrec( );

    if(M_data->dataFluid()->useShapeDerivatives())
    {
        *M_shapeDerivativesBlock *= 0.;
        M_shapeDerivativesBlock->openCrsMatrix( );
        shapeDerivatives( M_shapeDerivativesBlock ,*M_uk/*subX*/, M_domainVelImplicit, M_convectiveTermDer );
        //*M_shapeDerivativesBlock += *M_monolithicMatrix->getMatrix();
        M_shapeDerivativesBlock->GlobalAssemble( );
        M_monolithicMatrix->addToGlobalMatrix( M_shapeDerivativesBlock );
    }

    //M_solidDerBlock = M_solidBlockPrec; // an inexact Newton approximation of the Jacobian

    //The following part accounts for a possibly nonlinear structure model, should not be run when linear
    //elasticity is used
    if( M_data->dataSolid()->useExactJacobian() )
    {
        M_solid->updateJacobian( *M_uk, M_solidDerBlock ); // computing the derivatives if nonlinear (comment this for inexact Newton);
        *M_monolithicMatrix->getMatrix() *= 0;
        // doing nothing if linear
        M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
        *M_solidBlockPrec += *M_solidDerBlock;
        M_solidBlockPrec->GlobalAssemble();
        M_precPtr->replace_matrix( M_solidBlockPrec, 0 );
        M_monolithicMatrix->blockAssembling();
        M_monolithicMatrix->applyBoundaryConditions( dataFluid()->dataTime()->getTime());
        M_monolithicMatrix->GlobalAssemble();
        *M_monolithicMatrix->getMatrix() *= M_data->dataFluid()->dataTime()->getTimeStep();
    }

    if( M_precPtr->getBlockVector( ).size( )<3 )
    {
        M_precPtr->push_back_matrix( M_meshBlock, false );
        M_precPtr->setConditions( M_BChs );
        M_precPtr->setSpaces( M_FESpaces );
        M_precPtr->setOffsets( 3, M_offset, 0,  M_solidAndFluidDim + nDimensions*M_interface );
        M_precPtr->coupler( M_monolithicMap, M_dofStructureToHarmonicExtension->locDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->getTimeStep(), 2 );

        if(M_data->dataFluid()->useShapeDerivatives())
        {
            M_precPtr->push_back_coupling( M_shapeDerivativesBlock );
        }
    }
     else
    {
        //M_precPtr->replace_matrix( M_solidBlockPrec, 0 );
        //M_precPtr->replace_matrix( M_fluidBlock, 1 );
        M_precPtr->replace_matrix( M_meshBlock, 2 );

        if(M_data->dataFluid()->useShapeDerivatives())
        {
            M_precPtr->replace_coupling( M_shapeDerivativesBlock, 2 );
        }
    }
}

void MonolithicGI::solveJac(vector_type       &_step,
                            const vector_type &_res,
                            const Real       _linearRelTol)
{

    setupBlockPrec( );

    M_precPtr->blockAssembling( );
    M_precPtr->applyBoundaryConditions( dataFluid()->dataTime()->getTime() );
    M_precPtr->GlobalAssemble( );

    //boost::dynamic_pointer_cast<BlockMatrix>(M_precPtr)->getMatrix()->spy("P");

    M_linearSolver->setMatrix(*M_monolithicMatrix->getMatrix());

    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _res.NormInf(), "\n");
    M_solid->getDisplayer().leaderPrint("  M-  Solving Jacobian system ...              \n" );

    this->iterateMonolithic(_res, _step);
    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _step.NormInf(), "\n");
}

void MonolithicGI::initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                               FSIOperator::solid_type::value_type::Function const& p0,
                               FSIOperator::solid_type::value_type::Function const& d0,
                               FSIOperator::solid_type::value_type::Function const& w0,
                               FSIOperator::solid_type::value_type::Function const& df0 )
{
    super::initialize(u0, p0, d0, w0, df0);

    vector_type df(M_mmFESpace->map());
    M_mmFESpace->interpolate(df0, df, M_data->dataSolid()->dataTime()->getTime());
    M_un->add(df, M_solidAndFluidDim+getDimInterface());
    M_meshMotion->setDisplacement(df);
}

void MonolithicGI::shapeDerivatives(matrix_ptrtype sdMatrix, const vector_type& sol, bool domainVelImplicit, bool convectiveTermDer)
{
    double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();
    vector_ptrtype rhsNew(new vector_type(*M_monolithicMap));
    vector_type un(M_uFESpace->map()/*+M_pFESpace->map()*/);
    vector_type uk(M_uFESpace->map()+M_pFESpace->map());

    //vector_type meshVel(M_meshMotion->dispDiff(), Repeated);
    vector_ptrtype meshVel(new vector_type(M_mmFESpace->map()));

    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    if(domainVelImplicit)
    {
        vector_type meshDispOld(M_mmFESpace->map());
        meshVel->subset(sol, offset); //if the conv. term is to be condidered implicitly
        meshDispOld.subset(*M_un, offset);
        *meshVel -= meshDispOld;
    }
    else
    {
        meshVel->subset(*M_un, offset); //if the conv. term is to be condidered partly explicitly
        *meshVel -= M_meshMotion->dispOld();
    }

    if(convectiveTermDer)
        un.subset(sol, 0);
    else
        un.subset(*this->M_un, 0);

    *meshVel *= alpha;
    vector_ptrtype meshVelRep(new vector_type(M_mmFESpace->map(), Repeated));
    *meshVelRep = *meshVel;

    uk.subset(sol, 0);
    vector_type dvfm(M_uFESpace->map(), Repeated);
    vector_type vfm(M_uFESpace->map(), Repeated);
    this->transferMeshMotionOnFluid(*meshVelRep, vfm);

    M_fluid->updateShapeDerivatives(*sdMatrix,
                                    alpha,
                                    un,//un if !domainVelImplicit, otherwise uk
                                    uk,//uk
                                    vfm, //(xk-xn)/dt (FI), or (xn-xn-1)/dt (CE)//Repeated
                                    M_solidAndFluidDim+M_interface*nDimensions,
                                    *M_uFESpace,
                                    domainVelImplicit,
                                    convectiveTermDer
                                    );
}

void
MonolithicGI::setUp( const GetPot& dataFile )
{
    super::setUp(dataFile);

    M_domainVelImplicit     = dataFile( "fluid/domainVelImplicit", true);
    M_convectiveTermDer     = dataFile( "fluid/convectiveTermDer", false);
}


void
MonolithicGI::assembleMeshBlock(UInt iter)
{
    M_meshBlock.reset(new matrix_type(*M_monolithicMap));
    M_meshMotion->setMatrix(M_meshBlock);
    M_meshBlock->GlobalAssemble();
    UInt offset(M_solidAndFluidDim+nDimensions*M_interface);
    std::map<ID, ID>::const_iterator ITrow;
    std::map<ID, ID> locdofmap(M_dofStructureToHarmonicExtension->locDofMap());

    /******************alternative way************************/
//     BCFunctionBase bcf(fZero);
//     fluid_bchandler_type BCh(new fluid_bchandler_raw_type() );
//     BCh->addBC("Interface", 1, Essential, Full,
//                bcf, 3);

//     BCh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);

//     if ( !BCh->bdUpdateDone() )
//         BCh->bdUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

//     bcManage( *M_meshBlock, *M_rhsFull, *M_mmFESpace->mesh(), M_mmFESpace->dof(), *BCh, M_mmFESpace->feBd(), 1., dataFluid().dataTime()->getTime());
/********************************************************/

    for ( ID dim=0; dim < nDimensions; ++dim )
        for( ITrow = locdofmap.begin(); ITrow != locdofmap.end(); ++ITrow )
        {
            UInt i = ITrow->first;
            M_meshBlock->diagonalize(i+offset+dim*M_mmFESpace->dof().numTotalDof()-1 , 1.);
        }
}

namespace
{

BlockMatrix*    createAdditiveSchwarzGI(){
return new BlockMatrix(31);
}

BlockMatrix*    createAdditiveSchwarzRNGI(){
return new BlockMatrixRN(31);
}

BlockInterface*    createComposedDNGI(){
    const ComposedBlockOper::Block order[] = {  ComposedBlockOper::solid, ComposedBlockOper::fluid, ComposedBlockOper::mesh };
    const Int couplingsDNGI[] = { 0, 7, 16 };
    const std::vector<Int> couplingVectorDNGI(couplingsDNGI, couplingsDNGI+3);
    const std::vector<ComposedBlockOper::Block> orderVector(order, order+3);
    return new ComposedDN( couplingVectorDNGI, orderVector );
}

BlockInterface*    createComposedDN2GI(){
    const ComposedBlockOper::Block order[] = { ComposedBlockOper::fluid, ComposedBlockOper::solid, ComposedBlockOper::mesh };
    const Int couplingsDN2GI[] = { 8, 6, 16 };
    const std::vector<Int> couplingVectorDN2GI(couplingsDN2GI, couplingsDN2GI+3);
    const std::vector<ComposedBlockOper::Block> orderVector(order, order+3);
    return new ComposedDN( couplingVectorDN2GI, orderVector );
}

BlockInterface*    createComposedDNDGI(){
    const ComposedBlockOper::Block order[] = {  ComposedBlockOper::mesh, ComposedBlockOper::solid, ComposedBlockOper::fluid };
    const Int couplingsDNGI2[] = { 0, 7, 0 };
    const std::vector<Int> couplingVectorDNGI2(couplingsDNGI2, couplingsDNGI2+3);
    const std::vector<ComposedBlockOper::Block> orderVector(order, order+3);
    return new ComposedDND( couplingVectorDNGI2, orderVector );
}

BlockInterface*    createComposedDND2GI(){
    const ComposedBlockOper::Block order[] = { ComposedBlockOper::mesh, ComposedBlockOper::fluid , ComposedBlockOper::solid};
const Int couplingsDN2GI2[] = { 8, 6, 0 };
const std::vector<Int> couplingVectorDN2GI2(couplingsDN2GI2, couplingsDN2GI2+3);
const std::vector<ComposedBlockOper::Block> orderVector(order, order+3);
return new ComposedDND( couplingVectorDN2GI2, orderVector );
}
FSIOperator*    createFM(){ return new MonolithicGI(); }
}

bool MonolithicGI::reg =  BlockPrecFactory::instance().registerProduct("AdditiveSchwarzGI"  , &createAdditiveSchwarzGI )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDNGI"  , &createComposedDNGI )
    &&
    BlockMatrix::Factory::instance().registerProduct( "AdditiveSchwarzGI", &createAdditiveSchwarzGI )
    &&
    BlockMatrix::Factory::instance().registerProduct( "AdditiveSchwarzRNGI", &createAdditiveSchwarzRNGI )
    &&
    FSIFactory::instance().registerProduct( "monolithicGI", &createFM )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDNDGI"  , &createComposedDNDGI )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDND2GI"  , &createComposedDND2GI )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDN2GI"  , &createComposedDN2GI );

}
