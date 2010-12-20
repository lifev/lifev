/* -*- mode: c++ -*- */
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

#include <life/lifecore/life.hpp>

#include <life/lifesolver/VenantKirchhoffSolver.hpp>

#include <lifemc/lifesolver/FSIMonolithicGI.hpp>
#include <lifemc/lifesolver/MonolithicBlockComposedDN.hpp>
#include <lifemc/lifesolver/MonolithicBlockComposedDND.hpp>
#include <lifemc/lifesolver/MonolithicBlockMatrixRN.hpp>


// ===================================================
//! Constructors and Descructor
// ===================================================


namespace LifeV
{

FSIMonolithicGI::FSIMonolithicGI():
        super_Type(),
        M_mapWithoutMesh(),
        M_uk(),
        M_domainVelImplicit(true),
        M_convectiveTermDer(true),
        M_interface(0),
        M_meshBlock(),
        M_shapeDerivativesBlock()
{}


// ===================================================
//! Public Methods
// ===================================================


void
FSIMonolithicGI::setUp( const GetPot& dataFile )
{
    super_Type::setUp(dataFile);

    M_domainVelImplicit     = dataFile( "fluid/domainVelImplicit", true);
    M_convectiveTermDer     = dataFile( "fluid/convectiveTermDer", false);
}

void
FSIMonolithicGI::setupFluidSolid( UInt const fluxes )
{
    super_Type::setupFluidSolid( fluxes );
    UInt offset = M_monolithicMap->map(Unique)->NumGlobalElements();
    M_mapWithoutMesh.reset(new EpetraMap(*M_monolithicMap));

    *this->M_monolithicMap += this->M_mmFESpace->map();
    /* OBSOLETE
       if(M_data->dataFluid()->useShapeDerivatives())
       {
       M_epetraOper.reset( new Epetra_FullMonolithic(this));
       M_solid->setOperator(*M_epetraOper);
       }*/
    //std::cout<<"map global elements : "<<M_monolithicMap->getMap(Unique)->NumGlobalElements()<<std::endl;
    M_interface=M_monolithicMatrix->interface();

    vector_Type u0(*this->M_monolithicMap);
    M_bdf.reset(new BdfT<vector_Type>());
    M_bdf->setup(M_data->dataFluid()->dataTime()->orderBDF());
    M_bdf->setInitialCondition(u0);
    this->M_rhs.reset(new vector_Type(*this->M_monolithicMap));
    this->M_rhsFull.reset(new vector_Type(*this->M_monolithicMap));
    if(M_data->dataFluid()->useShapeDerivatives())
        M_shapeDerivativesBlock.reset(new matrix_Type(*M_monolithicMap));
    M_uk.reset (new vector_Type(*this->M_monolithicMap));
    M_un.reset (new vector_Type(*this->M_monolithicMap));

    M_meshMotion.reset(new meshMotion_Type(*M_mmFESpace,
                                               M_epetraComm,
                                               *M_monolithicMap,
                                               offset));
    M_fluid.reset     (new fluid_Type(M_data->dataFluid(),
                                          *M_uFESpace,
                                          *M_pFESpace,
                                          *M_mmFESpace,
                                          M_epetraComm,
                                          *M_monolithicMap));
    M_solid.reset(solid_Type::StructureSolverFactory::instance().createObject( M_data->dataSolid()->getSolidType( ) ));

    M_solid->setup(M_data->dataSolid(),
                   M_dFESpace,
                   M_epetraComm,
                   M_monolithicMap,
                   M_offset
                  );
}

void
FSIMonolithicGI::updateSystem()
{
    //M_meshMotion->dispOld() is at time n-1 !!
    //M_meshMotion->updateSystem();

    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    vectorPtr_Type meshDispDiff(new vector_Type(M_mmFESpace->map()));
    meshDispDiff->subset(*M_uk, offset); //if the conv. term is to be condidered implicitly
    M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.`
    super_Type::updateSystem();
    M_un.reset(new vector_Type(*M_uk));
}

void
FSIMonolithicGI::buildSystem ()
{
    super_Type::buildSystem();
    M_meshMotion->computeMatrix();
}

void
FSIMonolithicGI::evalResidual( vector_Type&       res,
                            const vector_Type& disp,
                            const UInt          iter )
{
    setDispSolid(disp);
    if (iter > 0)
    {
        this->M_solid->updateVel();
    }
    M_uk.reset(new vector_Type( disp ));
    this->M_beta.reset( new vector_Type(M_uFESpace->map()) );
    UInt offset( M_solidAndFluidDim + nDimensions*M_interface );

    vectorPtr_Type meshDispDiff( new vector_Type(M_mmFESpace->map()) );
    vectorPtr_Type meshDispOld( new vector_Type(M_mmFESpace->map()) );

    meshDispDiff->subset(disp, offset); //if the conv. term is to be condidered implicitly

    meshDispOld->subset(*M_un, offset);

    //meshDispDiff->subset(*M_uk, offset); //if the mesh motion is at the previous nonlinear step (FP) in the convective term
    //meshDispDiff->subset(*M_un, offset); //if we linearize in a semi-implicit way
    M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.
    double alpha = 1/M_data->dataFluid()->dataTime()->timeStep();
    vector_Type mmRep(*meshDispDiff, Repeated);// just to repeat dispDiff. No way witout copying?
    this->moveMesh(mmRep);// re-initialize the mesh points
    *meshDispDiff -= *meshDispOld;//relative displacement
    if (!M_domainVelImplicit)
    {
        meshDispDiff=meshDispOld;// at time n /*->subset(*M_un, offset)*/; //if the mesh motion is at the previous time step in the convective term
        *meshDispDiff -= M_meshMotion->dispOld();//at time n-1
    }
    *meshDispDiff *= -alpha;// -w, mesh velocity
    mmRep = *meshDispDiff;

    this->interpolateVelocity(mmRep, *this->M_beta);
    //            *this->M_beta *= -alpha; //HE solution scaled!
    vectorPtr_Type fluid(new vector_Type(this->M_uFESpace->map()));
    if (!M_convectiveTermDer)
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

    super_Type::evalResidual( disp, M_rhsFull, res, false );
}

void
FSIMonolithicGI::applyBoundaryConditions()
{
    M_monolithicMatrix->setRobin( M_robinCoupling, M_rhsFull );
    M_precPtr->setRobin(M_robinCoupling, M_rhsFull);

    if (!M_monolithicMatrix->set())
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
        M_monolithicMatrix->coupler(M_monolithicMap, M_dofStructureToHarmonicExtension->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep());
    }
    else
    {
        M_monolithicMatrix->replace_matrix(M_solidBlockPrec, 0);
        M_monolithicMatrix->replace_matrix(M_fluidBlock, 1);
        M_monolithicMatrix->replace_matrix(M_meshBlock, 2);
    }

    M_monolithicMatrix->blockAssembling();

    if ( !M_BCh_u->bcUpdateDone() )
        M_BCh_u->bcUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    M_BCh_d->setOffset(M_offset);
    if ( !M_BCh_d->bcUpdateDone() )
        M_BCh_d->bcUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    M_BCh_mesh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);
    if ( !M_BCh_mesh->bcUpdateDone() )
        M_BCh_mesh->bcUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

    M_monolithicMatrix->applyBoundaryConditions(dataFluid()->dataTime()->time(), M_rhsFull);
    M_monolithicMatrix->GlobalAssemble();
    //M_monolithicMatrix->matrix()->spy("FM");
}


void FSIMonolithicGI::solveJac(vector_Type       &_step,
                            const vector_Type &_res,
                            const Real       /*_linearRelTol*/)
{

    setupBlockPrec( );

    M_precPtr->blockAssembling( );
    M_precPtr->applyBoundaryConditions( dataFluid()->dataTime()->time() );
    M_precPtr->GlobalAssemble( );

    //boost::dynamic_pointer_cast<MonolithicBlockMatrix>(M_precPtr)->matrix()->spy("P");

    M_linearSolver->setMatrix(*M_monolithicMatrix->matrix());

    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _res.normInf(), "\n");
    M_solid->getDisplayer().leaderPrint("  M-  Solving Jacobian system ...              \n" );

    this->iterateMonolithic(_res, _step);
    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _step.normInf(), "\n");
}

void FSIMonolithicGI::initialize( FSI::fluidPtr_Type::value_type::function_Type const& u0,
                               FSI::solidPtr_Type::value_type::Function const& p0,
                               FSI::solidPtr_Type::value_type::Function const& d0,
                               FSI::solidPtr_Type::value_type::Function const& w0,
                               FSI::solidPtr_Type::value_type::Function const& df0 )
{
    super_Type::initialize(u0, p0, d0, w0, df0);
    vector_Type df(M_mmFESpace->map());
    M_mmFESpace->interpolate(df0, df, M_data->dataSolid()->getDataTime()->time());

    M_un->add(df, M_solidAndFluidDim+getDimInterface());
    M_meshMotion->setDisplacement(df);
}


int FSIMonolithicGI::setupBlockPrec( )
{
    super_Type::setupBlockPrec( );

    if (M_data->dataFluid()->useShapeDerivatives())
    {
        *M_shapeDerivativesBlock *= 0.;
        M_shapeDerivativesBlock->openCrsMatrix( );
        shapeDerivatives( M_shapeDerivativesBlock ,*M_uk/*subX*/, M_domainVelImplicit, M_convectiveTermDer );

        M_shapeDerivativesBlock->globalAssemble( );
        M_monolithicMatrix->addToGlobalMatrix( M_shapeDerivativesBlock );
    }

    //M_solidDerBlock = M_solidBlockPrec; // an inexact Newton approximation of the Jacobian

    //The following part accounts for a possibly nonlinear structure model, should not be run when linear
    //elasticity is used
    if ( M_data->dataSolid()->getUseExactJacobian() )
    {
        M_solid->updateJacobian( *M_uk, M_solidDerBlock ); // computing the derivatives if nonlinear (comment this for inexact Newton);
        *M_monolithicMatrix->matrix() *= 0;
        // doing nothing if linear
        M_solidBlockPrec.reset(new matrix_Type(*M_monolithicMap, 1));
        *M_solidBlockPrec += *M_solidDerBlock;
        M_solidBlockPrec->globalAssemble();
        M_precPtr->replace_matrix( M_solidBlockPrec, 0 );
        M_monolithicMatrix->blockAssembling();
        M_monolithicMatrix->applyBoundaryConditions( dataFluid()->dataTime()->time());
        M_monolithicMatrix->GlobalAssemble();
        *M_monolithicMatrix->matrix() *= M_data->dataFluid()->dataTime()->timeStep();
    }

    if ( M_precPtr->blockVector( ).size( )<3 )
    {
        M_precPtr->push_back_matrix( M_meshBlock, false );
        M_precPtr->setConditions( M_BChs );
        M_precPtr->setSpaces( M_FESpaces );
        M_precPtr->setOffsets( 3, M_offset, 0,  M_solidAndFluidDim + nDimensions*M_interface );
        M_precPtr->coupler( M_monolithicMap, M_dofStructureToHarmonicExtension->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), 2 );

        if (M_data->dataFluid()->useShapeDerivatives())
        {
            M_precPtr->push_back_coupling( M_shapeDerivativesBlock );
        }
    }
    else
    {
        //M_precPtr->replace_matrix( M_solidBlockPrec, 0 );
        //M_precPtr->replace_matrix( M_fluidBlock, 1 );
        M_precPtr->replace_matrix( M_meshBlock, 2 );

        if (M_data->dataFluid()->useShapeDerivatives())
        {
            M_precPtr->replace_coupling( M_shapeDerivativesBlock, 2 );
        }
    }
}


void FSIMonolithicGI::shapeDerivatives(matrixPtr_Type sdMatrix, const vector_Type& sol, bool domainVelImplicit, bool convectiveTermDer)
{
    double alpha = 1./M_data->dataFluid()->dataTime()->timeStep();
    vectorPtr_Type rhsNew(new vector_Type(*M_monolithicMap));
    vector_Type un(M_uFESpace->map()/*+M_pFESpace->map()*/);
    vector_Type uk(M_uFESpace->map()+M_pFESpace->map());

    //vector_Type meshVel(M_meshMotion->dispDiff(), Repeated);
    vectorPtr_Type meshVel(new vector_Type(M_mmFESpace->map()));

    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    if (domainVelImplicit)
    {
        vector_Type meshDispOld(M_mmFESpace->map());
        meshVel->subset(sol, offset); //if the conv. term is to be condidered implicitly
        meshDispOld.subset(*M_un, offset);
        *meshVel -= meshDispOld;
    }
    else
    {
        meshVel->subset(*M_un, offset); //if the conv. term is to be condidered partly explicitly
        *meshVel -= M_meshMotion->dispOld();
    }

    if (convectiveTermDer)
        un.subset(sol, 0);
    else
        un.subset(*this->M_un, 0);

    *meshVel *= alpha;
    vectorPtr_Type meshVelRep(new vector_Type(M_mmFESpace->map(), Repeated));
    *meshVelRep = *meshVel;

    uk.subset(sol, 0);
    vector_Type dvfm(M_uFESpace->map(), Repeated);
    vector_Type vfm(M_uFESpace->map(), Repeated);
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
FSIMonolithicGI::assembleMeshBlock(UInt /*iter*/)
{
    M_meshBlock.reset(new matrix_Type(*M_monolithicMap));
    M_meshMotion->setMatrix(M_meshBlock);
    M_meshBlock->globalAssemble();
    UInt offset(M_solidAndFluidDim+nDimensions*M_interface);
    std::map<ID, ID>::const_iterator ITrow;
    std::map<ID, ID> locdofmap(M_dofStructureToHarmonicExtension->localDofMap());

    /******************alternative way************************/
//     BCFunctionBase bcf(fZero);
//     fluidBchandlerPtr_Type BCh(new fluid_bchandler_raw_type() );
//     BCh->addBC("Interface", 1, Essential, Full,
//                bcf, 3);

//     BCh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);

//     if ( !BCh->bcUpdateDone() )
//         BCh->bcUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

//     bcManage( *M_meshBlock, *M_rhsFull, *M_mmFESpace->mesh(), M_mmFESpace->dof(), *BCh, M_mmFESpace->feBd(), 1., dataFluid().dataTime()->time());
    /********************************************************/

    for ( ID dim=0; dim < nDimensions; ++dim )
        for ( ITrow = locdofmap.begin(); ITrow != locdofmap.end(); ++ITrow )
        {
            UInt i = ITrow->first;
            M_meshBlock->diagonalize(i+offset+dim*M_mmFESpace->dof().numTotalDof()-1 , 1.);
        }
}



// ===================================================
//! Factory methods
// ===================================================

namespace
{

MonolithicBlockMatrix*    createAdditiveSchwarzGI()
{
    return new MonolithicBlockMatrix(31);
}

MonolithicBlockMatrix*    createAdditiveSchwarzRNGI()
{
    return new MonolithicBlockMatrixRN(31);
}

MonolithicBlock*    createComposedDNGI()
{
    const MonolithicBlockComposed::Block order[] = {  MonolithicBlockComposed::solid, MonolithicBlockComposed::fluid, MonolithicBlockComposed::mesh };
    const Int couplingsDNGI[] = { 0, 7, 16 };
    const std::vector<Int> couplingVectorDNGI(couplingsDNGI, couplingsDNGI+3);
    const std::vector<MonolithicBlockComposed::Block> orderVector(order, order+3);
    return new MonolithicBlockComposedDN( couplingVectorDNGI, orderVector );
}

MonolithicBlock*    createComposedDN2GI()
{
    const MonolithicBlockComposed::Block order[] = { MonolithicBlockComposed::fluid, MonolithicBlockComposed::solid, MonolithicBlockComposed::mesh };
    const Int couplingsDN2GI[] = { 8, 6, 16 };
    const std::vector<Int> couplingVectorDN2GI(couplingsDN2GI, couplingsDN2GI+3);
    const std::vector<MonolithicBlockComposed::Block> orderVector(order, order+3);
    return new MonolithicBlockComposedDN( couplingVectorDN2GI, orderVector );
}

MonolithicBlock*    createComposedDNDGI()
{
    const MonolithicBlockComposed::Block order[] = {  MonolithicBlockComposed::mesh, MonolithicBlockComposed::solid, MonolithicBlockComposed::fluid };
    const Int couplingsDNGI2[] = { 0, 7, 0 };
    const std::vector<Int> couplingVectorDNGI2(couplingsDNGI2, couplingsDNGI2+3);
    const std::vector<MonolithicBlockComposed::Block> orderVector(order, order+3);
    return new MonolithicBlockComposedDND( couplingVectorDNGI2, orderVector );
}

MonolithicBlock*    createComposedDND2GI()
{
    const MonolithicBlockComposed::Block order[] = { MonolithicBlockComposed::mesh, MonolithicBlockComposed::fluid , MonolithicBlockComposed::solid};
    const Int couplingsDN2GI2[] = { 8, 6, 0 };
    const std::vector<Int> couplingVectorDN2GI2(couplingsDN2GI2, couplingsDN2GI2+3);
    const std::vector<MonolithicBlockComposed::Block> orderVector(order, order+3);
    return new MonolithicBlockComposedDND( couplingVectorDN2GI2, orderVector );
}
FSI*    createFM() { return new FSIMonolithicGI(); }
}

// ===================================================
//! Products registration
// ===================================================

bool FSIMonolithicGI::reg =  BlockPrecFactory::instance().registerProduct("AdditiveSchwarzGI"  , &createAdditiveSchwarzGI )
                          &&
                          BlockPrecFactory::instance().registerProduct("ComposedDNGI"  , &createComposedDNGI )
                          &&
                          MonolithicBlockMatrix::Factory::instance().registerProduct( "AdditiveSchwarzGI", &createAdditiveSchwarzGI )
                          &&
                          MonolithicBlockMatrix::Factory::instance().registerProduct( "AdditiveSchwarzRNGI", &createAdditiveSchwarzRNGI )
                          &&
                          FSIFactory_Type::instance().registerProduct( "monolithicGI", &createFM )
                          &&
                          BlockPrecFactory::instance().registerProduct("ComposedDNDGI"  , &createComposedDNDGI )
                          &&
                          BlockPrecFactory::instance().registerProduct("ComposedDND2GI"  , &createComposedDND2GI )
                          &&
                          BlockPrecFactory::instance().registerProduct("ComposedDN2GI"  , &createComposedDN2GI );

}
