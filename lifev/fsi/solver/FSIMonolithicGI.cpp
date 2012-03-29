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
#include <lifev/structure/solver/VenantKirchhoffSolver.hpp>

#include <lifev/fsi/solver/FSIMonolithicGI.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDND.hpp>
#include <lifev/fsi/solver/MonolithicBlockMatrixRN.hpp>

namespace LifeV
{

  // ===================================================
  //  Constructors and Descructor
  // ===================================================
  FSIMonolithicGI::FSIMonolithicGI():
    super_Type              (),
    M_mapWithoutMesh        (),
    M_uk                    (),
    M_domainVelImplicit     (true),
    M_convectiveTermDer     (true),
    M_interface             (0),
    M_meshBlock             (),
    M_shapeDerivativesBlock (),
    M_solidDerBlock         ()
  {}

  // ===================================================
  //  Public Methods
  // ===================================================
  void
  FSIMonolithicGI::setUp( const GetPot& dataFile )
  {
    super_Type::setUp( dataFile );

    M_domainVelImplicit = M_data->dataFluid()->domainVelImplicit(); //dataFile( "fluid/domainVelImplicit", true );
    M_convectiveTermDer = M_data->dataFluid()->convectiveImplicit(); //dataFile( "fluid/convectiveTermDer", false );
  }

  void
  FSIMonolithicGI::setupFluidSolid( UInt const fluxes )
  {
    super_Type::setupFluidSolid( fluxes );
    UInt offset = M_monolithicMap->map(Unique)->NumGlobalElements();
    M_mapWithoutMesh.reset(new MapEpetra(*M_monolithicMap));

    *M_monolithicMap += M_mmFESpace->map();
    /* OBSOLETE
       if(M_data->dataFluid()->useShapeDerivatives())
       {
       M_epetraOper.reset( new Epetra_FullMonolithic(this));
       M_solid->setOperator(*M_epetraOper);
       }*/

    M_interface=M_monolithicMatrix->interface();

    M_beta.reset( new vector_Type(M_uFESpace->map()) );
    M_rhs.reset(new vector_Type(*M_monolithicMap));
    M_rhsFull.reset(new vector_Type(*M_monolithicMap));
    if(M_data->dataFluid()->useShapeDerivatives())
      M_shapeDerivativesBlock.reset(new FSIOperator::fluidPtr_Type::value_type::matrix_Type(*M_monolithicMap));
    M_uk.reset (new vector_Type(*M_monolithicMap));

    M_meshMotion.reset(new meshMotion_Type(*M_mmFESpace,
                                           M_epetraComm,
                                           *M_monolithicMap,
                                           offset));

    M_fluid.reset     (new fluid_Type(M_data->dataFluid(),
                                      *M_uFESpace,
                                      *M_pFESpace,
                                      *M_mmFESpace,
                                      M_epetraComm,
                                      *M_monolithicMap,
                                      fluxes));
    M_solid.reset(new solid_Type());

    M_solid->setup(M_data->dataSolid(),
                   M_dFESpace,
                   M_epetraComm,
                   M_monolithicMap,
                   M_offset
		   );
  }

  void
  FSIMonolithicGI::buildSystem()
  {
    super_Type::buildSystem();
    M_meshMotion->computeMatrix();
  }

  void
  FSIMonolithicGI::updateSystem()
  {
    super_Type::updateSystem();
  }

  void
  FSIMonolithicGI::evalResidual( vector_Type&       res,
				 const vector_Type& disp,
				 const UInt          iter )
  {
    res *= 0.;//this is important. Don't remove it!
    if ((iter==0)|| !M_data->dataFluid()->isSemiImplicit())
      {

	Real alpha( 1./M_data->dataFluid()->dataTime()->timeStep() );
	if(M_restarts)
	  {
            alpha = 1/M_data->restartTimeStep();
            M_restarts = false;
	  }

	M_uk.reset(new vector_Type( disp ));//M_uk should point to the same vector as disp, indeed, not be copied.
	//in that case this line is useless

	UInt offset( M_solidAndFluidDim + nDimensions*M_interface );

	vectorPtr_Type meshDisp( new vector_Type(M_mmFESpace->map()) );
	vectorPtr_Type meshVel( new vector_Type(M_mmFESpace->map()) );
	vectorPtr_Type mmRep( new vector_Type(M_mmFESpace->map(), Repeated ));
	meshDisp->subset(disp, offset); //if the conv. term is to be condidered implicitly

	if (!M_domainVelImplicit)//if the mesh motion is at the previous time step in the convective term
	  {
	    *meshVel = M_ALETimeAdvance->velocity( );
	    M_ALETimeAdvance->extrapolation(*mmRep);
	    moveMesh(*mmRep);// re-initialize the mesh points
	    if( iter==0 )
	      {
		M_ALETimeAdvance->updateRHSFirstDerivative(M_data->dataFluid()->dataTime()->timeStep());
		M_ALETimeAdvance->shiftRight(*meshDisp);
	      }
	  }
        else
	  {
            if( iter==0 )
	      {
                M_ALETimeAdvance->updateRHSFirstDerivative(M_data->dataFluid()->dataTime()->timeStep());
                M_ALETimeAdvance->shiftRight(*meshDisp);
                M_ALETimeAdvance->extrapolation(*mmRep);
	      }
            else
	      {
                M_ALETimeAdvance->setSolution(*meshDisp);
                *mmRep = *meshDisp;
	      }
            *meshVel = M_ALETimeAdvance->velocity();
            moveMesh(*mmRep);// re-initialize the mesh points
	  }

        *mmRep = *meshVel*(-1.);
        interpolateVelocity(*mmRep, *M_beta);

        vectorPtr_Type fluid(new vector_Type(M_uFESpace->map()));
        if (!M_convectiveTermDer)
	  {
            M_fluidTimeAdvance->extrapolation(*fluid );
	  }
        else
	  {
            fluid->subset(disp, 0);
	  }
        *M_beta += *fluid; /*M_un or disp, it could be also M_uk in a FP strategy*/

        assembleSolidBlock( iter, *M_uk );
        assembleFluidBlock( iter, *M_uk );
        assembleMeshBlock ( iter );
        *M_rhsFull = *M_rhs;

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
	    M_monolithicMatrix->coupler(M_monolithicMap, M_dofStructureToFluid->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), M_solidTimeAdvance->coefficientFirstDerivative( 0 ), M_solid->rescaleFactor());
	    M_monolithicMatrix->coupler( M_monolithicMap, M_dofStructureToHarmonicExtension->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep(), M_solidTimeAdvance->coefficientFirstDerivative( 0 ), M_solid->rescaleFactor(), 2);
	  }
	else
	  {
	    M_monolithicMatrix->replace_matrix(M_solidBlockPrec, 0);
	    M_monolithicMatrix->replace_matrix(M_fluidBlock, 1);
	    M_monolithicMatrix->replace_matrix(M_meshBlock, 2);
	  }

	M_monolithicMatrix->blockAssembling();
	super_Type::checkIfChangedFluxBC( M_monolithicMatrix );



	if( (M_data->dataSolid()->solidType().compare("exponential") && M_data->dataSolid()->solidType().compare("neoHookean")) )
	  applyBoundaryConditions();
      }
    M_monolithicMatrix->GlobalAssemble();
    //M_monolithicMatrix->matrix()->spy("FMFI");
    super_Type::evalResidual( disp, M_rhsFull, res, false );

    if(!(M_data->dataSolid()->solidType().compare("exponential") && M_data->dataSolid()->solidType().compare("neoHookean")) )
      {
	res += *M_meshBlock*disp;

	if ( !M_BCh_u->bcUpdateDone() )
	  M_BCh_u->bcUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
	M_BCh_d->setOffset(M_offset);
	if ( !M_BCh_d->bcUpdateDone() )
	  M_BCh_d->bcUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
	M_BCh_mesh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);
	if ( !M_BCh_mesh->bcUpdateDone() )
	  M_BCh_mesh->bcUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

	M_monolithicMatrix->applyBoundaryConditions(dataFluid()->dataTime()->time()/*, M_rhsFull*/);
	M_monolithicMatrix->GlobalAssemble();


	bcManageResidual(res, *M_rhsFull, disp, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_u, M_uFESpace->feBd(),  M_data->dataFluid()->dataTime()->time(), 1.);

	// below sol is repeated by BCManageResidual
	bcManageResidual(res, *M_rhsFull, disp, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), M_data->dataSolid()->dataTime()->time(), 1.);
	bcManageResidual(res, *M_rhsFull, disp, *M_mmFESpace->mesh(), M_mmFESpace->dof(), *M_BCh_mesh, M_mmFESpace->feBd(),  M_data->dataFluid()->dataTime()->time(), 1.);
	res -= *M_rhsFull;
      }
  }

  void
  FSIMonolithicGI::applyBoundaryConditions()
  {
    if ( !M_BCh_u->bcUpdateDone() )
      M_BCh_u->bcUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    M_BCh_d->setOffset(M_offset);
    if ( !M_BCh_d->bcUpdateDone() )
      M_BCh_d->bcUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    M_BCh_mesh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);
    if ( !M_BCh_mesh->bcUpdateDone() )
      M_BCh_mesh->bcUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

    M_monolithicMatrix->applyBoundaryConditions(dataFluid()->dataTime()->time(), M_rhsFull);
    // M_monolithicMatrix->GlobalAssemble();
    // M_monolithicMatrix->matrix()->spy("FMFI");
  }


  //============ Protected Methods ===================

  void FSIMonolithicGI::setupBlockPrec()
  {

    //M_solidDerBlock = M_solidBlockPrec; // an inexact Newton approximation of the Jacobian

    //The following part accounts for a possibly nonlinear structure model, should not be run when linear
    //elasticity is used

    if ( M_data->dataSolid()->getUseExactJacobian() && (M_data->dataSolid()->solidType().compare("exponential") && M_data->dataSolid()->solidType().compare("neoHookean")))
      {
        M_solid->material()->updateJacobianMatrix( *M_uk*M_solid->rescaleFactor()/**M_data->dataFluid()->dataTime()->timeStep()*/, dataSolid(), M_solid->displayerPtr() ); // computing the derivatives if nonlinear (comment this for inexact Newton);
        M_solidBlockPrec.reset(new matrix_Type(*M_monolithicMap, 1));
        *M_solidBlockPrec += *M_solid->Mass();
        *M_solidBlockPrec += *M_solid->material()->jacobian(); //stiffMatrix();
        M_solidBlockPrec->globalAssemble();
        //        *M_solidBlockPrec *= M_data->dataSolid()->dataTime()->timeStep();
        *M_solidBlockPrec *= M_solid->rescaleFactor();

	//         *M_solidBlockPrec += *M_solidDerBlock;
        //M_precPtr->replace_matrix( M_solidBlockPrec, 0 );
        M_monolithicMatrix->replace_matrix( M_solidBlockPrec, 0 );
      }

    // *M_monolithicMatrix->matrix() *= 0;
    M_monolithicMatrix->blockAssembling();
    //*M_monolithicMatrix->matrix() *= M_data->dataFluid()->dataTime()->timeStep();

    if (M_data->dataFluid()->useShapeDerivatives())
      {
        *M_shapeDerivativesBlock *= 0.;
        M_shapeDerivativesBlock->openCrsMatrix( );
        shapeDerivatives( M_shapeDerivativesBlock );
        M_shapeDerivativesBlock->globalAssemble( );
        M_monolithicMatrix->addToGlobalMatrix( M_shapeDerivativesBlock );
      }

    if ( M_data->dataFluid()->useShapeDerivatives() || M_data->dataSolid()->getUseExactJacobian() )
      {
        if ( !M_BCh_u->bcUpdateDone() )
	  M_BCh_u->bcUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
        M_BCh_d->setOffset(M_offset);
        if ( !M_BCh_d->bcUpdateDone() )
	  M_BCh_d->bcUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
        M_BCh_mesh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);
        if ( !M_BCh_mesh->bcUpdateDone() )
	  M_BCh_mesh->bcUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );

        M_monolithicMatrix->applyBoundaryConditions(dataFluid()->dataTime()->time());
        M_monolithicMatrix->GlobalAssemble();
	//M_monolithicMatrix->matrix()->spy("jacobian");
      }


    super_Type::setupBlockPrec( );

    if ( M_precPtr->blockVector( ).size( )<3 )
      {
        M_precPtr->push_back_matrix( M_meshBlock, false );
        M_precPtr->setConditions( M_BChs );
        M_precPtr->setSpaces( M_FESpaces );
        M_precPtr->setOffsets( 3, M_offset, 0,  M_solidAndFluidDim + nDimensions*M_interface );
        M_precPtr->coupler( M_monolithicMap, M_dofStructureToFluid/*HarmonicExtension*/->localDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->timeStep() ,M_solidTimeAdvance->coefficientFirstDerivative( 0 ), M_solid->rescaleFactor(), 2);

        if (M_data->dataFluid()->useShapeDerivatives())
	  {
            boost::shared_ptr<MatrixEpetra<Real> > staticCast=boost::static_pointer_cast<MatrixEpetra<Real> >(M_shapeDerivativesBlock);
            M_precPtr->push_back_coupling( staticCast );
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

  void FSIMonolithicGI::shapeDerivatives( FSIOperator::fluidPtr_Type::value_type::matrixPtr_Type sdMatrix )
  {
    Real alpha = 1./M_data->dataFluid()->dataTime()->timeStep();
    vectorPtr_Type rhsNew(new vector_Type(*M_monolithicMap));
    vector_Type un(M_uFESpace->map());
    vector_Type uk(M_uFESpace->map()+M_pFESpace->map());

    vectorPtr_Type meshVelRep(new vector_Type(M_mmFESpace->map(), Repeated));

    *meshVelRep=M_ALETimeAdvance->velocity();

    if ( M_convectiveTermDer )
      {
        un.subset(*M_uk, 0);
      }
    else
      M_fluidTimeAdvance->extrapolation(un);

    uk.subset(*M_uk, 0);
    vector_Type veloFluidMesh(M_uFESpace->map(), Repeated);
    this->transferMeshMotionOnFluid(*meshVelRep, veloFluidMesh);

    M_fluid->updateShapeDerivatives(*sdMatrix,
                                    alpha,
                                    un,//un if !M_convectiveTermDer, otherwise uk
                                    uk,//uk
                                    veloFluidMesh, //(xk-xn)/dt (FI), or (xn-xn-1)/dt (CE)//Repeated
                                    M_solidAndFluidDim+M_interface*nDimensions,
                                    *M_uFESpace,
                                    M_domainVelImplicit,
                                    M_convectiveTermDer
				    );
  }

  void
  FSIMonolithicGI::assembleMeshBlock(UInt /*iter*/)
  {
    M_meshBlock.reset(new matrix_Type(*M_monolithicMap));
    M_meshMotion->addSystemMatrixTo(M_meshBlock);
    M_meshBlock->globalAssemble();
    UInt offset(M_solidAndFluidDim+nDimensions*M_interface);
    std::map<ID, ID>::const_iterator ITrow;
    std::map<ID, ID> locdofmap(M_dofStructureToHarmonicExtension->localDofMap());

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

    for ( ID dim=0; dim < nDimensions; ++dim )
      for ( ITrow = locdofmap.begin(); ITrow != locdofmap.end(); ++ITrow )
        {
	  UInt i = ITrow->first;
	  M_meshBlock->diagonalize(i+offset+dim*M_mmFESpace->dof().numTotalDof() , 1.);
        }
  }

  // ===================================================
  //  Products registration
  // ===================================================
  bool FSIMonolithicGI::S_register =  BlockPrecFactory::instance().registerProduct("AdditiveSchwarzGI"  , &MonolithicBlockMatrix::createAdditiveSchwarz )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDNGI"  , &MonolithicBlockComposedDN::createComposedDNGI )
    &&
    MonolithicBlockMatrix::Factory_Type::instance().registerProduct( "AdditiveSchwarzGI", &MonolithicBlockMatrix::createAdditiveSchwarz )
    &&
    MonolithicBlockMatrix::Factory_Type::instance().registerProduct( "AdditiveSchwarzRNGI", &MonolithicBlockMatrixRN::createAdditiveSchwarzRN )
    &&
    FSIOperator::FSIFactory_Type::instance().registerProduct( "monolithicGI", &FSIMonolithicGI::instantiate )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDNDGI"  , &MonolithicBlockComposedDND::createComposedDNDGI )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDND2GI"  , &MonolithicBlockComposedDND::createComposedDND2GI )
    &&
    BlockPrecFactory::instance().registerProduct("ComposedDN2GI"  , &MonolithicBlockComposedDN::createComposedDN2GI );

}
