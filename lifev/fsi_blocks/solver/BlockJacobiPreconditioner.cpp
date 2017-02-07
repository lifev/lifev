
#include <lifev/fsi_blocks/solver/BlockJacobiPreconditioner.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>
#include <lifev/core/linear_algebra/AztecooOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BelosOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>

#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/LifeChrono.hpp>


namespace LifeV
{
namespace Operators
{
//===========================================================================//
// Constructors
//===========================================================================//

BlockJacobiPreconditioner::BlockJacobiPreconditioner():
    M_label("BlockJacobiPreconditioner"),
    M_useTranspose(false),
    M_approximatedStructureMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedGeometryOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedFluidMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedSchurComplementOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedSchurComplementCouplingOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_shapeDerivatives ( false ),
    M_subiterateFluidDirichlet ( false ),
    M_useStabilization ( false ),
    M_nonconforming ( false ),
    M_useBDFStructure ( false )
{

}

BlockJacobiPreconditioner::~BlockJacobiPreconditioner()
{

}

// show information about the class
void BlockJacobiPreconditioner::showMe(){
    //std::cout<<"Dimension u: "<< M_Bt->NumGlobalRows()<<
    //", Dimension p: "<<M_Bt->NumGlobalCols()<<std::endl;
    //std::cout<<"Pressure correction order "<<M_p<<std::endl;
}

void
BlockJacobiPreconditioner::setFluidPreconditioner(const std::string& type)
{
	M_FluidPrec.reset ( Operators::NSPreconditionerFactory::instance().createObject (type));
}

void
BlockJacobiPreconditioner::setMonolithicMap(const mapEpetraPtr_Type& monolithicMap)
{
	M_monolithicMap = monolithicMap;
}

void
BlockJacobiPreconditioner::setStructureBlock(const matrixEpetraPtr_Type & S)
{
	M_S = S;
	M_X_displacement.reset( new VectorEpetra_Type( M_S->map(), Unique) );
	M_Y_displacement.reset( new VectorEpetra_Type( M_S->map(), Unique) );
	M_structure = M_S->map().mapSize();
}

void
BlockJacobiPreconditioner::setGeometryBlock(const matrixEpetraPtr_Type & G)
{
	M_G = G;
	M_X_geometry.reset( new VectorEpetra_Type (M_G->map(), Unique) );
	M_Y_geometry.reset( new VectorEpetra_Type (M_G->map(), Unique) );
}

void
BlockJacobiPreconditioner::setFluidBlocks(const matrixEpetraPtr_Type & F,
		  	  	  	  	  	  	   const matrixEpetraPtr_Type & Btranspose,
		  	  	  	  	  	  	   const matrixEpetraPtr_Type & B)
{
	M_F = F;
	M_Btranspose = Btranspose;
	M_B = B;
	M_X_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique) );
	M_X_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
	M_Y_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique) );
	M_Y_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
	M_fluidVelocity = M_F->map().mapSize();
	M_fluid         = M_fluidVelocity + M_B->map().mapSize();
	M_useStabilization = false;
}

void
BlockJacobiPreconditioner::setFluidBlocks(const matrixEpetraPtr_Type & F,
		  	  	  	  	  	  	   	   	   	   const matrixEpetraPtr_Type & Btranspose,
		  	  	  	  	  	  	   	   	   	   const matrixEpetraPtr_Type & B,
		  	  	  	  	  	  	   	   	   	   const matrixEpetraPtr_Type & D)
{
	M_F = F;
	M_Btranspose = Btranspose;
	M_B = B;
	M_D.reset( new matrixEpetra_Type ( *D ) );
	M_X_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique) );
	M_X_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
	M_Y_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique) );
	M_Y_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
	M_fluidVelocity = M_F->map().mapSize();
	M_fluid         = M_fluidVelocity + M_B->map().mapSize();
	M_useStabilization = true;
}

void
BlockJacobiPreconditioner::setCouplingBlocks ( const matrixEpetraPtr_Type & C1transpose,
				   	   	   	   	  	    const matrixEpetraPtr_Type & C2transpose,
				   	   	   	   	  	    const matrixEpetraPtr_Type & C2,
				   	   	   	   	  	    const matrixEpetraPtr_Type & C1,
				   	   	   	   	  	    const matrixEpetraPtr_Type & C3)
{
	M_C1transpose = C1transpose;
	M_C2transpose = C2transpose;
	M_C2 = C2;
	M_C1 = C1;
	M_C3 = C3;
	M_X_lambda.reset(  new VectorEpetra_Type ( M_C1->map(), Unique) );
	M_Y_lambda.reset(  new VectorEpetra_Type ( M_C1->map(), Unique) );
	M_lambda = M_C1->map().mapSize();
}

void
BlockJacobiPreconditioner::setCouplingOperators_nonconforming( interpolationPtr_Type fluidToStructure,
																    interpolationPtr_Type structureToFluid,
																    mapEpetraPtr_Type lagrangeMap)
{
	M_FluidToStructureInterpolant = fluidToStructure;
	M_StructureToFluidInterpolant = structureToFluid;

	M_lagrangeMap = lagrangeMap;
	M_X_lambda.reset(  new VectorEpetra_Type ( *M_lagrangeMap, Unique) );
	M_Y_lambda.reset(  new VectorEpetra_Type ( *M_lagrangeMap, Unique) );
	M_lambda = M_lagrangeMap->map(Unique)->NumGlobalElements();
    M_nonconforming = true;
}

void
BlockJacobiPreconditioner::setUseShapeDerivatives(const bool & useShapeDerivatives)
{
	M_shapeDerivatives = useShapeDerivatives;
}

void
BlockJacobiPreconditioner::setSubiterateFluidDirichlet(const bool & subiterateFluidDirichlet)
{
	M_subiterateFluidDirichlet = subiterateFluidDirichlet;
}

void
BlockJacobiPreconditioner::setShapeDerivativesBlocks( const matrixEpetraPtr_Type & shapeVelocity,
   											   const matrixEpetraPtr_Type & shapePressure)
{
	M_shapeVelocity = shapeVelocity;
	M_shapePressure = shapePressure;
}

void
BlockJacobiPreconditioner::setPCDBlocks(const matrixEpetraPtr_Type & Fp,
		  	  	  	  	  	  	   	   	   	 const matrixEpetraPtr_Type & Mp,
		  	  	  	  	  	  	   	   	   	 const matrixEpetraPtr_Type & Mu)
{
	M_Fp = Fp;
	M_Mp = Mp;
	M_Mu = Mu;
}

void
BlockJacobiPreconditioner::setOptions(const Teuchos::ParameterList& solversOptions)
{

	std::shared_ptr<Teuchos::ParameterList> structureMomentumOptions;
	structureMomentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("StructureMomentumOperator")) );
	setStructureMomentumOptions(structureMomentumOptions);

	std::shared_ptr<Teuchos::ParameterList> geometryMomentumOptions;
	geometryMomentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("ALEOperator")) );
	setGeometryOptions(geometryMomentumOptions);

	std::shared_ptr<Teuchos::ParameterList> fluidMomentumOptions;
	fluidMomentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("FluidMomentumOperator")) );
	setFluidMomentumOptions(fluidMomentumOptions);

	std::shared_ptr<Teuchos::ParameterList> schurFluidOptions;
	schurFluidOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("ApproximatedSchurOperatorFluid")) );
	setSchurOptions(schurFluidOptions);

	if ( std::strcmp(M_FluidPrec->Label(),"aPCDOperator") == 0 )
	{
		std::shared_ptr<Teuchos::ParameterList> pressureMassOptions;
		pressureMassOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("PressureMassOperator")) );
		setPressureMassOptions(pressureMassOptions);
	}

}

void
BlockJacobiPreconditioner::setStructureMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_structureMomentumOptions = _oList;
}

void
BlockJacobiPreconditioner::setGeometryOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_geometryOptions = _oList;
}

void
BlockJacobiPreconditioner::setFluidMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_fluidMomentumOptions = _oList;
    M_FluidPrec->setMomentumOptions( M_fluidMomentumOptions );
}

void
BlockJacobiPreconditioner::setSchurOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurOptions = _oList;
    M_FluidPrec->setSchurOptions( M_schurOptions );
}

void
BlockJacobiPreconditioner::setPressureMassOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_pressureMassOptions = _oList;
    M_FluidPrec->setPressureMassOptions( M_pressureMassOptions );
}

void
BlockJacobiPreconditioner::updateApproximatedStructureMomentumOperator( )
{
	M_approximatedStructureMomentumOperator->SetRowMatrix(M_S->matrixPtr());
	M_approximatedStructureMomentumOperator->SetParameterList(*M_structureMomentumOptions);
	M_approximatedStructureMomentumOperator->Compute();

}

void
BlockJacobiPreconditioner::updateApproximatedGeometryOperator( )
{
    M_approximatedGeometryOperator->SetRowMatrix(M_G->matrixPtr());
    M_approximatedGeometryOperator->SetParameterList(*M_geometryOptions);
    M_approximatedGeometryOperator->Compute();
}

void
BlockJacobiPreconditioner::updateApproximatedFluidOperator( )
{
    matrixEpetraPtr_Type Fprecc;
    matrixEpetraPtr_Type Btprecc;
    Fprecc.reset ( new matrixEpetra_Type(*M_F) );
    Btprecc.reset ( new matrixEpetra_Type(*M_Btranspose) );

    if ( !M_myBC->bcUpdateDone() )
    	M_myBC->bcUpdate ( *M_velocityFESpace->mesh(), M_velocityFESpace->feBd(), M_velocityFESpace->dof() );

    bcManageMatrix ( *Fprecc, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *M_myBC, M_velocityFESpace->feBd(), 1.0, 0.0 );

    bcManageMatrix ( *Btprecc, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *M_myBC, M_velocityFESpace->feBd(), 0.0, 0.0 );

    if ( std::strcmp(M_FluidPrec->Label(),"aSIMPLEOperator")==0 )
    {
    	if ( !M_useStabilization )
    		M_FluidPrec->setUp( Fprecc, M_B, Btprecc);
    	else
    		M_FluidPrec->setUp( Fprecc, M_B, Btprecc, M_D);
    }
    else if ( std::strcmp(M_FluidPrec->Label(),"aPCDOperator")==0 )
    {
    	M_FluidPrec->setUp( Fprecc, M_B, Btprecc, M_Fp, M_Mp, M_Mu);
    }

    Displayer M_displayer( M_F->map().commPtr() );

    M_displayer.leaderPrint( "\tNS operator - set up the block operator...");
    LifeChrono chrono;
    chrono.start();

    Operators::NavierStokesOperator::operatorPtrContainer_Type operData(2,2);
    operData(0,0) = Fprecc->matrixPtr();
    operData(0,1) = Btprecc->matrixPtr();
    operData(1,0) = M_B->matrixPtr();
    if ( M_useStabilization )
    	operData(1,1) = M_D->matrixPtr();

    M_oper.reset( new Operators::NavierStokesOperator );
    M_oper->setUp(operData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    M_displayer.leaderPrint( "\tPreconditioner operator - set up the block operator...");
    chrono.reset();
    chrono.start();

    M_FluidPrec->setDomainMap(M_oper->OperatorDomainBlockMapPtr());
    M_FluidPrec->setRangeMap(M_oper->OperatorRangeBlockMapPtr());
    M_FluidPrec->updateApproximatedMomentumOperator();
    M_FluidPrec->updateApproximatedSchurComplementOperator();
    if ( std::strcmp(M_FluidPrec->Label(),"aPCDOperator")==0 )
    {
    	M_displayer.leaderPrint( "\tNS operator - UPDATING PRESSURE MASS...");
    	M_FluidPrec->updateApproximatedPressureMassOperator();
    }

    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    M_displayer.leaderPrint( "\tset up the Trilinos solver...");
    chrono.start();

    std::string solverType(M_fluidMomentumOptions->sublist("FluidSolver").get<std::string>("Linear Solver Type"));
    M_invOper.reset(Operators::InvertibleOperatorFactory::instance().createObject(solverType));
    M_invOper->setParameterList(M_fluidMomentumOptions->sublist("FluidSolver").sublist(solverType));

    M_invOper->setOperator(M_oper);
    M_invOper->setPreconditioner(M_FluidPrec);

    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

}

inline int
BlockJacobiPreconditioner::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    //! Input vector
    const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);

    //! Extract each component of the input vector
    M_X_velocity->zero();
    M_X_velocity->subset(X_vectorEpetra, M_F->map(), 0, 0);

    M_X_pressure->zero();
    M_X_pressure->subset(X_vectorEpetra, M_B->map(), M_fluidVelocity, 0 );

    M_X_displacement->zero();
    M_X_displacement->subset(X_vectorEpetra, M_S->map(), M_fluid, 0 );

    if ( M_nonconforming )
    {
    	M_X_lambda->zero();
    	M_X_lambda->subset(X_vectorEpetra, *M_lagrangeMap, M_fluid + M_structure, 0 );
    }
    else
    {
    	M_X_lambda->zero();
    	M_X_lambda->subset(X_vectorEpetra, M_C1->map(), M_fluid + M_structure, 0 );
    }

    M_X_geometry->zero();
    M_X_geometry->subset(X_vectorEpetra, M_G->map(), M_fluid + M_structure + M_lambda , 0 );

    M_Y_displacement->zero();
    M_approximatedStructureMomentumOperator->ApplyInverse(M_X_displacement->epetraVector(), M_Y_displacement->epetraVector() );

    if ( M_nonconforming )
    {
    	M_Y_geometry->zero();
    	VectorEpetraPtr_Type tmp_geo ( new VectorEpetra_Type ( M_X_geometry->map() ) );
    	tmp_geo->zero();
    	M_StructureToFluidInterpolant->updateRhs(M_Y_displacement);
    	M_StructureToFluidInterpolant->interpolate();
    	M_StructureToFluidInterpolant->solution(tmp_geo);

    	M_approximatedGeometryOperator->ApplyInverse(( *M_X_geometry ).epetraVector(), M_Y_geometry->epetraVector() );
    }
    else
    {
        M_approximatedGeometryOperator->ApplyInverse(( *M_X_geometry ).epetraVector(), M_Y_geometry->epetraVector() );
    }

    VectorEpetraPtr_Type Zlambda ( new VectorEpetra_Type ( M_X_lambda->map() ) );
    Zlambda->zero();

    if ( M_nonconforming )
    {
    	VectorEpetraPtr_Type tmp_z_lambda ( new VectorEpetra_Type ( M_X_lambda->map() ) );
    	VectorEpetraPtr_Type tmp_z_lambda_omega ( new VectorEpetra_Type ( M_X_velocity->map() ) );
    	tmp_z_lambda->zero();
    	tmp_z_lambda_omega->zero();

    	M_StructureToFluidInterpolant->updateRhs(M_Y_displacement);
    	M_StructureToFluidInterpolant->interpolate();
    	M_StructureToFluidInterpolant->solution(tmp_z_lambda_omega);

    	*tmp_z_lambda_omega /= M_timeStep;

        if ( M_useBDFStructure )
        {
            *tmp_z_lambda_omega *= M_bdfCoef;
        }
        else
        {
            *tmp_z_lambda_omega *= M_gamma;
            *tmp_z_lambda_omega /= M_beta;
        }

    	M_FluidToStructureInterpolant->restrictOmegaToGamma_Known(tmp_z_lambda_omega, tmp_z_lambda);

    	*Zlambda = ( *M_X_lambda - *tmp_z_lambda ) ;
    }
    else
    {
        *Zlambda = ( *M_X_lambda + *M_C2* ( *M_Y_displacement ) ) ;
    }


    VectorEpetra_Type Zf_velocity ( *M_X_velocity );
    VectorEpetra_Type Zf_pressure ( *M_X_pressure );
    VectorEpetra_Type Wf_velocity ( Zf_velocity.map() );
    Wf_velocity.zero();

    if ( M_nonconforming )
    {
        Wf_velocity += Zf_velocity;

    	if ( !M_myBC->bcUpdateDone() )
        {
            M_myBC->bcUpdate ( *M_velocityFESpace->mesh(), M_velocityFESpace->feBd(), M_velocityFESpace->dof() );
        }

    	bcManageRhs ( Wf_velocity, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *M_myBC, M_velocityFESpace->feBd(), 0.0, 0.0 );

    	VectorEpetraPtr_Type Zlambda_omega ( new VectorEpetra_Type ( Wf_velocity.map() ) );

    	M_FluidToStructureInterpolant->expandGammaToOmega_Known( Zlambda, Zlambda_omega );

    	Wf_velocity += *Zlambda_omega;
    }
    else
    {
    	Wf_velocity = Zf_velocity - ( *M_C1transpose * (*M_C1 *  Zf_velocity) - ( *M_C1transpose * ( *Zlambda ) ) );
    }

    M_Y_velocity->zero();
    M_Y_pressure->zero();

    M_FluidPrec->ApplyInverse( Wf_velocity, Zf_pressure, *M_Y_velocity, *M_Y_pressure);

    M_Y_lambda->zero();

    if ( M_nonconforming )
    {
    	VectorEpetraPtr_Type tmp_omega ( new VectorEpetra_Type ( M_Y_velocity->map() ) );
    	tmp_omega->zero();

    	VectorEpetraPtr_Type tmp_gamma ( new VectorEpetra_Type ( M_Y_lambda->map() ) );
    	tmp_gamma->zero();

    	*tmp_omega = Zf_velocity + ( *M_F * ( *M_Y_velocity ) + *M_Btranspose * ( *M_Y_pressure ) );

    	M_FluidToStructureInterpolant->restrictOmegaToGamma_Known( tmp_omega, tmp_gamma);

    	*M_Y_lambda += *tmp_gamma;
    }
    else
    {
        *M_Y_lambda = *M_C1 * ( Zf_velocity + ( *M_F * ( *M_Y_velocity )  ) );
    }

    VectorEpetra_Type Y_vectorEpetra(Y, M_monolithicMap, Unique);
    Y_vectorEpetra.zero();

    Y_vectorEpetra.subset(*M_Y_velocity, M_Y_velocity->map(), 0, 0 );
    Y_vectorEpetra.subset(*M_Y_pressure, M_Y_pressure->map(), 0, M_fluidVelocity );
    Y_vectorEpetra.subset(*M_Y_displacement, M_Y_displacement->map(), 0, M_fluid );
    Y_vectorEpetra.subset(*M_Y_lambda, M_Y_lambda->map(), 0, M_fluid + M_structure );
    Y_vectorEpetra.subset(*M_Y_geometry, M_Y_geometry->map(), 0, M_fluid + M_structure + M_lambda );
    Y = dynamic_cast<Epetra_MultiVector &>( Y_vectorEpetra.epetraVector() );

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
