
#include <lifev/fsi_blocks/solver/aSIMPLEFSIOperator.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>
#include <lifev/core/linear_algebra/AztecooOperator.hpp>
#include <lifev/core/linear_algebra/BelosOperator.hpp>
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

aSIMPLEFSIOperator::aSIMPLEFSIOperator():
    M_label("aSIMPLEFSIOperator"),
    M_useTranspose(false),
    M_approximatedStructureMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedGeometryOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedFluidMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedSchurComplementOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedSchurComplementCouplingOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_shapeDerivatives ( false )
{

}

aSIMPLEFSIOperator::~aSIMPLEFSIOperator()
{

}

// show information about the class
void aSIMPLEFSIOperator::showMe(){
    //std::cout<<"Dimension u: "<< M_Bt->NumGlobalRows()<<
    //", Dimension p: "<<M_Bt->NumGlobalCols()<<std::endl;
    //std::cout<<"Pressure correction order "<<M_p<<std::endl;
}

void
aSIMPLEFSIOperator::setMonolithicMap(const mapEpetraPtr_Type& monolithicMap)
{
	M_monolithicMap = monolithicMap;
}

void
aSIMPLEFSIOperator::setStructureBlock(const matrixEpetraPtr_Type & S)
{
	M_S = S;
	M_X_displacement.reset( new VectorEpetra_Type( M_S->map(), Unique) );
	M_structure = M_S->map().mapSize();
}

void
aSIMPLEFSIOperator::setGeometryBlock(const matrixEpetraPtr_Type & G)
{
	M_G = G;
	M_X_geometry.reset( new VectorEpetra_Type (M_G->map(), Unique) );
}

void
aSIMPLEFSIOperator::setFluidBlocks(const matrixEpetraPtr_Type & F,
		  	  	  	  	  	  	   const matrixEpetraPtr_Type & Btranspose,
		  	  	  	  	  	  	   const matrixEpetraPtr_Type & B)
{
	M_F = F;
	M_Btranspose = Btranspose;
	M_B = B;
	M_X_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique) );
	M_X_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
	M_fluidVelocity = M_F->map().mapSize();
	M_fluid         = M_fluidVelocity + M_B->map().mapSize();
}

void
aSIMPLEFSIOperator::setCouplingBlocks ( const matrixEpetraPtr_Type & C1transpose,
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
	M_X_lambda.reset(  new VectorEpetra_Type ( M_C1->map(), Unique) );\
	M_lambda = M_C1->map().mapSize();
}

void
aSIMPLEFSIOperator::setUseShapeDerivatives(const bool & useShapeDerivatives)
{
	M_shapeDerivatives = useShapeDerivatives;
}

void
aSIMPLEFSIOperator::setShapeDerivativesBlocks( const matrixEpetraPtr_Type & shapeVelocity,
   											   const matrixEpetraPtr_Type & shapePressure)
{
	M_shapeVelocity = shapeVelocity;
	M_shapePressure = shapePressure;
}

void
aSIMPLEFSIOperator::setOptions(const Teuchos::ParameterList& solversOptions)
{

	boost::shared_ptr<Teuchos::ParameterList> structureMomentumOptions;
	structureMomentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("StructureMomentumOperator")) );
	setStructureMomentumOptions(structureMomentumOptions);

	boost::shared_ptr<Teuchos::ParameterList> geometryMomentumOptions;
	geometryMomentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("ALEOperator")) );
	setGeometryOptions(geometryMomentumOptions);

	boost::shared_ptr<Teuchos::ParameterList> fluidMomentumOptions;
	fluidMomentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("FluidMomentumOperator")) );
	setFluidMomentumOptions(fluidMomentumOptions);

	boost::shared_ptr<Teuchos::ParameterList> schurFluidOptions;
	schurFluidOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("ApproximatedSchurOperatorFluid")) );
	setSchurOptions(schurFluidOptions);

}

void
aSIMPLEFSIOperator::setStructureMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_structureMomentumOptions = _oList;
}

void
aSIMPLEFSIOperator::setGeometryOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_geometryOptions = _oList;
}

void
aSIMPLEFSIOperator::setFluidMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_fluidMomentumOptions = _oList;
    M_FluidDirichlet.reset( new Operators::aSIMPLEOperator);
    M_FluidDirichlet->setMomentumOptions( M_fluidMomentumOptions );
}

void
aSIMPLEFSIOperator::setSchurOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurOptions = _oList;
    M_FluidDirichlet->setSchurOptions( M_schurOptions );
}

void
aSIMPLEFSIOperator::updateApproximatedStructureMomentumOperator( )
{
	M_approximatedStructureMomentumOperator->SetRowMatrix(M_S->matrixPtr());
	M_approximatedStructureMomentumOperator->SetParameterList(*M_structureMomentumOptions);
	M_approximatedStructureMomentumOperator->Compute();

}

void
aSIMPLEFSIOperator::updateApproximatedGeometryOperator( )
{
    M_approximatedGeometryOperator->SetRowMatrix(M_G->matrixPtr());
    M_approximatedGeometryOperator->SetParameterList(*M_geometryOptions);
    M_approximatedGeometryOperator->Compute();
}

void
aSIMPLEFSIOperator::updateApproximatedFluidOperator( )
{
    // Creating copies for simple on Navier-Stokes with Dirichlet BC on the interface
    matrixEpetraPtr_Type FDirichlet;
    matrixEpetraPtr_Type BtDirichlet;
    FDirichlet.reset ( new matrixEpetra_Type(*M_F) );
    BtDirichlet.reset ( new matrixEpetra_Type(*M_Btranspose) );

    // Apply Dirichlet BC on the interface
    matrixEpetra_Type C1tC1 (M_F->map(), 5);
    M_C1transpose->multiply( false, *M_C1, false, C1tC1, true);

    matrixEpetra_Type C1tC1F (M_F->map());
    C1tC1.multiply( false, *M_F, false, C1tC1F, true);

    *FDirichlet -= C1tC1F;
    *FDirichlet += C1tC1;

    matrixEpetra_Type C1tC1Bt (M_Btranspose->map());
    C1tC1.multiply( false, *M_Btranspose, false, C1tC1Bt, false);
    C1tC1Bt.globalAssemble( M_Btranspose->domainMapPtr(), M_C1transpose->rangeMapPtr());

    *BtDirichlet -= C1tC1Bt;

    M_FluidDirichlet->setUp( FDirichlet, M_B, BtDirichlet);

    // Setting M_oper
    Displayer M_displayer( M_F->map().commPtr() );

    //(1) Set up the OseenOperator
    M_displayer.leaderPrint( "\tNS operator - set up the block operator...");
    LifeChrono chrono;
    chrono.start();

    Operators::NavierStokesOperator::operatorPtrContainer_Type operData(2,2);
    operData(0,0) = M_FluidDirichlet->F()->matrixPtr();
    operData(0,1) = M_FluidDirichlet->Btranspose()->matrixPtr();
    operData(1,0) = M_FluidDirichlet->B()->matrixPtr();

    M_oper.reset( new Operators::NavierStokesOperator );
    M_oper->setUp(operData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(2) Set the data for the preconditioner
    M_displayer.leaderPrint( "\tPreconditioner operator - set up the block operator...");
    chrono.reset();
    chrono.start();

    //M_FluidDirichlet->setUp(M_F, M_B, M_Btranspose);
    M_FluidDirichlet->setDomainMap(M_oper->OperatorDomainBlockMapPtr());
    M_FluidDirichlet->setRangeMap(M_oper->OperatorRangeBlockMapPtr());
    M_FluidDirichlet->updateApproximatedMomentumOperator();
    M_FluidDirichlet->updateApproximatedSchurComplementOperator();
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(3) Set the solver for the linear system
    M_displayer.leaderPrint( "\tset up the Trilinos solver...");
    chrono.start();

    std::string solverType(M_fluidMomentumOptions->sublist("FluidSolver").get<std::string>("Linear Solver Type"));
    M_invOper.reset(Operators::InvertibleOperatorFactory::instance().createObject(solverType));
    M_invOper->setParameterList(M_fluidMomentumOptions->sublist("FluidSolver").sublist(solverType));

    M_invOper->setOperator(M_oper);
    M_invOper->setPreconditioner(M_FluidDirichlet);

    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

}

int
aSIMPLEFSIOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    //! Input vector
    const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);

    //! Extract each component of the input vector
    M_X_velocity->subset(X_vectorEpetra, M_F->map(), 0, 0);

    M_X_pressure->subset(X_vectorEpetra, M_B->map(), M_fluidVelocity, 0 );

    M_X_displacement->subset(X_vectorEpetra, M_S->map(), M_fluid, 0 );

    M_X_lambda->subset(X_vectorEpetra, M_C1->map(), M_fluid + M_structure, 0 );

    M_X_geometry->subset(X_vectorEpetra, M_G->map(), M_fluid + M_structure + M_lambda , 0 );

    //--------------------------------------------------//
    // First: apply the preconditioner of the structure //
    //--------------------------------------------------//

    VectorEpetra_Type Y_displacement(M_S->map(), Unique);
    M_approximatedStructureMomentumOperator->ApplyInverse(M_X_displacement->epetraVector(), Y_displacement.epetraVector() );

    //--------------------------------------------------//
    // Second: apply the preconditioner of the geometry //
    //--------------------------------------------------//

    VectorEpetra_Type Zg ( *M_X_geometry );
    Zg -= *M_C3*Y_displacement;

    VectorEpetra_Type Y_geometry ( M_X_geometry->map(), Unique);
    M_approximatedGeometryOperator->ApplyInverse(Zg.epetraVector(), Y_geometry.epetraVector() );

    //----------------------------------------------//
    // Third: Update shape derivatives and coupling with the solid (velocity)
    //----------------------------------------------//

    VectorEpetra_Type Zlambda ( *M_X_lambda );
    Zlambda -= *M_C2*Y_displacement;

    // Shape derivatives
    VectorEpetra_Type Zf_velocity ( *M_X_velocity );
    VectorEpetra_Type Zf_pressure ( *M_X_pressure );

    if ( M_shapeDerivatives )
    {
    	Zf_velocity -= *M_shapeVelocity * Y_geometry;
    	Zf_pressure	-= *M_shapePressure * Y_geometry;
    }

    //----------------------------------------------//
    // Forth: Fluid by condensation
    //----------------------------------------------//
    VectorEpetra_Type Y_velocity ( M_X_velocity->map() );
    VectorEpetra_Type Y_pressure ( M_X_pressure->map() );

    VectorEpetra_Type Wf_velocity ( Zf_velocity );

    Wf_velocity -= *M_C1transpose * (*M_C1 *  Zf_velocity);
    Wf_velocity += *M_C1transpose * ( Zlambda );

    // Solver of fluid with dirichlet BC on interface
    //----------------------------------------------//

    if (true)
    {
        M_FluidDirichlet->ApplyInverse( Wf_velocity, Zf_pressure, Y_velocity, Y_pressure);
    }
    else
    {

        // Solving the system
        BlockEpetra_Map upMap;
        upMap.setUp ( M_X_velocity->map().map(Unique), M_X_pressure->map().map(Unique));

        // Todo: optimize with pushing shared pointer into blockepetra.
        BlockEpetra_MultiVector up(upMap, 1), rhs(upMap, 1);
        rhs.block(0).Update(1.0, Wf_velocity.epetraVector(), 0.);
        rhs.block(1).Update(1.0, Zf_pressure.epetraVector(), 0.);
        // Solving the linear system
        M_invOper->ApplyInverse(rhs,up);

        Y_velocity.epetraVector().Update(1.0,up.block(0),0.0);
        Y_pressure.epetraVector().Update(1.0,up.block(1),0.0);


    }

    // Update coupling
    VectorEpetra_Type Y_lambda ( M_X_lambda->map(), Unique );

    VectorEpetra_Type Kf_velocity ( Zf_velocity );
    Kf_velocity -= (*M_F * Y_velocity + *M_Btranspose * Y_pressure);
    Y_lambda = *M_C1 * Kf_velocity;

    // output vector
    VectorEpetra_Type Y_vectorEpetra(Y, M_monolithicMap, Unique);

    //! Copy the single contributions into the optput vector
    Y_vectorEpetra.subset(Y_velocity, Y_velocity.map(), 0, 0 );
    Y_vectorEpetra.subset(Y_pressure, Y_pressure.map(), 0, M_fluidVelocity );
    Y_vectorEpetra.subset(Y_displacement, Y_displacement.map(), 0, M_fluid );
    Y_vectorEpetra.subset(Y_lambda, Y_lambda.map(), 0, M_fluid + M_structure );
    Y_vectorEpetra.subset(Y_geometry, Y_geometry.map(), 0, M_fluid + M_structure + M_lambda );

    Y = dynamic_cast<Epetra_MultiVector &>( Y_vectorEpetra.epetraVector() );


    return 0;
}

} /* end namespace Operators */

} /*end namespace */
