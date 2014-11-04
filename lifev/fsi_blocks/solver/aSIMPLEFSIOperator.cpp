
#include <lifev/fsi_blocks/solver/aSIMPLEFSIOperator.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>
#include <lifev/core/linear_algebra/AztecooOperator.hpp>
#include <lifev/core/linear_algebra/BelosOperator.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>

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
    M_approximatedSchurComplementCouplingOperator(new Operators::ApproximatedInvertibleRowMatrix)
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
}

void
aSIMPLEFSIOperator::setGeometryBlock(const matrixEpetraPtr_Type & G)
{
	M_G = G;
}

void
aSIMPLEFSIOperator::setFluidBlocks(const matrixEpetraPtr_Type & F,
		  	  	  	  	  	  	   const matrixEpetraPtr_Type & Btranspose,
		  	  	  	  	  	  	   const matrixEpetraPtr_Type & B)
{
	M_F = F;
	M_Btranspose = Btranspose;
	M_B = B;
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

	boost::shared_ptr<Teuchos::ParameterList> schurCouplingOptions;
	schurCouplingOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("ApproximatedSchurCouplingOperator")) );
	setSchurCouplingOptions(schurCouplingOptions);
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
}

void
aSIMPLEFSIOperator::setSchurOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurOptions = _oList;
}

void
aSIMPLEFSIOperator::setSchurCouplingOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurCouplingOptions = _oList;
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
aSIMPLEFSIOperator::updateApproximatedFluidMomentumOperator( )
{
	M_approximatedFluidMomentumOperator->SetRowMatrix( M_F->matrixPtr() );
	M_approximatedFluidMomentumOperator->SetParameterList(*M_fluidMomentumOptions);
	M_approximatedFluidMomentumOperator->Compute();

}

void
aSIMPLEFSIOperator::updateApproximatedSchurComplementOperator( )
{
	buildShurComplement();
	M_approximatedSchurComplementOperator->SetRowMatrix( M_schurComplement->matrixPtr() );
    M_approximatedSchurComplementOperator->SetParameterList(*M_schurOptions);
    M_approximatedSchurComplementOperator->Compute();
}

void
aSIMPLEFSIOperator::updateApproximatedSchurComplementCouplingOperator( )
{
	buildShurComplementCoupling();
	M_approximatedSchurComplementCouplingOperator->SetRowMatrix( M_schurComplementCoupling->matrixPtr() );
	M_approximatedSchurComplementCouplingOperator->SetParameterList(*M_schurCouplingOptions);
	M_approximatedSchurComplementCouplingOperator->Compute();
}

void
aSIMPLEFSIOperator::buildShurComplement( )
{
    Epetra_Vector diag( M_F->matrixPtr()->OperatorRangeMap() );
    M_invD.reset(new Epetra_Vector( M_F->matrixPtr()->OperatorRangeMap() ) );

    // extracting diag(F)
    M_F->matrixPtr()->ExtractDiagonalCopy(diag);

    // computing diag(F)^{-1}
    M_invD->Reciprocal(diag);

    // computing diag(F)^{-1}*M_Btranspose
    matrixEpetra_Type FBT (*M_Btranspose);
    FBT.matrixPtr()->LeftScale(*M_invD);

    M_schurComplement.reset ( new matrixEpetra_Type ( M_B->map() ) );

    // computing M_B*(diag(F)^{-1}*M_Btranspose)
    M_B->multiply (false, FBT, false, *M_schurComplement, false);
    M_schurComplement->globalAssemble();
}

void
aSIMPLEFSIOperator::buildShurComplementCoupling( )
{
	// computing diag(F)^{-1}*M_C1transpose
    matrixEpetra_Type invDC1transpose (*M_C1transpose);
    invDC1transpose.matrixPtr()->LeftScale(*M_invD);

    M_schurComplementCoupling.reset ( new matrixEpetra_Type ( M_C1->map() ) );

    // computing M_C1*(diag(F)^{-1}*M_C1transpose)
    M_C1->multiply (false, invDC1transpose, false, *M_schurComplementCoupling, false);
    M_schurComplementCoupling->globalAssemble();
}

int
aSIMPLEFSIOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    //! Input vector
    const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);

    //! Extract each component of the input vector
    VectorEpetra_Type X_velocity(M_F->map(), Unique);
    X_velocity.subset(X_vectorEpetra, M_F->map(), 0, 0);

    VectorEpetra_Type X_pressure(M_B->map(), Unique);
    X_pressure.subset(X_vectorEpetra, M_B->map(), M_F->map().mapSize(), 0 );

    VectorEpetra_Type X_displacement(M_S->map(), Unique);
    X_displacement.subset(X_vectorEpetra, M_S->map(), M_F->map().mapSize() + M_B->map().mapSize(), 0 );

    VectorEpetra_Type X_lambda(M_C1->map(), Unique);
    X_lambda.subset(X_vectorEpetra, M_C1->map(), M_F->map().mapSize() + M_B->map().mapSize() + M_S->map().mapSize(), 0 );

    VectorEpetra_Type X_geometry(M_G->map(), Unique);
    X_geometry.subset(X_vectorEpetra, M_G->map(), M_F->map().mapSize() + M_B->map().mapSize() + M_S->map().mapSize() + M_C1->map().mapSize(), 0 );

    //--------------------------------------------------//
    // First: apply the preconditioner of the structure //
    //--------------------------------------------------//

    VectorEpetra_Type Y_displacement(M_S->map(), Unique);
    M_approximatedStructureMomentumOperator->ApplyInverse(X_displacement.epetraVector(), Y_displacement.epetraVector() );


    // output vector
    VectorEpetra_Type Y_vectorEpetra(M_monolithicMap, Unique);


    //! Copy the single contributions into the optput vector
    Y_vectorEpetra.subset(Y_displacement, Y_displacement.map(), 0, M_F->map().mapSize() + M_B->map().mapSize() );

    Y = dynamic_cast<Epetra_MultiVector &>( Y_vectorEpetra.epetraVector() );

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
