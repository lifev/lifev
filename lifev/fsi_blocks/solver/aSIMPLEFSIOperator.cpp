
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
				   	   	   	   	  	    const matrixEpetraPtr_Type & C1)
{
	M_C1transpose = C1transpose;
	M_C2transpose = C2transpose;
	M_C2 = C2;
	M_C1 = C1;
}

void aSIMPLEFSIOperator::setOptions(const Teuchos::ParameterList& solversOptions)
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

void aSIMPLEFSIOperator::setStructureMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_structureMomentumOptions = _oList;
}

void aSIMPLEFSIOperator::setGeometryOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_geometryOptions = _oList;
}

void aSIMPLEFSIOperator::setFluidMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_fluidMomentumOptions = _oList;
}

void aSIMPLEFSIOperator::setSchurOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurOptions = _oList;
}

void aSIMPLEFSIOperator::setSchurCouplingOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurCouplingOptions = _oList;
}

void aSIMPLEFSIOperator::updateApproximatedStructureMomentumOperator( )
{
	M_approximatedStructureMomentumOperator->SetRowMatrix(M_S->matrixPtr());
	M_approximatedStructureMomentumOperator->SetParameterList(*M_structureMomentumOptions);
	M_approximatedStructureMomentumOperator->Compute();

}

void aSIMPLEFSIOperator::updateApproximatedGeometryOperator( )
{
    M_approximatedGeometryOperator->SetRowMatrix(M_G->matrixPtr());
    M_approximatedGeometryOperator->SetParameterList(*M_geometryOptions);
    M_approximatedGeometryOperator->Compute();
}

void aSIMPLEFSIOperator::updateApproximatedFluidMomentumOperator( )
{
	//M_approximatedFluidMomentumOperator->SetRowMatrix();
	M_approximatedFluidMomentumOperator->SetParameterList(*M_fluidMomentumOptions);
	M_approximatedFluidMomentumOperator->Compute();

}

void aSIMPLEFSIOperator::updateApproximatedSchurComplementOperator( )
{
	//buildShurComplement();
	//M_approximatedSchurComplementOperator->SetRowMatrix();
    M_approximatedSchurComplementOperator->SetParameterList(*M_schurOptions);
    M_approximatedSchurComplementOperator->Compute();
}

void aSIMPLEFSIOperator::updateApproximatedSchurComplementCouplingOperator( )
{
	//M_approximatedSchurComplementCouplingOperator->SetRowMatrix();
	M_approximatedSchurComplementCouplingOperator->SetParameterList(*M_schurCouplingOptions);
	M_approximatedSchurComplementCouplingOperator->Compute();
}

int aSIMPLEFSIOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");
    return 0;
}

} /* end namespace Operators */

} /*end namespace */
