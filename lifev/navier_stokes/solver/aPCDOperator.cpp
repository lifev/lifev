#include <lifev/navier_stokes/solver/aPCDOperator.hpp>

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

aPCDOperator::aPCDOperator():
    M_label("aPCDOperator"),
    M_useTranspose(false),
    M_approximatedMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedSchurComplementOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedPressureMassOperator(new Operators::ApproximatedInvertibleRowMatrix)
{

}

aPCDOperator::~aPCDOperator()
{

}

// show information about the class
void aPCDOperator::showMe(){
    //std::cout<<"Dimension u: "<< M_Bt->NumGlobalRows()<<
    //", Dimension p: "<<M_Bt->NumGlobalCols()<<std::endl;
    //std::cout<<"Pressure correction order "<<M_p<<std::endl;
}

void aPCDOperator::setUp(const matrixEpetraPtr_Type & F,
                         const matrixEpetraPtr_Type & B,
                         const matrixEpetraPtr_Type & Btranspose,
                         const matrixEpetraPtr_Type & Fp,
                         const matrixEpetraPtr_Type & Mp,
                         const matrixEpetraPtr_Type & Mu)
{
    M_F = F;
    M_B = B;
    M_Btranspose = Btranspose;
    M_Fp = Fp;
    M_Mp = Mp;
    M_Mu = Mu;

    M_comm = F->map().commPtr();
    M_monolithicMap.reset( new mapEpetra_Type ( M_F->map() ) );
    *M_monolithicMap += M_B->map();

    M_Zu.reset(new VectorEpetra_Type( F->map(), Unique ) );
    M_Zp.reset(new VectorEpetra_Type( B->map(), Unique ) );
    M_X_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique) );
    M_X_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
    M_Y_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique ) );
    M_Y_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
}

void aPCDOperator::setOptions(const Teuchos::ParameterList& solversOptions)
{
    boost::shared_ptr<Teuchos::ParameterList> schurOptions;
    schurOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("ApproximatedSchurOperator")) );
    setSchurOptions(schurOptions);

    boost::shared_ptr<Teuchos::ParameterList> momentumOptions;
    momentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("MomentumOperator")) );
    setMomentumOptions(momentumOptions);

    boost::shared_ptr<Teuchos::ParameterList> pressureMassOptions;
    pressureMassOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("PressureMassOperator")) );
    setPressureMassOptions(pressureMassOptions);
}

void aPCDOperator::setMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_momentumOptions = _oList;
}


void aPCDOperator::setSchurOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurOptions = _oList;
}

void aPCDOperator::setPressureMassOptions(const parameterListPtr_Type & _oList)
{
	ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
	M_pressureMassOptions = _oList;
}

void aPCDOperator::updateApproximatedMomentumOperator( )
{
    M_approximatedMomentumOperator->SetRowMatrix(M_F->matrixPtr());
    M_approximatedMomentumOperator->SetParameterList(*M_momentumOptions);
    M_approximatedMomentumOperator->Compute();

}

void aPCDOperator::updateApproximatedPressureMassOperator( )
{
	M_approximatedPressureMassOperator->SetRowMatrix(M_Mp->matrixPtr());
    M_approximatedPressureMassOperator->SetParameterList(*M_pressureMassOptions);
    M_approximatedPressureMassOperator->Compute();

}

void aPCDOperator::updateApproximatedSchurComplementOperator( )
{
    buildShurComplement();
    M_approximatedSchurComplementOperator->SetRowMatrix(M_schurComplement->matrixPtr());
    M_approximatedSchurComplementOperator->SetParameterList(*M_schurOptions);
    M_approximatedSchurComplementOperator->Compute();
}

void aPCDOperator::buildShurComplement( )
{
    Epetra_Vector diag( M_Mu->matrixPtr()->OperatorRangeMap() );
    M_invDiagMassVelocity.reset(new Epetra_Vector( M_Mu->matrixPtr()->OperatorRangeMap() ) );

    // extracting diag(Mu)
    M_Mu->matrixPtr()->ExtractDiagonalCopy(diag);

    // computing diag(Mu)^{-1}
    M_invDiagMassVelocity->Reciprocal(diag);

    // computing diag(Mu)^{-1}*M_Btranspose
    matrixEpetra_Type MuBT (*M_Btranspose);
    MuBT.matrixPtr()->LeftScale(*M_invDiagMassVelocity);

    M_schurComplement.reset ( new matrixEpetra_Type ( M_B->map() ) );

    // computing M_B*(diag(Mu)^{-1}*M_Btranspose)
    M_B->multiply (false, MuBT, false, *M_schurComplement, false);
    M_schurComplement->globalAssemble();
}

inline int aPCDOperator::ApplyInverse( VectorEpetra_Type const& X_velocity,
                                   VectorEpetra_Type const& X_pressure,
                                   VectorEpetra_Type & Y_velocity,
                                   VectorEpetra_Type & Y_pressure) const
{
	M_approximatedPressureMassOperator->ApplyInverse(X_pressure.epetraVector(), M_Zp->epetraVector() );

	M_approximatedSchurComplementOperator->ApplyInverse((*M_Fp*(*M_Zp)).epetraVector(), Y_pressure.epetraVector() );

	Y_pressure *= -1.0;

	*M_Zu = X_velocity - (*M_Btranspose*Y_pressure);

	M_approximatedMomentumOperator->ApplyInverse( M_Zu->epetraVector(), Y_velocity.epetraVector() );

    return 0;

}


int aPCDOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);

    // gather input values
    M_X_velocity->subset(X_vectorEpetra, M_F->map(), 0, 0);
    M_X_pressure->subset(X_vectorEpetra, M_B->map(), M_F->map().mapSize(), 0);

    ApplyInverse( *M_X_velocity, *M_X_pressure, *M_Y_velocity, *M_Y_pressure);

    // output vector
    VectorEpetra_Type Y_vectorEpetra(M_monolithicMap, Unique);

    // Copy the individual parts inside
    Y_vectorEpetra.subset(*M_Y_velocity, M_X_velocity->map(), 0, 0);
    Y_vectorEpetra.subset(*M_Y_pressure, M_X_pressure->map(), 0, M_X_velocity->map().mapSize() );

    Y = dynamic_cast<Epetra_MultiVector&>( Y_vectorEpetra.epetraVector() );

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
