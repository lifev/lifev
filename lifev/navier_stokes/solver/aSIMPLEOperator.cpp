
#include <lifev/navier_stokes/solver/aSIMPLEOperator.hpp>

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

aSIMPLEOperator::aSIMPLEOperator():
    M_label("aSIMPLEOperator"),
    M_useTranspose(false),
    M_approximatedMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
    M_approximatedSchurComplementOperator(new Operators::ApproximatedInvertibleRowMatrix)
{

}

aSIMPLEOperator::~aSIMPLEOperator()
{

}

// show information about the class
void aSIMPLEOperator::showMe(){
    //std::cout<<"Dimension u: "<< M_Bt->NumGlobalRows()<<
    //", Dimension p: "<<M_Bt->NumGlobalCols()<<std::endl;
    //std::cout<<"Pressure correction order "<<M_p<<std::endl;
}

void aSIMPLEOperator::setUp(const matrixEpetraPtr_Type & F,
                            const matrixEpetraPtr_Type & B,
                            const matrixEpetraPtr_Type & Btranspose)
{
    M_F = F;
    M_B = B;
    M_Btranspose = Btranspose;
    M_comm = F->map().commPtr();
    M_monolithicMap.reset( new mapEpetra_Type ( M_F->map() ) );
    *M_monolithicMap += M_B->map();
}

void aSIMPLEOperator::setOptions(const Teuchos::ParameterList& solversOptions)
{
    boost::shared_ptr<Teuchos::ParameterList> schurOptions;
    schurOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("ApproximatedSchurOperator")) );
    setSchurOptions(schurOptions);

    boost::shared_ptr<Teuchos::ParameterList> momentumOptions;
    momentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("MomentumOperator")) );
    setMomentumOptions(momentumOptions);
}

void aSIMPLEOperator::setMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_momentumOptions = _oList;
}


void aSIMPLEOperator::setSchurOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_schurOptions = _oList;
}

void aSIMPLEOperator::updateApproximatedMomentumOperator( )
{
    M_approximatedMomentumOperator->SetRowMatrix(M_F->matrixPtr());
    M_approximatedMomentumOperator->SetParameterList(*M_momentumOptions);
    M_approximatedMomentumOperator->Compute();

}

void aSIMPLEOperator::updateApproximatedSchurComplementOperator( )
{
    buildShurComplement();
    M_approximatedSchurComplementOperator->SetRowMatrix(M_schurComplement->matrixPtr());
    M_approximatedSchurComplementOperator->SetParameterList(*M_schurOptions);
    M_approximatedSchurComplementOperator->Compute();
}

void aSIMPLEOperator::buildShurComplement( )
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

    M_DBT.reset ( new matrixEpetra_Type( *M_Btranspose ) );
    M_DBT->matrixPtr()->LeftScale(*M_invD);
}

int aSIMPLEOperator::ApplyInverse( VectorEpetra_Type const& X_velocity,
                                   VectorEpetra_Type const& X_pressure,
                                   VectorEpetra_Type & Y_velocity,
                                   VectorEpetra_Type & Y_pressure) const
{
    VectorEpetra_Type Z ( X_velocity.map(), Unique );
    M_approximatedMomentumOperator->ApplyInverse(X_velocity.epetraVector(), Z.epetraVector() );

    VectorEpetra_Type K ( X_pressure, Unique );
    K -= *M_B*Z;

    VectorEpetra_Type W ( X_pressure.map(), Unique );
    M_approximatedSchurComplementOperator->ApplyInverse(K.epetraVector(), W.epetraVector());
    W *= (-1.0);

    Y_velocity = Z;

    Y_velocity -= *M_DBT*W;

    Y_pressure = W;

    return 0;

}


int aSIMPLEOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);

    // split the input vector into the velocity and pressure components
    VectorEpetra_Type X_velocity(M_F->map(), Unique);
    VectorEpetra_Type X_pressure(M_B->map(), Unique);
    VectorEpetra_Type Y_velocity( X_velocity.map(), Unique );
    VectorEpetra_Type Y_pressure ( X_pressure.map(), Unique);

    // gather input values
    X_velocity.subset(X_vectorEpetra, M_F->map(), 0, 0);
    X_pressure.subset(X_vectorEpetra, M_B->map(), M_F->map().mapSize(), 0);

    ApplyInverse( X_velocity, X_pressure, Y_velocity, Y_pressure);

    // output vector
    VectorEpetra_Type Y_vectorEpetra(M_monolithicMap, Unique);

    // Copy the individual parts inside
    Y_vectorEpetra.subset(Y_velocity, X_velocity.map(), 0, 0);
    Y_vectorEpetra.subset(Y_pressure, X_pressure.map(), 0, X_velocity.map().mapSize() );

    Y = dynamic_cast<Epetra_MultiVector&>( Y_vectorEpetra.epetraVector() );

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
