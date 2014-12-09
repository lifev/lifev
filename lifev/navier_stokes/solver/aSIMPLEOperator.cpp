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
    M_Z.reset(new VectorEpetra_Type( F->map(), Unique ) );
    M_X_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique) );
    M_X_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
    M_Y_velocity.reset( new VectorEpetra_Type (M_F->map(), Unique ) );
    M_Y_pressure.reset( new VectorEpetra_Type (M_B->map(), Unique) );
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

inline int aSIMPLEOperator::ApplyInverse( VectorEpetra_Type const& X_velocity,
                                   VectorEpetra_Type const& X_pressure,
                                   VectorEpetra_Type & Y_velocity,
                                   VectorEpetra_Type & Y_pressure) const
{
    M_approximatedMomentumOperator->ApplyInverse(X_velocity.epetraVector(), M_Z->epetraVector() );

    M_approximatedSchurComplementOperator->ApplyInverse( (*M_B*(*M_Z) - X_pressure ).epetraVector(), Y_pressure.epetraVector());

    Y_velocity = (*M_Z - *M_DBT*Y_pressure);

    return 0;

}


int aSIMPLEOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
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
