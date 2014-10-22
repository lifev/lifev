
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

int aSIMPLEOperator::ApplyInverse(const vector_Type &X, vector_Type &Y) const
{
    return 0;
}

// show information about the class
void aSIMPLEOperator::showMe(){
    //std::cout<<"Dimension u: "<< M_Bt->NumGlobalRows()<<
    //", Dimension p: "<<M_Bt->NumGlobalCols()<<std::endl;
    //std::cout<<"Pressure correction order "<<M_p<<std::endl;
}
    
void aSIMPLEOperator::setUp(const matrixEpetraPtr_Type & F,
                            const matrixEpetraPtr_Type & B,
                            const matrixEpetraPtr_Type & Btranspose,
                            const commPtr_Type & comm)
{
    M_F = F;
    M_B = B;
    M_Btranspose = Btranspose;
    M_comm = comm;
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
    M_approximatedSchurComplementOperator->SetRowMatrix(M_F->matrixPtr());
    M_approximatedSchurComplementOperator->SetParameterList(*M_schurOptions);
    M_approximatedSchurComplementOperator->Compute();
}

} /* end namespace Operators */
} /*end namespace */
