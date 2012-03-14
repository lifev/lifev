#include <lifev/operator/linear_algebra/AztecooOperator.hpp>

namespace LifeV
{

namespace Operators
{

AztecooOperator::AztecooOperator():
		InvertibleOperator(),
		M_linSolver(new SolverType)
{
	M_name = "AztecooOperator";
}

int AztecooOperator::doApplyInverse(const vector_Type& X, vector_Type& Y) const
{
	vector_Type Xcopy(X);
	Y.PutScalar(0.0);
	M_linSolver->SetUserOperator( M_oper.get());
	M_linSolver->SetParameters(*M_pList);
	M_linSolver->SetRHS(&Xcopy);
	M_linSolver->SetLHS(&Y);
	if(M_prec.get() != 0)
		EPETRA_CHK_ERR(M_linSolver->SetPrecOperator( (Epetra_Operator*) M_prec.get()) );

	int maxIter( M_pList->get<int>("max_iter"));
	double tol(  M_pList->get<double>("tol"));

	return M_linSolver->Iterate(maxIter, tol);
}

} // Namespace Operators

} // Namespace LifeV
