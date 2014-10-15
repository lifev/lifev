/*
 * ApproximatedInvertibleRowMatrix.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: uvilla
 */

#include "ApproximatedInvertibleRowMatrix.hpp"

namespace LifeV
{

namespace Operators
{

ApproximatedInvertibleRowMatrix::ApproximatedInvertibleRowMatrix():
	LinearOperator(),
	usePreconditionerAsApproximatedInverse(false)
{ }

ApproximatedInvertibleRowMatrix::~ApproximatedInvertibleRowMatrix()
{ }

void ApproximatedInvertibleRowMatrix::SetRowMatrix(const rowMatrixPtr_Type & rowMatrix)
{
	ASSERT_PRE(rowMatrix.get() != 0,
			"[ApproximatedInvertibleRowMatrix::SetRowMatrix] rowMatrix should be a valid pointer");
	ASSERT_PRE(rowMatrix->OperatorRangeMap().SameAs(rowMatrix->OperatorDomainMap()),
			"[ApproximatedInvertibleRowMatrix::SetRowMatrix] should be a square matrix");

	M_rowMatrix = rowMatrix;
}

void ApproximatedInvertibleRowMatrix::SetParameterList(const pList_Type pList)
{
	ASSERT_PRE(pList.isParameter("use preconditioner as approximated inverse"),
			"[ApproximatedInvertibleRowMatrix::SetParameterList] Not a valid list");
	ASSERT_PRE(pList.isParameter("preconditioner type"),
			"[ApproximatedInvertibleRowMatrix::SetParameterList] Not a valid list");
	ASSERT_PRE(pList.isSublist("solver"),
			"[ApproximatedInvertibleRowMatrix::SetParameterList] Not a valid list");
	ASSERT_PRE(pList.isSublist("preconditioner"),
			"[ApproximatedInvertibleRowMatrix::SetParameterList] Not a valid list");

	M_pList = pList;

}

int ApproximatedInvertibleRowMatrix::SetUseTranspose(bool /*UseTranspose*/)
{
	return -1;
}

int ApproximatedInvertibleRowMatrix::Compute()
{
	// Allocate the preconditioner and after allocate the solver
	std::string precType( M_pList.get<std::string>("preconditioner type") );

	M_prec.reset( RowMatrixPreconditionerFactory::instance().createObject(precType));
	M_prec->SetRowMatrix(M_rowMatrix);
	M_prec->SetParameterList(M_pList.sublist("preconditioner").sublist(precType));

	EPETRA_CHK_ERR( M_prec->Compute() );

	std::string solverType(M_pList.sublist("solver").get<std::string>("Linear Solver Type"));
	M_linSolver.reset( InvertibleOperatorFactory::instance().createObject(solverType));
	M_linSolver->setOperator(M_rowMatrix);
	M_linSolver->setPreconditioner(M_prec);
	M_linSolver->setParameterList(M_pList.sublist("solver").sublist(solverType));

	usePreconditionerAsApproximatedInverse = M_pList.get<bool>("use preconditioner as approximated inverse");

	return 0;

}

int ApproximatedInvertibleRowMatrix::Apply(const vector_Type& X, vector_Type& Y) const
{
	return M_rowMatrix->Apply(X,Y);
}

int ApproximatedInvertibleRowMatrix::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
	int ierr;
	if(usePreconditionerAsApproximatedInverse)
		ierr = M_prec->ApplyInverse(X,Y);
	else
		ierr = M_linSolver->ApplyInverse(X,Y);

	return ierr;
}

double ApproximatedInvertibleRowMatrix::NormInf() const
{
	return M_rowMatrix->NormInf();
}

const char * ApproximatedInvertibleRowMatrix::Label() const
{
	return "ApproximatedInvertibleRowMatrix";
}

bool ApproximatedInvertibleRowMatrix::UseTranspose() const
{
	return -1;
}

bool ApproximatedInvertibleRowMatrix::HasNormInf() const
{
	return M_rowMatrix->HasNormInf();
}

const LinearOperator::comm_Type & ApproximatedInvertibleRowMatrix::Comm() const
{
	return M_rowMatrix->Comm();
}

const LinearOperator::map_Type & ApproximatedInvertibleRowMatrix::OperatorDomainMap() const
{
	return M_rowMatrix->OperatorDomainMap();
}

const LinearOperator::map_Type & ApproximatedInvertibleRowMatrix::OperatorRangeMap() const
{
	return M_rowMatrix->OperatorRangeMap();
}

} // Namespace Operators

} // Namespace LifeV
