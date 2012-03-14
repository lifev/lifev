#include <lifev/operator/linear_algebra/InvertibleOperator.hpp>
#include <Teuchos_RCPBoostSharedPtrConversions.hpp>

namespace LifeV
{

namespace Operators
{

InvertibleOperator::InvertibleOperator():
	M_name("InvertibleOperator"),
	M_useTranspose(false)
{ }

int InvertibleOperator::SetUseTranspose(bool useTranspose)
{
	M_useTranspose = useTranspose;

	int ierr(0);
	if(M_useTranspose)
		ierr = -1;

	return ierr;
}

void InvertibleOperator::setOperator(const operatorPtr_Type& _oper)
{
	ASSERT_PRE(_oper.get() != this, "Can't self assign");
	ASSERT_PRE(_oper.get() != 0, "Can't assign a null pointer");
	M_oper = Teuchos::rcp(_oper);
	doSetOperator();
}

void InvertibleOperator::setPreconditioner(const operatorPtr_Type& _prec)
{
	ASSERT_PRE(_prec.get() != this, "Self Assignment is forbidden");
	ASSERT_PRE(_prec.get() != 0, "Can't assign a null pointer");
	M_prec = Teuchos::rcp(_prec);
	doSetPreconditioner();
}

void InvertibleOperator::setParameterList(const Teuchos::ParameterList& _pList)
{
	M_pList = Teuchos::rcp(new Teuchos::ParameterList(_pList), true);
	doSetParameterList();
}

int InvertibleOperator::Apply(const vector_Type& X, vector_Type& Y) const
{
	ASSERT_PRE(M_oper.assert_valid_ptr().get() != 0, "M_oper must be assigned");
	ASSERT_PRE(X.Map().SameAs(M_oper->OperatorDomainMap()), "X and domain map do no coincide \n");
	ASSERT_PRE(Y.Map().SameAs(M_oper->OperatorRangeMap()) , "Y and range map do no coincide \n");

	return M_oper->Apply(X,Y);
}

int InvertibleOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
	ASSERT_PRE(M_oper.assert_valid_ptr().get() != 0, "M_oper must be assigned \n");
	ASSERT_PRE(Y.Map().SameAs(M_oper->OperatorDomainMap()), "Y and domain map do no coincide \n");
	ASSERT_PRE(X.Map().SameAs(M_oper->OperatorRangeMap()) , "X and range map do no coincide \n");

	if (M_useTranspose)
		return -1;

	return doApplyInverse(X,Y);
}

} // Namespace Operators

} // Namespace LifeV
