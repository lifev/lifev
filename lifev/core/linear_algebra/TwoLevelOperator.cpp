/*
 * TwoLevelOperator.cpp
 *
 *  Created on: Oct 9, 2011
 *      Author: uvilla
 */

#include "TwoLevelOperator.hpp"

namespace LifeV
{

namespace Operators
{

TwoLevelOperator::TwoLevelOperator():
		LinearOperatorAlgebra()
{ }

TwoLevelOperator::~TwoLevelOperator()
{ }

void TwoLevelOperator::SetFineLevelOperator(const operatorPtr_Type & fineLevelOper)
{
	ASSERT_PRE(fineLevelOper.get() != 0, "TwoLevelOperator::SetFineLevelOperator Invalid Pointer");
	M_fineLevelOper = fineLevelOper;
}

void TwoLevelOperator::SetSmootherOperator(const operatorPtr_Type & smootherOper)
{
	ASSERT_PRE(smootherOper.get() != 0, "TwoLevelOperator::SetFineLevelOperator Invalid Pointer");
	M_smootherOper = smootherOper;
}


void TwoLevelOperator::SetCoarseLevelOperator(const operatorPtr_Type & coarseLevelOper )
{
	ASSERT_PRE(coarseLevelOper.get() != 0, "TwoLevelOperator::SetCoarseLevelOperator Invalid Pointer");
	M_coarseLevelOper = coarseLevelOper;
}

void TwoLevelOperator::SetRestrictionOperator(const operatorPtr_Type & restrictionOper )
{
	ASSERT_PRE(restrictionOper.get() !=0, "TwoLevelOperator::SetRestrictionOperator Invalid Pointer");
	M_restrictionOper = restrictionOper;
}

void TwoLevelOperator::SetEstensionOperator(const operatorPtr_Type & estensionOper)
{
	ASSERT_PRE(estensionOper.get() !=0, "TwoLevelOperator::SetEstensionOperator Invalid Pointer");
	M_estensionOper = estensionOper;
}

int TwoLevelOperator::checkConsistency()
{
	int returnValue(0);
	bool verbose( 0 == M_fineLevelOper->Comm().MyPID() );

	//Check that the Comm is the same:
	if( &(M_fineLevelOper->Comm()) != &(M_coarseLevelOper->Comm()))/* == &(M_smootherOper->Comm())
			== &(M_restrictionOper->Comm()) == &(M_estensionOper->Comm())) )*/
	{
		if(verbose)
			std::cout<< "[TwoLevelOperator::checkConsistency] Comm must be the same for all operators! \n";
		--returnValue;
	}

	//Check that maps are correct
	// inv(this) = inv(S) (I - A E inv(Ac) R) (I - A*inv(S))

	if(! M_fineLevelOper->OperatorRangeMap().SameAs(M_restrictionOper->OperatorDomainMap()))
	{
		if(verbose)
			std::cout<< "[TwoLevelOperator::checkConsistency] Not compatible maps! \n";
		--returnValue;
	}

	if( ! M_restrictionOper->OperatorRangeMap().SameAs(M_coarseLevelOper->OperatorRangeMap()))
	{
		if(verbose)
			std::cout<< "[TwoLevelOperator::checkConsistency] Not compatible maps! \n";
		--returnValue;
	}

	if( ! M_coarseLevelOper->OperatorDomainMap().SameAs(M_estensionOper->OperatorDomainMap()))
	{
		if(verbose)
			std::cout<< "[TwoLevelOperator::checkConsistency] Not compatible maps! \n";
		--returnValue;
	}

	if( ! M_estensionOper->OperatorRangeMap().SameAs(M_fineLevelOper->OperatorDomainMap()))
	{
		if(verbose)
			std::cout<< "[TwoLevelOperator::checkConsistency] Not compatible maps! \n";
		--returnValue;
	}

	return returnValue;

}

int TwoLevelOperator::SetUseTranspose(bool /*UseTranspose*/)
{
	return -1.0;
}
//@}

//! @name Mathematical functions
//@{
int TwoLevelOperator::Apply(const vector_Type& X, vector_Type& Y) const
{
	return M_fineLevelOper->Apply(X,Y);
}

// inv(this) = inv(S) (I - A E inv(Ac) R) (I - A*inv(S))
int TwoLevelOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{

	ASSERT_PRE(X.NumVectors() == Y.NumVectors(),
			"[TwoLevelOperator::ApplyInverse] X and Y should have the same number of vectors.");
	ASSERT_PRE(X.Map().SameAs(M_smootherOper->OperatorRangeMap()),
			"[TwoLevelOperator::ApplyInverse] The map of X is not compatible with this operator.");
	ASSERT_PRE(Y.Map().SameAs(M_smootherOper->OperatorDomainMap()),
			"[TwoLevelOperator::ApplyInverse] The map of X is not compatible with this operator.");

	UInt nVectors(X.NumVectors());

	vector_Type Z(Y.Map(), nVectors), Rfine(Y.Map(), nVectors), Yfine(Y.Map(), nVectors);
	vector_Type Rcoarse(M_coarseLevelOper->OperatorDomainMap(), nVectors);
	vector_Type Ycoarse(M_coarseLevelOper->OperatorDomainMap(), nVectors);

	// y_fine = inv(S)*x
	// r_fine = x - A*y_fine
	EPETRA_CHK_ERR(M_smootherOper->ApplyInverse(X,Yfine));
	EPETRA_CHK_ERR(M_fineLevelOper->Apply(Yfine, Z));
	EPETRA_CHK_ERR(Rfine.Update(1.0, X, -1.0, Z, 0.0));

	// y_coarse = inv(Ac)*(R * r_fine)
	EPETRA_CHK_ERR(M_restrictionOper->Apply(Rfine,Rcoarse));
	EPETRA_CHK_ERR(M_coarseLevelOper->ApplyInverse(Rcoarse, Ycoarse));
	// y_fine += E*y_coarse
	EPETRA_CHK_ERR(M_estensionOper->Apply(Ycoarse, Y));
	EPETRA_CHK_ERR(Yfine.Update(1.0, Y, 1.0));
	// r_fine -= A*(E*y_coarse)
	EPETRA_CHK_ERR(M_fineLevelOper->Apply(Y, Z));
	EPETRA_CHK_ERR(Rfine.Update(-1.0, Z, 1.0));

	// y = inv(S)*r_fine
	// y = y_fine + y
	EPETRA_CHK_ERR(M_smootherOper->ApplyInverse(Rfine,Y));
	EPETRA_CHK_ERR(Y.Update(1.0, Yfine, 1.0));

	return 0;
}

double TwoLevelOperator::NormInf() const
{
	return M_fineLevelOper->NormInf();
}
//@}

const char * TwoLevelOperator::Label() const
{
	return "TwoLevelOperator";
}

bool TwoLevelOperator::UseTranspose() const
{
	return false;
}

bool TwoLevelOperator::HasNormInf() const
{
	return M_fineLevelOper->HasNormInf();
}

const LinearOperatorAlgebra::comm_Type & TwoLevelOperator::Comm() const
{
	return M_fineLevelOper->Comm();
}

const LinearOperatorAlgebra::map_Type & TwoLevelOperator::OperatorDomainMap() const
{
	return M_fineLevelOper->OperatorDomainMap();
}

const LinearOperatorAlgebra::map_Type & TwoLevelOperator::OperatorRangeMap() const
{
	return M_fineLevelOper->OperatorRangeMap();
}


} /* end namespace Operators */

} /* end namespace LifeV */
