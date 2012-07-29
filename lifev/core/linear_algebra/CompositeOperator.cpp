/*
 * CompositeOperator.cpp
 *
 *  Created on: Jul 28, 2012
 *      Author: uvilla
 */

#include "CompositeOperator.hpp"

namespace LifeV {
namespace Operators {

CompositeOperator::CompositeOperator():
		nAllocatedVectorsInMultVect(-1),
		isAlreadyInverted(false)
{
	vects.push_back( static_cast<vector_Type *>(NULL) );
}

CompositeOperator::~CompositeOperator()
{
	for( UInt i(1); i < vects.size()-1; ++i)
		delete vects[i];
}

int CompositeOperator::pushBack( const operatorPtr_Type & op, const bool inverted_)
{
	if(inverted_ && !op->OperatorDomainMap().SameAs( op->OperatorRangeMap() ) )
		return -1;

	if(ops.size() != 0 && !op->OperatorDomainMap().SameAs( OperatorRangeMap()))
		return -1;

	ops.push_back(op);
	inverted.push_back(inverted_);
	vects.push_back( static_cast<vector_Type *>(NULL) );

	return 0;
}

int CompositeOperator::Apply(const vector_Type& X, vector_Type& Y) const
{

	if(!X.Map().PointSameAs(OperatorDomainMap()))
		return -2;
	if(!Y.Map().PointSameAs(OperatorRangeMap()))
		return -3;

	if(isAlreadyInverted)
		return -1;

	EPETRA_CHK_ERR( allocateTmpVects(X, Y) );
	// Now vects[0] points to X and vects[end] points to Y

	for(UInt i(0); i < ops.size(); ++i)
	{
		if(!inverted[i])
		{
			EPETRA_CHK_ERR( ops[i]->Apply(*vects[i], *vects[i+1]) );
		}
		else
		{
			EPETRA_CHK_ERR( ops[i]->ApplyInverse(*vects[i], *vects[i+1]) );
		}
	}

	// since vects[end] points to Y there is nothing else to be done :)

	return 0;
}

int CompositeOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
	if(!X.Map().PointSameAs(OperatorDomainMap()))
		return -2;
	if(!Y.Map().PointSameAs(OperatorRangeMap()))
		return -3;

	if(!isAlreadyInverted)
		return -1;

	EPETRA_CHK_ERR( allocateTmpVects(X, Y) );
	// Now vects[0] points to X and vects[end] points to Y

	for(UInt i(0); i < ops.size(); ++i)
	{
		if(!inverted[i])
		{
			EPETRA_CHK_ERR( ops[i]->Apply(*vects[i], *vects[i+1]) );
		}
		else
		{
			EPETRA_CHK_ERR( ops[i]->ApplyInverse(*vects[i], *vects[i+1]) );
		}
	}

	return 0;
}

int CompositeOperator::allocateTmpVects(const vector_Type& X, vector_Type& Y) const
{
	if(X.NumVectors() != Y.NumVectors())
		return -1;

	int nVect(X.NumVectors());

	if( nAllocatedVectorsInMultVect != nVect)
	{
		deleteTmpVects();
		for(UInt i(1); i<vects.size()-1; ++i)
			vects[i] = new vector_Type(ops[i]->OperatorDomainMap(), nVect);
	}

	vects.front() = const_cast<vector_Type *>(&X);
	vects.back() = &Y;

	return 0;

}

void CompositeOperator::deleteTmpVects() const
{
for( UInt i(1); i < vects.size()-1; ++i)
	delete vects[i];
}

} /* namespace Operators */
} /* namespace LifeV */
