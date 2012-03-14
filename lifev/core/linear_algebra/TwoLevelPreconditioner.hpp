/*
 * TwoLevelPreconditioner.hpp
 *
 *  Created on: Oct 12, 2011
 *      Author: uvilla
 */

#ifndef TWOLEVELPRECONDITIONER_HPP_
#define TWOLEVELPRECONDITIONER_HPP_

#include <lifev/operator/linear_algebra/RowMatrixPreconditioner.hpp>

namespace LifeV
{

namespace Operators
{

class TwoLevelPreconditioner : public RowMatrixPreconditioner
{
public:
	TwoLevelPreconditioner();
	virtual ~TwoLevelPreconditioner();

protected:
    virtual int myCompute();

};

inline RowMatrixPreconditioner* createTwoLevel() { return new TwoLevelPreconditioner(); }
namespace
{
	static bool registerTL = RowMatrixPreconditionerFactory::instance().registerProduct( "TwoLevel", &createTwoLevel );
}


}

}

#endif /* TWOLEVELPRECONDITIONER_HPP_ */
