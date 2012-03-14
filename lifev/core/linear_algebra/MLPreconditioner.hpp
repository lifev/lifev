/*
 * MLPreconditioner.hpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */

#ifndef MLPRECONDITIONER_HPP_
#define MLPRECONDITIONER_HPP_

#include <lifev/operator/linear_algebra/RowMatrixPreconditioner.hpp>

namespace LifeV
{

namespace Operators
{

class MLPreconditioner : public RowMatrixPreconditioner
{
public:
	MLPreconditioner();
	virtual ~MLPreconditioner();

protected:
	int myCompute();

};

inline RowMatrixPreconditioner* createML() { return new MLPreconditioner(); }
namespace
{
	static bool registerML = RowMatrixPreconditionerFactory::instance().registerProduct( "ML", &createML );
}


}

}

#endif /* MLPRECONDITIONER_HPP_ */
