/*
 * IfpackPreconditioner.hpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */

#ifndef IFPACKPRECONDITIONER_HPP_
#define IFPACKPRECONDITIONER_HPP_

#include <lifev/operator/linear_algebra/RowMatrixPreconditioner.hpp>

namespace LifeV
{

namespace Operators
{

class IfpackPreconditioner : public RowMatrixPreconditioner
{

public:
	IfpackPreconditioner();
	virtual ~IfpackPreconditioner();

protected:
	int myCompute();

};

inline RowMatrixPreconditioner* createIfpack() { return new IfpackPreconditioner(); }
namespace
{
	static bool registerIF = RowMatrixPreconditionerFactory::instance().registerProduct( "Ifpack", &createIfpack );
}

} /* end Operators namespace */

} /* end LifeV namespace */



#endif /* IFPACKPRECONDITIONER_HPP_ */
