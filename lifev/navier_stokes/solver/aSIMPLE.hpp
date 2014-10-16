/*
 * OseenOperatorManagerLSC.hpp
 *
 *  Created on: Sep 29, 2010
 *      Author: uvilla
 */

#ifndef aSIMPLE_HPP_
#define aSIMPLE_HPP_

#include <lifev/navier_stokes/solver/SolverManager.hpp>

#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <lifev/core/linear_algebra/LumpedOperator.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>

namespace LifeV
{

class aSIMPLE : public SolverManager
{
public:

	aSIMPLE();
	virtual ~aSIMPLE(){};

protected:

    void getMatrices(const matrix_Type& F, const matrix_Type& B, const matrix_Type& Btranspose);
    
	void updateApproximatedMomentumOperator( );

	void updateApproximatedSchurComplementOperator( );

	Operators::ApproximatedInvertibleRowMatrix * M_aSIMPLEapproximatedMomentumOperator;
    
    Operators::ApproximatedInvertibleRowMatrix * M_aSIMPLEapproximatedSchurComplementOperator;
};

inline SolverManager* create_aSIMPLE() { return new aSIMPLE(); }
namespace
{
	static bool register_aSIMPLE = SolverFactory::instance().registerProduct( "aSIMPLE", &create_aSIMPLE );
}


} /*end namespace */



#endif /* aSIMPLE_HPP_ */
