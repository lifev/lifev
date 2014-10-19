/*
 * OseenOperatorManagerLSC.hpp
 *
 *  Created on: Sep 29, 2010
 *      Author: uvilla
 */

#ifndef solverManager_aSIMPLE_HPP_
#define solverManager_aSIMPLE_HPP_

#include <lifev/navier_stokes/solver/SolverManager.hpp>

#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <lifev/navier_stokes/solver/aSIMPLE.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>

namespace LifeV
{

class solverManager_aSIMPLE : public SolverManager
{
public:

	solverManager_aSIMPLE();
	virtual ~solverManager_aSIMPLE(){};

protected:

    void getMatrices(const matrixPtr_Type& F, const matrixPtr_Type& B, const matrixPtr_Type& Btranspose);
    
	void updateApproximatedMomentumOperator( );

	void updateApproximatedSchurComplementOperator( );

	Operators::ApproximatedInvertibleRowMatrix * M_aSIMPLEapproximatedMomentumOperator;
    
    Operators::ApproximatedInvertibleRowMatrix * M_aSIMPLEapproximatedSchurComplementOperator;
};

inline SolverManager* create_aSIMPLE() { return new solverManager_aSIMPLE(); }
namespace
{
	static bool register_aSIMPLE = SolverFactory::instance().registerProduct( "aSIMPLE", &create_aSIMPLE );
}


} /*end namespace */



#endif /* solverManager_aSIMPLE_HPP_ */
