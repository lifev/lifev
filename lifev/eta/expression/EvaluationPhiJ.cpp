#include <lifev/eta/expression/EvaluationPhiJ.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{
const flag_Type EvaluationPhiJ<1>::S_globalUpdateFlag = ET_UPDATE_NONE;

const flag_Type EvaluationPhiJ<1>::S_testUpdateFlag = ET_UPDATE_NONE;

const flag_Type EvaluationPhiJ<1>::S_solutionUpdateFlag = ET_UPDATE_PHI;

} // Namespace ExpressionAssembly

} // Namespace LifeV
