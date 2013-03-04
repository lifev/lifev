#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

//! Outer product
/*!
@param vector second operand
*/
namespace LifeV
{

MatrixSmall<3, 3>
VectorSmall<3>::outerProduct ( VectorSmall<3> const& vector ) const
{
    MatrixSmall<3, 3> result;

    for ( UInt i = 0; i < 3; i++ )
        for ( UInt j = 0; j < 3; j++ )
        {
            result[i][j] = M_coords[ i ] * vector.M_coords[ j ];
        }
    return result;
}

} //namespace LifeV

