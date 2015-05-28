/*
 * laplacianFunctor.hpp
 *
 *  Created on: Jan 20, 2015
 *      Author: dalsanto
 */

#ifndef _LAPLACIANFUNCTOR_HPP_
#define _LAPLACIANFUNCTOR_HPP_

#include <lifev/core/LifeV.hpp>


namespace LifeV
{
template < typename Return_Type>
class laplacianFunctor
{

public:
    typedef Return_Type return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

        Real x = spaceCoordinates[0];
        Real y = spaceCoordinates[1];
        Real z = spaceCoordinates[2];

        return myFunction( 0, x, y, z, 0 );

    }

//    laplacianFunctor() {}
//    laplacianFunctor (const laplacianFunctor&) {}
    laplacianFunctor ( return_Type (*newFunction)(const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/) )
        : myFunction(newFunction) {  }
    ~laplacianFunctor() {}

private:

    return_Type (*myFunction)( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/ );

};

} // end namespace LifeV


#endif /* LIFEV_ETA_EXAMPLES_LAPLACIAN_LAPLACIANFUNCTOR_HPP_ */
