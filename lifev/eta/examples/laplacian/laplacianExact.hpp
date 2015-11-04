/*
 * laplacianExact.hpp
 *
 *  Created on: Jan 19, 2015
 *      Author: dalsanto
 */

#ifndef _LAPLACIANEXACT_HPP_
#define _LAPLACIANEXACT_HPP_

#include <lifev/core/LifeV.hpp>


namespace LifeV
{

class laplacianExact
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

        Real x = spaceCoordinates[0];
        Real y = spaceCoordinates[1];
        Real z = spaceCoordinates[2];

        return sin( M_PI * x ) * sin( M_PI * y ) * sin( M_PI * z );
    }

    laplacianExact() {}
    laplacianExact (const laplacianExact&) {}
    ~laplacianExact() {}
};

} // end namespace LifeV




#endif /* _LAPLACIANEXACT_HPP_ */
