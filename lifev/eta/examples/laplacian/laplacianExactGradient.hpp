/*
 * laplacianExactGradient.hpp
 *
 *  Created on: Jan 20, 2015
 *      Author: dalsanto
 */

#ifndef _LAPLACIANEXACTGRADIENT_HPP_
#define _LAPLACIANEXACTGRADIENT_HPP_

#include <lifev/core/LifeV.hpp>


namespace LifeV
{

class laplacianExactGradient
{
public:
    typedef VectorSmall<3> return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

        Real x = spaceCoordinates[0];
        Real y = spaceCoordinates[1];
        Real z = spaceCoordinates[2];

        return_Type v;

        v[0] = M_PI * cos( M_PI * x ) * sin( M_PI * y ) * sin( M_PI * z );
        v[1] = M_PI * sin( M_PI * x ) * cos( M_PI * y ) * sin( M_PI * z );
        v[2] = M_PI * sin( M_PI * x ) * sin( M_PI * y ) * cos( M_PI * z );

        return v;
    }

    laplacianExactGradient() {}
    laplacianExactGradient (const laplacianExactGradient&) {}
    ~laplacianExactGradient() {}
};

} // end namespace LifeV





#endif /* _LAPLACIANEXACTGRADIENT_HPP_ */
