/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Copyright (C) 2015 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/*!
 @file laplacianSource.hpp
 @author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
 @date 2015-01-19
 */

#ifndef LAPLACIANSOURCE_HPP
#define LAPLACIANSOURCE_HPP

#include <lifev/core/LifeV.hpp>


namespace LifeV
{

class laplacianSource
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

        Real x = spaceCoordinates[0];
        Real y = spaceCoordinates[1];
        Real z = spaceCoordinates[2];
        
        return 3 * M_PI * M_PI * sin( M_PI * x ) * sin( M_PI * y ) * sin( M_PI * z );
    }

    laplacianSource() {}
    laplacianSource (const laplacianSource&) {}
    ~laplacianSource() {}
};

} // end namespace LifeV

#endif // LAPLACIANSOURCE_HPP
