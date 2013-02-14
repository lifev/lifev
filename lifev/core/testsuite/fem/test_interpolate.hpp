//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Interpolate test

    @author Mauro Perego <mperego@fsu.edu>
    @contributor
    @maintainer Mauro Perego <mperego@fsu.edu>

    @date 07-01-2010

The program tests the interpolation methods between different finite elements (mostly between scalar continuous finite elements).
Also it test the interpolation of an analytical function into a finite element space.
 */

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <string>
#include <fstream>
#include <boost/bind.hpp>

using namespace LifeV;

Real linearFunction (const Real& t, const Real& x, const Real& y, const Real& z, const ID& ic);

Real linearBubbleFunction (const Real& t, const Real& x, const Real& y, const Real& z, const ID& ic);

Real quadraticFunction (const Real& t, const Real& x, const Real& y, const Real& z, const ID& ic);

Real bilinearFunction (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& ic);


//! Interpolate the analytical function into the FE spaces specified in originalFeSpaceVecPtr, obtaining the set of FE vectors.
//! then, interpolate this set of FE vectors into vectors having FE spaces specified in finalFeSpaceVecPtr (all the possible combination).
//! Compute the errors between these FE vectors and the analytical solution, and compare them with the ones provided in the array errorArray.
//! return true if all the errors are equal to the errors in errorArray, within a tolerance eps. return false otherwise.
template<typename MeshType, typename MapType, typename Fct>
bool check_interpolate ( std::vector< boost::shared_ptr < FESpace<MeshType, MapType> > >& originalFeSpaceVecPtr,
                         std::vector< boost::shared_ptr < FESpace<MeshType, MapType> > >& finalFeSpaceVecPtr,
                         const MapEpetraType& outputMapType,  Fct& function,
                         const Real errorArray [], const string stringArray [], Real eps, Real time, UInt verbose)
{
    std::vector< boost::shared_ptr <VectorEpetra> > interpVecPtr (originalFeSpaceVecPtr.size() );
    bool check (true);

    for (UInt i = 0; i < originalFeSpaceVecPtr.size(); i++)
    {
        boost::shared_ptr <VectorEpetra> tmp (new VectorEpetra (originalFeSpaceVecPtr[i]->map(), outputMapType) );
        originalFeSpaceVecPtr[i]->interpolate ( static_cast<typename FESpace<MeshType, MapType>::function_Type> ( function ), *tmp, time);
        interpVecPtr[i] = tmp;
    }

    Real err_rel (0);
    for (UInt i = 0; i < originalFeSpaceVecPtr.size(); i++)
        for (UInt j = 0; j < finalFeSpaceVecPtr.size(); j++)
        {
            VectorEpetra interpolated = finalFeSpaceVecPtr[j]->feToFEInterpolate (*originalFeSpaceVecPtr[i], *interpVecPtr[i], outputMapType);
            finalFeSpaceVecPtr[j]->l2Error ( function, VectorEpetra (interpolated, Repeated), time, &err_rel );
            check &= fabs (err_rel - errorArray[finalFeSpaceVecPtr.size() * i + j]) < eps;

            if (verbose)
            {
                UInt index = finalFeSpaceVecPtr.size() * i + j;
                cout.precision (7);
                std::cout << stringArray[index] << ": " << std::setw (15) << std::setprecision (10) << err_rel << " (expected " << errorArray[index] << ")\t";
                std::cout << "\n";
            }
        }
    if (verbose)
    {
        std::cout << std::endl;
    }
    return check;
}
