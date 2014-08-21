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

#ifndef ANALYTICALSOLUTION_H_
#define ANALYTICALSOLUTION_H_

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

class AnalyticalSolution
{
public:
    // Give to the class a "functor" interface
    Real operator() (const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic) const
    {
        return u_ex ( t, x, y, z, ic );
    }
    // This method is required by the error calculation routines
    Real grad (UInt icoor, const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic) const
    {
        return grad_ex ( icoor, t, x, y, z, ic );
    }
    // The exact solution to the problem
    static Real u_ex (const Real& t, const Real& x, const Real& y, const Real& z, const UInt& i)
    {
    	Real a = 1.57;

    	Real d = 0.78;

    	Real nu = 0.01;

    	Real e = std::exp (-d*d*t*nu);

    	switch (i)
    	{
    	case 0:
    		return -a * e * ( std::exp (a * x) * std::sin (a * y + d * z) + std::exp (a * z) * std::cos (a * x + d * y) );
    	case 1:
    		return -a * e * ( std::exp (a * y) * std::sin (a * z + d * x) + std::exp (a * x) * std::cos (a * y + d * z) );
    	case 2:
    		return -a * e * ( std::exp (a * z) * std::sin (a * x + d * y) + std::exp (a * y) * std::cos (a * z + d * x) );
    	default:
    		exit (1);
    	}
    }
    // The gradient of the exact solution
    static Real grad_ex (const UInt& icoor, const Real& t, const Real& x, const Real& y, const Real& z, const UInt& i)
    {
    	Real a = 1.57;

    	Real d = 0.78;

    	Real nu = 0.01;

    	Real e = std::exp (-d* d* t * nu);

    	switch (icoor)
    	{
    	case 0:
    		switch (i)
    		{
    		case 0:
    			return -a* e * ( a* std::exp (a* x) * std::sin (a* y + d* z) - a* std::exp (a* z) * std::sin (a* x + d* y) );
    		case 1:
    			return -a* e * ( d* std::exp (a* y) * std::cos (a* z + d* x) + a* std::exp (a* x) * std::cos (a* y + d* z) );
    		case 2:
    			return -a* e * ( a* std::exp (a* z) * std::cos (a* x + d* y) - d* std::exp (a* y) * std::sin (a* z + d* x) );
    		default:
    			exit (1);
    		}
    		case 1:   // u_y
    		switch (i)
    		{
    		case 0:
    			return -a* e * ( a* std::exp (a* x) * std::cos (a* y + d* z) - d* std::exp (a* z) * std::sin (a* x + d* y) );
    		case 1:
    			return -a* e * ( a* std::exp (a* y) * std::sin (a* z + d* x) - a* std::exp (a* x) * std::sin (a* y + d* z) );
    		case 2:
    			return -a* e * ( d* std::exp (a* z) * std::cos (a* x + d* y) + a* std::exp (a* y) * std::cos (a* z + d* x) );
    		default:
    			exit (1);
    		}
    		case 2:
    			switch (i)
    			{
    			case 0:
    				return -a* e * ( d* std::exp (a* x) * std::cos (a* y + d* z) + a* std::exp (a* z) * std::cos (a* x + d* y) );
    			case 1:
    				return -a* e * ( a* std::exp (a* y) * std::cos (a* z + d* x) - d* std::exp (a* x) * std::sin (a* y + d* z) );
    			case 2:
    				return -a* e * ( a* std::exp (a* z) * std::sin (a* x + d* y) - a* std::exp (a* y) * std::sin (a* z + d* x) );
    			default:
    				exit (1);
    			}
    			default:
    				exit (1);
    	}
    	return 1.;
    }

};

} //namespace LifeV

#endif
