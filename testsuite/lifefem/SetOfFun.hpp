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
    @brief Quadrature Rule test

	@author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @author Umberto Villa <uvilla@emory.edu>
    @contributor
    @maintainer Umberto Villa <uvilla@emory.edu>

    @date 02-03-2010

Definition of all the polynomials we need to integrate to check the degree of exactness.
The last function is an exponential function to check the convergence rate


*/

#ifndef SETOFFUN_HPP_
#include <life/lifecore/life.hpp>

namespace LifeV
{
class SetofFun
{
public:
    SetofFun();
    Real val(int fun, Real& x, Real& y, Real& z);
    UInt degree(UInt fun);
    Real ex_int(UInt fun);
    std::string name(UInt fun);
    UInt nfun();
private:
    int _i;
    std::vector<Real> _integral;
    std::vector<std::string> _name;
    std::vector<UInt> _degree;
};
}
#define SETOFFUN_HPP_


#endif /* SETOFFUN_HPP_ */
