/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/* ========================================================

Definition of all the polynomials we need to integrate to check the degree of exactness.
The last function is an exponential function to check the convergence rate


\author Umberto Villa <uvilla@emory.edu>
\date 02/03/2010
*/

#ifndef SETOFFUN_HPP_
#include <life/lifecore/life.hpp>

namespace LifeV
{
class SetofFun{
	public:
	SetofFun();
	Real val(int fun, Real& x, Real& y, Real& z);
	int degree(int fun);
	Real ex_int(int fun);
	std::string name(int fun);
	int nfun();
	private:
	int _i;
	std::vector<Real> _integral;
	std::vector<std::string> _name;
	std::vector<int> _degree;
};
}
#define SETOFFUN_HPP_


#endif /* SETOFFUN_HPP_ */
