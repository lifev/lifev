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

#include "SetOfFun.hpp"

//      IMPLEMENTATION
namespace LifeV
{
SetofFun::SetofFun():
	_i(31),
	_integral(_i+1),
	_name(_i+1),
	_degree(_i+1)
	{
_integral[0]= 1.;			_name[0] = "p(x) = 1          ";	_degree[0] = 0;
_integral[1]= .5;			_name[1] = "p(x) = x          ";	_degree[1] = 1;
_integral[2]= 1./3.;        _name[2] = "p(x) = x^2        ";	_degree[2] = 2;
_integral[3]= .25;			_name[3] = "p(x) = x*y        ";	_degree[3] = 2;
_integral[4]= .25;			_name[4] = "p(x) = x^3        ";	_degree[4] = 3;
_integral[5]= 1./6.;		_name[5] = "p(x) = x^2*y      ";	_degree[5] = 3;
_integral[6]= .125;			_name[6] = "p(x) = x*y*z      ";	_degree[6] = 3;
_integral[7]= 0.2;			_name[7] = "p(x) = x^4        ";	_degree[7] = 4;
_integral[8]= 1./8.;		_name[8] = "p(x) = x^3*y      ";	_degree[8] = 4;
_integral[9]= 1./9.;		_name[9] = "p(x) = x^2*y^2    ";	_degree[9] = 4;
_integral[10]= 1./12.;		_name[10] ="p(x) = x^2*y*z    ";	_degree[10]= 4;
_integral[11]= 1./6.;		_name[11] ="p(x) = x^5        ";	_degree[11]= 5;
_integral[12]= 1./10.;		_name[12] ="p(x) = x^4*y      ";	_degree[12]= 5;
_integral[13]= 1./12.;		_name[13] ="p(x) = x^3*y^2    ";	_degree[13]= 5;
_integral[14]= 1./16.;		_name[14] ="p(x) = x^3*y*z    ";	_degree[14]= 5;
_integral[15]= 1./18.;		_name[15] ="p(x) = x^2*y^2*z  ";	_degree[15]= 5;
_integral[16]= 1./7.;		_name[16] ="p(x) = x^6        ";	_degree[16]= 6;
_integral[17]= 1./12.;		_name[17] ="p(x) = x^5*y      ";	_degree[17]= 6;
_integral[18]= 1./15.;		_name[18] ="p(x) = x^4*y^2    ";	_degree[18]= 6;
_integral[19]= 1./16.;		_name[19] ="p(x) = x^3*y^3    ";	_degree[19]= 6;
_integral[20]= 1./20.;		_name[20] ="p(x) = x^4*y*z    ";	_degree[20]= 6;
_integral[21]= 1./24.;		_name[21] ="p(x) = x^3*y^2*z  ";	_degree[21]= 6;
_integral[22]= 1./27.;		_name[22] ="p(x) = x^2*y^2*z^2";	_degree[22]= 6;
_integral[23]= 1./8.;		_name[23] ="p(x) = x^7        ";	_degree[23]= 7;
_integral[24]= 1./14.;		_name[24] ="p(x) = x^6*y      ";	_degree[24]= 7;
_integral[25]= 1./18.;		_name[25] ="p(x) = x^5*y^2    ";	_degree[25]= 7;
_integral[26]= 1./20.;		_name[26] ="p(x) = x^4*y^3    ";	_degree[26]= 7;
_integral[27]= 1./24.;		_name[27] ="p(x) = x^5*y*z    ";	_degree[27]= 7;
_integral[28]= 1./30.;		_name[28] ="p(x) = x^4*y^2*z  ";	_degree[28]= 7;
_integral[29]= 1./32.;		_name[29] ="p(x) = x^3*y^3*z  ";	_degree[29]= 7;
_integral[30]= 1./36.;		_name[30] ="p(x) = x^3*y^2*z^2";	_degree[30]= 7;
_integral[31]= 3.12912420246652;	_name[31] ="f(x) = exp(x^2+y^2+z^2)"; _degree[31]=100;
}

Real SetofFun::val(int fun, Real& x,Real& y, Real& z){
switch(fun){
case 0:
	return 1.;
case 1:
	return x;
case 2:
	return x*x;
case 3:
	return x*y;
case 4:
	return x*x*x;
case 5:
	return x*x*y;
case 6:
	return x*y*z;
case 7:
	return x*x*x*x;
case 8:
	return x*x*x*y;
case 9:
	return x*x*y*y;
case 10:
	return x*x*y*z;
case 11:
	return x*x*x*x*x;
case 12:
	return x*x*x*x*y;
case 13:
	return x*x*x*y*y;
case 14:
	return x*x*x*y*z;
case 15:
	return x*x*y*y*z;
case 16:
	return x*x*x*x*x*x;
case 17:
	return x*x*x*x*x*y;
case 18:
	return x*x*x*x*y*y;
case 19:
	return x*x*x*y*y*y;
case 20:
	return x*x*x*x*y*z;
case 21:
	return x*x*x*y*y*z;
case 22:
	return x*x*y*y*z*z;
case 23:
	return x*x*x*x*x*x*x;
case 24:
	return x*x*x*x*x*x*y;
case 25:
	return x*x*x*x*x*y*y;
case 26:
	return x*x*x*x*y*y*y;
case 27:
	return x*x*x*x*x*y*z;
case 28:
	return x*x*x*x*y*y*z;
case 29:
	return x*x*x*y*y*y*z;
case 30:
	return x*x*x*y*y*z*z;
case 31:
	return exp(x*x+y*y+z*z);
default:
	return 0.;
}
}

std::string SetofFun::name(UInt fun){
	return _name[fun];
}

UInt SetofFun::degree(UInt fun){
	return _degree[fun];
}

Real SetofFun::ex_int(UInt fun){
	return _integral[fun];
}

UInt SetofFun::nfun(){return _i;}

} /*end namespace */
