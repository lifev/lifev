/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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
/*!
  \file oneDFunctions.cpp
  \author Tiziano Passerini
  \author Lucia Mirabella
  \date 04/2007

  \brief This file contains a class for the boundary conditions for 1D tubes.

*/

#include <life/lifefem/oneDBCFunctions.hpp>

namespace LifeV
{
  Reimann::Reimann( const BasicOneDMesh& mesh,
		    const NonLinearFluxFun1D& fluxFun,
		    const ScalVec& U1_thistime, const ScalVec& U2_thistime,  
		    const ScalVec& W1_thistime, const ScalVec& W2_thistime,  
		    const std::string & border, const std::string & var ) :
    OneDBCFunctionBase( ),
    _M_mesh(mesh),  
    _M_fluxFun(fluxFun),
    _M_U1_thistime(U1_thistime),
    _M_U2_thistime(U2_thistime),
    _M_W1_thistime(W1_thistime),
    _M_W2_thistime(W2_thistime),
    _M_dof1D(_M_mesh.numVertices()),
    _M_dimDof(_M_dof1D.numTotalDof()),
    _M_border(border),
    _M_var(var)
  {
    _M_oneDBCFunctionsMapStringValues["left"] = OneDBCLeftBoundary;
    _M_oneDBCFunctionsMapStringValues["right"] = OneDBCRightBoundary;

    switch( _M_oneDBCFunctionsMapStringValues[_M_border] )
      {
      case OneDBCLeftBoundary:
	_M_boundaryDof=0;
	break;
      case OneDBCRightBoundary:	
	_M_boundaryDof=_M_dimDof - 1;
	break;
      default:
	std::cout << "\n[Reimann::Reimann] incorrect boundary identifier: " << _M_border << std::endl;
      }    	
  }

void Reimann::update_U_boundary()
{
    _M_U_boundary.first=_M_U1_thistime(_M_boundaryDof);
    _M_U_boundary.second=_M_U2_thistime(_M_boundaryDof);
    _M_W_boundary.first=_M_W1_thistime(_M_boundaryDof);
    _M_W_boundary.second=_M_W2_thistime(_M_boundaryDof);
} 
	

  Real Reimann::evaluate( const Real& /*time*/ )
  {
	update_U_boundary();

    return ( ( _M_var == "W1" ) ? _M_W_boundary.first : _M_W_boundary.second );

  }



  Compatibility::Compatibility( const BasicOneDMesh& mesh,
				const NonLinearFluxFun1D& fluxFun, const NonLinearSourceFun1D& sourceFun,
				const ScalVec& U1_thistime, const ScalVec& U2_thistime,  
				const ScalVec& W1_thistime, const ScalVec& W2_thistime,  
				const Real& dt, const std::string & border, const std::string & var ) :
    Reimann( mesh, fluxFun, U1_thistime, U2_thistime, W1_thistime, W2_thistime, border, var),
	_M_sourceFun(sourceFun),
    _M_nb_elem(_M_dimDof-1),
    _M_time_step (dt)//data_file("time/timestep",0.1)),
  {
    _M_oneDBCFunctionsMapStringValues["W1"] = OneDBCW1;
    _M_oneDBCFunctionsMapStringValues["W2"] = OneDBCW2;

    switch( _M_oneDBCFunctionsMapStringValues[_M_border] )
      {
      case OneDBCLeftBoundary:
    _M_internalBoundaryDof=_M_boundaryDof + 1;
    _M_boundaryEdge=_M_mesh.edgeList(1);
    _M_boundaryPoint=_M_boundaryEdge.pt1().x();
    _M_internalBdPoint=_M_boundaryEdge.pt2().x();
	break;
      case OneDBCRightBoundary:	
    _M_internalBoundaryDof=_M_boundaryDof - 1;
    _M_boundaryEdge=_M_mesh.edgeList(_M_nb_elem);
    _M_boundaryPoint=_M_boundaryEdge.pt2().x();
    _M_internalBdPoint=_M_boundaryEdge.pt1().x();
	break;
      default:
	std::cout << "\n[Compatibility::Compatibility] incorrect boundary identifier: " << _M_border << std::endl;
      }    	
  }

	
 		
void Compatibility::computeEigenValuesVectors()
{
    _M_fluxFun.jacobian_EigenValues_Vectors( _M_U1_thistime(_M_boundaryDof),
    					_M_U2_thistime(_M_boundaryDof),
					    _M_eigval1, _M_eigval2,
					    _M_left_eigvec1.first, _M_left_eigvec1.second,
					    _M_left_eigvec2.first, _M_left_eigvec2.second,
					    _M_boundaryDof);
}


void Compatibility::update_U_internalBd()
{
    _M_U_internalBd.first=_M_U1_thistime(_M_internalBoundaryDof);
    _M_U_internalBd.second=_M_U2_thistime(_M_internalBoundaryDof);
} 
	

Real Compatibility::extrapolate_L_dot_U( Real const& eigval, Vec2D const& eigvec )
{
	Real L_dot_U_extrap;
    //! Quasi linear source term
    Vec2D qlSource;

    Vec2D U_charact_pt=_interpolLinear( _M_boundaryPoint, _M_internalBdPoint,
					_M_time_step, eigval, _M_U_boundary, _M_U_internalBd);
  
    L_dot_U_extrap=dot( eigvec, U_charact_pt);

    qlSource.first=_M_sourceFun.QuasiLinearSource(U_charact_pt.first,
    												U_charact_pt.second, 1,
    												_M_boundaryDof);
  
    qlSource.second=_M_sourceFun.QuasiLinearSource(U_charact_pt.first,
    												U_charact_pt.second, 2,
												    _M_boundaryDof);
  
    L_dot_U_extrap-=_M_time_step * dot(eigvec, qlSource);

	return L_dot_U_extrap;

}



Real Compatibility::extrapolate_W( OneDBCStringValue const& _W )
{
	Real W_out;
	
	update_U_boundary();
	update_U_internalBd();

	computeEigenValuesVectors();

    switch( _W )
      {
      case OneDBCW1:
	W_out = extrapolate_L_dot_U(_M_eigval1, _M_left_eigvec1) 
		- dot( _M_left_eigvec1, _M_U_boundary ) + _M_W_boundary.first;

	break;	
      case OneDBCW2:
	W_out = extrapolate_L_dot_U(_M_eigval2, _M_left_eigvec2) 
		- dot( _M_left_eigvec2, _M_U_boundary ) + _M_W_boundary.second;

	break;	
      default:
	std::cout << "\n[Compatibility::extrapolate_W] incorrect variable identifier: " << _W << std::endl;
      }

    return W_out;
}	


  Real Compatibility::evaluate( Real const& /*time*/ )
  {
    Debug( 6030 ) << "[Compatibility::evaluate] variable "
    << _M_var << ", code " << _M_oneDBCFunctionsMapStringValues[_M_var] << "\n";
 
  	return extrapolate_W( _M_oneDBCFunctionsMapStringValues[_M_var] );
  }



  Vec2D Compatibility::_interpolLinear(const Real& point_bound, const Real& point_internal,
				       const Real& deltaT, const Real& eigenvalue,
				       const Vec2D& U_bound, const Vec2D& U_intern) const
  {
    Real deltaX = std::abs(point_bound - point_internal);
  
    Real cfl =  eigenvalue * deltaT / deltaX;

    Real weight;   //! weight in the linear approximation

    Debug( 6030 ) << "[Compatibility::_interpolLinear] point_bound "
    << point_bound << ", point_internal " << point_internal
    << ", deltaT " << deltaT << ", deltaX " << deltaX 
    << ", eigenvalue " << eigenvalue << ", cfl " << cfl << "\n";
 
     if ( point_bound < point_internal ) //! the edge is on the left of the domain
      { 
	ASSERT( -1. < cfl && cfl < 0. ,
		"This characteristics is wrong!\nEither it is not outcoming (eigenvalue>0 at the left of the domain),\n or CFL is too high.");

	weight = - cfl;
      } 
    else   //! the edge is on the right of the domain
      {
	ASSERT( 0. < cfl && cfl < 1. ,
		"This characteristics is wrong!\nEither it is not outcoming (eigenvalue<0 at the right of the domain),\n or CFL is too high.");

	weight = cfl;
      }

    Vec2D u_interp( ( 1 - weight ) * U_bound.first  + weight * U_intern.first ,
		    ( 1 - weight ) * U_bound.second + weight * U_intern.second );
    return u_interp;

  }
	
  Real dot(const Vec2D& vec1, const Vec2D& vec2)
  {
    return vec1.first * vec2.first + vec1.second * vec2.second;
  }


}
