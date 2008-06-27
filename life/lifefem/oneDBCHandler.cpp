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
  \file oneDBC.hpp
  \author Lucia Mirabella
  \date 08/2006

  \brief This file contains a class for the boundary conditions for 1D tubes.

*/

#include <life/lifefem/oneDBCHandler.hpp>

namespace LifeV
{

  OneDBC::OneDBC( const ScalVec& U1_thistime, const ScalVec& U2_thistime, 
                  const ScalVec& W1_thistime, const ScalVec& W2_thistime, 
                  const NonLinearFluxFun1D& fluxFun,
                  const std::string& side ):
    _M_U1_thistime(U1_thistime),
    _M_U2_thistime(U2_thistime),
    _M_W1_thistime(W1_thistime),
    _M_W2_thistime(W2_thistime),
    _M_fluxFun(fluxFun)
  {
	std::string prefix = "boundcond/" + side;
	Debug( 6030 ) << "[OneDBC::OneDBC] data file prefix = " << prefix << "\n";

  	(side == "left") ? _M_boundaryDof=0 : _M_boundaryDof = U1_thistime.size() - 1; //data_file("discretization/nb_elem",10);
    
//    std::string prefix_internal = prefix + "/internal";
    _M_isInternal=false; //data_file(prefix_internal.c_str(),0);
        	
//    if ( _M_isInternal ){
      	_M_variable_at_RHS["first"]="not set";
    	_M_variable_at_RHS["second"]="not set";
    	//}
//    else
//	{
//		std::string prefix_var = prefix + "/first";
//		_M_variable_at_RHS["first"] = data_file(prefix_var.c_str(),"W1");
//		prefix_var = prefix + "/second";
//    	_M_variable_at_RHS["second"] = data_file(prefix_var.c_str(),"W2");
//	}
	
    _M_oneDBCMapStringValues["W1"] = OneDBCW1;
    _M_oneDBCMapStringValues["W2"] = OneDBCW2;
    _M_oneDBCMapStringValues["A"] = OneDBCA;
    _M_oneDBCMapStringValues["Q"] = OneDBCQ;
    //    	std::cout << "\n[OneDBC::OneDBC], " <<_M_time_step << "\t" << Pa0 << std::endl;

  }
   


  Vec2D OneDBC::Uboundary(const ScalVec& U1, const ScalVec& U2) const
  {
    Vec2D Ubound(U1(_M_boundaryDof),U2(_M_boundaryDof));
    return Ubound;
  }
  


  void OneDBC::compute_resBC(const Real& time_val) //, OneDBCHandler::OneDBCFunctionPointer f1, OneDBCHandler::OneDBCFunctionPointer f2 )
  {
    Vec2D rhsBC;
    
    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    Real  eigval1, eigval2;
    //! Left eigen vectors for the eigen values eigval1 and eigval2
    Vec2D left_eigvec1, left_eigvec2;
    
	Vec2D U_boundary, W_boundary;
	
    U_boundary.first = _M_U1_thistime(_M_boundaryDof);
    U_boundary.second = _M_U2_thistime(_M_boundaryDof);
    W_boundary.first = _M_W1_thistime(_M_boundaryDof);
    W_boundary.second = _M_W2_thistime(_M_boundaryDof);

    _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
					    eigval1, eigval2,
					    left_eigvec1.first, left_eigvec1.second,
					    left_eigvec2.first, left_eigvec2.second,
					    _M_boundaryDof);  

    Real LnUn, add;

    rhsBC.first=_M_f1->evaluate(time_val);
    rhsBC.second=_M_f2->evaluate(time_val);

	Debug(6030) << "[OneDBC::compute_resBC] at node " << _M_boundaryDof
			<< " f1(time_val) = " << rhsBC.first
			<< ", f2(time_val) = " << rhsBC.second << "\n";  

	switch( _M_oneDBCMapStringValues[ _M_variable_at_RHS["first"] ] ) 
      {
      case OneDBCW1: //"W1"
	//       ASSERT(eigval1<0. && eigval2<0., "The eigenvalues do no have the expected signs (lam1<0 and lam2<0).");
	_M_line1 = left_eigvec1;
	LnUn = dot( left_eigvec1, U_boundary );
	add = LnUn - W_boundary.first;
	rhsBC.first += add;

	break;
      case OneDBCW2: //"W2"
	_M_line1 = left_eigvec2;  
	LnUn = dot( left_eigvec2, U_boundary );
	add = LnUn - W_boundary.second;
	rhsBC.first += add;
	break;
      case OneDBCA: //"A"
      	_M_line1.first = 1.; _M_line1.second = 0.;  
	break;
      case OneDBCQ: //"Q"
      	_M_line1.first = 0.; _M_line1.second = 1.;  
    	break;	
      default: std::cout << "\n[OneDBC::compute_resBC] Wrong boundary variable as first condition at node " << _M_boundaryDof;
      }
    
    switch( _M_oneDBCMapStringValues[ _M_variable_at_RHS["second"] ] ) 
      {
      case OneDBCW1: //"W1"
	//       ASSERT(eigval1<0. && eigval2<0., "The eigenvalues do no have the expected signs (lam1<0 and lam2<0).");
	_M_line2 = left_eigvec1;  
	LnUn = dot( left_eigvec1, U_boundary );
	add = LnUn - W_boundary.first;
	rhsBC.second += add;
	break;
      case OneDBCW2: //"W2"
	_M_line2 = left_eigvec2;  
	LnUn = dot( left_eigvec2, U_boundary );
	add = LnUn - W_boundary.second;
	rhsBC.second += add;
	break;
      case OneDBCA: //"A"
      	_M_line2.first = 1.;  _M_line2.second = 0.;  
	break;
      case OneDBCQ: //"Q"
      	_M_line2.first = 0.;  _M_line2.second = 1.;  
    	break;	
      default: std::cout << "\n[OneDBC::compute_resBC] Wrong boundary variable as second condition at node " << _M_boundaryDof;
      }

    Debug(6030) << "[OneDBC::compute_resBC] add term for second line = " << add << "\n";

    _M_resBC=_solveLinearSyst2x2(_M_line1, _M_line2, rhsBC); 
   
  }
   


  void OneDBC::applyBC(const Real& time_val,
//		       OneDBCHandler::OneDBCFunctionPointer f1,
//		       OneDBCHandler::OneDBCFunctionPointer f2,
		       Vec2D& BC_dir)
  { 
	if( _M_isInternal )
	{
    	Debug(6030) << "[OneDBC::compute_resBC] found internal boundary\n";
	}
	else
	{
    compute_resBC(time_val);
    BC_dir.first=_M_resBC.first;
    BC_dir.second=_M_resBC.second;
	}
	
	Debug(6030) << "[OneDBC::applyBC] at node " << _M_boundaryDof
			<< " imposing [ A, Q ] = [ " << BC_dir.first
			<< ", " << BC_dir.second << " ]\n";  
  }  


  /*! Matrix A is given by two pairs corresponding to the 2 _M_lines.
    A = [_M_line1;
    _M_line2 ]
    return A^{-1} * rhs2d
  */
  Vec2D OneDBC::_solveLinearSyst2x2(const Vec2D& line1, const Vec2D& line2,
				    const Vec2D& rhs2d) const
  {
    Real aa11, aa12, aa21, aa22;

    aa11 = line1.first;  aa12 = line1.second;
    aa21 = line2.first;  aa22 = line2.second;

    Real determinant = aa11 * aa22 - aa12 * aa21;

    ASSERT( determinant != 0,
            "Error: the 2x2 system on the boundary is not invertible. \nCheck the boundary conditions.");
    Vec2D res( (   aa22 * rhs2d.first - aa12 * rhs2d.second ) / determinant,
               ( - aa21 * rhs2d.first + aa11 * rhs2d.second ) / determinant );

    return res;
  }


	
  OneDBCHandler::OneDBCHandler( const ScalVec& U1_thistime, const ScalVec& U2_thistime, 
                                const ScalVec& W1_thistime, const ScalVec& W2_thistime, 
                                const NonLinearFluxFun1D& fluxFun ) :
    _M_U1_thistime(U1_thistime),
    _M_U2_thistime(U2_thistime),
    _M_W1_thistime(W1_thistime),
    _M_W2_thistime(W2_thistime),
    _M_fluxFun(fluxFun),
    leftBoundary( U1_thistime, U2_thistime, W1_thistime, W2_thistime, fluxFun, "left" ),
    rightBoundary( U1_thistime, U2_thistime, W1_thistime, W2_thistime, fluxFun, "right" )
  {  	
    leftBCmap.insert( make_pair("first", false));
    leftBCmap.insert( make_pair("second", false));
    rightBCmap.insert( make_pair("first", false));
    rightBCmap.insert( make_pair("second", false));

    _M_oneDBCHandlerMapStringValues["left"] = OneDBCLeftBoundary;
    _M_oneDBCHandlerMapStringValues["right"] = OneDBCRightBoundary;
    _M_oneDBCHandlerMapStringValues["first"] = OneDBCFirstRHS;
    _M_oneDBCHandlerMapStringValues["second"] = OneDBCSecondRHS;
  }
  


  void OneDBCHandler::setDefaultBC( const BasicOneDMesh& mesh,
  				 const NonLinearSourceFun1D& sourceFun,
			     const Real& dt )
  {
    if(!leftBCmap["first"])
      {
	OneDBCFunctionPointer point ( new Reimann( mesh, _M_fluxFun,
							_M_U1_thistime, _M_U2_thistime,  
							_M_W1_thistime, _M_W2_thistime,  
						  "left", "W1" ) );
	leftBCmap["first"] = true;
	leftBoundary.variable("first")="W1";
	leftBoundary.f1() = point;
      Debug( 6030 ) << "[OneDBCHandler::setDefaultBC] imposing constant W1 at left boundary (first line).\n";
      }
    if(!leftBCmap["second"])
      {		OneDBCFunctionPointer point ( new Compatibility( mesh,
								_M_fluxFun, sourceFun, _M_U1_thistime, _M_U2_thistime,  
							_M_W1_thistime, _M_W2_thistime,  
								dt, "left", "W2" ) );
	leftBCmap["second"] = true;
	leftBoundary.variable("second")="W2";
      leftBoundary.f2() = point;
      Debug( 6030 ) << "[OneDBCHandler::setDefaultBC] imposing compatibility condition for W2 at left boundary (second line).\n";
      }
    if(!rightBCmap["first"])
      {		
      	OneDBCFunctionPointer point ( new Reimann( mesh, _M_fluxFun, _M_U1_thistime, _M_U2_thistime,  
							_M_W1_thistime, _M_W2_thistime,  
						    "right", "W2" ) );
	rightBCmap["first"] = true;
	rightBoundary.variable("first")="W2";
      rightBoundary.f1() = point;
      Debug( 6030 ) << "[OneDBCHandler::setDefaultBC] imposing constant W2 at right boundary (first line).\n";
      }
    if(!rightBCmap["second"])
	{
	  OneDBCFunctionPointer point ( new Compatibility( mesh,
								_M_fluxFun, sourceFun, _M_U1_thistime, _M_U2_thistime,  
							_M_W1_thistime, _M_W2_thistime,  
								dt, "right", "W1" ) );
	rightBCmap["second"] = true;
	rightBoundary.variable("second")="W1";
	  rightBoundary.f2() = point;
      Debug( 6030 ) << "[OneDBCHandler::setDefaultBC] imposing compatibility condition for W1 at right boundary (second line).\n";
	}

  }

void OneDBCHandler::applyBC(const Real& time_val,
		 Vec2D& left_BC_dir,
		 Vec2D& right_BC_dir)
{
	leftBoundary.applyBC( time_val, left_BC_dir );
	rightBoundary.applyBC( time_val, right_BC_dir );
};
		 
		 
		 
  void OneDBCHandler::setBC( OneDBCFunctionPointer funptr,
			     std::string const& border,
			     std::string const& line,
			     std::string const& var )
  {
    switch( _M_oneDBCHandlerMapStringValues[border] ){
    case OneDBCLeftBoundary:
	leftBCmap[line] = true;
	leftBoundary.variable(line)=var;
      switch( _M_oneDBCHandlerMapStringValues[line] ){
      case OneDBCFirstRHS:
	leftBoundary.f1() = funptr;
      Debug( 6030 ) << "[OneDBCHandler::setBC] imposing user defined function at left boundary (first line).\n";
	break;
      case OneDBCSecondRHS:
	leftBoundary.f2() = funptr; 
      Debug( 6030 ) << "[OneDBCHandler::setBC] imposing user defined function at left boundary (second line).\n";
	break;
      default:
      std::cout << "\n[OneDBCHandler::setBC] error while tring to impose "
      			<< line << " boundary condition at right boundary" << std::endl;
      }
      break;	
    case OneDBCRightBoundary:
	rightBCmap[line] = true;
	rightBoundary.variable(line)=var;
      switch( _M_oneDBCHandlerMapStringValues[line] ){
      case OneDBCFirstRHS:
	rightBoundary.f1() = funptr;
      Debug( 6030 ) << "[OneDBCHandler::setBC] imposing user defined function at right boundary (first line).\n";
	break;
      case OneDBCSecondRHS:
	rightBoundary.f2() = funptr;
      Debug( 6030 ) << "[OneDBCHandler::setBC] imposing user defined function at right boundary (second line).\n";
	break;
      default:
      std::cout << "\n[OneDBCHandler::setBC] error while tring to impose "
      			<< line << " boundary condition at right boundary" << std::endl;
      }
      break;	
    default:
      std::cout << "\n[OneDBCHandler::setBC] incorrect boundary identifier: "
      			<< border << std::endl;
    }
  }
  
  
  
  OneDBCFunctionPointer&
  OneDBCHandler::leftBCFunction( const std::string& line )
  {
  	switch( _M_oneDBCHandlerMapStringValues[line] )
  	{
  		case OneDBCFirstRHS:
  		return leftBoundary.f1();
  		break;
  		case OneDBCSecondRHS:
  		return leftBoundary.f2();
  		break;
  		default:
  		std::cout << "\n[OneDBCHandler::leftBC] wrong line identifier: " << line << std::endl;
  		abort();
  	}
  }
		


	OneDBCFunctionPointer&
  OneDBCHandler::rightBCFunction( const std::string& line )
  {
  	switch( _M_oneDBCHandlerMapStringValues[line] )
  	{
  		case OneDBCFirstRHS:
  		return rightBoundary.f1();
  		break;
  		case OneDBCSecondRHS:
  		return rightBoundary.f2();
  		break;
  		default:
  		std::cout << "\n[OneDBCHandler::rightBC] wrong line identifier: " << line << std::endl;
  		abort();
  	}
  }
	
}
