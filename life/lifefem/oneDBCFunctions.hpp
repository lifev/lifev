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

#ifndef _ONED_FUNCTIONS_HPP_
#define _ONED_FUNCTIONS_HPP_

#include <life/lifearray/vecUnknown.hpp>
#include <life/lifemesh/basicOneDMesh.hpp>
#include <life/lifesolver/vectorFunction1D.hpp>
#include <life/lifefem/dofOneD.hpp>
#include <boost/shared_ptr.hpp>

namespace LifeV
{
  /** @name Typedefs
   */
  //@{
  //! 2D vector       
  typedef std::pair< Real, Real > Vec2D; 
  //! Vector containing the nodal values of the scalar unknown
  typedef ScalUnknown<Vector> ScalVec;   
  //! 2D dot product
  Real dot(const Vec2D& vec1, const Vec2D& vec2);
  // Value-Defintions of the different String values
  enum OneDBCStringValue { OneDBCLeftBoundary, OneDBCRightBoundary,
  							OneDBCW1, OneDBCW2, OneDBCA, OneDBCQ,
  							OneDBCFirstRHS, OneDBCSecondRHS };



  //! Forward declaration
  class OneDBCFunctionBase;
  //! Pointer to BC function
	typedef boost::shared_ptr<OneDBCFunctionBase> OneDBCFunctionPointer;
  //@}

 
	


	class OneDBCFunctionBase
  {
  public:
    //		OneDBCFunctionBase( const OneDModelSolver& solver, const OneDNonLinModelParam& param ){};
		
    virtual Real evaluate( const Real& /*time*/ ){ return 0.;}
    
    virtual ~OneDBCFunctionBase() {}
    
  };



  class Const : public OneDBCFunctionBase
  {
  public:
    Const(const Real val) : _M_val(val) {}
    Real evaluate( const Real& /*time*/ ) {return _M_val;}
    ~Const() {}
  private:
    Real _M_val;
  };



  class Reimann : public OneDBCFunctionBase
  {
  public:

    Reimann( const BasicOneDMesh& mesh,
	     const NonLinearFluxFun1D& fluxFun,
	     const ScalVec& U1_thistime, const ScalVec& U2_thistime,  
	     const ScalVec& W1_thistime, const ScalVec& W2_thistime,  
	     const std::string & border, const std::string & var );
	
    virtual Real evaluate( const Real& time );
    
    ~Reimann() {}
	
  protected:
	void update_U_boundary();
    // 	virtual Real function( const Real& time );
    //! Reference to the mesh
    const BasicOneDMesh& _M_mesh;
    //! Reference to the solver non linear flux functions
    const NonLinearFluxFun1D& _M_fluxFun;
    //! Reference to the solver current unknowns (A) 
    const ScalVec& _M_U1_thistime;
    //! Reference to the solver current unknowns (Q)
    const ScalVec& _M_U2_thistime;
    //! Reference to the solver current unknowns (W1) 
    const ScalVec& _M_W1_thistime;
    //! Reference to the solver current unknowns (W2)
    const ScalVec& _M_W2_thistime;
    //! Mesh vertices
    DofOneD _M_dof1D;
    //! Number of Dof
    UInt _M_dimDof;
    //! Boundary Dof (right or left)
    UInt _M_boundaryDof;
    //! Value of U at the boundary
    Vec2D _M_U_boundary;
    //! Value of W at the boundary
    Vec2D _M_W_boundary;
//! boundary
    std::string _M_border;	
 //! variable
    std::string _M_var;
    
    std::map<std::string, OneDBCStringValue> _M_oneDBCFunctionsMapStringValues;
  };



  class Compatibility : public Reimann
  {
  public:
    Compatibility( const BasicOneDMesh& mesh,
		   const NonLinearFluxFun1D& fluxFun, const NonLinearSourceFun1D& sourceFun,
		   const ScalVec& U1_thistime, const ScalVec& U2_thistime,  
		   const ScalVec& W1_thistime, const ScalVec& W2_thistime,  
		   const Real& dt, const std::string & border, const std::string & var);
	
    virtual Real evaluate( const Real& time );
    
    ~Compatibility() {}
	
  protected:
	void update_U_internalBd();

	void computeEigenValuesVectors();

	Real extrapolate_L_dot_U( Real const& eigval, Vec2D const& eigvec );

	Real extrapolate_W( OneDBCStringValue const& _W );
	
    Vec2D _interpolLinear(const Real& point_bound, const Real& point_internal,
			  const Real& deltaT, const Real& eigenvalue,
			  const Vec2D& U_bound, const Vec2D& U_intern) const;
    //! Reference to the solver non linear source functions
    const NonLinearSourceFun1D& _M_sourceFun;
    //! Number of elements
    UInt _M_nb_elem;  
    //! Time step
    const Real& _M_time_step;   
    //! Dof of the internal node adjacent to the boundary
    UInt _M_internalBoundaryDof;
    //! Boundary Edge
    Edge1D _M_boundaryEdge;    
    //! Boundary point and internal boundary point
    Real _M_boundaryPoint, _M_internalBdPoint;
    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    Real _M_eigval1, _M_eigval2;
    //! Left eigen vectors for the eigen values eigval1 and eigval2
    Vec2D _M_left_eigvec1, _M_left_eigvec2;
    //! Value of U at the neighboring internal node
    Vec2D _M_U_internalBd;
  };

}

#endif
