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

#ifndef _ONED_BCHANDLER_HPP_
#define _ONED_BCHANDLER_HPP_

#include <life/lifesolver/oneDModelHandler.hpp>
#include <life/lifesolver/oneDNonLinModelParam.hpp>
#include <life/lifesolver/vectorFunction1D.hpp>
#include <life/lifefem/oneDBCFunctions.hpp>
#include <boost/shared_ptr.hpp>

namespace LifeV
{
	
  /*
    class OneDBCHandler {

    public:
    void addFunction( std::string, const OneDBCFunctionBase& );
	
    Real evaluate( std::string, const Real& );
	
    private:
    std::map<std::string, const OneDBCFunctionBase&> BCList;
    };
  */



  /*!
   *\class OneDBC
   * Abstract class which contains members and methods common to OneDBCLeft and OneDBCRight classes.
   */
  class OneDBC 
  {
  public:
    /** @name Constructor
     */
    //@{
    //! Constructor taking data_file, onedparam, mesh, fluxFun, U1_thistime, U_2thistime
    OneDBC( const ScalVec& U1_thistime, const ScalVec& U2_thistime, 
            const ScalVec& W1_thistime, const ScalVec& W2_thistime, 
            const NonLinearFluxFun1D& fluxFun,
            const std::string& side ); 
    //@}
    /** @name Destructor
     */
    //@{  
    ~OneDBC(){}
    //@}
    /** @name Public Methods
     */
    //@{
    //! Compute [A,Q] at the boundary
    Vec2D Uboundary(const ScalVec& U1, const ScalVec& U2) const;
    //! Return the boundary Dof
    inline UInt giveBoundaryDof() const{return _M_boundaryDof;}
    //! Apply boundary conditions
    void applyBC(const Real& time_val,
//		 OneDBCHandler::OneDBCFunctionPointer f1,
//		 OneDBCHandler::OneDBCFunctionPointer f2,
		 Vec2D& BC_dir);
	
	OneDBCFunctionPointer& f1(){ return _M_f1; }

	OneDBCFunctionPointer& f2(){ return _M_f2; }
		 
	    inline bool& isInternal() { return _M_isInternal; }

    inline std::string& variable( const std::string& line ) 
    {
		return _M_variable_at_RHS[line];
    }
    //@}   
  protected:
    /** @name Protected Methods
     */
    //@{
    //! Impose the chosen boundary condition
    void compute_resBC( const Real& time_val ); //, OneDBCHandler::OneDBCFunctionPointer f1, OneDBCHandler::OneDBCFunctionPointer f2 );    
  
    //! solve a 2x2 linear system by the Cramer method (for the boundary systems)
    Vec2D _solveLinearSyst2x2(const Vec2D& line1, const Vec2D& line2,
                              const Vec2D& rhs2d) const;
    //@}

    OneDBCFunctionPointer _M_f1, _M_f2;
    
    bool _M_isInternal;
	
	std::map<std::string, std::string> _M_variable_at_RHS;
	
	std::map<std::string, OneDBCStringValue> _M_oneDBCMapStringValues;

    //! Right hand side for the 2x2 linear system to be solved at each side
//    Vec2D _M_rhsBC;
    //! Result of the 2x2 linear system to be solved at each side
    Vec2D _M_resBC;
    
//    Vec2D _M_U_boundary;
    
    //! Reference to the solver current unknowns (A) 
    const ScalVec& _M_U1_thistime;
    //! Reference to the solver current unknowns (Q)
    const ScalVec& _M_U2_thistime;
    //! Reference to the solver current unknowns (W1) 
    const ScalVec& _M_W1_thistime;
    //! Reference to the solver current unknowns (W2)
    const ScalVec& _M_W2_thistime;
    //! Value of U at the boundary, at the neighboring internal node
    Vec2D _M_line1, _M_line2;
    //! Boundary Dof (right or left)
    UInt _M_boundaryDof;
    //! Reference to the solver non linear flux functions
    const NonLinearFluxFun1D& _M_fluxFun;
  };



  class OneDBCHandler
  {
  public:
		
    OneDBCHandler( const ScalVec& U1_thistime, const ScalVec& U2_thistime, 
  				const ScalVec& W1_thistime, const ScalVec& W2_thistime, 
		  		const NonLinearFluxFun1D& fluxFun );
		
    void setBC( OneDBCFunctionPointer funptr, std::string const& border, std::string const& line, std::string const& var );
    
    void setDefaultBC(const BasicOneDMesh& mesh,
	       const NonLinearSourceFun1D& sourceFun,
	       const Real& dt);
    //		OneDBCFunctionBase giveFunc(std::string const& border, std::string const& var );

    //! Apply boundary conditions
    void applyBC(const Real& time_val,
		 Vec2D& left_BC_dir,
		 Vec2D& right_BC_dir);
	
    OneDBCFunctionPointer& leftBCFunction( const std::string & line );	
    OneDBCFunctionPointer& rightBCFunction( const std::string & line );
    
    inline bool& leftBCReady( const std::string & line ){ return leftBCmap[line]; }	
    inline bool& rightBCReady( const std::string & line ){ return leftBCmap[line]; }	
    	
	inline void setBCLeft_internalnode(){ leftBoundary.isInternal()=true; }
	inline void setBCRight_internalnode(){ rightBoundary.isInternal()=true; }

  private:
//    std::vector< OneDBCFunctionPointer > leftBCvec, rightBCvec;
    std::map< std::string, bool > leftBCmap, rightBCmap;		

	std::map<std::string, OneDBCStringValue> _M_oneDBCHandlerMapStringValues;
    //! Reference to the solver current unknowns (A) 
    const ScalVec& _M_U1_thistime;
    //! Reference to the solver current unknowns (Q)
    const ScalVec& _M_U2_thistime;
    //! Reference to the solver current unknowns (W1) 
    const ScalVec& _M_W1_thistime;
    //! Reference to the solver current unknowns (W2)
    const ScalVec& _M_W2_thistime;
    //! Reference to the solver non linear flux functions
    const NonLinearFluxFun1D& _M_fluxFun;

    OneDBC leftBoundary, rightBoundary;
  };

}
#endif
