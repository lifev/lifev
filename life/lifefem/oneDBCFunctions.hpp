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
#include <life/lifesolver/oneDModelHandler.hpp>
#include <boost/shared_ptr.hpp>

namespace ublas = boost::numeric::ublas;

namespace LifeV
{
/** @name Typedefs
 */
//@{
// Value-Defintions of the different String values
enum OneDBCStringValue { OneDBCLeftBoundary, OneDBCRightBoundary,
                         OneDBCW1, OneDBCW2, OneDBCA, OneDBCQ, OneDBCFUN,
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


template< class FLUX >
class Reimann : public OneDBCFunctionBase
{
public:

    Reimann( const BasicOneDMesh& mesh,
             const FLUX& fluxFun,
             const OneDModelHandler::ScalVec_vector& U_thistime,
             const std::string & border, const std::string & var );

    virtual Real evaluate( const Real& time );

    ~Reimann() {}

protected:
    void update_U_boundary();
    // 	virtual Real function( const Real& time );
    //! Reference to the mesh
    const BasicOneDMesh& _M_mesh;
    //! Reference to the solver non linear flux functions
    const FLUX& _M_fluxFun;
    //! Reference to the solver current unknowns (U)
    const OneDModelHandler::ScalVec_vector& _M_U_thistime;
    //! Mesh vertices
    DofOneD _M_dof1D;
    //! Number of Dof
    UInt _M_dimDof;
    //! Boundary Dof (right or left)
    UInt _M_boundaryDof;
    //! Value of U at the boundary
    OneDModelHandler::Vec2D _M_U_boundary;
    //! Value of W at the boundary
    OneDModelHandler::Vec2D _M_W_boundary;
    //! boundary
    std::string _M_border;
    //! variable
    std::string _M_var;

    std::map<std::string, OneDBCStringValue> _M_oneDBCFunctionsMapStringValues;
};


template< class FLUX >
Reimann<FLUX>::Reimann( const BasicOneDMesh& mesh,
                        const FLUX& fluxFun,
                        const OneDModelHandler::ScalVec_vector& U_thistime,
                        const std::string & border, const std::string & var )
    :
    OneDBCFunctionBase( ),
    _M_mesh(mesh),
    _M_fluxFun(fluxFun),
    _M_U_thistime(U_thistime),
    _M_dof1D(_M_mesh.numVertices()),
    _M_dimDof(_M_dof1D.numTotalDof()),
    _M_U_boundary(2),
    _M_W_boundary(2),
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

template< class FLUX >
void
Reimann<FLUX>::update_U_boundary()
{
    _M_U_boundary[0]=_M_U_thistime[0](_M_boundaryDof);
    _M_U_boundary[1]=_M_U_thistime[1](_M_boundaryDof);
    _M_W_boundary[0]=_M_U_thistime[2](_M_boundaryDof);
    _M_W_boundary[1]=_M_U_thistime[3](_M_boundaryDof);
}


template< class FLUX >
Real
Reimann<FLUX>::evaluate( const Real& /*time*/ )
{
    update_U_boundary();

    return ( ( _M_var == "W1" ) ? _M_W_boundary[0] : _M_W_boundary[1] );

}



template< class FLUX, class SOURCE >
class Compatibility : public Reimann<FLUX>
{
public:
    Compatibility( const BasicOneDMesh& mesh,
                   const FLUX& fluxFun, const SOURCE& sourceFun,
                   const OneDModelHandler::ScalVec_vector& U_thistime,
                   const Real& dt, const std::string & border, const std::string & var);

    virtual Real evaluate( const Real& time );

    ~Compatibility() {}

protected:
    void update_U_internalBd();

    void computeEigenValuesVectors();

    Real extrapolate_L_dot_U( Real const& eigval, OneDModelHandler::Vec2D const& eigvec );

    Real extrapolate_W( OneDBCStringValue const& _W );

    OneDModelHandler::Vec2D _interpolLinear(const Real& point_bound, const Real& point_internal,
                                            const Real& deltaT, const Real& eigenvalue,
                                            const OneDModelHandler::Vec2D& U_bound, const OneDModelHandler::Vec2D& U_intern) const;
    //! Reference to the solver non linear source functions
    const SOURCE& _M_sourceFun;
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
    OneDModelHandler::Vec2D _M_left_eigvec1, _M_left_eigvec2;
    //! Value of U at the neighboring internal node
    OneDModelHandler::Vec2D _M_U_internalBd;
};

template< class FLUX, class SOURCE >
Compatibility<FLUX, SOURCE>::Compatibility( const BasicOneDMesh& mesh,
                                            const FLUX& fluxFun, const SOURCE& sourceFun,
                                            const OneDModelHandler::ScalVec_vector& U_thistime,
                                            const Real& dt, const std::string & border, const std::string & var )
    :
    Reimann<FLUX>( mesh, fluxFun, U_thistime, /*W_thistime,*/ border, var),
    _M_sourceFun(sourceFun),
    _M_nb_elem(this->_M_dimDof-1),
    _M_time_step (dt), //data_file("time/timestep",0.1)),
    _M_left_eigvec1(2),
    _M_left_eigvec2(2),
    _M_U_internalBd(2)
{
    this->_M_oneDBCFunctionsMapStringValues["W1"] = OneDBCW1;
    this->_M_oneDBCFunctionsMapStringValues["W2"] = OneDBCW2;

    switch( this->_M_oneDBCFunctionsMapStringValues[this->_M_border] )
        {
        case OneDBCLeftBoundary:
            _M_internalBoundaryDof=this->_M_boundaryDof + 1;
            _M_boundaryEdge=this->_M_mesh.edgeList(1);
            _M_boundaryPoint=_M_boundaryEdge.pt1().x();
            _M_internalBdPoint=_M_boundaryEdge.pt2().x();
            break;
        case OneDBCRightBoundary:
            _M_internalBoundaryDof=this->_M_boundaryDof - 1;
            _M_boundaryEdge=this->_M_mesh.edgeList(_M_nb_elem);
            _M_boundaryPoint=_M_boundaryEdge.pt2().x();
            _M_internalBdPoint=_M_boundaryEdge.pt1().x();
            break;
        default:
            std::cout << "\n[Compatibility::Compatibility] incorrect boundary identifier: " << this->_M_border << std::endl;
        }
}



template< class FLUX, class SOURCE >
void
Compatibility<FLUX, SOURCE>::computeEigenValuesVectors()
{
    this->_M_fluxFun.jacobian_EigenValues_Vectors( this->_M_U_thistime[0](this->_M_boundaryDof),
                                                   this->_M_U_thistime[1](this->_M_boundaryDof),
                                                   _M_eigval1, _M_eigval2,
                                                   _M_left_eigvec1[0], _M_left_eigvec1[1],
                                                   _M_left_eigvec2[0], _M_left_eigvec2[1],
                                                   this->_M_boundaryDof );
}


template< class FLUX, class SOURCE >
void
Compatibility<FLUX, SOURCE>::update_U_internalBd()
{
    _M_U_internalBd[0]=this->_M_U_thistime[0](_M_internalBoundaryDof);
    _M_U_internalBd[1]=this->_M_U_thistime[1](_M_internalBoundaryDof);
}


template< class FLUX, class SOURCE >
Real
Compatibility<FLUX, SOURCE>::extrapolate_L_dot_U( Real const& eigval,
                                                  OneDModelHandler::Vec2D const& eigvec )
{
    ASSERT_PRE( eigvec.size() == 2, "extrapolate_L_dot_U work only for 2D vectors");

    Real L_dot_U_extrap;
    //! Quasi linear source term
    OneDModelHandler::Vec2D qlSource(2);

    OneDModelHandler::Vec2D U_charact_pt=_interpolLinear( _M_boundaryPoint, _M_internalBdPoint,
                                                          _M_time_step, eigval, this->_M_U_boundary,
                                                          _M_U_internalBd);

    L_dot_U_extrap=dot( eigvec, U_charact_pt);

    Debug( 6315 ) << "[extrapolate_L_dot_U] eigvec.size() = " << eigvec.size()
                  << ", U_charact_pt.size() = " << U_charact_pt.size() << "\n";

    qlSource[0]=_M_sourceFun.QuasiLinearSource(U_charact_pt[0],
                                               U_charact_pt[1], 1,
                                               this->_M_boundaryDof);

    qlSource[1]=_M_sourceFun.QuasiLinearSource(U_charact_pt[0],
                                               U_charact_pt[1], 2,
                                               this->_M_boundaryDof);

    L_dot_U_extrap-=_M_time_step * dot(eigvec, qlSource);

    return L_dot_U_extrap;

}


template< class FLUX, class SOURCE >
Real
Compatibility<FLUX, SOURCE>::extrapolate_W( OneDBCStringValue const& _W )
{
    Real W_out(0.);

    this->update_U_boundary();
    this->update_U_internalBd();

    computeEigenValuesVectors();

    switch( _W )
        {
        case OneDBCW1:
            W_out = extrapolate_L_dot_U(_M_eigval1, _M_left_eigvec1)
                - dot( _M_left_eigvec1, this->_M_U_boundary ) + this->_M_W_boundary[0];

            break;
        case OneDBCW2:
            W_out = extrapolate_L_dot_U(_M_eigval2, _M_left_eigvec2)
                - dot( _M_left_eigvec2, this->_M_U_boundary ) + this->_M_W_boundary[1];

            break;
        default:
            std::cout << "\n[Compatibility::extrapolate_W] incorrect variable identifier: " << _W << std::endl;
        }

    return W_out;
}


template< class FLUX, class SOURCE >
Real
Compatibility<FLUX, SOURCE>::evaluate( Real const& /*time*/ )
{
    Debug( 6315 ) << "[Compatibility::evaluate] variable "
                  << this->_M_var << ", code "
                  << this->_M_oneDBCFunctionsMapStringValues[this->_M_var] << "\n";

    return extrapolate_W( this->_M_oneDBCFunctionsMapStringValues[this->_M_var] );
}


template< class FLUX, class SOURCE >
OneDModelHandler::Vec2D
Compatibility<FLUX, SOURCE>::_interpolLinear(const Real& point_bound, const Real& point_internal,
                                             const Real& deltaT, const Real& eigenvalue,
                                             const OneDModelHandler::Vec2D& U_bound, const OneDModelHandler::Vec2D& U_intern) const
{
    ASSERT_PRE( U_bound.size() == 2 && U_intern.size() == 2,
                "_interpolLinear work only for 2D vectors");

    Real deltaX = std::abs(point_bound - point_internal);

    Real cfl =  eigenvalue * deltaT / deltaX;

    Real weight;   //! weight in the linear approximation

    Debug( 6315 ) << "[Compatibility::_interpolLinear] point_bound "
                  << point_bound << ", point_internal " << point_internal
                  << ", deltaT " << deltaT << ", deltaX " << deltaX
                  << ", eigenvalue " << eigenvalue
                  << ", A boundary " << U_bound[0] << ", Q boundary " << U_bound[1]
                  << ", A internal " << U_intern[0] << ", Q internal " << U_intern[1]
                  << ", cfl " << cfl << "\n";

    if ( point_bound < point_internal ) //! the edge is on the left of the domain
        {
            ASSERT( -1. < cfl && cfl < 0. ,
                    "This characteristics is wrong!\nEither it is not outcoming " \
                    "(eigenvalue>0 at the left of the domain),\n or CFL is too high.");

            weight = - cfl;
        }
    else   //! the edge is on the right of the domain
        {
            ASSERT( 0. < cfl && cfl < 1. ,
                    "This characteristics is wrong!\nEither it is not outcoming " \
                    "(eigenvalue<0 at the right of the domain),\n or CFL is too high.");

            weight = cfl;
        }

    OneDModelHandler::Vec2D u_interp(2);
    for( UInt i=0; i<2; ++i )
        u_interp[i] = ( 1 - weight ) * U_bound[i]  + weight * U_intern[i];
    return u_interp;

}

}

#endif
