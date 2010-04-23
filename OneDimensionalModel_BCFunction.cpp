//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for the boundary conditions for 1D tubes.
 *
 *  @version 1.0
 *  @author Lucia Mirabella
 *  @date 01-08-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 *
 */

#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
Riemann::Riemann( const FESpace_Type&                   fespace,
                  const Flux_PtrType                    fluxFun,
                  const std::vector<Vector_Type>&       U_thistime,
                  const std::string&                    border,
                  const std::string&                    var ):
    OneDimensionalModel_BCFunction            ( ),
    M_FESpace                     ( fespace ),
    M_fluxFun                     ( fluxFun ),
    M_U_thistime                  ( U_thistime ),
    M_dimDof                      ( fespace.dof().numTotalDof() ),
    M_U_boundary                  ( 2 ),
    M_W_boundary                  ( 2 ),
    M_border                      ( border ),
    M_var                         ( var )
{
    Debug( 6315 ) << "[OneDBCFUnction] Riemann \n";

    M_oneDBCFunctionsMapStringValues["left" ] = OneDBCLeftBoundary;
    M_oneDBCFunctionsMapStringValues["right"] = OneDBCRightBoundary;

    switch( M_oneDBCFunctionsMapStringValues[M_border] )
    {
        case OneDBCLeftBoundary:
            M_boundaryDof = 1;
        break;

        case OneDBCRightBoundary:
            M_boundaryDof = M_dimDof;
        break;

        default:
            std::cout << "\n[Riemann::Riemann] incorrect boundary identifier: " << M_border << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
Real
Riemann::evaluate( const Real& /*time*/ )
{
    update_U_boundary();

    return ( ( M_var == "W1" ) ? M_W_boundary[0] : M_W_boundary[1] );
}

// ===================================================
// Protected Methods
// ===================================================
void
Riemann::update_U_boundary()
{
    M_U_boundary[0] = M_U_thistime[0](M_boundaryDof);
    M_U_boundary[1] = M_U_thistime[1](M_boundaryDof);
    M_W_boundary[0] = M_U_thistime[2](M_boundaryDof);
    M_W_boundary[1] = M_U_thistime[3](M_boundaryDof);
}



// ===================================================
// Constructors & Destructor
// ===================================================
Compatibility::Compatibility( const FESpace_Type& fespace,
                              const Flux_PtrType                    fluxFun,
                              const Source_PtrType                  sourceFun,
                              const std::vector<Vector_Type>&       U_thistime,
                              const Real&                           dt,
                              const std::string&                    border,
                              const std::string&                    var ):
    Riemann        ( fespace, fluxFun, U_thistime, /*W_thistime,*/ border, var),
    M_sourceFun    ( sourceFun ),
    M_nb_elem      ( this->M_dimDof - 1 ),
    M_time_step    ( dt ),
    M_left_eigvec1 ( 2 ),
    M_left_eigvec2 ( 2 ),
    M_U_internalBd ( 2 )
{
    this->M_oneDBCFunctionsMapStringValues["W1"] = OneDBCW1;
    this->M_oneDBCFunctionsMapStringValues["W2"] = OneDBCW2;

    switch( this->M_oneDBCFunctionsMapStringValues[this->M_border] )
    {
        case OneDBCLeftBoundary:
            M_internalBoundaryDof = this->M_boundaryDof + 1;
            M_boundaryEdge        = this->M_FESpace.mesh()->edgeList(1);
            M_boundaryPoint       = M_boundaryEdge.point(1).x();
            M_internalBdPoint     = M_boundaryEdge.point(2).x();
        break;

        case OneDBCRightBoundary:
            M_internalBoundaryDof = this->M_boundaryDof - 1;
            M_boundaryEdge        = this->M_FESpace.mesh()->edgeList(M_nb_elem);
            M_boundaryPoint       = M_boundaryEdge.point(2).x();
            M_internalBdPoint     = M_boundaryEdge.point(1).x();
        break;

        default:
            std::cout << "\n[Compatibility::Compatibility] incorrect boundary identifier: " << this->M_border << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
Real
Compatibility::evaluate( Real const& /*time*/ )
{
    Debug( 6315 ) << "[Compatibility::evaluate] variable "
                  << this->M_var << ", code "
                  << this->M_oneDBCFunctionsMapStringValues[this->M_var] << "\n";

    return extrapolate_W( this->M_oneDBCFunctionsMapStringValues[this->M_var] );
}

// ===================================================
// Protected Methods
// ===================================================
void
Compatibility::computeEigenValuesVectors()
{
    this->M_fluxFun->jacobian_EigenValues_Vectors( this->M_U_thistime[0](this->M_boundaryDof),
                                                   this->M_U_thistime[1](this->M_boundaryDof),
                                                   M_eigval1, M_eigval2,
                                                   M_left_eigvec1[0], M_left_eigvec1[1],
                                                   M_left_eigvec2[0], M_left_eigvec2[1],
                                                   this->M_boundaryDof );
}

void
Compatibility::update_U_internalBd()
{
    M_U_internalBd[0]=this->M_U_thistime[0](M_internalBoundaryDof);
    M_U_internalBd[1]=this->M_U_thistime[1](M_internalBoundaryDof);
}

Real
Compatibility::extrapolate_L_dot_U( Real const& eigval,
                                                  Vec2D const& eigvec )
{
    ASSERT_PRE( eigvec.size() == 2, "extrapolate_L_dot_U work only for 2D vectors");

    Real L_dot_U_extrap;
    Vec2D qlSource(2); // Quasi linear source term
    Vec2D U_charact_pt=_interpolLinear( M_boundaryPoint, M_internalBdPoint,
                                        M_time_step, eigval, this->M_U_boundary, M_U_internalBd );

    L_dot_U_extrap = dot( eigvec, U_charact_pt);

    Debug( 6315 ) << "[extrapolate_L_dot_U] eigvec.size() = " << eigvec.size()
                  << ", U_charact_pt.size() = " << U_charact_pt.size() << "\n";

    qlSource[0]=M_sourceFun->QuasiLinearSource( U_charact_pt[0],
                                               U_charact_pt[1], 1,
                                               this->M_boundaryDof - 1);

    qlSource[1]=M_sourceFun->QuasiLinearSource( U_charact_pt[0],
                                               U_charact_pt[1], 2,
                                               this->M_boundaryDof - 1);

    L_dot_U_extrap-=M_time_step * dot(eigvec, qlSource);

    return L_dot_U_extrap;
}

Real
Compatibility::extrapolate_W( OneDBCStringValue const& _W )
{
    Real W_out(0.);

    this->update_U_boundary();
    this->update_U_internalBd();

    computeEigenValuesVectors();

    switch( _W )
    {
        case OneDBCW1:
            W_out = extrapolate_L_dot_U(M_eigval1, M_left_eigvec1)
            - dot( M_left_eigvec1, this->M_U_boundary ) + this->M_W_boundary[0];

//             std::cout << extrapolate_L_dot_U(M_eigval1, M_left_eigvec1) << " "
//                       << M_left_eigvec1[0] << " " << M_left_eigvec1[1] << " "
//                       << this->M_U_boundary[0]  << " " << this->M_U_boundary[1]  << " "
//                       << this->M_W_boundary[1] << std::endl;
        break;

        case OneDBCW2:
            W_out = extrapolate_L_dot_U(M_eigval2, M_left_eigvec2)
            - dot( M_left_eigvec2, this->M_U_boundary ) + this->M_W_boundary[1];

//             std::cout << extrapolate_L_dot_U(M_eigval2, M_left_eigvec2) << " "
//                       << M_left_eigvec2[0] << " " << M_left_eigvec2[1] << " "
//                       << this->M_U_boundary[0]  << " " << this->M_U_boundary[1]  << " "
//                       << this->M_W_boundary[1] << std::endl;
        break;

        default:
        std::cout << "\n[Compatibility::extrapolate_W] incorrect variable identifier: " << _W << std::endl;
    }
    return W_out;
}

Vec2D
Compatibility::_interpolLinear( const Real& point_bound, const Real& point_internal,
                                              const Real& deltaT,      const Real& eigenvalue,
                                              const Vec2D& U_bound,    const Vec2D& U_intern) const
{
    ASSERT_PRE( U_bound.size() == 2 && U_intern.size() == 2, "_interpolLinear work only for 2D vectors");

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

    Vec2D u_interp(2);
    for( UInt i=0; i<2; ++i )
        u_interp[i] = ( 1 - weight ) * U_bound[i]  + weight * U_intern[i];

    return u_interp;
}

}
