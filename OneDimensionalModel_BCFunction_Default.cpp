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
 *  @brief File containing some functions for the boundary conditions of 1D models.
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

#include <lifemc/lifefem/OneDimensionalModel_BCFunction_Default.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Default::OneDimensionalModel_BCFunction_Default( const Flux_PtrType     flux,
                                                                                const Source_PtrType   source,
                                                                                const OneD_BCSide&     side,
                                                                                const OneD_BC&         bcType ):
    M_Flux                          ( flux ),
    M_Source                        ( source ),
    M_boundaryDof                   (),
    M_bcType                        ( bcType )
{
    switch( side )
    {
        case OneD_left:
            M_boundaryDof = 1;
        break;

        case OneD_right:
            M_boundaryDof = M_Flux->Physics()->Data()->NumberOfElements() + 1;
        break;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Default::setSolution( const Solution_PtrType solution )
{
    M_Solution = solution;
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Riemann::OneDimensionalModel_BCFunction_Riemann( const Flux_PtrType     flux,
                                                                                const Source_PtrType   source,
                                                                                const OneD_BCSide&     side,
                                                                                const OneD_BC&         bcType ):
    super                           ( flux, source, side, bcType ),
    M_U_boundary                    (),
    M_W_boundary                    ()
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Riemann::operator()( const Real& /*time*/, const Real& /*timeStep*/ )
{
    update_U_boundary();

    return ( ( M_bcType == OneD_W1 ) ? M_W_boundary[0] : M_W_boundary[1] );
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Riemann::update_U_boundary()
{
    M_U_boundary[0] = (*(*this->M_Solution)["A"])(M_boundaryDof);
    M_U_boundary[1] = (*(*this->M_Solution)["Q"])(M_boundaryDof);
    M_W_boundary[0] = (*(*this->M_Solution)["W1"])(M_boundaryDof);
    M_W_boundary[1] = (*(*this->M_Solution)["W2"])(M_boundaryDof);
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Compatibility::OneDimensionalModel_BCFunction_Compatibility( const Flux_PtrType     flux,
                                                                                            const Source_PtrType   source,
                                                                                            const OneD_BCSide&     side,
                                                                                            const OneD_BC&         bcType ):
    super                           ( flux, source, side, bcType ),
    M_internalBoundaryDof           (),
    M_boundaryPoint                 (),
    M_internalBdPoint               (),
    M_eigval1                       (),
    M_eigval2                       (),
    M_left_eigvec1                  (),
    M_left_eigvec2                  ()
{
    Mesh_Type::EdgeType boundaryEdge;
    switch( side )
    {
        case OneD_left:
            M_internalBoundaryDof = this->M_boundaryDof + 1;
            boundaryEdge          = this->M_Flux->Physics()->Data()->mesh()->edgeList(1);
            M_boundaryPoint[0]    = boundaryEdge.point(1).x();
            M_boundaryPoint[1]    = boundaryEdge.point(1).y();
            M_boundaryPoint[2]    = boundaryEdge.point(1).z();
            M_internalBdPoint[0]  = boundaryEdge.point(2).x();
            M_internalBdPoint[1]  = boundaryEdge.point(2).y();
            M_internalBdPoint[2]  = boundaryEdge.point(2).z();
        break;

        case OneD_right:
            M_internalBoundaryDof = this->M_boundaryDof - 1;
            boundaryEdge          = this->M_Flux->Physics()->Data()->mesh()->edgeList(this->M_boundaryDof - 1);
            M_boundaryPoint[0]    = boundaryEdge.point(2).x();
            M_boundaryPoint[1]    = boundaryEdge.point(2).y();
            M_boundaryPoint[2]    = boundaryEdge.point(2).z();
            M_internalBdPoint[0]  = boundaryEdge.point(1).x();
            M_internalBdPoint[1]  = boundaryEdge.point(1).y();
            M_internalBdPoint[2]  = boundaryEdge.point(1).z();
        break;

        default:
            std::cout << "\n[Compatibility::Compatibility] incorrect boundary identifier: " << side << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Compatibility::operator()( const Real& /*time*/, const Real& timeStep )
{
    Debug( 6315 ) << "[Compatibility::evaluate] variable " << this->M_bcType << "\n";
    return extrapolate_W( this->M_bcType, timeStep );
}

// ===================================================
// Protected Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Compatibility::extrapolate_W( const OneD_BC& W, const Real& timeStep )
{
    Real W_out(0.);

    this->update_U_boundary();
    computeEigenValuesVectors();

    switch( W )
    {
        case OneD_W1:
            W_out = extrapolate_L_dot_U( M_eigval1, M_left_eigvec1, timeStep )
                  - dot( M_left_eigvec1, this->M_U_boundary ) + this->M_W_boundary[0];

//             std::cout << extrapolate_L_dot_U(M_eigval1, M_left_eigvec1) << " "
//                       << M_left_eigvec1[0] << " " << M_left_eigvec1[1] << " "
//                       << this->M_U_boundary[0]  << " " << this->M_U_boundary[1]  << " "
//                       << this->M_W_boundary[1] << std::endl;
        break;

        case OneD_W2:
            W_out = extrapolate_L_dot_U( M_eigval2, M_left_eigvec2, timeStep )
                  - dot( M_left_eigvec2, this->M_U_boundary ) + this->M_W_boundary[1];

//             std::cout << extrapolate_L_dot_U(M_eigval2, M_left_eigvec2) << " "
//                       << M_left_eigvec2[0] << " " << M_left_eigvec2[1] << " "
//                       << this->M_U_boundary[0]  << " " << this->M_U_boundary[1]  << " "
//                       << this->M_W_boundary[1] << std::endl;
        break;

        default:
        std::cout << "\n[Compatibility::extrapolate_W] incorrect variable identifier: " << W << std::endl;
    }
    return W_out;
}

void
OneDimensionalModel_BCFunction_Compatibility::computeEigenValuesVectors()
{
    this->M_Flux->jacobian_EigenValues_Vectors( (*(*this->M_Solution)["A"])(this->M_boundaryDof),
                                                (*(*this->M_Solution)["Q"])(this->M_boundaryDof),
                                                M_eigval1, M_eigval2,
                                                M_left_eigvec1[0], M_left_eigvec1[1],
                                                M_left_eigvec2[0], M_left_eigvec2[1],
                                                this->M_boundaryDof );
}

Real
OneDimensionalModel_BCFunction_Compatibility::extrapolate_L_dot_U( const Real& eigval, const Container2D_Type& eigvec, const Real& timeStep )
{
    ASSERT_PRE( eigvec.size() == 2, "extrapolate_L_dot_U work only for 2D vectors");

    Real L_dot_U_extrap;
    Container2D_Type qlSource; // Quasi linear source term
    Container2D_Type U_charact_pt=_interpolLinear( eigval, this->M_U_boundary, timeStep );

    L_dot_U_extrap = dot( eigvec, U_charact_pt );

    Debug( 6315 ) << "[extrapolate_L_dot_U] eigvec.size() = " << eigvec.size()
                  << ", U_charact_pt.size() = " << U_charact_pt.size() << "\n";

    qlSource[0]=this->M_Source->QuasiLinearSource( U_charact_pt[0],
                                                   U_charact_pt[1], 1,
                                                   this->M_boundaryDof - 1);

    qlSource[1]=this->M_Source->QuasiLinearSource( U_charact_pt[0],
                                                   U_charact_pt[1], 2,
                                                   this->M_boundaryDof - 1);

    L_dot_U_extrap -= timeStep * dot(eigvec, qlSource);

    return L_dot_U_extrap;
}

Container2D_Type
OneDimensionalModel_BCFunction_Compatibility::_interpolLinear( const Real& eigenvalue, const Container2D_Type& U_bound, const Real& timeStep ) const
{
    ASSERT_PRE( U_bound.size() == 2, "_interpolLinear work only for 2D vectors");

    Real deltaX = std::sqrt( std::pow(M_boundaryPoint[0] - M_internalBdPoint[0], 2) +
                             std::pow(M_boundaryPoint[1] - M_internalBdPoint[1], 2) +
                             std::pow(M_boundaryPoint[2] - M_internalBdPoint[2], 2) );

    Real cfl =  eigenvalue * timeStep / deltaX;

    Real weight;   //! weight in the linear approximation

    Debug( 6315 ) << "[Compatibility::_interpolLinear] point_bound ("
                  << M_boundaryPoint[0] << "," << M_boundaryPoint[1] << "," << M_boundaryPoint[2] << "), point_internal ("
                  << M_internalBdPoint[0] << "," << M_internalBdPoint[1] << "," << M_internalBdPoint[2]
                  << "), deltaT " << timeStep << ", deltaX " << deltaX
                  << ", eigenvalue " << eigenvalue
                  << ", A boundary " << U_bound[0] << ", Q boundary " << U_bound[1]
                  << ", A internal " << (*(*this->M_Solution)["A"])(M_internalBoundaryDof) << ", Q internal " << (*(*this->M_Solution)["Q"])(M_internalBoundaryDof)
                  << ", cfl " << cfl << "\n";

    if ( M_internalBoundaryDof == 2 ) //! the edge is on the left of the domain
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

    Container2D_Type u_interp;

    u_interp[0] = ( 1 - weight ) * U_bound[0]  + weight * (*(*this->M_Solution)["A"])(M_internalBoundaryDof);
    u_interp[1] = ( 1 - weight ) * U_bound[1]  + weight * (*(*this->M_Solution)["Q"])(M_internalBoundaryDof);

    return u_interp;
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Absorbing::OneDimensionalModel_BCFunction_Absorbing( const Flux_PtrType     flux,
                                                                                    const Source_PtrType   source,
                                                                                    const OneD_BCSide&     side,
                                                                                    const OneD_BC&         bcType ):
    super                           ( flux, source, side, bcType )
{}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_BCFunction_Absorbing::operator()( Real const& /*time*/, const Real& timeStep )
{
    Real W_out(0.);
    Real a1, a2, a11, a22, b1, b2, c1, c2;

    this->update_U_boundary();

    Debug( 6030 ) << "[Absorbing::Absorbing] at node " << this->M_boundaryDof
                  << ", A  = " << this->M_U_boundary[0] << "( " << (*(*this->M_Solution)["A"])[1] << " ) "
                  << ", Q  = " << this->M_U_boundary[1]
                  << ", W1 = " << this->M_W_boundary[0]
                  << ", W2 = " << this->M_W_boundary[1]
                  << "\n";

    this->computeEigenValuesVectors();

    a1 = M_Flux->Physics()->elasticPressure(this->M_U_boundary[0], this->M_boundaryDof - 1); // elastic pressure at previous time step

    a2 = this->M_U_boundary[1]; // flux at previous time step

    b1 = M_Flux->Physics()->pressure_WDiff( this->M_W_boundary[0], this->M_W_boundary[1], 1, this->M_boundaryDof - 1);  // dP / dW1

    b2 = this->M_U_boundary[0] / 2; // dQ / dW1

    c1 = M_Flux->Physics()->pressure_WDiff( this->M_W_boundary[0], this->M_W_boundary[1], 2, this->M_boundaryDof - 1);  // dP / dW2

    c2 = b2; // dQ / dW2

    Debug( 6030 ) << "[Absorbing::operator()] P(A) = " << a1 << "\n";
    Debug( 6030 ) << "[Absorbing::operator()] P(W1,W2) = "
                  << M_Flux->Physics()->pressure_W(this->M_W_boundary[0], this->M_W_boundary[1], this->M_boundaryDof - 1) << "\n";

    a11 = a1 - b1*this->M_W_boundary[0] - c1*this->M_W_boundary[1];
    a22 = a2 - b2*this->M_W_boundary[0] - c2*this->M_W_boundary[1];

    switch( this->M_bcType )
    {
    case OneD_W1:
        W_out = this->extrapolate_L_dot_U(this->M_eigval2, this->M_left_eigvec2, timeStep )
              - dot( this->M_left_eigvec2, this->M_U_boundary ) + this->M_W_boundary[1];

        break;
    case OneD_W2:
        W_out = this->extrapolate_L_dot_U(this->M_eigval1, this->M_left_eigvec1, timeStep )
              - dot( this->M_left_eigvec1, this->M_U_boundary ) + this->M_W_boundary[0];

        break;
    default:
        std::cout << "\n[Absorbing::operator()] incorrect variable identifier: " << this->M_bcType << std::endl;
    }

    Debug( 6030 ) << "[Absorbing::operator()] extrapolated exiting characteristic = " << W_out << "\n";

    Real resistance(b1 / b2);

    this->resistance( resistance );
    Debug( 6030 ) << "[Absorbing::operator()] imposing condition, R = " << resistance << "\n";

    Debug( 6030 ) << "[Absorbing::operator()] a1 = " << a1 << "\n";
    Debug( 6030 ) << "[Absorbing::operator()] b1 = " << b1 << "\n";
    Debug( 6030 ) << "[Absorbing::operator()] c1 = " << c1 << "\n";
    Debug( 6030 ) << "[Absorbing::operator()] a2 = " << a2 << "\n";
    Debug( 6030 ) << "[Absorbing::operator()] b2 = " << b2 << "\n";
    Debug( 6030 ) << "[Absorbing::operator()] c2 = " << c2 << "\n";

    return W_out * ((b2*resistance-b1)/(c1-c2*resistance)) + ((a22*resistance-a11)/(c1-c2*resistance));
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Absorbing::resistance( Real& /*resistance*/ )
{
    //Do nothing = absorbing!
}



// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction_Resistance::OneDimensionalModel_BCFunction_Resistance( const Flux_PtrType     flux,
                                                                                      const Source_PtrType   source,
                                                                                      const OneD_BCSide&     side,
                                                                                      const OneD_BC&         bcType,
                                                                                      const Real&            resistance ):
    super                           ( flux, source, side, bcType ),
    M_resistance                    ( resistance )
{}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalModel_BCFunction_Resistance::resistance( Real& resistance )
{
    resistance = M_resistance;
}

}
