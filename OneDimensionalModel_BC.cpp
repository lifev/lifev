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
 *  @brief File containing a class for the boundary conditions of the 1D model.
 *
 *  @version 1.0
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *  @date 01-28-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 */


#include <lifemc/lifefem/OneDimensionalModel_BC.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BC::OneDimensionalModel_BC( const std::vector<Vector_Type>&  U_thistime,
                                                const Flux_PtrType               fluxFun,
                                                const Real&                      dimDof,
                                                const std::string&               side ) :
    M_isInternal                (),
    M_variable_at_line          (),
    M_matrixrow_at_line         (),
    M_rhs_at_line               (),
    M_OneDimensionalModel_BCMapStringValues     (),
    M_resBC                     ( 2 ),
    M_U_thistime                ( U_thistime ),
    M_boundaryDof               (),
    M_fluxFun                   ( fluxFun )
{
    (side == "left") ? M_boundaryDof = 1 : M_boundaryDof = dimDof;

    M_isInternal=false;

    M_variable_at_line["first"]="not set";
    M_variable_at_line["second"]="not set";

    M_matrixrow_at_line["first"]=Vec2D(2);
    M_matrixrow_at_line["second"]=Vec2D(2);

    M_OneDimensionalModel_BCMapStringValues["W1"]  = OneDBCW1;
    M_OneDimensionalModel_BCMapStringValues["W2"]  = OneDBCW2;
    M_OneDimensionalModel_BCMapStringValues["A"]   = OneDBCA;
    M_OneDimensionalModel_BCMapStringValues["Q"]   = OneDBCQ;
    M_OneDimensionalModel_BCMapStringValues["fun"] = OneDBCFUN;
}

// ===================================================
// Methods
// ===================================================
Vec2D
OneDimensionalModel_BC::Uboundary( const ScalVec& U1, const ScalVec& U2 ) const
{
    Vec2D Ubound(2);
    Ubound[0] = U1(M_boundaryDof); Ubound[1] = U2(M_boundaryDof);
    return Ubound;
}

void
OneDimensionalModel_BC::applyBC( const Real& time_val, Vec2D& BC_dir )
{
    //std::cout << "BC_dir = " << BC_dir[0] << " " << BC_dir[1] << std::endl;
    ASSERT_PRE( BC_dir.size() == 2,
                "applyBC works only for 2D vectors");

    if( M_isInternal )
        Debug(6311) << "[OneDimensionalModel_BC::compute_resBC] found internal boundary\n";
    else
    {
        compute_resBC(time_val);

        for( UInt i = 0; i < 2; ++i )
            BC_dir[i]=M_resBC[i];
    }

    Debug(6311) << "[OneDimensionalModel_BC::applyBC] at node " << M_boundaryDof
                << " imposing [ A, Q ] = [ " << BC_dir[0]
                << ", " << BC_dir[1] << " ]\n";
}

// ===================================================
// Get Methods
// ===================================================
UInt
OneDimensionalModel_BC::boundaryDof() const
{
    return M_boundaryDof;
}

OneDimensionalModel_BC::OneDimensionalModel_BCFunction_PtrType&
OneDimensionalModel_BC::rhs( const std::string& line )
{
    return M_rhs_at_line[line];
}

std::string&
OneDimensionalModel_BC::variable( const std::string& line )
{
    return M_variable_at_line[line];
}

Vec2D&
OneDimensionalModel_BC::matrixrow( const std::string& line )
{
    return M_matrixrow_at_line[line];
}

bool&
OneDimensionalModel_BC::isInternal()
{
    return M_isInternal;
}

// ===================================================
// Protected Methods
// ===================================================
void
OneDimensionalModel_BC::compute_resBC( const Real& time_val )
{
    Vec2D rhsBC(2);

    // Eigen values of the jacobian diffFlux (= dF/dU = H)
    Real  eigval1, eigval2;

    // Left eigen vectors for the eigen values eigval1 and eigval2
    Vec2D left_eigvec1(2), left_eigvec2(2);
    Vec2D left_eigvec_first(2), left_eigvec_second(2);

    left_eigvec1[0] = 0.;
    left_eigvec1[1] = 0.;

    left_eigvec2[0] = 0.;
    left_eigvec2[1] = 0.;

    Vec2D U_boundary(2), W_boundary(2);

    for( UInt i = 0; i < 2; ++i )
    {
        U_boundary[i] = M_U_thistime[i    ](M_boundaryDof);
        W_boundary[i] = M_U_thistime[2 + i](M_boundaryDof) ;
        //std::cout << "bdof " << M_boundaryDof << " : " << M_U_thistime[2 + i](M_boundaryDof) << std::endl;
    }

    Real Aboundary = U_boundary[0];
    Real Qboundary = U_boundary[1];

    M_fluxFun->jacobian_EigenValues_Vectors( Aboundary, Qboundary,
                                             eigval1, eigval2,
                                             left_eigvec1[0], left_eigvec1[1],
                                             left_eigvec2[0], left_eigvec2[1],
                                             M_boundaryDof );

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] 1\n";

    rhsBC[0] = M_rhs_at_line["first"]->evaluate(time_val);
    rhsBC[1] = M_rhs_at_line["second"]->evaluate(time_val);

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] rhsBC[0] = " << rhsBC[0] << "\n";;
    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] rhsBC[1] = " << rhsBC[1] << "\n";;

    // this is not general
    // PRECONDITION (typical situation):
    //   on first line,  left boundary, I impose W1
    //   on second line, left boundary, I impose W2
    //   on first line,  right boundary, I impose W2
    //   on second line, right boundary, I impose W1
    // the code does not check for coherence (you cannot impose
    // the same variable on both lines!)

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] 2\n";
    M_OneDimensionalModel_BCMapStringValues[ M_variable_at_line["first"] ] == OneDBCW1 ? //"W1"
        left_eigvec_first = left_eigvec1 :
        left_eigvec_first = left_eigvec2;

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] 3\n";
    M_OneDimensionalModel_BCMapStringValues[ M_variable_at_line["second"] ] == OneDBCW1 ? //"W1"
        left_eigvec_second = left_eigvec1 :
        left_eigvec_second = left_eigvec2;

    compute_resBC_line("first",  left_eigvec_first, U_boundary, W_boundary, rhsBC[0]);
    compute_resBC_line("second", left_eigvec_second, U_boundary, W_boundary, rhsBC[1]);

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC] solving linear system with "
                << "\n\tfirst line = " << M_matrixrow_at_line["first"][0]
                << ", " << M_matrixrow_at_line["first"][1]
                << "\n\tsecond line = " << M_matrixrow_at_line["second"][0]
                << ", " << M_matrixrow_at_line["second"][1]
                << "\n\trhs = " << rhsBC[0]
                << ", " << rhsBC[1]
                << "\n";

    M_resBC=_solveLinearSyst2x2(M_matrixrow_at_line["first"], M_matrixrow_at_line["second"], rhsBC);
}

void
OneDimensionalModel_BC::compute_resBC_line( std::string line, Vec2D left_eigvec, Vec2D U, Vec2D W, Real& rhs )
{
    ASSERT_PRE( left_eigvec.size() == 2 && U.size() == 2 && W.size() == 2,
                "compute_resBC_line works only for 2D vectors");

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] "
                << "M_matrixrow_at_line.size() = "
                << M_matrixrow_at_line.size() << ", "
                << "line = " << line << ", "
                << "M_matrixrow_at_line[line].size() = "
                << M_matrixrow_at_line[line].size() << ", "
                << "\n";

    Real LnUn, add;

    switch( M_OneDimensionalModel_BCMapStringValues[ M_variable_at_line[line] ] )
    {
        case OneDBCW1: //"W1"
            //       ASSERT(eigval1<0. && eigval2<0.,
            // "The eigenvalues do no have the expected signs (lam1<0 and lam2<0).");
            M_matrixrow_at_line[line] = left_eigvec;
            LnUn = dot( left_eigvec, U );
            add = LnUn - W[0];
            rhs += add;

        break;
        case OneDBCW2: //"W2"
            M_matrixrow_at_line[line] = left_eigvec;
            LnUn = dot( left_eigvec, U );
            add = LnUn - W[1];
            rhs += add;
        break;
        case OneDBCA: //"A"
            M_matrixrow_at_line[line][0] = 1.; M_matrixrow_at_line[line][1] = 0.;
        break;
        case OneDBCQ: //"Q"
            M_matrixrow_at_line[line][0] = 0.; M_matrixrow_at_line[line][1] = 1.;
        break;
        case OneDBCFUN: //linear combination of A, Q
        break;
        default: std::cout << "\n[OneDimensionalModel_BC::compute_resBC] Wrong boundary variable as " << line
                           << " condition at node " << M_boundaryDof;
    }

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] to impose variable "
                << M_variable_at_line[line]
                << ", " << line << " line = "
                << M_matrixrow_at_line[line][0] << ", "
                << M_matrixrow_at_line[line][1] << "\n";
}

Vec2D
OneDimensionalModel_BC::_solveLinearSyst2x2( const Vec2D& line1,
                                             const Vec2D& line2,
                                             const Vec2D& rhs2d ) const
{
    ASSERT_PRE( line1.size() == 2 && line2.size() == 2 && rhs2d.size() == 2,
                "_solveLinearSyst2x2 works only for 2D vectors");

    long double aa11, aa12, aa21, aa22;

    aa11 = line1[0];  aa12 = line1[1];
    aa21 = line2[0];  aa22 = line2[1];

    long double determinant = aa11 * aa22 - aa12 * aa21;

    Debug(6311) << "[OneDimensionalModel_BC::_solveLinearSyst2x2] solving linear system with "
                << "\n\tfirst line = " << aa11
                << ", " << aa12
                << "\n\tsecond line = " << aa21
                << ", " << aa22
                << "\n\trhs = " << rhs2d[0]
                << ", " << rhs2d[1]
                << "\n\tdet1 = " << aa11 * aa22
                << ", det2 = " << aa12 * aa21
                << "\n\tdet = " << determinant
                << "\n";

    ASSERT( determinant != 0,
            "Error: the 2x2 system on the boundary is not invertible."
            "\nCheck the boundary conditions.");

    Vec2D res(2);
    res[0] = ( aa22 * rhs2d[0] - aa12 * rhs2d[1] ) / determinant;
    res[1] = ( - aa21 * rhs2d[0] + aa11 * rhs2d[1] ) / determinant;
    return res;
}

}
