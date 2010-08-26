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
OneDimensionalModel_BC::OneDimensionalModel_BC( const OneD_BCSide& side ) :
    M_bcType                    (),
    M_bcSide                    ( side ),
    M_bcFunction                (),
    M_isInternal                ( false ),
    M_matrixrow_at_line         (),
    M_resBC                     ()
{
    M_matrixrow_at_line[ OneD_first ]  = Container2D_Type();
    M_matrixrow_at_line[ OneD_second ] = Container2D_Type();
}

OneDimensionalModel_BC::OneDimensionalModel_BC( const OneDimensionalModel_BC& BC ) :
    M_bcType                    ( BC.M_bcType ),
    M_bcSide                    ( BC.M_bcSide ),
    M_bcFunction                ( BC.M_bcFunction ),
    M_isInternal                ( BC.M_isInternal ),
    M_matrixrow_at_line         ( BC.M_matrixrow_at_line ),
    M_resBC                     ( BC.M_resBC )
{}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_BC::applyBC( const Real&             time,
                                 const Real&             timeStep,
                                 const Solution_Type&    solution,
                                 const Flux_PtrType&     flux,
                                       Container2D_Type& BC_dir )
{
    ASSERT_PRE( BC_dir.size() == 2, "applyBC works only for 2D vectors");

    if( M_isInternal )
        Debug(6311) << "[OneDimensionalModel_BC::compute_resBC] found internal boundary\n";
    else
    {
        // Update the time step (Compatibility condition need it)
        compute_resBC( time, timeStep, solution, flux );

        for( UInt i(0) ; i < 2 ; ++i )
            BC_dir[i]=M_resBC[i];
    }

    Debug(6311) << "[OneDimensionalModel_BC::applyBC] on side " << M_bcSide
                << " imposing [ A, Q ] = [ " << BC_dir[0] << ", " << BC_dir[1] << " ]\n";
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_BC::setType( const OneD_BCLine& line, const OneD_BC& bc )
{
    M_bcType[line] = bc;
}

void
OneDimensionalModel_BC::setBCFunction( const OneD_BCLine& line, const BCFunction_Type& BCFunction )
{
    M_bcFunction[line] = BCFunction; //FactoryClone_OneDimensionalModel_BCFunction::instance().createObject( &BCFunction );
}

void
OneDimensionalModel_BC::setInternalFlag( const bool& flag )
{
    M_isInternal = flag;
}

void
OneDimensionalModel_BC::setMatrixRow( const OneD_BCLine& line, const Container2D_Type& matrixrow )
{
    M_matrixrow_at_line[line] = matrixrow;
}

// ===================================================
// Get Methods
// ===================================================
const OneD_BC&
OneDimensionalModel_BC::type( const OneD_BCLine& line )
{
    return M_bcType[line];
}

OneDimensionalModel_BC::BCFunction_Type&
OneDimensionalModel_BC::BCFunction( const OneD_BCLine& line )
{
    return M_bcFunction[line]; //FactoryClone_OneDimensionalModel_BCFunction::instance().createObject( &rhs );
}

const bool&
OneDimensionalModel_BC::isInternal()
{
    return M_isInternal;
}

// ===================================================
// Private Methods
// ===================================================
void
OneDimensionalModel_BC::compute_resBC( const Real&             time,
                                       const Real&             timeStep,
                                       const Solution_Type&    solution,
                                       const Flux_PtrType&     flux )
{
    Container2D_Type rhsBC;

    // Eigen values of the jacobian diffFlux (= dF/dU = H)
    Real  eigval1, eigval2;

    // Left eigen vectors for the eigen values eigval1 and eigval2
    Container2D_Type left_eigvec1, left_eigvec2;
    Container2D_Type left_eigvec_first, left_eigvec_second;

    left_eigvec1[0] = 0.;
    left_eigvec1[1] = 0.;

    left_eigvec2[0] = 0.;
    left_eigvec2[1] = 0.;

    UInt dof;
    ( M_bcSide == OneD_left ) ? dof = 1 : dof = flux->Physics()->Data()->NumberOfElements() + 1;

    Container2D_Type U_boundary, W_boundary;

    U_boundary[0] = (*solution.find("A")->second)(dof);
    U_boundary[1] = (*solution.find("Q")->second)(dof);
    W_boundary[0] = (*solution.find("W1")->second)(dof);
    W_boundary[1] = (*solution.find("W2")->second)(dof);

    flux->jacobian_EigenValues_Vectors( U_boundary[0], U_boundary[1],
                                        eigval1, eigval2,
                                        left_eigvec1[0], left_eigvec1[1],
                                        left_eigvec2[0], left_eigvec2[1],
                                        dof );

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] 1\n";

    // Compute RHS
    rhsBC[0] = M_bcFunction[OneD_first](time, timeStep);
    rhsBC[1] = M_bcFunction[OneD_second](time, timeStep);

    // The flow rate  and the pressure are given positive with respect to
    // the normal direction: (left) <=  --------  => (right)
    if ( M_bcSide == OneD_left )
    {
        if ( M_bcType[OneD_first] == OneD_Q || M_bcType[OneD_first] == OneD_P )
            rhsBC[0] *= -1;
        if ( M_bcType[OneD_second] == OneD_Q || M_bcType[OneD_second] == OneD_P )
            rhsBC[1] *= -1;
    }

    // If is a pressure BC, convert it into an area
    if ( M_bcType[OneD_first] == OneD_P )
        rhsBC[0] = flux->Physics()->A_from_P( rhsBC[0], dof - 1 ); // Index start from 0

    if ( M_bcType[OneD_second] == OneD_P )
        rhsBC[1] = flux->Physics()->A_from_P( rhsBC[1], dof - 1 ); // Index start from 0

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] rhsBC[0] = " << rhsBC[0] << "\n";
    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] rhsBC[1] = " << rhsBC[1] << "\n";

    // this is not general
    // PRECONDITION (typical situation):
    //   on first line,  left boundary, I impose W1
    //   on second line, left boundary, I impose W2
    //   on first line,  right boundary, I impose W2
    //   on second line, right boundary, I impose W1
    // the code does not check for coherence (you cannot impose
    // the same variable on both lines!)

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] 2\n";
    M_bcType[OneD_first] == OneD_W1 ? //"W1"
        left_eigvec_first = left_eigvec1 :
        left_eigvec_first = left_eigvec2;

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] 3\n";
    M_bcType[OneD_second] == OneD_W1 ? //"W1"
        left_eigvec_second = left_eigvec1 :
        left_eigvec_second = left_eigvec2;

    compute_resBC_line(OneD_first,  left_eigvec_first, U_boundary, W_boundary, rhsBC[0]);
    compute_resBC_line(OneD_second, left_eigvec_second, U_boundary, W_boundary, rhsBC[1]);

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC] solving linear system with "
                << "\n\tfirst line = " << M_matrixrow_at_line[OneD_first][0]
                << ", " << M_matrixrow_at_line[OneD_first][1]
                << "\n\tsecond line = " << M_matrixrow_at_line[OneD_second][0]
                << ", " << M_matrixrow_at_line[OneD_second][1]
                << "\n\trhs = " << rhsBC[0]
                << ", " << rhsBC[1]
                << "\n";

    M_resBC=_solveLinearSyst2x2(M_matrixrow_at_line[OneD_first], M_matrixrow_at_line[OneD_second], rhsBC);
}

void
OneDimensionalModel_BC::compute_resBC_line( OneD_BCLine line, Container2D_Type left_eigvec, Container2D_Type U, Container2D_Type W, Real& rhs )
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

    switch( M_bcType[line] )
    {
        case OneD_W1: //"W1"
            //       ASSERT(eigval1<0. && eigval2<0.,
            // "The eigenvalues do no have the expected signs (lam1<0 and lam2<0).");
            M_matrixrow_at_line[line] = left_eigvec;
            LnUn = dot( left_eigvec, U );
            add = LnUn - W[0];
            rhs += add;

        break;
        case OneD_W2: //"W2"
            M_matrixrow_at_line[line] = left_eigvec;
            LnUn = dot( left_eigvec, U );
            add = LnUn - W[1];
            rhs += add;
        break;
        case OneD_A: //"A"
        case OneD_P: //"P"
            M_matrixrow_at_line[line][0] = 1.; M_matrixrow_at_line[line][1] = 0.;
        break;
        case OneD_Q: //"Q"
            M_matrixrow_at_line[line][0] = 0.; M_matrixrow_at_line[line][1] = 1.;
        break;
        default: std::cout << "\n[OneDimensionalModel_BC::compute_resBC] Wrong boundary variable as " << line
                           << " condition on side " << M_bcSide;
    }

    Debug(6311) << "[OneDimensionalModel_BC::compute_resBC_line] to impose variable "
                << M_bcType[line]
                << ", " << line << " line = "
                << M_matrixrow_at_line[line][0] << ", "
                << M_matrixrow_at_line[line][1] << "\n";
}

Container2D_Type
OneDimensionalModel_BC::_solveLinearSyst2x2( const Container2D_Type& line1,
                                             const Container2D_Type& line2,
                                             const Container2D_Type& rhs2d ) const
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

    Container2D_Type res;
    res[0] = ( aa22 * rhs2d[0] - aa12 * rhs2d[1] ) / determinant;
    res[1] = ( - aa21 * rhs2d[0] + aa11 * rhs2d[1] ) / determinant;
    return res;
}

}
