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
 *  @file
 *  @brief File containing a class for the boundary conditions of the 1D model.
 *
 *  @version 1.0
 *  @date 01-28-2006
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *
 *  @version 2.0
 *  @date 20-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifefem/OneDimensionalModel_BC.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BC::OneDimensionalModel_BC( const bcSide_Type& bcSide ) :
        M_bcType                    (),
        M_bcSide                    ( bcSide ),
        M_bcFunction                (),
        M_isInternal                ( false ),
        M_bcMatrix                  (),
        M_bcRHS                     ()
{
    M_bcMatrix[ OneDimensional::first ]  = container2D_Type();
    M_bcMatrix[ OneDimensional::second ] = container2D_Type();
}

OneDimensionalModel_BC::OneDimensionalModel_BC( const OneDimensionalModel_BC& bc ) :
        M_bcType                    ( bc.M_bcType ),
        M_bcSide                    ( bc.M_bcSide ),
        M_bcFunction                ( bc.M_bcFunction ),
        M_isInternal                ( bc.M_isInternal ),
        M_bcMatrix                  ( bc.M_bcMatrix ),
        M_bcRHS                     ( bc.M_bcRHS )
{}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_BC::applyBC( const Real& time, const Real& timeStep, const solution_Type& solution,
                                 const fluxPtr_Type& flux, container2D_Type& bc )
{

#ifdef HAVE_LIFEV_DEBUG
    ASSERT_PRE( bc.size() == 2, "applyBC works only for 2D vectors");

    if ( M_isInternal )
        Debug(6311) << "[OneDimensionalModel_BC::compute_resBC] found internal boundary\n";
    else
#endif
    {
        UInt dof;
        ( M_bcSide == OneDimensional::left ) ? dof = 0 : dof = flux->physics()->data()->numberOfNodes() - 1;

        container2D_Type boundaryU;
        boundaryU[0] = (*solution.find("A")->second)(dof + 1);
        boundaryU[1] = (*solution.find("Q")->second)(dof + 1);

        // Eigenvalues and eigenvectors of the jacobian diffFlux (= dF/dU = H)
        container2D_Type eigenvalues;
        container2D_Type leftEigenvector1, leftEigenvector2;

        flux->eigenValuesEigenVectors( boundaryU[0], boundaryU[1],
                                       eigenvalues, leftEigenvector1, leftEigenvector2, dof );

        computeMatrixAndRHS( time, timeStep, flux, OneDimensional::first,
                             leftEigenvector1, leftEigenvector2, dof, M_bcRHS[0]);

        computeMatrixAndRHS( time, timeStep, flux, OneDimensional::second,
                             leftEigenvector1, leftEigenvector2, dof, M_bcRHS[1]);

        bc = solveLinearSystem( M_bcMatrix[OneDimensional::first], M_bcMatrix[OneDimensional::second], M_bcRHS );
    }

#ifdef HAVE_LIFEV_DEBUG
    Debug(6311) << "[OneDimensionalModel_BC::applyBC] on bcSide " << M_bcSide << " imposing [ A, Q ] = [ " << bc[0] << ", " << bc[1] << " ]\n";
#endif
}

// ===================================================
// Private Methods
// ===================================================
void
OneDimensionalModel_BC::computeMatrixAndRHS( const Real& time, const Real& timeStep, const fluxPtr_Type& flux, const bcLine_Type& line,
                                             const container2D_Type& leftEigenvector1, const container2D_Type& leftEigenvector2,
                                             const UInt& dof, Real& rhs )
{
    // This is not general (typical situation):
    //     on first line,  left boundary,  I impose W1
    //     on second line, left boundary,  I impose W2
    //     on first line,  right boundary, I impose W2
    //     on second line, right boundary, I impose W1
    // The code does not check for coherence (you cannot impose the same variable on both lines!)

    // Compute Matrix & RHS
    rhs = M_bcFunction[ line ](time, timeStep);
    switch ( M_bcType[line] )
    {
    case OneDimensional::W1:
        M_bcMatrix[line] = leftEigenvector1;
        break;
    case OneDimensional::W2:
        M_bcMatrix[line] = leftEigenvector2;
        break;
    case OneDimensional::A:
        M_bcMatrix[line][0] = 1.;
        M_bcMatrix[line][1] = 0.;
        break;
    case OneDimensional::P:
        rhs = flux->physics()->fromPToA( rhs, dof );
        M_bcMatrix[line][0] = 1.;
        M_bcMatrix[line][1] = 0.;
        break;
    case OneDimensional::Q:
        // Flow rate is positive with respect to the outgoing normal
        if ( M_bcSide == OneDimensional::left )
            rhs *= -1;
        M_bcMatrix[line][0] = 0.;
        M_bcMatrix[line][1] = 1.;
        break;
    default:
        std::cout << "\n[OneDimensionalModel_BC::compute_resBC] Wrong boundary variable as " << line
                  << " condition on bcSide " << M_bcSide;
    }

#ifdef HAVE_LIFEV_DEBUG
    Debug(6311) << "[OneDimensionalModel_BC::compute_MatrixAndRHS] to impose variable "
    << M_bcType[line] << ", " << line << " line = "
    << M_bcMatrix[line][0] << ", "
    << M_bcMatrix[line][1] << "\n";
#endif
}

OneDimensionalModel_BC::container2D_Type
OneDimensionalModel_BC::solveLinearSystem( const container2D_Type& line1,
                                           const container2D_Type& line2,
                                           const container2D_Type& rhs ) const
{
#ifdef HAVE_LIFEV_DEBUG
    ASSERT_PRE( line1.size() == 2 && line2.size() == 2 && rhs.size() == 2,
                "_solveLinearSyst2x2 works only for 2D vectors");
#endif

    Real determinant = line1[0] * line2[1] - line1[1] * line2[0];

#ifdef HAVE_LIFEV_DEBUG
    ASSERT( determinant != 0,
            "Error: the 2x2 system on the boundary is not invertible."
            "\nCheck the boundary conditions.");
#endif

    container2D_Type solution;

    solution[0] = (   line2[1] * rhs[0] - line1[1] * rhs[1] ) / determinant;
    solution[1] = ( - line2[0] * rhs[0] + line1[0] * rhs[1] ) / determinant;

    return solution;
}

}
