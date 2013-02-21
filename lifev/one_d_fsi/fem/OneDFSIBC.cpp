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

#include <lifev/one_d_fsi/fem/OneDFSIBC.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSIBC::OneDFSIBC ( const bcSide_Type& bcSide ) :
    M_bcType                    (),
    M_bcSide                    ( bcSide ),
    M_bcFunction                ()
{
}

OneDFSIBC::OneDFSIBC ( const OneDFSIBC& bc ) :
    M_bcType                    ( bc.M_bcType ),
    M_bcSide                    ( bc.M_bcSide ),
    M_bcFunction                ( bc.M_bcFunction )
{}

// ===================================================
// Methods
// ===================================================
void
OneDFSIBC::applyBC ( const Real& time, const Real& timeStep, const solution_Type& solution,
                     const fluxPtr_Type& fluxPtr, vectorPtrContainer_Type& rhs )
{
    UInt iNode;
    ( M_bcSide == OneDFSI::left ) ? iNode = 0 : iNode = fluxPtr->physics()->data()->numberOfNodes() - 1;

    container2D_Type boundaryU;
    boundaryU[0] = (*solution.find ("A")->second) (iNode);
    boundaryU[1] = (*solution.find ("Q")->second) (iNode);

    // Eigenvalues and eigenvectors of the jacobian diffFlux (= dF/dU = H)
    container2D_Type eigenvalues;
    container2D_Type leftEigenvector1, leftEigenvector2;

    fluxPtr->eigenValuesEigenVectors ( boundaryU[0], boundaryU[1],
                                       eigenvalues, leftEigenvector1, leftEigenvector2, iNode );

    std::map<bcLine_Type, container2D_Type> bcMatrix;
    bcMatrix[ OneDFSI::first ]  = container2D_Type();
    bcMatrix[ OneDFSI::second ] = container2D_Type();

    container2D_Type bcRHS;

    // First line of Matrix and RHS
    computeMatrixAndRHS ( time, timeStep, fluxPtr, OneDFSI::first,
                          leftEigenvector1, leftEigenvector2, iNode, bcMatrix, bcRHS[0] );

    // Second line of Matrix and RHS
    computeMatrixAndRHS ( time, timeStep, fluxPtr, OneDFSI::second,
                          leftEigenvector1, leftEigenvector2, iNode, bcMatrix, bcRHS[1] );


    container2D_Type bc = solveLinearSystem ( bcMatrix[OneDFSI::first], bcMatrix[OneDFSI::second], bcRHS );

    // Set the BC in the RHS
    (*rhs[0]) ( iNode ) = bc[0];
    (*rhs[1]) ( iNode ) = bc[1];

#ifdef HAVE_LIFEV_DEBUG
    debugStream (6311) << "[OneDFSIBC::applyBC] on bcSide " << M_bcSide << " imposing [ A, Q ] = [ " << bc[0] << ", " << bc[1] << " ]\n";
#endif
}

void
OneDFSIBC::applyViscoelasticBC ( const fluxPtr_Type& fluxPtr, matrix_Type& matrix, vector_Type& rhs )
{
    UInt iNode;
    ( M_bcSide == OneDFSI::left ) ? iNode = 0 : iNode = fluxPtr->physics()->data()->numberOfNodes() - 1;

    switch ( M_bcType.find ( OneDFSI::first )->second )
    {
        case OneDFSI::A:
        case OneDFSI::P:
        case OneDFSI::S:

#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
            break;
#endif

        case OneDFSI::Q:
        case OneDFSI::W1:
        case OneDFSI::W2:

            matrix.diagonalize ( iNode, 1 );
            rhs ( iNode ) = 0;

            break;

        default:

            std::cout << "Warning: bcType \"" << M_bcType.find ( OneDFSI::first )->second << "\"not available!" << std::endl;

            break;
    }
}

// ===================================================
// Private Methods
// ===================================================
void
OneDFSIBC::computeMatrixAndRHS ( const Real& time, const Real& timeStep, const fluxPtr_Type& fluxPtr, const bcLine_Type& line,
                                 const container2D_Type& leftEigenvector1, const container2D_Type& leftEigenvector2,
                                 const UInt& iNode, std::map<bcLine_Type, container2D_Type>& bcMatrix, Real& bcRHS )
{
    // This is not general (typical situation):
    //     on first line,  left boundary,  I impose W1
    //     on second line, left boundary,  I impose W2
    //     on first line,  right boundary, I impose W2
    //     on second line, right boundary, I impose W1
    // The code does not check for coherence (you cannot impose the same variable on both lines!)

    // Compute Matrix & RHS
    bcRHS = M_bcFunction[ line ] (time, timeStep);
    switch ( M_bcType[line] )
    {
        case OneDFSI::W1:
            bcMatrix[line] = leftEigenvector1;
            break;
        case OneDFSI::W2:
            bcMatrix[line] = leftEigenvector2;
            break;
        case OneDFSI::A:
            bcMatrix[line][0] = 1.;
            bcMatrix[line][1] = 0.;
            break;
        case OneDFSI::S:
            // The normal stress has opposite sign with respect to the pressure
            bcRHS *= -1;
            // The break here is missing on purpose!
        case OneDFSI::P:
            bcRHS = fluxPtr->physics()->fromPToA ( bcRHS, timeStep, iNode );
            bcMatrix[line][0] = 1.;
            bcMatrix[line][1] = 0.;
            break;
        case OneDFSI::Q:
            // Flow rate is positive with respect to the outgoing normal
            if ( M_bcSide == OneDFSI::left )
            {
                bcRHS *= -1;
            }
            bcMatrix[line][0] = 0.;
            bcMatrix[line][1] = 1.;
            break;
        default:
            std::cout << "\n[OneDFSIBC::computeMatrixAndRHS] Wrong boundary variable as " << line << " condition on bcSide " << M_bcSide;
            break;
    }

#ifdef HAVE_LIFEV_DEBUG
    debugStream (6311) << "[OneDFSIBC::computeMatrixAndRHS] to impose variable "
                       << M_bcType[line] << ", " << line << " line = " << bcMatrix[line][0] << ", " << bcMatrix[line][1] << "\n";
#endif
}

OneDFSIBC::container2D_Type
OneDFSIBC::solveLinearSystem ( const container2D_Type& line1,
                               const container2D_Type& line2,
                               const container2D_Type& rhs ) const
{
#ifdef HAVE_LIFEV_DEBUG
    ASSERT_PRE ( line1.size() == 2 && line2.size() == 2 && rhs.size() == 2,
                 "OneDFSIBC::solveLinearSystem works only for 2D vectors");
#endif

    Real determinant = line1[0] * line2[1] - line1[1] * line2[0];

#ifdef HAVE_LIFEV_DEBUG
    ASSERT ( determinant != 0,
             "Error: the 2x2 system on the boundary is not invertible."
             "\nCheck the boundary conditions.");
#endif

    container2D_Type solution;

    solution[0] = (   line2[1] * rhs[0] - line1[1] * rhs[1] ) / determinant;
    solution[1] = ( - line2[0] * rhs[0] + line1[0] * rhs[1] ) / determinant;

    return solution;
}

}
