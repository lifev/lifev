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

#include <lifemc/lifefem/OneDimensionalBC.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalBC::OneDimensionalBC( const bcSide_Type& bcSide ) :
        M_bcType                    (),
        M_bcSide                    ( bcSide ),
        M_bcFunction                ()
{
}

OneDimensionalBC::OneDimensionalBC( const OneDimensionalBC& bc ) :
        M_bcType                    ( bc.M_bcType ),
        M_bcSide                    ( bc.M_bcSide ),
        M_bcFunction                ( bc.M_bcFunction )
{}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalBC::applyBC( const Real& time, const Real& timeStep, const solution_Type& solution,
                           const fluxPtr_Type& flux, vectorPtrContainer_Type& rhs )
{
    UInt iNode;

    ( M_bcSide == OneDimensional::left ) ? iNode = 0 : iNode = flux->physics()->data()->numberOfNodes() - 1;

    container2D_Type boundaryU;
    boundaryU[0] = (*solution.find("A")->second)(iNode);
    boundaryU[1] = (*solution.find("Q")->second)(iNode);

    // Eigenvalues and eigenvectors of the jacobian diffFlux (= dF/dU = H)
    container2D_Type eigenvalues;
    container2D_Type leftEigenvector1, leftEigenvector2;

    flux->eigenValuesEigenVectors( boundaryU[0], boundaryU[1],
                                   eigenvalues, leftEigenvector1, leftEigenvector2, iNode );

    std::map<bcLine_Type, container2D_Type> bcMatrix;
    bcMatrix[ OneDimensional::first ]  = container2D_Type();
    bcMatrix[ OneDimensional::second ] = container2D_Type();

    container2D_Type bcRHS;

    // First line of Matrix and RHS
    computeMatrixAndRHS( time, timeStep, flux, OneDimensional::first,
                         leftEigenvector1, leftEigenvector2, iNode, bcMatrix, bcRHS[0] );

    // Second line of Matrix and RHS
    computeMatrixAndRHS( time, timeStep, flux, OneDimensional::second,
                         leftEigenvector1, leftEigenvector2, iNode, bcMatrix, bcRHS[1] );


    container2D_Type bc = solveLinearSystem( bcMatrix[OneDimensional::first], bcMatrix[OneDimensional::second], bcRHS );

    // Set the BC in the RHS
    (*rhs[0])( iNode ) = bc[0];
    (*rhs[1])( iNode ) = bc[1];

#ifdef HAVE_LIFEV_DEBUG
    Debug(6311) << "[OneDimensionalBC::applyBC] on bcSide " << M_bcSide << " imposing [ A, Q ] = [ " << bc[0] << ", " << bc[1] << " ]\n";
#endif
}

void
OneDimensionalBC::applyViscoelasticBC( const Real& timeStep, const vector_Type& area, const vector_Type& flowRate, const fluxPtr_Type& flux, matrix_Type& matrix, vector_Type& rhs )
{
    UInt iNode;
    UInt iNodeInternal;

    if ( M_bcSide == OneDimensional::left )
    {
        iNode = 0;
        iNodeInternal = iNode + 1;
    }
    else
    {
        iNode = flux->physics()->data()->numberOfNodes() - 1;
        iNodeInternal = iNode - 1;
    }

    matrix.globalAssemble();
    switch ( M_bcType.find( OneDimensional::first )->second )
    {
    case OneDimensional::Q:
    case OneDimensional::W1:
    case OneDimensional::W2:

        matrix.diagonalize( iNode, 1 );
        rhs( iNode ) = 0;

        break;

    case OneDimensional::A:
    case OneDimensional::P:
    {
        Real massCoefficient = 1 / ( 0.5 * ( area[ iNode ] + area[ iNodeInternal ] ) );
        Real stiffnessCoefficient  = timeStep * 0.5* ( flux->physics()->data()->viscoelasticCoefficient( iNode ) + flux->physics()->data()->viscoelasticCoefficient( iNodeInternal ) )
                                   / flux->physics()->data()->densityRho() * massCoefficient * std::sqrt( massCoefficient );
        if ( M_bcSide == OneDimensional::left )
            rhs( iNode ) -= stiffnessCoefficient * flux->physics()->data()->computeSpatialDerivativeAtNode( flowRate, iNode, 1 );
        else
            rhs( iNode ) += stiffnessCoefficient * flux->physics()->data()->computeSpatialDerivativeAtNode( flowRate, iNode, 1 );

        break;
    }
    default:

        std::cout << "Warning: bcType \"" << M_bcType.find( OneDimensional::first )->second << "\"not available!" << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================
void
OneDimensionalBC::computeMatrixAndRHS( const Real& time, const Real& timeStep, const fluxPtr_Type& flux, const bcLine_Type& line,
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
    bcRHS = M_bcFunction[ line ](time, timeStep);
    switch ( M_bcType[line] )
    {
    case OneDimensional::W1:
        bcMatrix[line] = leftEigenvector1;
        break;
    case OneDimensional::W2:
        bcMatrix[line] = leftEigenvector2;
        break;
    case OneDimensional::A:
        bcMatrix[line][0] = 1.;
        bcMatrix[line][1] = 0.;
        break;
    case OneDimensional::P:
        bcRHS = flux->physics()->fromPToA( bcRHS, timeStep, iNode );
        bcMatrix[line][0] = 1.;
        bcMatrix[line][1] = 0.;
        break;
    case OneDimensional::Q:
        // Flow rate is positive with respect to the outgoing normal
        if ( M_bcSide == OneDimensional::left )
            bcRHS *= -1;
        bcMatrix[line][0] = 0.;
        bcMatrix[line][1] = 1.;
        break;
    default:
        std::cout << "\n[OneDimensionalBC::computeMatrixAndRHS] Wrong boundary variable as " << line << " condition on bcSide " << M_bcSide;
    }

#ifdef HAVE_LIFEV_DEBUG
    Debug(6311) << "[OneDimensionalBC::computeMatrixAndRHS] to impose variable "
    << M_bcType[line] << ", " << line << " line = " << bcMatrix[line][0] << ", " << bcMatrix[line][1] << "\n";
#endif
}

OneDimensionalBC::container2D_Type
OneDimensionalBC::solveLinearSystem( const container2D_Type& line1,
                                     const container2D_Type& line2,
                                     const container2D_Type& rhs ) const
{
#ifdef HAVE_LIFEV_DEBUG
    ASSERT_PRE( line1.size() == 2 && line2.size() == 2 && rhs.size() == 2,
                "OneDimensionalBC::solveLinearSystem works only for 2D vectors");
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
