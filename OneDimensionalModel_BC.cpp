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
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *  @date 01-28-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 */

#include <lifemc/lifefem/OneDimensionalModel_BC.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BC::OneDimensionalModel_BC( const OneD_BCSide& side ) :
        M_bcType                    (),
        M_bcSide                    ( side ),
        M_bcFunction                (),
        M_isInternal                ( false ),
        M_bcMatrix                  (),
        M_bcRHS                     ()
{
    M_bcMatrix[ OneD_first ]  = container2D_Type();
    M_bcMatrix[ OneD_second ] = container2D_Type();
}

OneDimensionalModel_BC::OneDimensionalModel_BC( const OneDimensionalModel_BC& BC ) :
        M_bcType                    ( BC.M_bcType ),
        M_bcSide                    ( BC.M_bcSide ),
        M_bcFunction                ( BC.M_bcFunction ),
        M_isInternal                ( BC.M_isInternal ),
        M_bcMatrix                  ( BC.M_bcMatrix ),
        M_bcRHS                     ( BC.M_bcRHS )
{}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_BC::applyBC( const Real&             time,
                                 const Real&             timeStep,
                                 const Solution_Type&    solution,
                                 const Flux_PtrType&     flux,
                                 container2D_Type& BC )
{
    ASSERT_PRE( BC.size() == 2, "applyBC works only for 2D vectors");

    if ( M_isInternal )
        Debug(6311) << "[OneDimensionalModel_BC::compute_resBC] found internal boundary\n";
    else
    {
        //Container2D_Type leftEigenvector_first, leftEigenvector_second;

        UInt dof;
        ( M_bcSide == OneD_left ) ? dof = 0 : dof = flux->physics()->Data()->NumberOfNodes() - 1;

        container2D_Type U_boundary;
        U_boundary[0] = (*solution.find("A")->second)(dof + 1);
        U_boundary[1] = (*solution.find("Q")->second)(dof + 1);

        // Eigenvalues and eigenvectors of the jacobian diffFlux (= dF/dU = H)
        container2D_Type eigenvalues;
        container2D_Type leftEigenvector1, leftEigenvector2;

        flux->eigenValuesEigenVectors( U_boundary[0], U_boundary[1],
                                       eigenvalues, leftEigenvector1, leftEigenvector2, dof );

        computeMatrixAndRHS( time, timeStep, flux, OneD_first,
                             leftEigenvector1, leftEigenvector2, dof, M_bcRHS[0]);

        computeMatrixAndRHS( time, timeStep, flux, OneD_second,
                             leftEigenvector1, leftEigenvector2, dof, M_bcRHS[1]);

        BC = solveLinearSystem( M_bcMatrix[OneD_first], M_bcMatrix[OneD_second], M_bcRHS );
    }

    Debug(6311) << "[OneDimensionalModel_BC::applyBC] on side " << M_bcSide
    << " imposing [ A, Q ] = [ " << BC[0] << ", " << BC[1] << " ]\n";
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

/*
void
OneDimensionalModel_BC::setMatrixRow( const OneD_BCLine& line, const Container2D_Type& matrixrow )
{
    M_bcMatrix[line] = matrixrow;
}
*/

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
OneDimensionalModel_BC::computeMatrixAndRHS( const Real& time, const Real& timeStep, const Flux_PtrType& flux, const OneD_BCLine& line,
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
    case OneD_W1:
        M_bcMatrix[line] = leftEigenvector1;
        break;
    case OneD_W2:
        M_bcMatrix[line] = leftEigenvector2;
        break;
    case OneD_A:
        M_bcMatrix[line][0] = 1.;
        M_bcMatrix[line][1] = 0.;
        break;
    case OneD_P:
        rhs = flux->physics()->A_from_P( rhs, dof );
        M_bcMatrix[line][0] = 1.;
        M_bcMatrix[line][1] = 0.;
        break;
    case OneD_Q:
        // Flow rate is positive with respect to the outgoing normal
        if ( M_bcSide == OneD_left )
            rhs *= -1;
        M_bcMatrix[line][0] = 0.;
        M_bcMatrix[line][1] = 1.;
        break;
    default:
        std::cout << "\n[OneDimensionalModel_BC::compute_resBC] Wrong boundary variable as " << line
                  << " condition on side " << M_bcSide;
    }

    Debug(6311) << "[OneDimensionalModel_BC::compute_MatrixAndRHS] to impose variable "
    << M_bcType[line] << ", " << line << " line = "
    << M_bcMatrix[line][0] << ", "
    << M_bcMatrix[line][1] << "\n";
}

container2D_Type
OneDimensionalModel_BC::solveLinearSystem( const container2D_Type& line1,
                                           const container2D_Type& line2,
                                           const container2D_Type& rhs ) const
{
    ASSERT_PRE( line1.size() == 2 && line2.size() == 2 && rhs.size() == 2,
                "_solveLinearSyst2x2 works only for 2D vectors");

    Real determinant = line1[0] * line2[1] - line1[1] * line2[0];

    ASSERT( determinant != 0,
            "Error: the 2x2 system on the boundary is not invertible."
            "\nCheck the boundary conditions.");

    container2D_Type solution;

    solution[0] = (   line2[1] * rhs[0] - line1[1] * rhs[1] ) / determinant;
    solution[1] = ( - line2[0] * rhs[0] + line1[0] * rhs[1] ) / determinant;

    return solution;
}

}
