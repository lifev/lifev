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
 *  @brief Rythmos Model Interface
 *
 *  @date 21-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/zero_dimensional/solver/ZeroDimensionalRythmosModelInterface.hpp>

namespace LifeV
{

#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
// ===================================================
// Constructors
// ===================================================
RythmosModelInterface::RythmosModelInterface ( Int numCircuitElements,
                                               Epetra_Comm* comm,
                                               zeroDimensionalCircuitDataPtr_Type circuitData ) :
    M_numCircuitElements ( numCircuitElements ), M_numMyElements ( 0 ), M_myPID ( comm->MyPID() ), M_numProc ( comm->NumProc() ),

    M_standardMap ( 0 ), M_initialSolutionY ( 0 ), M_initialSolutionYp ( 0 ), M_circuitData ( circuitData )
{
    M_comm = comm->Clone();
    ( *M_comm ).PrintInfo ( std::cout );
    M_commSharedPtr.reset ( M_comm );
    M_mapEpetraPtr.reset ( new MapEpetra ( M_numCircuitElements,
                                           M_commSharedPtr ) );
    M_standardMap = ( M_mapEpetraPtr->map ( Unique ) ).get();
    M_numMyElements = M_standardMap->NumMyElements();
    M_A.reset ( new matrix_Type ( *M_mapEpetraPtr,
                                  M_numCircuitElements ) );
    M_B.reset ( new matrix_Type ( *M_mapEpetraPtr,
                                  M_numCircuitElements ) );

    for ( Int i = 0; i < M_numCircuitElements; i++ )
    {
        for ( Int j = 0; j < M_numCircuitElements; j++ )
        {
            M_A->addToCoefficient ( i, j, 1.0 );
            M_B->addToCoefficient ( i, j, 1.0 );
        }
    }
    M_A->globalAssemble();
    M_B->globalAssemble();

    M_graph = new Epetra_CrsGraph ( M_A->matrixPtr()->Graph() );
    M_graphSharedPtr.reset ( M_graph );
    M_C.reset ( new vector_Type ( *M_mapEpetraPtr ) );

    M_fA.reset ( new vectorEpetra_Type ( *M_standardMap ) );
    M_fB.reset ( new vectorEpetra_Type ( *M_standardMap ) );

    M_Y0.reset ( new vectorEpetra_Type ( *M_standardMap ) );
    M_Yp0.reset ( new vectorEpetra_Type ( *M_standardMap ) );

    M_initialSolutionY = M_Y0.get();
    M_initialSolutionYp = M_Yp0.get();
    ;
    initializeSolnY();
    initializeSolnYp();
}

// Destructor
RythmosModelInterface::~RythmosModelInterface()
{

    M_circuitData.reset();
    M_mapEpetraPtr.reset();
    M_graphSharedPtr.reset();
    M_A.reset();
    M_B.reset();
    M_C.reset();
    M_Y0.reset();
    M_Yp0.reset();
}

bool RythmosModelInterface::computeF ( __attribute__ ( (unused) )  const Epetra_Vector& x,
                                       __attribute__ ( (unused) )  Epetra_Vector& FVec,
                                       __attribute__ ( (unused) )  NOX::Epetra::Interface::Required::FillType fillType )
{
    std::cout << "Warning: I should not be here, 0D, Rythmos Model Interface";
    return false;
}

bool RythmosModelInterface::computeJacobian ( __attribute__ ( (unused) )  const Epetra_Vector& x,
                                              __attribute__ ( (unused) )  Epetra_Operator& my_Jac )
{
    std::cout << "Warning: I should not be here, 0D, Rythmos Model Interface";
    return false;
}

bool RythmosModelInterface::computePrecMatrix ( __attribute__ ( (unused) )  const Epetra_Vector& x )
{
    std::cout << "Warning: I should not be here, 0D, Rythmos Model Interface";
    return false;
}
bool RythmosModelInterface::computePreconditioner ( __attribute__ ( (unused) )   const Epetra_Vector& x,
                                                    __attribute__ ( (unused) )   Epetra_Operator& my_Prec,
                                                    __attribute__ ( (unused) )   Teuchos::ParameterList* precParams )
{
    cout << "ERROR: Interface::preconditionVector()  " << endl;
    std::cout << "Error: I should not be here, 0D, Rythmos Model Interface";
    throw "Interface Error";
}

bool RythmosModelInterface::evaluate ( __attribute__ ( (unused) )   Real t,
                                       __attribute__ ( (unused) )  const Epetra_Vector* x,
                                       __attribute__ ( (unused) )   Epetra_Vector* f )
{
    std::cout << "Warning: I should not be here. 0D, Rythmos Inreface, ::evaluate" << std::endl;
    return true;
}

Epetra_Vector& RythmosModelInterface::getSolutionY()
{
    return *M_initialSolutionY;
}

Epetra_Vector& RythmosModelInterface::getSolutionYp()
{
    return *M_initialSolutionYp;
}

Epetra_Map& RythmosModelInterface::getMap()
{
    return *M_standardMap;
}

bool RythmosModelInterface::initializeSolnY()
{
    M_initialSolutionY->PutScalar ( 0 );
    return true;
}

bool RythmosModelInterface::initializeSolnY ( const vectorEpetra_Type& y )
{
    ( *M_initialSolutionY ) = y;
    return true;
}

bool RythmosModelInterface::initializeSolnYp ( const vectorEpetra_Type& yp )
{
    ( *M_initialSolutionYp ) = yp;
    return true;
}

bool RythmosModelInterface::initializeSolnYp()
{
    M_initialSolutionYp->PutScalar ( 0 );
    return true;
}
Epetra_CrsGraph& RythmosModelInterface::getGraph()
{
    return *M_graph;
}

bool RythmosModelInterface::evaluateFImplicit ( const Real& t,
                                                const Epetra_Vector* x,
                                                const Epetra_Vector* x_dot,
                                                Epetra_Vector* f )
{

    //updateCircuitData from Y and Yp
#ifdef HAVE_LIFEV_DEBUG
    x->Print (std::cout);
    x_dot->Print (std::cout);
#endif
    M_circuitData->updateCircuitDataFromY ( t,
                                            x,
                                            x_dot );
    M_circuitData->updateABC ( *M_A,
                               *M_B,
                               *M_C );
    //calculate the sum
#ifdef HAVE_LIFEV_DEBUG
    M_A->matrixPtr()->Print (std::cout);
    M_B->matrixPtr()->Print (std::cout);
    M_C->epetraVector().Print (std::cout);
#endif

    M_A->matrixPtr()->Multiply1 ( false,
                                  *x_dot,
                                  *M_fA );
#ifdef HAVE_LIFEV_DEBUG
    M_A->matrixPtr()->Print (std::cout);
#endif

    M_B->matrixPtr()->Multiply1 ( false,
                                  *x,
                                  *M_fB );
#ifdef HAVE_LIFEV_DEBUG
    M_B->matrixPtr()->Print (std::cout);
#endif
    for ( Int i = 0; i < M_numCircuitElements; i++ )
    {
        ( *f ) [i] = ( *M_fA ) [i] + ( *M_fB ) [i] + ( *M_C ) [i];
    }
#ifdef HAVE_LIFEV_DEBUG
    f->Print (std::cout);
#endif
    return true;
}

bool RythmosModelInterface::evaluateWImplicit ( const Real& t,
                                                const Real& alpha,
                                                const Real& beta,
                                                const Epetra_Vector* x,
                                                const Epetra_Vector* x_dot,
                                                Epetra_CrsMatrix* W )
{
    //updateCircuitData from Y and Yp
#ifdef HAVE_LIFEV_DEBUG
    x->Print (std::cout);
    x_dot->Print (std::cout);
#endif
    M_circuitData->updateCircuitDataFromY ( t,
                                            x,
                                            x_dot );
    M_circuitData->updateABC ( *M_A,
                               *M_B,
                               *M_C );
#ifdef HAVE_LIFEV_DEBUG
    M_A->matrixPtr()->Print (std::cout);
    M_B->matrixPtr()->Print (std::cout);
    M_C->epetraVector().Print (std::cout);
#endif
    M_A->operator*= ( alpha );
    M_B->operator*= ( beta );
    M_A->operator+= ( *M_B );
    Epetra_CrsMatrix* tmp = M_A->matrixPtr().get();
    ( *W ) = ( *tmp );
#ifdef HAVE_LIFEV_DEBUG
    W->Print (std::cout);
#endif
    return true;
}

void RythmosModelInterface::extractSolution ( const Real& t,
                                              const vectorEpetra_Type& y,
                                              const vectorEpetra_Type& yp )
{
    M_circuitData->extractSolutionFromY ( t, y, yp );
}

#endif /* HAVE_NOX_THYRA && HAVE_TRILINOS_RYTHMOS */

} // LifeV namespace
