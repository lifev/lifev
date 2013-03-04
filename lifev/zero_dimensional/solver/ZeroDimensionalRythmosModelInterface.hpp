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
 *  @version alpha (experimental)
 *
 *  @date 21-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalRythmosModelInterface_H
#define ZeroDimensionalRythmosModelInterface_H 1

// Include definitions
#include <lifev/zero_dimensional/solver/ZeroDimensionalDefinitions.hpp>

#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H> // base class
#include <NOX_Epetra_Interface_Jacobian.H> // base class
#include <NOX_Epetra_Interface_Preconditioner.H> // base class
#include <NOX_Thyra.H>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LIFEV
#include <lifev/zero_dimensional/solver/ZeroDimensionalCircuitData.hpp>

namespace NOX
{
namespace Parameter
{
class List;
}
}
namespace LifeV
{

//! Rhytmos model interface.
/*!
 * Rhytmos solver interface will communicate with this class.
 * this class have access to circuit data. The main task of this class is to
 * provide the residual and jacobian at every step to Rhythmos solver interface.
 */
class RythmosModelInterface : public NOX::Epetra::Interface::Required,
    public NOX::Epetra::Interface::Jacobian,
    public NOX::Epetra::Interface::Preconditioner
{
public:

    //! Constructor
    RythmosModelInterface (Int NumGlobalElements, Epetra_Comm* comm, zeroDimensionalCircuitDataPtr_Type circuitData);

    //! Destructor
    virtual ~RythmosModelInterface();

    //! This method is empty.
    virtual bool computeF (const Epetra_Vector& x, Epetra_Vector& FVec, FillType fillType = Residual);

    //! This method is empty.
    virtual bool computeJacobian (const Epetra_Vector& x, Epetra_Operator& Jac);

    //! This method is empty
    virtual bool computePrecMatrix (const Epetra_Vector& x);

    //! This method is empty
    virtual bool computePreconditioner (const Epetra_Vector& x, Epetra_Operator& Prec, Teuchos::ParameterList* precParams = 0);

    //! get solution vector
    Epetra_Vector& getSolutionY();

    //! get derivative of solution vector respect to time
    Epetra_Vector& getSolutionYp();

    Epetra_Map& getMap();

    Epetra_CrsGraph& getGraph();

    //! this method empty.
    virtual bool evaluate (Real t, const Epetra_Vector* x, Epetra_Vector* f );

    //! compute Implicit residual.
    virtual bool evaluateFImplicit (const Real& t, const Epetra_Vector* x, const Epetra_Vector* x_dot, Epetra_Vector* f );

    //! compute jacobian.
    virtual bool evaluateWImplicit (const Real& t, const Real& alpha, const Real& beta, const Epetra_Vector* x, const Epetra_Vector* x_dot, Epetra_CrsMatrix* W );

    virtual bool initializeSolnY();

    bool initializeSolnY (const vectorEpetra_Type& y);

    virtual bool initializeSolnYp();

    bool initializeSolnYp (const vectorEpetra_Type& yp);

    //! hafter complete Rythmos step, this method will update circuit data.
    void extractSolution (const Real& t1, const vectorEpetra_Type& y , const vectorEpetra_Type& yp );

    Int numCircuitElements()
    {
        return M_numCircuitElements;
    }

protected:
    Int                                             M_numCircuitElements; // Total Number of elements
    Int                                             M_numMyElements; // Number of elements owned by this process
    Int                                             M_myPID; // Process number
    Int                                             M_numProc; // Total number of processes

    Epetra_CrsGraph*                                M_graph;
    boost::shared_ptr<Epetra_CrsGraph>              M_graphSharedPtr;
    Epetra_Comm*                                    M_comm;
    boost::shared_ptr<Epetra_Comm>                  M_commSharedPtr;
    Epetra_Map*                                     M_standardMap;
    Epetra_Vector*                                  M_initialSolutionY;
    Epetra_Vector*                                  M_initialSolutionYp;

    zeroDimensionalCircuitDataPtr_Type              M_circuitData;
    matrixPtr_Type                                  M_A;
    matrixPtr_Type                                  M_B;
    vectorPtr_Type                                  M_C;
    vectorEpetraPtr_Type                            M_Y0;
    vectorEpetraPtr_Type                            M_Yp0;
    boost::shared_ptr<MapEpetra>                    M_mapEpetraPtr;
    vectorEpetraPtr_Type                            M_fA;//dummy variable
    vectorEpetraPtr_Type                            M_fB;//dummy variable
};

typedef boost::shared_ptr< RythmosModelInterface >  rythmosModelInterfacePtr_Type;
typedef Teuchos::RCP< RythmosModelInterface >       rythmosModelInterfacePtrRCP_Type;
} // LifeV namespace

#endif /* HAVE_NOX_THYRA && HAVE_TRILINOS_RYTHMOS */

#endif //ZeroDimensionalRythmosModelInterface_H
