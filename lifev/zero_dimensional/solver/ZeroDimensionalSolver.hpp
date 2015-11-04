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
 *  @brief Rythmos solver
 *  @version alpha (experimental)
 *
 *  @date 16-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalSolver_H
#define ZeroDimensionalSolver_H 1


// Include definitions
#include <lifev/zero_dimensional/solver/ZeroDimensionalDefinitions.hpp>

// Includes for Rythmos:
#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
#include <Rythmos_StepperBase.hpp>
#include <Rythmos_ForwardEulerStepper.hpp>
#include <Rythmos_BackwardEulerStepper.hpp>
#include <Rythmos_ExplicitRKStepper.hpp>
#include <Rythmos_ImplicitBDFStepper.hpp>
#include <Rythmos_ImplicitRKStepper.hpp>
#include <Rythmos_RKButcherTableau.hpp>
#include <Rythmos_RKButcherTableauBuilder.hpp>
#include <Rythmos_TimeStepNonlinearSolver.hpp>

// Includes for Thyra:
#include <Thyra_DiagonalEpetraLinearOpWithSolveFactory.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraModelEvaluator.hpp>
#include <Thyra_NonlinearSolver_NOX.hpp>

// Includes for Stratimikos:
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#endif /* HAVE_NOX_THYRA && HAVE_TRILINOS_RYTHMOS */

// Includes for Teuchos:
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCPBoostSharedPtrConversions.hpp>


// LIFEV includes
#include <lifev/zero_dimensional/solver/ZeroDimensionalRythmosSolverInterface.hpp>
#include <lifev/zero_dimensional/solver/ZeroDimensionalData.hpp>

namespace LifeV
{

//! Rhytmos methods
enum EMethod { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF, METHOD_IRK };

//! time step method
enum STEP_METHOD { STEP_METHOD_FIXED, STEP_METHOD_VARIABLE };

//! ZeroDimensional Solver
#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
class ZeroDimensionalSolver
{
public:

    //! Constructor
    explicit ZeroDimensionalSolver ( Int numCircuitElements,
                                     std::shared_ptr< Epetra_Comm> comm,
                                     zeroDimensionalCircuitDataPtr_Type circuitData );
    //! Destructor
    virtual ~ZeroDimensionalSolver() {}

    //! setup solver
    void setup ( const ZeroDimensionalData::solverData_Type& data );

    //! integrate the system between t1 and t2
    void takeStep (Real t1, Real t2);

private:

    rythmosSolverInterfacePtr_Type                M_solverInterface;
    rythmosModelInterfacePtr_Type                 M_modelInterface;
    rythmosSolverInterfacePtrRCP_Type             M_solverInterfaceRCP;
    rythmosModelInterfacePtrRCP_Type              M_modelInterfaceRCP;
    std::shared_ptr< Epetra_Comm>               M_comm;
    Teuchos::RCP< Epetra_Comm>                    M_commRCP;
    Teuchos::RCP<Rythmos::StepperBase<Real> >     M_stepperPtr;
    Teuchos::RCP<Teuchos::FancyOStream>           M_out;
    STEP_METHOD                                   M_stepMethod;
    Real                                          M_finalTime;
    Real                                          M_startTime;
    Int                                           M_numberTimeStep;
    Int                                           M_outputLevel;
    Teuchos::EVerbosityLevel                      M_outputLevelTeuchos;
    std::string                                   M_method;

};

#else

class ZeroDimensionalSolver
{
public:

    //! Constructor
    explicit ZeroDimensionalSolver (Int /*numCircuitElements*/,
                                    std::shared_ptr<Epetra_Comm> /*comm*/,
                                    zeroDimensionalCircuitDataPtr_Type /*circuitData*/) {}
    //! Destructor
    virtual ~ZeroDimensionalSolver() {}

    //! setup solver
    void setup (const ZeroDimensionalData::solverData_Type& /*data*/) {}

    //! integrate the system between t1 and t2
    void takeStep (Real /*t1*/, Real /*t2*/) {}
};

#endif /* HAVE_NOX_THYRA && HAVE_TRILINOS_RYTHMOS */

} // LifeV namespace

#endif //ZeroDimensionalSolver_H
