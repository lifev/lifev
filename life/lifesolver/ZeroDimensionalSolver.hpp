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
 *  @Rythmos solver
 *
 *  @version 1.0
 *  @date 16-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalSolver_H
#define ZeroDimensionalSolver_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// Includes for Rythmos:
#include <Rythmos_ConfigDefs.h>
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
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraModelEvaluator.hpp>
#include <Thyra_NonlinearSolver_NOX.hpp>

#include <Thyra_DiagonalEpetraLinearOpWithSolveFactory.hpp>

// Includes for Stratimikos:
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

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

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LIFEV - MATHCARD
#include <lifemc/lifesolver/ZeroDimensionalRythmosSolverInterface.hpp>
#include <lifemc/lifesolver/ZeroDimensionalData.hpp>

//! Rhytmos methods
enum EMethod { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF, METHOD_IRK };

//! time step method
enum STEP_METHOD { STEP_METHOD_FIXED, STEP_METHOD_VARIABLE };

//! ZeroDimentional Solver
class ZeroDimensionalSolver
{
public:

    //! Constructor
    explicit ZeroDimensionalSolver(int numGlobalElements,
                                   boost::shared_ptr< Epetra_Comm> comm,
                                   LifeV::zeroDimensionalCircuitDataPtr_Type circuitData);
    //! Destructor
    virtual ~ZeroDimensionalSolver() {}

    //! setup solver
    void setup(const LifeV::ZeroDimensionalData::solverData_Type&  data);

    //! integrate the system between t1 and t2
    void takeStep(double t1, double t2);
private:

    LifeV::rythmosSolverInterfacePtr_Type      M_solverInterface;
    LifeV::rythmosModelInterfacePtr_Type       M_modelInterface;
    LifeV::rythmosSolverInterfacePtrRCP_Type   M_solverInterfaceRCP;
    LifeV::rythmosModelInterfacePtrRCP_Type    M_modelInterfaceRCP;
    boost::shared_ptr< Epetra_Comm>            M_comm;
    Teuchos::RCP< Epetra_Comm>                 M_commRCP;
    Teuchos::RCP<Rythmos::StepperBase<double> >M_stepperPtr;
    Teuchos::RCP<Teuchos::FancyOStream>        M_out;
    STEP_METHOD                                M_step_method;
    double                                     M_finalTime;
    double                                     M_startTime;
    int                                        M_numberTimeStep;
    int                                        M_outputLevel;
    Teuchos::EVerbosityLevel                   M_outputLevelTeuchos;
    std::string                                M_method;

};

#endif //ZeroDimensionalSolver_H
