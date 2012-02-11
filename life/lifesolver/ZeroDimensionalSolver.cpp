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
 *  @Rythmos solver setup
 *
 *  @version 1.0
 *  @date 16-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <life/lifesolver/ZeroDimensionalSolver.hpp>

ZeroDimensionalSolver::ZeroDimensionalSolver(int numGlobalElements,
                                             boost::shared_ptr< Epetra_Comm> comm,
                                             LifeV::zeroDimensionalCircuitDataPtr_Type circuitData){
    M_comm.swap(comm);
    M_commRCP.reset();
    M_commRCP = Teuchos::rcp(M_comm);
    M_modelInterface.reset(new LifeV::RythmosModelInterface(numGlobalElements,
                                                     M_comm.operator ->(),
                                                     circuitData));
    M_modelInterfaceRCP = Teuchos::rcp(M_modelInterface);

    M_solverInterface.reset(new LifeV::RythmosSolverInterface(numGlobalElements,
                                                       M_commRCP,
                                                       M_modelInterfaceRCP));

    M_solverInterfaceRCP.reset();
    M_solverInterfaceRCP = Teuchos::rcp(M_solverInterface);

    M_out                = Teuchos::VerboseObjectBase::getDefaultOStream();

}
void
ZeroDimensionalSolver::setup(const LifeV::ZeroDimensionalData::solverData_Type&  data){
    std::string commandLine = "--linear-solver-params-used-file=";
    commandLine.append(data.linearSolverParamsFile);
    char * argv[1];
    argv[0] =  new char[commandLine.size()];
    for (int i = 0; i < ((int)commandLine.size()); ++i) {
        (argv[0])[i] = commandLine[i];
    }
    int argc = 1;
    bool verbose;
    bool  success = true; // determine if the run was successfull
    try { // catch exceptions
      if (data.fixTimeStep){
         M_step_method = STEP_METHOD_FIXED;
      }else{
     M_step_method = STEP_METHOD_VARIABLE;
      }


    int numElements = M_modelInterface->numGlobalElements(); // number of elements in vector
    EMethod method_val = METHOD_BE;
    double maxError ;
    double reltol;
    double abstol;
    int maxOrder ;
    bool useNOX ;

    if (!data.method.compare("BE"))
        method_val = METHOD_BE;
    if (!data.method.compare("BDF"))
        method_val = METHOD_BDF;
    if (!data.method.compare("IRK"))
        method_val = METHOD_IRK;

    M_numberTimeStep = data.numberTimeStep;
    maxError = data.maxError;
    reltol = data.reltol;
    abstol = data.abstol;
    maxOrder = data.maxOrder;
    verbose = data.verbose;
    M_outputLevel = data.verbosLevel;
    M_outputLevel = min(max(M_outputLevel,-1),4);
    useNOX = data.useNOX;
    switch (M_outputLevel) {
        case -1:
            M_outputLevelTeuchos = Teuchos::VERB_DEFAULT;
            break;
        case 0:
            M_outputLevelTeuchos = Teuchos::VERB_NONE;
                break;
        case 1:
            M_outputLevelTeuchos = Teuchos::VERB_LOW;
                break;
        case 2:
            M_outputLevelTeuchos = Teuchos::VERB_MEDIUM;
                break;
        case 3:
            M_outputLevelTeuchos = Teuchos::VERB_HIGH;
                break;
        case 4:
            M_outputLevelTeuchos = Teuchos::VERB_EXTREME;
                break;
        default:
            break;
    }
    std::string extraLSParamsFile = data.extraLSParamsFile;
    // Parse the command-line options:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
    Stratimikos::DefaultLinearSolverBuilder lowsfCreator;
    lowsfCreator.setupCLP(&clp);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    delete  argv [0];
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return ;
    lowsfCreator.readParameters(M_out.get());
    if(extraLSParamsFile.length())
      Teuchos::updateParametersFromXmlFile( "./"+extraLSParamsFile, &*lowsfCreator.getNonconstParameterList() );
    *M_out << "\nThe parameter list after being read in:\n";
    lowsfCreator.getParameterList()->print(*M_out,2,true,false);

    // Set up the parameter list for the application:
    Teuchos::ParameterList params;
    params.set( "NumElements", numElements );
    // Create the factory for the LinearOpWithSolveBase object
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      W_factory;
      W_factory = lowsfCreator.createLinearSolveStrategy("");
      *M_out
        << "\nCreated a LinearOpWithSolveFactory described as:\n"
        << Teuchos::describe(*W_factory,Teuchos::VERB_MEDIUM);

    // create interface to problem
    Teuchos::RCP<Thyra::ModelEvaluator<double> >
      model = Teuchos::rcp(new Thyra::EpetraModelEvaluator(M_solverInterfaceRCP,W_factory));

   Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();

    std::string method;
    if (method_val == METHOD_BE) {
      Teuchos::RCP<Thyra::NonlinearSolverBase<double> > nonlinearSolver;
      if (useNOX) {
        Teuchos::RCP<Thyra::NOXNonlinearSolver> _nonlinearSolver = Teuchos::rcp(new Thyra::NOXNonlinearSolver);
        nonlinearSolver = _nonlinearSolver;
      }
      else {
        Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<double> >
          _nonlinearSolver = Teuchos::rcp(new Rythmos::TimeStepNonlinearSolver<double>());
        Teuchos::RCP<Teuchos::ParameterList>
          nonlinearSolverPL = Teuchos::parameterList();
        nonlinearSolverPL->set("Default Tol",double(maxError));
        _nonlinearSolver->setParameterList(nonlinearSolverPL);
        nonlinearSolver = _nonlinearSolver;
      }
      model->setDefaultVerbLevel(M_outputLevelTeuchos);
      nonlinearSolver->setVerbLevel(M_outputLevelTeuchos);
      M_stepperPtr = Teuchos::rcp(new Rythmos::BackwardEulerStepper<double>(model,nonlinearSolver));
      M_stepperPtr->setVerbLevel(M_outputLevelTeuchos);
      model->describe(*M_out,M_outputLevelTeuchos);
      model->description();
      nonlinearSolver->describe(*M_out,M_outputLevelTeuchos);
      nonlinearSolver->description();
      M_stepperPtr->describe(*M_out,M_outputLevelTeuchos);
      M_stepperPtr->description();
      method = "Backward Euler";
    } else if (method_val == METHOD_BDF) {
      Teuchos::RCP<Thyra::NonlinearSolverBase<double> > nonlinearSolver;
      if (useNOX) {
        Teuchos::RCP<Thyra::NOXNonlinearSolver> _nonlinearSolver = Teuchos::rcp(new Thyra::NOXNonlinearSolver);
        nonlinearSolver = _nonlinearSolver;
      }
      else {
        Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<double> >
          _nonlinearSolver = Teuchos::rcp(new Rythmos::TimeStepNonlinearSolver<double>());
        Teuchos::RCP<Teuchos::ParameterList>
          nonlinearSolverPL = Teuchos::parameterList();
        nonlinearSolverPL->set("Default Tol",double(maxError));
        _nonlinearSolver->setParameterList(nonlinearSolverPL);
        nonlinearSolver = _nonlinearSolver;
      }
      nonlinearSolver->setVerbLevel(M_outputLevelTeuchos);
      Teuchos::RCP<Teuchos::ParameterList> BDFparams = Teuchos::rcp(new Teuchos::ParameterList);
      Teuchos::RCP<Teuchos::ParameterList> BDFStepControlPL = Teuchos::sublist(BDFparams, "Step Control Settings");
      BDFStepControlPL->set( "maxOrder", maxOrder );
      BDFStepControlPL->set( "relErrTol", reltol );
      BDFStepControlPL->set( "absErrTol", abstol );
      M_stepperPtr = Teuchos::rcp(new Rythmos::ImplicitBDFStepper<double>(model,nonlinearSolver,BDFparams));
      method = "Implicit BDF";
    } else if (method_val == METHOD_IRK) {
      Teuchos::RCP<Thyra::NonlinearSolverBase<double> > nonlinearSolver;
      if (useNOX) {
        Teuchos::RCP<Thyra::NOXNonlinearSolver> _nonlinearSolver = Teuchos::rcp(new Thyra::NOXNonlinearSolver);
        nonlinearSolver = _nonlinearSolver;
      }
      else {
        Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<double> >
          _nonlinearSolver = Teuchos::rcp(new Rythmos::TimeStepNonlinearSolver<double>());
        Teuchos::RCP<Teuchos::ParameterList>
          nonlinearSolverPL = Teuchos::parameterList();
        nonlinearSolverPL->set("Default Tol",double(maxError));
        _nonlinearSolver->setParameterList(nonlinearSolverPL);
        nonlinearSolver = _nonlinearSolver;
      }
      nonlinearSolver->setVerbLevel(M_outputLevelTeuchos);
      Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory
        = lowsfCreator.createLinearSolveStrategy("");
      Teuchos::RCP<Rythmos::RKButcherTableauBase<double> > rkbt = Rythmos::createRKBT<double>("Implicit 3 Stage 6th order Gauss");
      M_stepperPtr = Rythmos::implicitRKStepper<double>(model,nonlinearSolver,irk_W_factory,rkbt);
      method = "Implicit RK";
    } else {
      TEST_FOR_EXCEPT(true);
    }

    M_stepperPtr->setInitialCondition(model_ic);

    M_method = method;
   }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*M_out,success)
    return ;

}
void
ZeroDimensionalSolver::takeStep(double t0,double t1){

    M_finalTime = t1;
    M_startTime = t0;
    double dt = (t1-t0)/M_numberTimeStep;
    double time = t0;
    int numSteps = 0;
    Thyra::ModelEvaluatorBase::InArgs< double >  myIn = M_stepperPtr->getInitialCondition();
    myIn.describe(*M_out, M_outputLevelTeuchos);

    std::cout<<"numSteps = " <<numSteps << std::endl;
    if (M_step_method == STEP_METHOD_FIXED)
    {
      // Integrate forward with fixed step sizes:
      for (int i=1 ; i<=M_numberTimeStep ; ++i)
      {
        double dt_taken = M_stepperPtr->takeStep(dt,Rythmos::STEP_TYPE_FIXED);
        numSteps++;
        if (dt_taken != dt)
        {
          cerr << "Error, M_stepper took step of dt = " << dt_taken << " when asked to take step of dt = " << dt << std::endl;
          break;
        }
        time += dt_taken;
      }
    }
    else // M_step_method == STEP_METHOD_VARIABLE
    {
      while (time < M_finalTime)
      {
        double dt_taken = M_stepperPtr->takeStep(0.0,Rythmos::STEP_TYPE_VARIABLE);
        numSteps++;
        if (dt_taken < 0)
        {
          cerr << "Error, M_stepper failed for some reason with step taken = " << dt_taken << endl;
          break;
        }
        time += dt_taken;
        *M_out << "Took stepsize of: " << dt_taken << " time = " << time << endl;
      }
    }
    *M_out << "Integrated to time = " << time << endl;
    // Get solution M_out of M_stepper:
    Rythmos::TimeRange< double > timeRange =M_stepperPtr->getTimeRange();
    Rythmos::Array< double > time_vec;
    double timeToEvaluate = timeRange.upper();
    time_vec.push_back(timeToEvaluate);
    Rythmos::Array< Teuchos::RCP< const Thyra::VectorBase< double > > >  x_vec;
    Rythmos::Array< Teuchos::RCP< const Thyra::VectorBase< double > > >  xdot_vec;
    Rythmos::Array< Teuchos::ScalarTraits<double>::magnitudeType >  accuracy_vec;

    M_stepperPtr->getPoints(time_vec,
              &x_vec,
              &xdot_vec,
              &accuracy_vec
                );
    const Rythmos::StepStatus<double> stepStatus = M_stepperPtr->getStepStatus();
    Teuchos::RCP<const Thyra::VectorBase<double> > &x_computed_thyra_ptr = x_vec[0];
    Teuchos::RCP<const Thyra::VectorBase<double> > &x_dot_computed_thyra_ptr = xdot_vec[0];

    // Convert Thyra::VectorBase to Epetra_Vector
    Teuchos::RCP<const Epetra_Vector>   x_computed_ptr = Thyra::get_Epetra_Vector(*(M_solverInterfaceRCP->get_x_map()),x_computed_thyra_ptr);
    const Epetra_Vector &x_computed = *x_computed_ptr;

    Teuchos::RCP<const Epetra_Vector>   x_dot_computed_ptr = Thyra::get_Epetra_Vector(*(M_solverInterfaceRCP->get_x_map()),x_dot_computed_thyra_ptr);
    const Epetra_Vector &x_dot_computed = *x_dot_computed_ptr;


//---------------Replace the solution as initial solution for the next iteration
#ifdef HAVE_LIFEV_DEBUG
    *M_out << "X_computed" << endl;
    x_computed.Print(*M_out);
    *M_out << "X_dot_computed" << endl;
    x_dot_computed.Print(*M_out);
#endif
    M_modelInterface->initializeSolnY(x_computed);
    M_modelInterface->initializeSolnYp(x_dot_computed);
    M_modelInterface->deepUpdate(time,x_computed , x_dot_computed);
//---------------update the Circuit data from the the final solution
    int MyPID = x_computed.Comm().MyPID();
    if (MyPID == 0)
    {
      *M_out << "Integrating 1DfemTransient t = " << t0
                << " to t = " << t1 << std::endl;
      *M_out << "using " << M_method << "." << std::endl;
      *M_out << "Took " << numSteps << " steps." << std::endl;
    }
     *M_out << "\nEnd Result: TEST PASSED" << endl;

}
