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
 *  @brief Rythmos solver Interface.
 *
 *  @date 21-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/zero_dimensional/solver/ZeroDimensionalRythmosSolverInterface.hpp>

namespace LifeV
{

#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
// ===================================================
// Constructors
// ===================================================
RythmosSolverInterface::RythmosSolverInterface ( Int numCircuitElements,
                                                 Teuchos::RCP< Epetra_Comm >& epetra_comm_ptr,
                                                 rythmosModelInterfacePtrRCP_Type theModel ) :
    M_epetraCommPtr ( epetra_comm_ptr ), M_numElements ( numCircuitElements ), M_problemInterfacePtr ( theModel ), M_comm ( epetra_comm_ptr )
{
    initialize();
}

void RythmosSolverInterface::initialize()
{
    M_epetraMapPtr = Teuchos::rcp ( new Epetra_Map ( M_problemInterfacePtr->getMap() ) );
    M_Wgraph = Teuchos::rcp ( new Epetra_CrsGraph ( M_problemInterfacePtr->getGraph() ) );
}

// Overridden from EpetraExt::ModelEvaluator
Teuchos::RCP< const Epetra_Map > RythmosSolverInterface::get_x_map() const
{
    return M_epetraMapPtr;
}

Teuchos::RCP< const Epetra_Map > RythmosSolverInterface::get_f_map() const
{
    return M_epetraMapPtr;
}

Teuchos::RCP< const Epetra_Vector > RythmosSolverInterface::get_x_init() const
{
    Epetra_Vector& soln = M_problemInterfacePtr->getSolutionY();
    Teuchos::RCP< Epetra_Vector > x_init = Teuchos::rcp ( new Epetra_Vector ( soln ) );
    return x_init;
}

Teuchos::RCP< const Epetra_Vector > RythmosSolverInterface::get_x_dot_init() const
{
    Epetra_Vector& soln = M_problemInterfacePtr->getSolutionYp();
    Teuchos::RCP< Epetra_Vector > x_dot_init = Teuchos::rcp ( new Epetra_Vector ( soln ) );
    return x_dot_init;
}

Teuchos::RCP< Epetra_Operator > RythmosSolverInterface::create_W() const
{
    Teuchos::RCP< Epetra_Operator > W = Teuchos::rcp ( new Epetra_CrsMatrix ( ::Copy,
                                                                              *M_Wgraph ) );
    return W;
}

EpetraExt::ModelEvaluator::InArgs RythmosSolverInterface::createInArgs() const
{
    InArgsSetup inArgs;
    inArgs.setSupports ( IN_ARG_x, true );
    inArgs.setSupports ( IN_ARG_x_dot, true );
    inArgs.setSupports ( IN_ARG_alpha, true );
    inArgs.setSupports ( IN_ARG_beta, true );
    inArgs.setSupports ( IN_ARG_t, true );
    return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs RythmosSolverInterface::createOutArgs() const
{
    OutArgsSetup outArgs;
    outArgs.setSupports ( OUT_ARG_f, true );
    outArgs.setSupports ( OUT_ARG_W, true );
    outArgs.setSupports ( OUT_ARG_W, true );
    outArgs.set_W_properties ( DerivativeProperties ( DERIV_LINEARITY_NONCONST, DERIV_RANK_UNKNOWN, true ) );
    return outArgs;
}

void RythmosSolverInterface::evalModel ( const InArgs& inArgs,
                                         const OutArgs& outArgs ) const
{
    Teuchos::RCP< const Epetra_Vector > x = inArgs.get_x();
    Teuchos::RCP< const Epetra_Vector > xdot = inArgs.get_x_dot();
    Real t = inArgs.get_t();

#ifdef HAVE_LIFEV_DEBUG
    std::cout << "RythmosSolverInterface::evalModel ---------------------------{" << std::endl;
    std::cout << "x = " << std::endl;
    x->Print (std::cout);
    std::cout << "xdot = " << std::endl;
    xdot->Print (std::cout);
#endif // HAVE_LIFEV_DEBUG

    Teuchos::RCP< Epetra_Vector > f;
    if ( ( f = outArgs.get_f() ).get() )
    {
        M_problemInterfacePtr->evaluateFImplicit ( t, &*x, &*xdot, &*f );
#ifdef HAVE_LIFEV_DEBUG
        std::cout << "f = " << std::endl;
        f->Print (std::cout);
#endif // HAVE_LIFEV_DEBUG

    }
    Teuchos::RCP< Epetra_Operator > W;
    if ( ( W = outArgs.get_W() ).get() )
    {
        const Real alpha = inArgs.get_alpha();
        const Real beta = inArgs.get_beta();
        Epetra_CrsMatrix& jac = Teuchos::dyn_cast< Epetra_CrsMatrix > ( *W );
        M_problemInterfacePtr->evaluateWImplicit ( t, alpha, beta, &*x, &*xdot, &jac );
#ifdef HAVE_LIFEV_DEBUG
        std::cout << "jac = " << std::endl;
        jac.Print (std::cout);
#endif // HAVE_LIFEV_DEBUG
    }
#ifdef HAVE_LIFEV_DEBUG
    std::cout << "RythmosSolverInterface::evalModel ---------------------------}" << std::endl;
#endif // HAVE_LIFEV_DEBUG
}

#endif /* HAVE_NOX_THYRA && HAVE_TRILINOS_RYTHMOS */

} // LifeV namespace
