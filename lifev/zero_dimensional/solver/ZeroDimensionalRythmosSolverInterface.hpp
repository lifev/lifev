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
 *  @version alpha (experimental)
 *
 *  @date 21-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalRythmosSolverInterface_H
#define ZeroDimensionalRythmosSolverInterface_H 1

#include <lifev/zero_dimensional/solver/ZeroDimensionalRythmosModelInterface.hpp>

#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
namespace LifeV
{
//! Rythmos solver interface.
/*!
 * This class will communicate with with Rythmos solver and model interface.
 */
class RythmosSolverInterface : public EpetraExt::ModelEvaluator
{
public:

    // Constructor
    RythmosSolverInterface (Int numCircuitElements,
                            Teuchos::RCP<Epetra_Comm>& epetra_comm_ptr,
                            rythmosModelInterfacePtrRCP_Type theModel);

    // Initialization
    void initialize();

    Teuchos::RCP<const Epetra_Map> get_x_map() const;

    Teuchos::RCP<const Epetra_Map> get_f_map() const;

    Teuchos::RCP<const Epetra_Vector> get_x_init() const;

    Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;

    Teuchos::RCP<Epetra_Operator> create_W() const;

    InArgs createInArgs() const;

    OutArgs createOutArgs() const;

    void evalModel ( const InArgs& inArgs, const OutArgs& outArgs ) const;

private:

    Teuchos::RCP<Epetra_Comm>            M_epetraCommPtr;
    Teuchos::RCP<const Epetra_Map>       M_epetraMapPtr;
    Int                                  M_numElements;
    Teuchos::RCP<Epetra_CrsGraph>        M_Wgraph;
    rythmosModelInterfacePtrRCP_Type     M_problemInterfacePtr;
    Teuchos::RCP<Epetra_Comm>            M_comm;
};

typedef boost::shared_ptr< RythmosSolverInterface > rythmosSolverInterfacePtr_Type;
typedef Teuchos::RCP< RythmosSolverInterface > rythmosSolverInterfacePtrRCP_Type;

} // LifeV namespace

#endif /* HAVE_NOX_THYRA && HAVE_TRILINOS_RYTHMOS */

#endif //ZeroDimensionalRythmosSolverInterface_H
