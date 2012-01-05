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
    @file
    @brief AztecooOperator

    @author Umberto Villa <umberto.villa@gmail.com>

    @date 03-09-2010
 */

#include <life/lifeoperators/AztecooOperator.hpp>

namespace LifeV
{

namespace Operators
{

AztecooOperator::AztecooOperator():
		SolverOperator(),
		M_linSolver(new SolverType)
{
	M_name = "AztecooOperator";
}

int AztecooOperator::doApplyInverse(const vector_Type& X, vector_Type& Y) const
{
	vector_Type Xcopy(X);
	Y.PutScalar(0.0);
	M_linSolver->SetUserOperator( M_oper.get());
	M_linSolver->SetParameters(*M_pList);
	M_linSolver->SetRHS(&Xcopy);
	M_linSolver->SetLHS(&Y);
	if(M_prec.get() != 0)
		EPETRA_CHK_ERR(M_linSolver->SetPrecOperator( (Epetra_Operator*) M_prec.get()) );

	int maxIter( M_pList->get<int>("max_iter"));
	double tol(  M_pList->get<double>("tol"));

	return M_linSolver->Iterate(maxIter, tol);
}

} // Namespace Operators

} // Namespace LifeV
