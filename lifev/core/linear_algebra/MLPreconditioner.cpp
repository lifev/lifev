/*
 * MlPreconditioner.cpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */

#include "MLPreconditioner.hpp"

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <ml_MultiLevelPreconditioner.h>
#include <ml_epetra.h>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"


namespace LifeV {

namespace Operators {

MLPreconditioner::MLPreconditioner() :
		RowMatrixPreconditioner()
	{ }

MLPreconditioner::~MLPreconditioner()
{ }

int MLPreconditioner::myCompute()
{
	ASSERT_PRE(M_pList.isSublist("options"), "RowMatrixPreconditioner::SetParameterList(): pList must have a options subList");

	ML_Epetra::MultiLevelPreconditioner * prec;

	prec = new ML_Epetra::MultiLevelPreconditioner(*M_rowMatrix, false);
	prec->SetParameterList(M_pList.sublist("options"));

	ML_CHK_ERR(prec->ComputePreconditioner());

	M_prec.reset(prec);

	return 0;
}


} /*end Operators namespace */

} /*end LifeV namespace */
