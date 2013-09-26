/*
 * MlPreconditioner.cpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */

#include <lifev/core/LifeV.hpp>
#include <lifev/operator/linear_algebra/MLPreconditioner.hpp>

// Tell the compiler to ignore specific kind of warnings
LIFEV_SUPPRESS_WARNINGS

#include <ml_MultiLevelPreconditioner.h>
#include <ml_epetra.h>

// Tell the compiler to restore the warnings
LIFEV_RESTORE_WARNINGS

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
