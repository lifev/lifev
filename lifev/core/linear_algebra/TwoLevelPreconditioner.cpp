/*
 * TwoLevelPreconditioner.cpp
 *
 *  Created on: Oct 12, 2011
 *      Author: uvilla
 */

#include "TwoLevelPreconditioner.hpp"

#include "TwoLevelOperator.hpp"
#include "ApproximatedInvertibleRowMatrix.hpp"

#include <Trilinos_version.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <EpetraExt_MatrixMatrix.h>

namespace LifeV
{

namespace Operators
{

TwoLevelPreconditioner::TwoLevelPreconditioner():
		RowMatrixPreconditioner()
{ }

TwoLevelPreconditioner::~TwoLevelPreconditioner()
{ }

int TwoLevelPreconditioner::myCompute()
{
	ASSERT_PRE(M_pList.isParameter("EstensionMatrix"), "[TwoLevelPreconditioner::myCompute] pList should contain an EstensionMatrix");
	ASSERT_PRE(M_pList.isSublist("CoarseLevel"), "[TwoLevelPreconditioner::myCompute] pList should contain a CoarseLevel sublist");
	ASSERT_PRE(M_pList.isSublist("FineLevel"), "[TwoLevelPreconditioner::myCompute] pList should contain a FineLevel sublist");
	ASSERT_PRE(M_rowMatrix.get() != 0, "[TwoLevelPreconditioner::myCompute] You need to SetRowMatrix first");

	rowMatrixPtr_Type E, R;
	E = M_pList.get<rowMatrixPtr_Type>("EstensionMatrix");

	if(M_pList.isParameter("RestrictionMatrix"))
	{
		R = M_pList.get<rowMatrixPtr_Type>("RestrictionMatrix");
	}
	else
	{
		bool MakeDataContiguous = true;
#if TRILINOS_MAJOR_VERSION < 11
	EpetraExt::RowMatrix_Transpose transposer ( MakeDataContiguous );
#else
	EpetraExt::RowMatrix_Transpose transposer ( 0, !MakeDataContiguous );
#endif
		R.reset( new rowMatrix_Type(dynamic_cast<rowMatrix_Type &>(transposer(*E))));
		M_pList.set("RestricitionMatrix", R);
	}

	rowMatrixPtr_Type AE(new rowMatrix_Type(Copy, E->RangeMap(), 0) );
	rowMatrixPtr_Type RAE(new rowMatrix_Type(Copy, E->DomainMap(), 0) );
    EPETRA_CHK_ERR(EpetraExt::MatrixMatrix::Multiply(*M_rowMatrix, false, *E, false, *AE));
    EPETRA_CHK_ERR(EpetraExt::MatrixMatrix::Multiply(*R,           false, *AE, false, *RAE));

    std::shared_ptr<RowMatrixPreconditioner> S(RowMatrixPreconditionerFactory::instance().createObject("Ifpack"));
    S->SetRowMatrix(M_rowMatrix);
    S->SetParameterList(M_pList.sublist("FineLevel"));
    EPETRA_CHK_ERR(S->Compute());

    std::shared_ptr<ApproximatedInvertibleRowMatrix> Ac(new ApproximatedInvertibleRowMatrix);
    Ac->SetRowMatrix(RAE);
    Ac->SetParameterList(M_pList.sublist("CoarseLevel"));
    EPETRA_CHK_ERR(Ac->Compute());

    TwoLevelOperator * prec(new TwoLevelOperator);
    prec->SetFineLevelOperator(M_rowMatrix);
    prec->SetSmootherOperator(S);
    prec->SetRestrictionOperator(R);
    prec->SetEstensionOperator(E);
    prec->SetCoarseLevelOperator(Ac);

    M_prec.reset(prec);

    return 0;
}

}

}
