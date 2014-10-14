/*
 * IfpackPreconditioner.cpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */

#include <Ifpack.h>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>

namespace LifeV {

namespace Operators {

IfpackPreconditioner::IfpackPreconditioner():
    RowMatrixPreconditioner()
{ }

IfpackPreconditioner::~IfpackPreconditioner() { }

int IfpackPreconditioner::myCompute()
{
    ASSERT_PRE(M_pList.isSublist("options"), "RowMatrixPreconditioner::SetParameterList(): pList must have a options subList");

    Ifpack FactoryIfpack;
    int returnCode(0);

    bool verbose(0 == M_rowMatrix->Comm().MyPID());
    std::string defaultPreconditioner("ILU");
    int defaultOverlap(1);

    if(!M_pList.isParameter("preconditioner") )
    {
        returnCode = 1;
        if(verbose)
            std::cout << "IfpackPreconditioner::myCompute: Using default preconditioner "
                        << defaultPreconditioner << "\n";

    }

    if(!M_pList.isParameter("overlap") )
    {
        returnCode = 1;
        if(verbose)
            std::cout << "IfpackPreconditioner::myCompute: Using default overlap " << defaultOverlap << "\n";
    }

    Ifpack_Preconditioner * prec;

    prec = FactoryIfpack.Create(M_pList.get("preconditioner", defaultPreconditioner), M_rowMatrix.get(),
                                        M_pList.get("overlap", defaultOverlap));

    prec->SetParameters(M_pList.sublist("options"));
    prec->Initialize();
    prec->Compute();

    M_prec.reset(prec);

    return returnCode;
}


} /* end Operators namespace */

} /* end LifeV namespace */
