/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-11-09

  Copyright (C) 2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file EpetraPreconditioner.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2006-11-09
 */

#include "MLPreconditioner.hpp"


namespace LifeV
{

MLPreconditioner::MLPreconditioner():
        super(),
        M_Prec()
{}

MLPreconditioner::~MLPreconditioner()
{}

// MLPreconditioner::MLPreconditioner(operator_type& oper):
//   M_Prec(),
//   M_Oper()
// {
//     buildPreconditioner(oper);
// }




void
MLPreconditioner::setDataFromGetPot( const GetPot& dataFile, const std::string& section )
{


    std::string defaultParameters = dataFile((section + "/default_parameters").data(), "SA");
    bool displayList = dataFile((section + "/displayList").data(),     false);

    ML_Epetra::SetDefaults(defaultParameters, M_List);





    // use Uncoupled scheme to create the aggregate
    //    M_List.set("aggregation: type", "Uncoupled");

    // fix the smoother
    //    M_List.set("smoother: type","symmetric Gauss-Seidel");

    //    M_List.set("output", 10);
    //    M_List.set("print unused", 1);

//     M_List.set("smoother: pre or post", "both");
//     M_List.set("PDE equations", 1);

    // fix the smoother to be IFPACK; can be set using (level X) syntax
    M_List.set("smoother: type", "Aztec");

    // now we have to specify which IFPACK preconditioner should be
    // built. Any value that is valid for the IFPACK factory. We also need
    // to define the overlap (>= 0).

    //M_List.set("smoother: type (level 0)","IFPACK");

    //    M_List.set("smoother: ifpack type", "Aztec");
//     M_List.set("smoother: ifpack overlap", M_overlapLevel);

    M_List.set("smoother: type (level 0)","IFPACK");
    M_List.set("smoother: type (level 1)","IFPACK");

    //extra parameters:
//     M_List.sublist("smoother: ifpack list").set("fact: drop tolerance",
//                                                 dropTolerance);
//     M_List.sublist("smoother: ifpack list").set("fact: ilut level-of-fill",
//                                                 levelOfFill);
//     M_List.sublist("smoother: ifpack list").set("fact: absolute threshold",
//                                                 athr);
//     M_List.sublist("smoother: ifpack list").set("fact: relative threshold",
//                                                 rthr);

//     M_List.sublist("smoother: ifpack list").set("partitioner: type",
//                                                 "user");

//     M_List.sublist("smoother: ifpack list").set("relaxation: type",
//                                                 "symmetric Gauss-Seidel");
//     M_List.sublist("smoother: ifpack list").set("relaxation: damping factor",
//                                                 omega);
//     M_List.sublist("smoother: ifpack list").set("relaxation: sweeps",
//                                                 1);
//     M_List.sublist("smoother: ifpack list").set("relaxation: zero starting solution",
//                                                 false);

    if (displayList) M_List.print(std::cout);

    //M_List.sublist("smoother: ifpack list").set("schwarz: filter singletons", true);

    //    M_List.sublist("smoother: ifpack list").set("amesos: solver type", "Amesos_Lapack");



}

int MLPreconditioner::buildPreconditioner(operator_type& oper)
{

    M_Oper = oper;

    /*
    //List.set("schwarz: combine mode", "Zero"); //
    M_List.set("schwarz: filter singletons", true);

//    List.set("amesos: solver type", "Amesos_Lapack");

//     M_List.set("PrintTiming", true);
//     M_List.set("PrintStatus", true);


    Ifpack factory;

    M_Prec.reset(factory.Create(M_precType, &M_Oper->getEpetraMatrix(), M_overlapLevel));
//    M_Prec.reset(new prec_type(&A.getEpetraMatrix(), OverlapLevel));
    if ( !M_Prec.get() )
    { //! if not filled, I do not know how to diagonalize.
      ERROR_MSG( "Preconditioner not setted, something went wrong in its computation\n" );
    }

    IFPACK_CHK_ERR(M_Prec->SetParameters(M_List));
    IFPACK_CHK_ERR(M_Prec->Initialize());
    IFPACK_CHK_ERR(M_Prec->Compute());
    */

    M_Prec.reset(new prec_raw_type(M_Oper->getEpetraMatrix(), M_List, true));

    //    ML_CHK_ERR(M_Prec->SetParameters(M_List));

//     M_Prec->Initialize();
//     M_Prec->Compute();

    return EXIT_SUCCESS;



//     Teuchos::ParameterList MLTestList;
//     ML_Epetra::SetDefaults("SA",MLTestList);
    //    MLTestList.set("test: Jacobi", false);
//     MLTestList.set("test: Gauss-Seidel", false);
//     MLTestList.set("test: block Gauss-Seidel", false);
//     MLTestList.set("test: symmetric Gauss-Seidel", false);
//     MLTestList.set("test: ParaSails", false);
//     MLTestList.set("test: IFPACK", false);
//     MLTestList.set("test: Aztec", false);
//     MLTestList.set("test: ML", false);
//     int NumPreCycles = 5;
//     int NumPostCycles = 1;
//     int NumMLCycles = 10;
//     M_Prec->AnalyzeHierarchy(true, NumPreCycles, NumPostCycles, NumMLCycles);
//    M_Prec->TestSmoothers(MLTestList);
}

double
MLPreconditioner::Condest()
{
    return 0.;
 }

EpetraPreconditioner::prec_raw_type*
MLPreconditioner::getPrec()
{
    return M_Prec.get();
}

void
MLPreconditioner::precReset()
{
    M_Oper.reset();
    M_Prec.reset();
}




} // namespace LifeV
