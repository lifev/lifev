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
//#ifdef HAVE_TRILINOS_ML


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
MLPreconditioner::setDataFromGetPot( const GetPot&          dataFile,
                                     const std::string&     section)
{
    createMLList(dataFile, section, M_List);

}

int
MLPreconditioner::buildPreconditioner(operator_type& oper)
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



void
createMLList( const GetPot&              dataFile,
              const std::string&         section,
              Teuchos::ParameterList&    list)
{

    std::string defList      = dataFile((section + "/ML/default_parameter_list").data(), "SA");

    if (defList != "none")
        ML_Epetra::SetDefaults(defList, list);

    int MLOutput             = dataFile((section + "/ML/MLOuput").data(),       0);
    int printUnused          = dataFile((section + "/ML/print_unused").data(), -2);
    int MLPrintParameterList = dataFile((section + "/ML/displayList").data(),      0);
    int PDEEquations         = dataFile((section + "/ML/pde_equations").data(),    1);

    int CycleApplications    = dataFile((section + "/ML/cycle_applications").data(), 1);
    int MaxLevels            = dataFile((section + "/ML/max_levels").data(),         2);
    std::string IncOrDec     = dataFile((section + "/ML/inc_or_dec").data(),         "increasing");
    std::string PrecType     = dataFile((section + "/ML/prec_type").data(),          "MGV");
    //    int NumProjectedModes    = dataFile((section + "/ML/number_of_prejected_modes").data(), 0);

    list.set("ML output",               MLOutput);
    list.set("print unused",            printUnused);
    // list.set("ML print parameter list", MLPrintParameterList);
    list.set("PDE equations",           PDEEquations);


    list.set("cycle applications", CycleApplications);
    list.set("max levels", MaxLevels);
    list.set("increasing or decreasing", IncOrDec);
    list.set("prec type", PrecType);

    std::string eigenAnalysisType = dataFile((section + "/ML/eigne-analysis/type").data(), "cg");
    int eigenAnalysisIterations   = dataFile((section + "/ML/eigne-analysis/iterations").data(), 10);

    list.set("eigen-analysis: type",       eigenAnalysisType);
    list.set("eigen-analysis: iterations", eigenAnalysisIterations);

    // Aggregation options

    std::string AggregationType                  = dataFile((section + "/ML/aggregation/type").data(), "Uncoupled");
    double      AggregationThreshold             = dataFile((section + "/ML/aggregation/threshold").data(), 0.0 );
    double      AggregationDampingFactor         = dataFile((section + "/ML/aggregation/damping_factor").data(), 4./3. );

    int AggregationSmoothingSweeps               = dataFile((section + "/ML/aggregation/smoothing_sweeps").data(),    1);
    int AggregationGlobalAggregates              = dataFile((section + "/ML/aggregation/global_aggregates").data(),   1);
    int AggregationLocalAggregates               = dataFile((section + "/ML/aggregation/local_aggregates").data(),    1);
    int AggregationNodesPerAggregate             = dataFile((section + "/ML/aggregation/nodes_per_aggregate").data(), 1);
    int AggregationNextLevelAggregatesPerProcess = dataFile((section + "/ML/aggregation/next-level_aggregates_per_process").data(), 128);

    bool AggregationUseTentativeRestriction      = dataFile((section + "/ML/aggregation/tentative_restriction").data(), false);
    bool AggregationSymmetrize                   = dataFile((section + "/ML/aggregation/symmetrize").data(),            false);

    list.set("aggregation: type",                              AggregationType);
    list.set("aggregation: threshold",                         AggregationThreshold);
    list.set("aggregation: damping factor",                    AggregationDampingFactor);
    list.set("aggregation: smoothing sweeps",                  AggregationSmoothingSweeps);
    list.set("aggregation: global aggregates",                 AggregationGlobalAggregates);
    list.set("aggregation: local aggregates",                  AggregationLocalAggregates);
    list.set("aggregation: nodes per aggregate",               AggregationNodesPerAggregate);
    list.set("aggregation: use tentative restriction",         AggregationUseTentativeRestriction);
    list.set("aggregation: symmetrize",                        AggregationSymmetrize);

    bool   EnergyMinimizationEnable              = dataFile((section + "/ML/energy_minimization/enable").data(),        false);
    int    EnergyMinimizationType                = dataFile((section + "/ML/energy_minimization/type").data(),          2);
    double EnergyMinimizationDropTol             = dataFile((section + "/ML/energy_minimization/droptol").data(),       0.);
    bool   EnergyMinimizationCheap               = dataFile((section + "/ML/energy_minimization/cheap").data(),         false);

    list.set("energy minimization: enable",                    EnergyMinimizationEnable);
    list.set("energy minimization: type",                      EnergyMinimizationType);
    list.set("energy minimization: droptol",                   EnergyMinimizationDropTol);
    list.set("energy minimization: cheap",                     EnergyMinimizationCheap);

    // Smoothing parameters

    std::string SmootherType                = dataFile((section + "/ML/smoothers/type").data(), "IFPACK");
    int SmootherSweeps                      = dataFile((section + "/ML/smoothers/sweeps").data(), 2);
    double SmootherDampingFactor            = dataFile((section + "/ML/smoothers/damping_factor").data(), 1.0);
    std::string SmootherPreOrPost           = dataFile((section + "/ML/smoothers/pre_or_post").data(), "both");

    double SmootherChebyshevAlpha           = dataFile((section + "/ML/smoothers/Chebyshev_alpha").data(), 20.);
    bool SmootherHiptmairEfficientSymmetric = dataFile((section + "/ML/smoothers/Hiptmair_efficient_symmetric").data(), true);

    std::string SubSmootherType             = dataFile((section + "/ML/subsmoothers/type").data(), "Chebyshev");
    double SubSmootherChebyshevAlpha        = dataFile((section + "/ML/subsmoothers/Chebyshev_alpha").data(), 20.);
    //    double SubSmootherSGSDampingFactor      = dataFile((section + "/ML/subsmoothers/SGS_damping_factor").data(), 1.);
    int SubSmootherEdgeSweeps               = dataFile((section + "/ML/subsmoothers/edge_sweeps").data(), 2);
    int SubSmootherNodeSweeps               = dataFile((section + "/ML/subsmoothers/node_sweeps").data(), 2);

    list.set("smoother: type",                         SmootherType);
    list.set("smoother: sweeps",                       SmootherSweeps);
    list.set("smoother: damping factor",               SmootherDampingFactor);
    list.set("smoother: pre or post",                  SmootherPreOrPost);
    list.set("smoother: Chebyshev alpha",              SmootherChebyshevAlpha);
    list.set("smoother: Hiptmair efficient symmetric", SmootherHiptmairEfficientSymmetric);

    list.set("subsmoother: type",                      SubSmootherType);
    list.set("subsmoother: Chebyshev alpha",           SubSmootherChebyshevAlpha);
    //    list.set("subsmoother: SGS damping factor",        SubSmootherSGSDampingFactor);
    list.set("subsmoother: edge sweeps",               SubSmootherEdgeSweeps);
    list.set("subsmoother: node sweeps",               SubSmootherNodeSweeps);

    // Coarsest Grid Parameters

    int CoarseMaxSize                 = dataFile((section + "/ML/coarse/max_size").data(), 128);
    std::string CoarseType            = dataFile((section + "/ML/coarse/type").data(), "Chebyshev");
    std::string CoarsePreOrPost       = dataFile((section + "/ML/coarse/pre_or_post").data(), "post");
    double CoarseDampingFactor        = dataFile((section + "/ML/coarse/damping_factor").data(), 1.0);
    std::string CoarseSubsmootherType = dataFile((section + "/ML/coarse/subsmoother_type").data(), "Chebyshev");
    int CoarseNodeSweeps              = dataFile((section + "/ML/coarse/node_sweeps").data(), 2);
    int CoarseEdgeSweeps              = dataFile((section + "/ML/coarse/edge_sweeps").data(), 2);
    double CoarseChebyshevAlpha       = dataFile((section + "/ML/coarse/Chebyshev_alpha").data(), 30.);
    int CoarseMaxProcesses            = dataFile((section + "/ML/coarse/max_processes").data(),-1);

    list.set("coarse: max size",  CoarseMaxSize);
    //    list.set("coarse: type ", CoarseType);
    list.set("coarse: pre or post", CoarsePreOrPost);
    list.set("coarse: damping factor", CoarseDampingFactor);
    list.set("coarse: subsmoother type", CoarseSubsmootherType);
    list.set("coarse: node sweeps", CoarseNodeSweeps);
    list.set("coarse: edge sweeps", CoarseEdgeSweeps);
    list.set("coarse: Chebyshev alpha", CoarseChebyshevAlpha);
    list.set("coarse: max processes", CoarseMaxProcesses);
    list.set("coarse: max size", CoarseMaxSize);

    // Load-balancing Options

    int RepartitionEnable              = dataFile((section + "/ML/repartition/enable").data(), 0);
    std::string RepartitionPartitioner = dataFile((section + "/ML/repartition/partitioner").data(), "ParMETIS");
    double RepartitionMaxMinRatio      = dataFile((section + "/ML/repartition/max_min_ratio").data(), 1.3);
    int RepartitionMinPerProc          = dataFile((section + "/ML/repartition/min_per_proc").data(), 512);
    double RepartitionNodeMaxMinRatio  = dataFile((section + "/ML/repartition/node_max_min_ratio").data(), 1.3);
    int RepartitionNodeMinPerProc      = dataFile((section + "/ML/repartition/node_min_per_proc").data(), 170);
    int RepartitionZoltanDimensions    = dataFile((section + "/ML/repartition/Zoltan_dimensions").data(), 2);


    list.set("repartition: enable",             RepartitionEnable);
    list.set("repartition: partitioner",        RepartitionPartitioner);
    list.set("repartition: max min ratio",      RepartitionMaxMinRatio);
    list.set("repartition: min per proc",       RepartitionMinPerProc);
    list.set("repartition: node max min ratio", RepartitionNodeMaxMinRatio);
    list.set("repartition: node min per proc",  RepartitionNodeMinPerProc);
    list.set("repartition: Zoltan dimensions",  RepartitionZoltanDimensions);

    // use Uncoupled scheme to create the aggregate
    //    list.set("aggregation: type", "Uncoupled");

    // fix the smoother
    //    list.set("smoother: type","symmetric Gauss-Seidel");

    //    list.set("output", 10);
    //    list.set("print unused", 1);

//     list.set("smoother: pre or post", "both");
//     list.set("PDE equations", 1);

    // fix the smoother to be IFPACK; can be set using (level X) syntax
    //    list.set("smoother: type", "Aztec");

    // now we have to specify which IFPACK preconditioner should be
    // built. Any value that is valid for the IFPACK factory. We also need
    // to define the overlap (>= 0).

    //list.set("smoother: type (level 0)","IFPACK");

    //    list.set("smoother: ifpack type", "Aztec");
//     list.set("smoother: ifpack overlap", M_overlapLevel);

//     list.set("smoother: type (level 0)","IFPACK");
//     list.set("smoother: type (level 1)","IFPACK");

    if (SmootherType == "IFPACK")
        {
            Teuchos::ParameterList& IFPACKList = list.sublist("smoother: ifpack list");

            createIfpackList(dataFile, section, IFPACKList);
        }

    //extra parameters:
//     list.sublist("smoother: ifpack list").set("fact: drop tolerance",
//                                                 dropTolerance);
//     list.sublist("smoother: ifpack list").set("fact: ilut level-of-fill",
//                                                 levelOfFill);
//     list.sublist("smoother: ifpack list").set("fact: absolute threshold",
//                                                 athr);
//     list.sublist("smoother: ifpack list").set("fact: relative threshold",
//                                                 rthr);

//     list.sublist("smoother: ifpack list").set("partitioner: type",
//                                                 "user");

//     list.sublist("smoother: ifpack list").set("relaxation: type",
//                                                 "symmetric Gauss-Seidel");
//     list.sublist("smoother: ifpack list").set("relaxation: damping factor",
//                                                 omega);
//     list.sublist("smoother: ifpack list").set("relaxation: sweeps",
//                                                 1);
//     list.sublist("smoother: ifpack list").set("relaxation: zero starting solution",
//                                                 false);


    if (MLPrintParameterList) list.print(std::cout);

    //list.sublist("smoother: ifpack list").set("schwarz: filter singletons", true);

    //    list.sublist("smoother: ifpack list").set("amesos: solver type", "Amesos_Lapack");

}


} // namespace LifeV
//#endif //#ifdef HAVE_TRILINOS_ML
