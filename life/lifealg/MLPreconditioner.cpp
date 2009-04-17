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
        M_Prec(),
        M_precType(),
        M_analyze(false)
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
    bool found;
    M_analyze = dataFile((section + "/ML/analyze_smoother").data(), false, found);

    createMLList(dataFile, section, M_List);

    std::string CT;

    Teuchos::ParameterList& SmootherIFSubList = M_List.sublist("smoother: ifpack list");
    createIfpackList(dataFile, section, SmootherIFSubList);


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

    M_Prec.reset(new prec_raw_type(M_Oper->getEpetraMatrix(), this->getList(), true));

    if (M_analyze)
        {
            ML_Epetra::MultiLevelPreconditioner* prec;
            prec = dynamic_cast<ML_Epetra::MultiLevelPreconditioner*> (M_Prec.get());
            int NumPreCycles = 5;
            int NumPostCycles = 1;
            int NumMLCycles = 10;
            prec->AnalyzeHierarchy(true, NumPreCycles, NumPostCycles, NumMLCycles);

            //prec->TestSmoothers();
        }


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




    list.setName("ML paramters list");

    std::string defList      = dataFile((section + "/ML/default_parameter_list").data(), "SA");
    if (defList != "none")
        ML_Epetra::SetDefaults(defList, list);

    bool found;

    int MLPrintParameterList = dataFile((section + "/displayList").data(),      0, found);

    int MLOutput             = dataFile((section + "/ML/MLOuput").data(),       0, found);
    if (found) list.set("ML output",               MLOutput);

    int printUnused          = dataFile((section + "/ML/print_unused").data(), -2, found);
    if (found) list.set("print unused",            printUnused);

    int PDEEquations         = dataFile((section + "/ML/pde_equations").data(),    1, found);
    if (found) list.set("PDE equations",           PDEEquations);

    int CycleApplications    = dataFile((section + "/ML/cycle_applications").data(), 1, found);
    if (found) list.set("cycle applications", CycleApplications);

    int MaxLevels            = dataFile((section + "/ML/max_levels").data(),         2, found);
    if (found) list.set("max levels", MaxLevels);

    std::string IncOrDec     = dataFile((section + "/ML/inc_or_dec").data(),         "increasing", found);
    if (found) list.set("increasing or decreasing", IncOrDec);

    std::string PrecType     = dataFile((section + "/ML/prec_type").data(),          "MGV", found);
    if (found) list.set("prec type", PrecType);

    //    int NumProjectedModes    = dataFile((section + "/ML/number_of_prejected_modes").data(), 0);

    // if (found) list.set("ML print parameter list", MLPrintParameterList);

    std::string eigenAnalysisType = dataFile((section + "/ML/eigne-analysis/type").data(), "cg", found);
    if (found) list.set("eigen-analysis: type",       eigenAnalysisType);

    int eigenAnalysisIterations   = dataFile((section + "/ML/eigne-analysis/iterations").data(), 10, found);
    if (found) list.set("eigen-analysis: iterations", eigenAnalysisIterations);

    // Aggregation options

    std::string AggregationType                  = dataFile((section + "/ML/aggregation/type").data(), "Uncoupled", found);
    if (found) list.set("aggregation: type",                              AggregationType);

    double      AggregationThreshold             = dataFile((section + "/ML/aggregation/threshold").data(), 0.0 , found);
    if (found) list.set("aggregation: threshold",                         AggregationThreshold);

    double      AggregationDampingFactor         = dataFile((section + "/ML/aggregation/damping_factor").data(), 4./3. , found);
    if (found) list.set("aggregation: damping factor",                    AggregationDampingFactor);

    int AggregationSmoothingSweeps               = dataFile((section + "/ML/aggregation/smoothing_sweeps").data(),    1, found);
    if (found) list.set("aggregation: smoothing sweeps",                 AggregationSmoothingSweeps);

    int AggregationGlobalAggregates              = dataFile((section + "/ML/aggregation/global_aggregates").data(),   1, found);
    if (found) list.set("aggregation: global aggregates",                 AggregationGlobalAggregates);

    int AggregationLocalAggregates               = dataFile((section + "/ML/aggregation/local_aggregates").data(),    1, found);
    if (found) list.set("aggregation: local aggregates",                  AggregationLocalAggregates);

    int AggregationNodesPerAggregate             = dataFile((section + "/ML/aggregation/nodes_per_aggregate").data(), 1, found);
    if (found) list.set("aggregation: nodes per aggregate",               AggregationNodesPerAggregate);

    int AggregationNextLevelAggregatesPerProcess = dataFile((section + "/ML/aggregation/next-level_aggregates_per_process").data(), 128, found);
    if (found) list.set("aggregation: next level aggregates per process", AggregationNextLevelAggregatesPerProcess);

    bool AggregationUseTentativeRestriction      = dataFile((section + "/ML/aggregation/tentative_restriction").data(), false, found);
    if (found) list.set("aggregation: use tentative restriction",         AggregationUseTentativeRestriction);

    bool AggregationSymmetrize                   = dataFile((section + "/ML/aggregation/symmetrize").data(),            false, found);
    if (found) list.set("aggregation: symmetrize",                        AggregationSymmetrize);

    bool   EnergyMinimizationEnable              = dataFile((section + "/ML/energy_minimization/enable").data(),        false, found);
    if (found) list.set("energy minimization: enable",                    EnergyMinimizationEnable);

    int    EnergyMinimizationType                = dataFile((section + "/ML/energy_minimization/type").data(),          2, found);
    if (found) list.set("energy minimization: type",                      EnergyMinimizationType);

    double EnergyMinimizationDropTol             = dataFile((section + "/ML/energy_minimization/droptol").data(),       0., found);
    if (found) list.set("energy minimization: droptol",                   EnergyMinimizationDropTol);

    bool   EnergyMinimizationCheap               = dataFile((section + "/ML/energy_minimization/cheap").data(),         false, found);
    if (found) list.set("energy minimization: cheap",                     EnergyMinimizationCheap);

    // Smoothing parameters


    std::string SmootherType                = dataFile((section + "/ML/smoother/type").data(), "IFPACK", found);
    if (found) list.set("smoother: type",                         SmootherType);

    int SmootherSweeps                      = dataFile((section + "/ML/smoother/sweeps").data(), 2, found);
    if (found) list.set("smoother: sweeps",                       SmootherSweeps);

    double SmootherDampingFactor            = dataFile((section + "/ML/smoother/damping_factor").data(), 1.0, found);
    if (found) list.set("smoother: damping factor",               SmootherDampingFactor);

    std::string SmootherPreOrPost           = dataFile((section + "/ML/smoother/pre_or_post").data(), "both", found);
    if (found) list.set("smoother: pre or post",                  SmootherPreOrPost);

    double SmootherChebyshevAlpha           = dataFile((section + "/ML/smoother/Chebyshev_alpha").data(), 20., found);
    if (found) list.set("smoother: Chebyshev alpha",              SmootherChebyshevAlpha);

    bool SmootherHiptmairEfficientSymmetric = dataFile((section + "/ML/smoother/Hiptmair_efficient_symmetric").data(), true, found);
    if (found) list.set("smoother: Hiptmair efficient symmetric", SmootherHiptmairEfficientSymmetric);


    // subsmoother parameter

    std::string SubSmootherType             = dataFile((section + "/ML/subsmoother/type").data(), "Chebyshev", found);
    if (found) list.set("subsmoother: type",                      SubSmootherType);

    double SubSmootherChebyshevAlpha        = dataFile((section + "/ML/subsmoother/Chebyshev_alpha").data(), 20., found);
    if (found) list.set("subsmoother: Chebyshev alpha",           SubSmootherChebyshevAlpha);

    //    double SubSmootherSGSDampingFactor      = dataFile((section + "/ML/subsmoothers/SGS_damping_factor").data(), 1., found);
    int SubSmootherEdgeSweeps               = dataFile((section + "/ML/subsmoother/edge_sweeps").data(), 2, found);
    if (found) list.set("subsmoother: edge sweeps",               SubSmootherEdgeSweeps);

    int SubSmootherNodeSweeps               = dataFile((section + "/ML/subsmoother/node_sweeps").data(), 2, found);
    if (found) list.set("subsmoother: node sweeps",               SubSmootherNodeSweeps);




//     if (found) list.set("smoother: type (level 0)","IFPACK");
//     if (found) list.set("smoother: type (level 1)","IFPACK");

    //    if (found) list.set("subsmoother: SGS damping factor",        SubSmootherSGSDampingFactor);

    // Coarsest Grid Parameters

    int CoarseMaxSize                 = dataFile((section + "/ML/coarse/max_size").data(), 128, found);
    if (found) list.set("coarse: max size",  CoarseMaxSize);

    std::string CoarseType            = dataFile((section + "/ML/coarse/type").data(), "Chebyshev", found);
    if (found) list.set("coarse: type", CoarseType);

    std::string CoarsePreOrPost       = dataFile((section + "/ML/coarse/pre_or_post").data(), "post", found);
    if (found) list.set("coarse: pre or post", CoarsePreOrPost);

    int CoarseSweeps                  = dataFile((section + "/ML/coarse/sweeps").data(), 2, found);
    if (found) list.set("coarse: sweeps", CoarseSweeps);

    double CoarseDampingFactor        = dataFile((section + "/ML/coarse/damping_factor").data(), 1.0, found);
    if (found) list.set("coarse: damping factor", CoarseDampingFactor);

    std::string CoarseSubsmootherType = dataFile((section + "/ML/coarse/subsmoother_type").data(), "Chebyshev", found);
    if (found) list.set("coarse: subsmoother type", CoarseSubsmootherType);

    int CoarseNodeSweeps              = dataFile((section + "/ML/coarse/node_sweeps").data(), 2, found);
    if (found) list.set("coarse: node sweeps", CoarseNodeSweeps);

    int CoarseEdgeSweeps              = dataFile((section + "/ML/coarse/edge_sweeps").data(), 2, found);
    if (found) list.set("coarse: edge sweeps", CoarseEdgeSweeps);

    double CoarseChebyshevAlpha       = dataFile((section + "/ML/coarse/Chebyshev_alpha").data(), 30., found);
    if (found) list.set("coarse: Chebyshev alpha", CoarseChebyshevAlpha);

    int CoarseMaxProcesses            = dataFile((section + "/ML/coarse/max_processes").data(),-1, found);
    if (found) list.set("coarse: max processes", CoarseMaxProcesses);


    // Load-balancing Options

    int RepartitionEnable              = dataFile((section + "/ML/repartition/enable").data(), 0, found);
    if (found) list.set("repartition: enable",             RepartitionEnable);

    std::string RepartitionPartitioner = dataFile((section + "/ML/repartition/partitioner").data(), "ParMETIS", found);
    if (found) list.set("repartition: partitioner",        RepartitionPartitioner);

    double RepartitionMaxMinRatio      = dataFile((section + "/ML/repartition/max_min_ratio").data(), 1.3, found);
    if (found) list.set("repartition: max min ratio",      RepartitionMaxMinRatio);

    int RepartitionMinPerProc          = dataFile((section + "/ML/repartition/min_per_proc").data(), 512, found);
    if (found) list.set("repartition: min per proc",       RepartitionMinPerProc);

    double RepartitionNodeMaxMinRatio  = dataFile((section + "/ML/repartition/node_max_min_ratio").data(), 1.3, found);
    if (found) list.set("repartition: node max min ratio", RepartitionNodeMaxMinRatio);

    int RepartitionNodeMinPerProc      = dataFile((section + "/ML/repartition/node_min_per_proc").data(), 170, found);
    if (found) list.set("repartition: node min per proc",  RepartitionNodeMinPerProc);

    int RepartitionZoltanDimensions    = dataFile((section + "/ML/repartition/Zoltan_dimensions").data(), 2, found);
    if (found) list.set("repartition: Zoltan dimensions",  RepartitionZoltanDimensions);


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



    //list.sublist("smoother: ifpack list").set("schwarz: filter singletons", true);

    //    list.sublist("smoother: ifpack list").set("amesos: solver type", "Amesos_Lapack");

    if (MLPrintParameterList) list.print(std::cout);
}


} // namespace LifeV
//#endif //#ifdef HAVE_TRILINOS_ML
