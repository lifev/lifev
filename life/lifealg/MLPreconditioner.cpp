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
    @brief ML preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-0006
 */

#include <lifeconfig.h>
#include "MLPreconditioner.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MLPreconditioner::MLPreconditioner():
        super(),
        M_Oper(),
        M_Prec(),
        M_analyze(false)
{}

MLPreconditioner::~MLPreconditioner()
{}


// ===================================================
// Methods
// ===================================================
int
MLPreconditioner::buildPreconditioner(operator_type& oper)
{

    //the Trilinos::MultiLevelPreconditioner unsafely access to the area of memory co-owned by M_Oper.
    //to avoid the risk of dandling pointers always deallocate M_Prec first and then M_Oper
    M_Prec.reset();
    M_Oper = oper->getMatrixPtr();

    M_precType = M_List.get("prec type", "undefined??");
    M_precType += "_ML";

    // <one-level-postsmoothing> / <two-level-additive>
    // <two-level-hybrid> / <two-level-hybrid2>

    M_Prec.reset(new prec_raw_type(*M_Oper, this->getList(), true));

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

//     ML_CHK_ERR(M_Prec->SetParameters(M_List));
//     M_Prec->Initialize();
//     M_Prec->Compute();

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

void
MLPreconditioner::precReset()
{
    //the Trilinos::MultiLevelPreconditioner unsafely access to the area of memory co-owned by M_Oper.
    //to avoid the risk of dandling pointers always deallocate M_Prec first and then M_Oper

    M_Prec.reset();
    M_Oper.reset();

    this->M_preconditionerCreated = false;
}

void
MLPreconditioner::createMLList(       list_Type&    list,
                                      const GetPot&       dataFile,
                                      const std::string&  section,
                                      const std::string&  subSection )
{
    list.setName("ML paramters list");

    std::string defList      = dataFile((section + "/" + subSection + "/default_parameter_list").data(), "SA");
    if (defList != "none")
        ML_Epetra::SetDefaults(defList, list);

    bool found;

    int MLPrintParameterList = dataFile((section + "/displayList").data(),      0, found);

    int MLOutput             = dataFile((section + "/" + subSection + "/MLOuput").data(),       0, found);
    if (found) list.set("ML output",               MLOutput);

    int printUnused          = dataFile((section + "/" + subSection + "/print_unused").data(), -2, found);
    if (found) list.set("print unused",            printUnused);

    int PDEEquations         = dataFile((section + "/" + subSection + "/pde_equations").data(),    1, found);
    if (found) list.set("PDE equations",           PDEEquations);

    int CycleApplications    = dataFile((section + "/" + subSection + "/cycle_applications").data(), 1, found);
    if (found) list.set("cycle applications", CycleApplications);

    int MaxLevels            = dataFile((section + "/" + subSection + "/max_levels").data(),         2, found);
    if (found) list.set("max levels", MaxLevels);

    std::string IncOrDec     = dataFile((section + "/" + subSection + "/inc_or_dec").data(),         "increasing", found);
    if (found) list.set("increasing or decreasing", IncOrDec);

    std::string precType     = dataFile((section + "/" + subSection + "/prec_type").data(),          "MGV", found);
    if (found) list.set("prec type", precType);

    //    int NumProjectedModes    = dataFile((section + "/" + subSection + "/number_of_prejected_modes").data(), 0);

    // if (found) list.set("ML print parameter list", MLPrintParameterList);

    std::string eigenAnalysisType = dataFile((section + "/" + subSection + "/eigne-analysis/type").data(), "cg", found);
    if (found) list.set("eigen-analysis: type",       eigenAnalysisType);

    int eigenAnalysisIterations   = dataFile((section + "/" + subSection + "/eigne-analysis/iterations").data(), 10, found);
    if (found) list.set("eigen-analysis: iterations", eigenAnalysisIterations);

    // Aggregation options

    std::string AggregationType                  = dataFile((section + "/" + subSection + "/aggregation/type").data(), "Uncoupled", found);
    if (found) list.set("aggregation: type",                              AggregationType);

    double      AggregationThreshold             = dataFile((section + "/" + subSection + "/aggregation/threshold").data(), 0.0 , found);
    if (found) list.set("aggregation: threshold",                         AggregationThreshold);

    double      AggregationDampingFactor         = dataFile((section + "/" + subSection + "/aggregation/damping_factor").data(), 4./3. , found);
    if (found) list.set("aggregation: damping factor",                    AggregationDampingFactor);

    int AggregationSmoothingSweeps               = dataFile((section + "/" + subSection + "/aggregation/smoothing_sweeps").data(),    1, found);
    if (found) list.set("aggregation: smoothing sweeps",                 AggregationSmoothingSweeps);

    int AggregationGlobalAggregates              = dataFile((section + "/" + subSection + "/aggregation/global_aggregates").data(),   1, found);
    if (found) list.set("aggregation: global aggregates",                 AggregationGlobalAggregates);

    int AggregationLocalAggregates               = dataFile((section + "/" + subSection + "/aggregation/local_aggregates").data(),    1, found);
    if (found) list.set("aggregation: local aggregates",                  AggregationLocalAggregates);

    int AggregationNodesPerAggregate             = dataFile((section + "/" + subSection + "/aggregation/nodes_per_aggregate").data(), 1, found);
    if (found) list.set("aggregation: nodes per aggregate",               AggregationNodesPerAggregate);

    int AggregationNextLevelAggregatesPerProcess = dataFile((section + "/" + subSection + "/aggregation/next-level_aggregates_per_process").data(), 128, found);
    if (found) list.set("aggregation: next level aggregates per process", AggregationNextLevelAggregatesPerProcess);

    bool AggregationUseTentativeRestriction      = dataFile((section + "/" + subSection + "/aggregation/tentative_restriction").data(), false, found);
    if (found) list.set("aggregation: use tentative restriction",         AggregationUseTentativeRestriction);

    bool AggregationSymmetrize                   = dataFile((section + "/" + subSection + "/aggregation/symmetrize").data(),            false, found);
    if (found) list.set("aggregation: symmetrize",                        AggregationSymmetrize);

    bool   EnergyMinimizationEnable              = dataFile((section + "/" + subSection + "/energy_minimization/enable").data(),        false, found);
    if (found) list.set("energy minimization: enable",                    EnergyMinimizationEnable);

    int    EnergyMinimizationType                = dataFile((section + "/" + subSection + "/energy_minimization/type").data(),          2, found);
    if (found) list.set("energy minimization: type",                      EnergyMinimizationType);

    double EnergyMinimizationDropTol             = dataFile((section + "/" + subSection + "/energy_minimization/droptol").data(),       0., found);
    if (found) list.set("energy minimization: droptol",                   EnergyMinimizationDropTol);

    bool   EnergyMinimizationCheap               = dataFile((section + "/" + subSection + "/energy_minimization/cheap").data(),         false, found);
    if (found) list.set("energy minimization: cheap",                     EnergyMinimizationCheap);

    // Smoothing parameters


    std::string SmootherType                = dataFile((section + "/" + subSection + "/smoother/type").data(), "IFPACK", found);
    if (found) list.set("smoother: type",                         SmootherType);

    int SmootherSweeps                      = dataFile((section + "/" + subSection + "/smoother/sweeps").data(), 2, found);
    if (found) list.set("smoother: sweeps",                       SmootherSweeps);

    double SmootherDampingFactor            = dataFile((section + "/" + subSection + "/smoother/damping_factor").data(), 1.0, found);
    if (found) list.set("smoother: damping factor",               SmootherDampingFactor);

    std::string SmootherPreOrPost           = dataFile((section + "/" + subSection + "/smoother/pre_or_post").data(), "both", found);
    if (found) list.set("smoother: pre or post",                  SmootherPreOrPost);

    double SmootherChebyshevAlpha           = dataFile((section + "/" + subSection + "/smoother/Chebyshev_alpha").data(), 20., found);
    if (found) list.set("smoother: Chebyshev alpha",              SmootherChebyshevAlpha);

    bool SmootherHiptmairEfficientSymmetric = dataFile((section + "/" + subSection + "/smoother/Hiptmair_efficient_symmetric").data(), true, found);
    if (found) list.set("smoother: Hiptmair efficient symmetric", SmootherHiptmairEfficientSymmetric);


    // subsmoother parameter

    std::string SubSmootherType             = dataFile((section + "/" + subSection + "/subsmoother/type").data(), "Chebyshev", found);
    if (found) list.set("subsmoother: type",                      SubSmootherType);

    double SubSmootherChebyshevAlpha        = dataFile((section + "/" + subSection + "/subsmoother/Chebyshev_alpha").data(), 20., found);
    if (found) list.set("subsmoother: Chebyshev alpha",           SubSmootherChebyshevAlpha);

    //    double SubSmootherSGSDampingFactor      = dataFile((section + "/" + subSection + "/subsmoothers/SGS_damping_factor").data(), 1., found);
    int SubSmootherEdgeSweeps               = dataFile((section + "/" + subSection + "/subsmoother/edge_sweeps").data(), 2, found);
    if (found) list.set("subsmoother: edge sweeps",               SubSmootherEdgeSweeps);

    int SubSmootherNodeSweeps               = dataFile((section + "/" + subSection + "/subsmoother/node_sweeps").data(), 2, found);
    if (found) list.set("subsmoother: node sweeps",               SubSmootherNodeSweeps);




//     if (found) list.set("smoother: type (level 0)","IFPACK");
//     if (found) list.set("smoother: type (level 1)","IFPACK");

    //    if (found) list.set("subsmoother: SGS damping factor",        SubSmootherSGSDampingFactor);

    // Coarsest Grid Parameters

    int CoarseMaxSize                 = dataFile((section + "/" + subSection + "/coarse/max_size").data(), 128, found);
    if (found) list.set("coarse: max size",  CoarseMaxSize);

    std::string CoarseType            = dataFile((section + "/" + subSection + "/coarse/type").data(), "Chebyshev", found);
    if (found) list.set("coarse: type", CoarseType);

    std::string CoarsePreOrPost       = dataFile((section + "/" + subSection + "/coarse/pre_or_post").data(), "post", found);
    if (found) list.set("coarse: pre or post", CoarsePreOrPost);

    int CoarseSweeps                  = dataFile((section + "/" + subSection + "/coarse/sweeps").data(), 2, found);
    if (found) list.set("coarse: sweeps", CoarseSweeps);

    double CoarseDampingFactor        = dataFile((section + "/" + subSection + "/coarse/damping_factor").data(), 1.0, found);
    if (found) list.set("coarse: damping factor", CoarseDampingFactor);

    std::string CoarseSubsmootherType = dataFile((section + "/" + subSection + "/coarse/subsmoother_type").data(), "Chebyshev", found);
    if (found) list.set("coarse: subsmoother type", CoarseSubsmootherType);

    int CoarseNodeSweeps              = dataFile((section + "/" + subSection + "/coarse/node_sweeps").data(), 2, found);
    if (found) list.set("coarse: node sweeps", CoarseNodeSweeps);

    int CoarseEdgeSweeps              = dataFile((section + "/" + subSection + "/coarse/edge_sweeps").data(), 2, found);
    if (found) list.set("coarse: edge sweeps", CoarseEdgeSweeps);

    double CoarseChebyshevAlpha       = dataFile((section + "/" + subSection + "/coarse/Chebyshev_alpha").data(), 30., found);
    if (found) list.set("coarse: Chebyshev alpha", CoarseChebyshevAlpha);

    int CoarseMaxProcesses            = dataFile((section + "/" + subSection + "/coarse/max_processes").data(),-1, found);
    if (found) list.set("coarse: max processes", CoarseMaxProcesses);


    // Load-balancing Options
    int RepartitionEnable              = dataFile((section + "/" + subSection + "/repartition/enable").data(), 0, found);
    if (found) list.set("repartition: enable",             RepartitionEnable);

    std::string RepartitionPartitioner = dataFile((section + "/" + subSection + "/repartition/partitioner").data(), "ParMETIS", found);
    if (found) list.set("repartition: partitioner",        RepartitionPartitioner);

    double RepartitionMaxMinRatio      = dataFile((section + "/" + subSection + "/repartition/max_min_ratio").data(), 1.3, found);
    if (found) list.set("repartition: max min ratio",      RepartitionMaxMinRatio);

    int RepartitionMinPerProc          = dataFile((section + "/" + subSection + "/repartition/min_per_proc").data(), 512, found);
    if (found) list.set("repartition: min per proc",       RepartitionMinPerProc);

    double RepartitionNodeMaxMinRatio  = dataFile((section + "/" + subSection + "/repartition/node_max_min_ratio").data(), 1.3, found);
    if (found) list.set("repartition: node max min ratio", RepartitionNodeMaxMinRatio);

    int RepartitionNodeMinPerProc      = dataFile((section + "/" + subSection + "/repartition/node_min_per_proc").data(), 170, found);
    if (found) list.set("repartition: node min per proc",  RepartitionNodeMinPerProc);

    int RepartitionZoltanDimensions    = dataFile((section + "/" + subSection + "/repartition/Zoltan_dimensions").data(), 2, found);
    if (found) list.set("repartition: Zoltan dimensions",  RepartitionZoltanDimensions);


    if (MLPrintParameterList)
    {
        std::cout << "  Parameters List: " << std::endl;
        list.print(std::cout);
        std::cout << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
MLPreconditioner::setDataFromGetPot( const GetPot&          dataFile,
                                     const std::string&     section )
{
    M_analyze = dataFile((section + "/" + "ML" + "/analyze_smoother").data(), false); // To be moved in createMLList

    // ML List
    createMLList(M_List, dataFile, section, "ML" );

    // IfPack list
    list_Type& SmootherIFSubList = M_List.sublist("smoother: ifpack list");
    IfpackPreconditioner::createIfpackList(SmootherIFSubList, dataFile, section, "ML");
}


// ===================================================
// Get Methods
// ===================================================
Real
MLPreconditioner::Condest()
{
    return 0.;
}

EpetraPreconditioner::prec_raw_type*
MLPreconditioner::getPrec()
{
    return M_Prec.get();
}

} // namespace LifeV
