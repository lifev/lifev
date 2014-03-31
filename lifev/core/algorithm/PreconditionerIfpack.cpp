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
    @brief Ifpack preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-2006
 */

#include <lifev/core/LifeV.hpp>
#include "PreconditionerIfpack.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
PreconditionerIfpack::PreconditionerIfpack ( boost::shared_ptr<Epetra_Comm> comm ) :
    super (),
    M_preconditioner(),
    M_comm ( comm ),
    M_overlapLevel (0),
    M_operator()
{

}

PreconditionerIfpack::~PreconditionerIfpack()
{

}


// ===================================================
// Methods
// ===================================================
Int
PreconditionerIfpack::buildPreconditioner ( operator_type& matrix )
{
    M_operator     = matrix->matrixPtr();

    M_overlapLevel = this->M_list.get ( "overlap level", -1 );

    M_precType     = this->M_list.get ( "prectype", "Amesos" );

#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory factory;
#else
    Ifpack factory;
#endif

    M_preconditioner.reset ( factory.Create ( M_precType, M_operator.get(), M_overlapLevel ) );

    M_precType += "_Ifpack";

    if ( !M_preconditioner.get() )
    {
        ERROR_MSG ( "Preconditioner not set, something went wrong in its computation\n" );
    }

    IFPACK_CHK_ERR ( M_preconditioner->SetParameters ( this->M_list ) );
    IFPACK_CHK_ERR ( M_preconditioner->Initialize() );
    IFPACK_CHK_ERR ( M_preconditioner->Compute() );

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

void
PreconditionerIfpack::resetPreconditioner()
{
    M_operator.reset();
    M_preconditioner.reset();

    this->M_preconditionerCreated = false;
}

void
PreconditionerIfpack::createParametersList ( list_Type&         list,
                                             const GetPot&      dataFile,
                                             const std::string& section,
                                             const std::string& subSection )
{
    createIfpackList ( list, dataFile, section, subSection, M_comm->MyPID() == 0 );
}

void
PreconditionerIfpack::createIfpackList ( list_Type&         list,
                                         const GetPot&      dataFile,
                                         const std::string& section,
                                         const std::string& subSection,
                                         const bool&        verbose )
{
    // See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    // for more informations on the parameters

    Int overlapLevel     = dataFile ( (section + "/" + subSection + "/overlap").data(), 0 );

    std::string precType = dataFile ( (section + "/" + subSection + "/prectype").data(), "Amesos" );

    list.set ( "prectype", precType );
    list.set ( "overlap level", overlapLevel );

    bool displayList = dataFile ( (section + "/displayList").data(), false );

    std::string relaxationType              = dataFile ( (section + "/" + subSection + "/relaxation/type").data(), "Jacobi" );
    Int         relaxationSweeps            = dataFile ( (section + "/" + subSection + "/relaxation/sweeps").data(), 1 );
    Real        relaxationDampingFactor     = dataFile ( (section + "/" + subSection + "/relaxation/damping_factor").data(), 1.0 );
    Real        relaxationMinDiagValue      = dataFile ( (section + "/" + subSection + "/relaxation/min_diagonal_value").data(), 0.0 );
    bool        relaxationZeroStartSolution = dataFile ( (section + "/" + subSection + "/relaxation/zero_starting_solution").data(), true );

    list.set ( "relaxation: type", relaxationType );
    list.set ( "relaxation: sweeps", relaxationSweeps );
    list.set ( "relaxation: damping factor", relaxationDampingFactor );
    list.set ( "relaxation: min diagonal value", relaxationMinDiagValue );
    list.set ( "relaxation: zero starting solution", relaxationZeroStartSolution );

    std::string partitionerType         = dataFile ( (section + "/" + subSection + "/partitioner/type").data(), "metis" );
    Int         partitionerOverlap      = dataFile ( (section + "/" + subSection + "/partitioner/overlap").data(), 0 );
    Int         partitionerLocalParts   = dataFile ( (section + "/" + subSection + "/partitioner/local_parts").data(), 1 );
    Int         partitionerRootNode     = dataFile ( (section + "/" + subSection + "/partitioner/root_node").data(), 0 );
    bool        partitionerUseSymmGraph = dataFile ( (section + "/" + subSection + "/partitioner/use_symmetric_graph").data(), true );

    list.set ( "partitioner: type",                partitionerType );
    list.set ( "partitioner: overlap",             partitionerOverlap );
    list.set ( "partitioner: local parts",         partitionerLocalParts );
    list.set ( "partitioner: root node",           partitionerRootNode );
    list.set ( "partitioner: use symmetric graph", partitionerUseSymmGraph );

    std::string amesosSolverType = dataFile ( (section + "/" + subSection + "/amesos/solvertype").data(), "Amesos_KLU" );

    list.set ( "amesos: solver type", amesosSolverType );

    Int    levelOfFill     = dataFile ( (section + "/" + subSection + "/fact/level-of-fill").data(),      4  );
    Real   ILUTlevelOfFill = dataFile ( (section + "/" + subSection + "/fact/ilut_level-of-fill").data(), 4. );
    Real   athr            = dataFile ( (section + "/" + subSection + "/fact/absolute_threshold").data(), 0. );
    Real   rthr            = dataFile ( (section + "/" + subSection + "/fact/relative_threshold").data(), 1. );
    Real   relaxValue      = dataFile ( (section + "/" + subSection + "/fact/relax_value").data(), 0. );
    Real   dropTolerance   = dataFile ( (section + "/" + subSection + "/fact/drop_tolerance").data(), 1e-5 );

    list.set ( "fact: drop tolerance", dropTolerance );
    list.set ( "fact: level-of-fill", levelOfFill );
    list.set ( "fact: ilut level-of-fill", ILUTlevelOfFill );
    list.set ( "fact: absolute threshold", athr );
    list.set ( "fact: relative threshold", rthr );
    list.set ( "fact: relax value", relaxValue );

    Int combineMode = dataFile ( (section + "/" + subSection + "/schwarz/combine_mode").data(), 0 );
    Epetra_CombineMode schwarzCombineMode;

    switch (combineMode)
    {
        case 0 :
            schwarzCombineMode = Add;
            break;
        case 1 :
            schwarzCombineMode = Zero;
            break;
        case 2 :
            schwarzCombineMode = Insert;
            break;
        case 3 :
            schwarzCombineMode = Average;
            break;
        case 4 :
            schwarzCombineMode = AbsMax;
            break;
        default:
            schwarzCombineMode = Zero;
    }

    bool schwarzComputeCondest        = dataFile ( (section + "/" + subSection + "/schwarz/compute_condest").data(), true );
    std::string schwarzReorderingType = dataFile ( (section + "/" + subSection + "/schwarz/reordering_type").data(), "none" );
    bool schwarzFilterSingletons      = dataFile ( (section + "/" + subSection + "/schwarz/filter_singletons").data(), true );

    list.set ( "schwarz: combine mode", schwarzCombineMode);
    list.set ( "schwarz: compute condest", schwarzComputeCondest);
    list.set ( "schwarz: reordering type", schwarzReorderingType);
    list.set ( "schwarz: filter singletons", schwarzFilterSingletons);

    // New parameters related to ShyLU and parallel subdomain problems
    Int subdomainSize = dataFile ( (section + "/" + subSection + "/subdomain/number_of_processors").data(), 1);
    list.set ("subdomain: number-of-processors", subdomainSize);

    // ShyLU parameters
    Teuchos::ParameterList shyluList;
    std::string outerSolverLibrary = dataFile ( (section + "/" + subSection + "/shylu/outer_solver_library").data(), "Belos");
    shyluList.set ("Outer Solver Library", outerSolverLibrary);
    std::string separatorType = dataFile ( (section + "/" + subSection + "/shylu/separator_type").data(), "Wide");
    shyluList.set ("Separator Type", separatorType);
    std::string schurApproxMethod = dataFile ( (section + "/" + subSection + "/shylu/schur_approx_method").data(), "A22AndBlockDiagonals");
    shyluList.set ("Schur Approximation Method", schurApproxMethod);
    double relativeThreshold = dataFile ( (section + "/" + subSection + "/shylu/relative_threshold").data(), 1e-3);
    shyluList.set ("Relative Threshold", relativeThreshold);
    double diagonalFactor = dataFile ( (section + "/" + subSection + "/shylu/diagonal_factor").data(), 0.02);
    shyluList.set ("Diagonal Factor", diagonalFactor);
    std::string schurComplementSolver = dataFile ( (section + "/" + subSection + "/shylu/schur_complement_solver").data(), "AztecOO-Exact");
    shyluList.set ("Schur Complement Solver", schurComplementSolver);
    std::string schurAmesosSolver = dataFile ( (section + "/" + subSection + "/shylu/schur_amesos_solver").data(), "Amesos_Klu");
    shyluList.set ("Schur Amesos Solver", schurAmesosSolver);
    std::string schurPrec = dataFile ( (section + "/" + subSection + "/shylu/schur_prec").data(), "ILU stand-alone");
    shyluList.set ("Schur Preconditioner", schurPrec);
    Int shyluSymmetry = dataFile ( (section + "/" + subSection + "/shylu/symmetry").data(), 1);
    shyluList.set ("Symmetry", 1);
    Int innerMaxIter = dataFile ( (section + "/" + subSection + "/shylu/inner_solver_iterations").data(), 5);
    shyluList.set ("Inner Solver MaxIters", innerMaxIter);
    double innerTol = dataFile ( (section + "/" + subSection + "/shylu/inner_solver_tolerance").data(), 1e-10);
    shyluList.set ("Inner Solver Tolerance", innerTol);
    bool silentSubiterations = dataFile ( (section + "/" + subSection + "/shylu/silent_subiterations").data(), true);
    shyluList.set ("Silent subiterations", silentSubiterations);
    std::string shyluDiagSolver = dataFile ( (section + "/" + subSection + "/shylu/diag_solver").data(), "Amesos_Klu");
    shyluList.set ("Diagonal Block Solver", shyluDiagSolver);
    double iqrKrylovDim = dataFile ( (section + "/" + subSection + "/shylu/iqr_krylov_dim").data(), 0.5);
    shyluList.set ("IQR Krylov Dim", iqrKrylovDim);
    Int iqrNumIter = dataFile ( (section + "/" + subSection + "/shylu/iqr_num_iter").data(), 0);
    shyluList.set ("IQR Number Iterations", iqrNumIter);
    bool iqrScaling = dataFile ( (section + "/" + subSection + "/shylu/iqr_scaling").data(), true);
    shyluList.set ("IQR Scaling", iqrScaling);
    std::string iqrInitialPrecType = dataFile ( (section + "/" + subSection + "/shylu/iqr_initial_prec_type").data(), "Amesos");
    shyluList.set ("IQR Initial Prec Type", iqrInitialPrecType);
    std::string iqrInitialPrecAmesosType = dataFile ( (section + "/" + subSection + "/shylu/iqr_initial_prec_amesos_type").data(), "Amesos_Klu");
    shyluList.set ("IQR Initial Prec Amesos Type", iqrInitialPrecAmesosType);

    list.set ("ShyLU list", shyluList);

    if ( displayList && verbose )
    {
        std::cout << "Ifpack parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        list.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

Int
PreconditionerIfpack::ApplyInverse ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
{
    return M_preconditioner->ApplyInverse ( vector1, vector2 );
}

Int
PreconditionerIfpack::Apply ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
{
    return M_preconditioner->Apply ( vector1, vector2 );
}

void
PreconditionerIfpack::showMe ( std::ostream& output ) const
{
    output << "showMe must be implemented for the PreconditionerIfpack class" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerIfpack::setDataFromGetPot ( const GetPot&      dataFile,
                                          const std::string& section )
{
    createIfpackList ( this->M_list, dataFile, section, "ifpack", M_comm->MyPID() == 0 );
}

Int
PreconditionerIfpack::SetUseTranspose ( bool useTranspose )
{
    return M_preconditioner->SetUseTranspose ( useTranspose );
}

// ===================================================
// Get Methods
// ===================================================
Real
PreconditionerIfpack::condest()
{
    return M_preconditioner->Condest();
}

Preconditioner::prec_raw_type*
PreconditionerIfpack::preconditioner()
{
    return M_preconditioner.get();
}

PreconditionerIfpack::super::prec_type
PreconditionerIfpack::preconditionerPtr()
{
    return M_preconditioner;
}

std::string
PreconditionerIfpack::preconditionerType()
{
    return M_precType;
}

const Int&
PreconditionerIfpack::getOverlapLevel() const
{
    return M_overlapLevel;
}

bool
PreconditionerIfpack::UseTranspose()
{
    return M_preconditioner->UseTranspose();
}

const Epetra_Map&
PreconditionerIfpack::OperatorRangeMap() const
{
    return M_preconditioner->OperatorRangeMap();
}

const Epetra_Map&
PreconditionerIfpack::OperatorDomainMap() const
{
    return M_preconditioner->OperatorDomainMap();
}

} // namespace LifeV
