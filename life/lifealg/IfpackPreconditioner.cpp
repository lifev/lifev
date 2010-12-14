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

#include <life/lifecore/life.hpp>
#include "IfpackPreconditioner.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
IfpackPreconditioner::IfpackPreconditioner():
        super (),
        M_preconditioner(),
        M_overlapLevel(0),
        M_operator()
{

}

IfpackPreconditioner::~IfpackPreconditioner()
{

}


// ===================================================
// Methods
// ===================================================
Int
IfpackPreconditioner::buildPreconditioner( operator_type& matrix )
{
    M_operator     = matrix->getMatrixPtr();

    M_overlapLevel = this->M_list.get( "overlap level", -1 );

    M_precType     = this->M_list.get( "prectype", "Amesos" );

    Ifpack factory;

    M_preconditioner.reset( factory.Create( M_precType, M_operator.get(), M_overlapLevel ) );

    M_precType += "_Ifpack";

    if ( !M_preconditioner.get() )
    {
        ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
    }

    IFPACK_CHK_ERR( M_preconditioner->SetParameters( this->M_list ) );
    IFPACK_CHK_ERR( M_preconditioner->Initialize() );
    IFPACK_CHK_ERR( M_preconditioner->Compute() );

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

void
IfpackPreconditioner::precReset()
{
    M_operator.reset();
    M_preconditioner.reset();

    this->M_preconditionerCreated = false;
}

void
IfpackPreconditioner::createList( list_Type&         list,
                                  const GetPot&      dataFile,
                                  const std::string& section,
                                  const std::string& subSection )
{
    createIfpackList( list, dataFile, section, subSection );
}

void
IfpackPreconditioner::createIfpackList( list_Type&         list,
                                        const GetPot&      dataFile,
                                        const std::string& section,
                                        const std::string& subSection )
{
    // See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    // for more informations on the parameters

    Int overlapLevel     = dataFile( (section + "/" + subSection + "/overlap").data(), 0 );

    std::string precType = dataFile( (section + "/" + subSection + "/prectype").data(), "Amesos" );

    list.set( "prectype", precType );
    list.set( "overlap level", overlapLevel );

    bool displayList = dataFile( (section + "/displayList").data(), false );

    std::string relaxationType              = dataFile( (section + "/" + subSection + "/relaxation/type").data(), "Jacobi" );
    Int         relaxationSweeps            = dataFile( (section + "/" + subSection + "/relaxation/sweeps").data(), 1 );
    Real        relaxationDampingFactor     = dataFile( (section + "/" + subSection + "/relaxation/damping_factor").data(), 1.0 );
    Real        relaxationMinDiagValue      = dataFile( (section + "/" + subSection + "/relaxation/min_diagonal_value").data(), 0.0 );
    bool        relaxationZeroStartSolution = dataFile( (section + "/" + subSection + "/relaxation/zero_starting_solution").data(), true );

    list.set( "relaxation: type", relaxationType );
    list.set( "relaxation: sweeps", relaxationSweeps );
    list.set( "relaxation: damping factor", relaxationDampingFactor );
    list.set( "relaxation: min diagonal value", relaxationMinDiagValue );
    list.set( "relaxation: zero starting solution", relaxationZeroStartSolution );

    std::string partitionerType         = dataFile( (section + "/" + subSection + "/partitioner/type").data(), "metis" );
    Int         partitionerOverlap      = dataFile( (section + "/" + subSection + "/partitioner/overlap").data(), 0 );
    Int         partitionerLocalParts   = dataFile( (section + "/" + subSection + "/partitioner/local_parts").data(), 1 );
    Int         partitionerRootNode     = dataFile( (section + "/" + subSection + "/partitioner/root_node").data(), 0 );
    bool        partitionerUseSymmGraph = dataFile( (section + "/" + subSection + "/partitioner/use_symmetric_graph").data(), true );

    list.set( "partitioner: type",                partitionerType );
    list.set( "partitioner: overlap",             partitionerOverlap );
    list.set( "partitioner: local parts",         partitionerLocalParts );
    list.set( "partitioner: root node",           partitionerRootNode );
    list.set( "partitioner: use symmetric graph", partitionerUseSymmGraph );

    std::string amesosSolverType = dataFile( (section + "/" + subSection + "/amesos/solvertype").data(), "Amesos_KLU" );

    list.set( "amesos: solver type", amesosSolverType );

    Real   levelOfFill     = dataFile( (section + "/" + subSection + "/fact/level-of-fill").data(),      4. );
    Real   ILUTlevelOfFill = dataFile( (section + "/" + subSection + "/fact/ilut_level-of-fill").data(), 4. );
    Real   athr            = dataFile( (section + "/" + subSection + "/fact/absolute_threshold").data(), 0. );
    Real   rthr            = dataFile( (section + "/" + subSection + "/fact/relative_threshold").data(), 1. );
    Real   relaxValue      = dataFile( (section + "/" + subSection + "/fact/relax_value").data(), 0. );
    Real   dropTolerance   = dataFile( (section + "/" + subSection + "/fact/drop_tolerance").data(), 1e-5 );

    list.set( "fact: drop tolerance", dropTolerance );
    list.set( "fact: level-of-fill", levelOfFill );
    list.set( "fact: ilut level-of-fill", ILUTlevelOfFill );
    list.set( "fact: absolute threshold", athr );
    list.set( "fact: relative threshold", rthr );
    list.set( "fact: relax value", relaxValue );

    Int combineMode = dataFile( (section + "/" + subSection + "/schwarz/combine_mode").data(), 0 );
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

    bool schwarzComputeCondest        = dataFile( (section + "/" + subSection + "/schwarz/compute_condest").data(), true );
    std::string schwarzReorderingType = dataFile( (section + "/" + subSection + "/schwarz/reordering_type").data(), "none" );
    bool schwarzFilterSingletons      = dataFile( (section + "/" + subSection + "/schwarz/filter_singletons").data(), true );

    list.set( "schwarz: combine mode", schwarzCombineMode);
    list.set( "schwarz: compute condest", schwarzComputeCondest);
    list.set( "schwarz: reordering type", schwarzReorderingType);
    list.set( "schwarz: filter singletons", schwarzFilterSingletons);

    if ( displayList ) list.print( std::cout );
}

Int
IfpackPreconditioner::ApplyInverse( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
{
    return M_preconditioner->ApplyInverse( vector1, vector2 );
}

Int
IfpackPreconditioner::Apply( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
{
    return M_preconditioner->Apply( vector1, vector2 );
}

void
IfpackPreconditioner::showMe( std::ostream& output ) const
{
    output << "showMe must be implemented for the IfpackPreconditioner class" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
IfpackPreconditioner::setDataFromGetPot( const GetPot&      dataFile,
                                         const std::string& section )
{
    createIfpackList( this->M_list, dataFile, section, "ifpack" );
}

Int
IfpackPreconditioner::SetUseTranspose( bool useTranspose )
{
    return M_preconditioner->SetUseTranspose( useTranspose );
}

// ===================================================
// Get Methods
// ===================================================
bool
IfpackPreconditioner::set() const
{
    return M_preconditioner;
}

Real
IfpackPreconditioner::Condest()
{
    return M_preconditioner->Condest();
}

EpetraPreconditioner::prec_raw_type*
IfpackPreconditioner::getPrec()
{
    return M_preconditioner.get();
}

IfpackPreconditioner::super::prec_type
IfpackPreconditioner::getPrecPtr()
{
    return M_preconditioner;
}

std::string
IfpackPreconditioner::precType()
{
    return M_precType;
}

const Int&
IfpackPreconditioner::getOverlapLevel() const
{
    return M_overlapLevel;
}

bool
IfpackPreconditioner::UseTranspose()
{
    return M_preconditioner->UseTranspose();
}

const Epetra_Map &
IfpackPreconditioner::OperatorRangeMap() const
{
    return M_preconditioner->OperatorRangeMap();
}

const Epetra_Map &
IfpackPreconditioner::OperatorDomainMap() const
{
    return M_preconditioner->OperatorDomainMap();
}

} // namespace LifeV
