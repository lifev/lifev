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

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-2006
 */

#include <Teuchos_RCPBoostSharedPtrConversions.hpp>

#include <lifev/core/algorithm/PreconditionerML2.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerMGOperators.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
PreconditionerML2::PreconditionerML2 ( boost::shared_ptr<Epetra_Comm> comm ) :
    super(),
    M_operator(),
    M_preconditioner(),
    M_restriction(),
    M_prolongation(),
    M_useP2ToP1Wrapper ( false ),
    M_uFESpace(),
    M_pFESpace(),
    M_comm ( comm )
{

}

PreconditionerML2::~PreconditionerML2()
{

}


// ===================================================
// Methods
// ===================================================
void
PreconditionerML2::createParametersList ( list_Type&         /*list*/,
                                          const GetPot&      /*dataFile*/,
                                          const std::string& /*section*/,
                                          const std::string& /*subSection*/ )
{
    // Not needed
}

Int
PreconditionerML2::buildPreconditioner ( operator_type& matrix )
{
    //the Trilinos::MultiLevelPreconditioner unsafely access to the area of memory co-owned by M_operator.
    //to avoid the risk of dandling pointers always deallocate M_preconditioner first and then M_operator
    M_preconditioner.reset();
    M_operator = matrix->matrixPtr();

    bool verbose = M_comm->MyPID() == 0;

    if ( M_useP2ToP1Wrapper )
    {
        if ( ( M_uFESpace.get() == NULL ) || ( M_pFESpace.get() == NULL ) )
        {
            if ( verbose)
            {
                std::cout << "You must specified manually the pointers to the FESpaces" << std::endl;
            }
            exit ( -1 );
        }

        if ( M_restriction.get() == 0 )
        {
            if ( verbose)
            {
                std::cout << "Building the restriction" << std::endl;
            }
            M_restriction.reset ( new matrix_Type ( matrix->map() ) );
            buildRestrictionP2ToP1 ( M_uFESpace, M_pFESpace, M_restriction );
        }
        if ( M_prolongation.get() == 0 )
        {
            if ( verbose)
            {
                std::cout << "Building the prolongation" << std::endl;
            }
            M_prolongation.reset ( new matrix_Type ( matrix->map() ) );
            buildProlongationP1ToP2 ( M_uFESpace, M_pFESpace, M_prolongation );
        }

        Teuchos::RCP<Epetra_CrsMatrix> rcpA = Teuchos::rcp<Epetra_CrsMatrix> ( M_operator.get(), false );
        Teuchos::RCP<Epetra_CrsMatrix> rcpP = Teuchos::rcp ( M_prolongation->matrixPtr() );
        Teuchos::RCP<Epetra_CrsMatrix> rcpR = Teuchos::rcp ( M_restriction->matrixPtr() );

        ML_Epetra::LevelWrap* MLPrec = new ML_Epetra::LevelWrap (rcpA, rcpP, rcpR, M_MLList, true);
        M_preconditioner.reset ( MLPrec );
    }
    else
    {
        M_preconditioner.reset ( new ML_Epetra::MultiLevelPreconditioner ( *M_operator, M_MLList, true ) );
    }

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

void
PreconditionerML2::resetPreconditioner()
{
    //the Trilinos::MultiLevelPreconditioner unsafely access to the area of memory co-owned by M_operator.
    //to avoid the risk of dandling pointers always deallocate M_preconditioner first and then M_operator

    M_preconditioner.reset();
    M_operator.reset();

    this->M_preconditionerCreated = false;
}

void PreconditionerML2::showMe ( std::ostream& output ) const
{
    output << "showMe must be implemented for the PreconditionerML2 class" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerML2::setDataFromGetPot ( const GetPot&      dataFile,
                                       const std::string& section )
{
    bool displayList = dataFile ( ( section + "/displayList" ).data(), false );
    const std::string MLParamFile = dataFile ( ( section + "/ML2/parameters_file" ).data(), "MLParamList.xml" );
    Teuchos::RCP< Teuchos::ParameterList > MLParameters = Teuchos::rcp ( new Teuchos::ParameterList );
    MLParameters = Teuchos::getParametersFromXmlFile ( MLParamFile );
    this->setParameters ( *MLParameters );

    bool verbose = M_comm->MyPID() == 0;

    if ( displayList && verbose )
    {
        std::cout << "ML2 parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        M_MLList.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

void
PreconditionerML2::setParameters ( Teuchos::ParameterList& list )
{
    std::string defaultListName = list.get ( "default list", "none" );
    if ( defaultListName == "wrap" )
    {
        ML_Epetra::SetDefaultsLevelWrap ( M_MLList );
    }
    else if ( defaultListName != "none" )
    {
        ML_Epetra::SetDefaults ( defaultListName, M_MLList );
    }
    M_MLList.setParameters ( list );
    M_MLList.remove ( "default list", false );

    M_useP2ToP1Wrapper = list.get ( "use P2 to P1 wrapper", false );
    M_MLList.remove ( "use P2 to P1 wrapper", false );
}

void
PreconditionerML2::setFESpace ( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace )
{
    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
}


// ===================================================
// Get Methods
// ===================================================
Real
PreconditionerML2::condest()
{
    return 0.0;
}

Preconditioner::prec_raw_type*
PreconditionerML2::preconditioner()
{
    return M_preconditioner.get();
}

} // namespace LifeV
