//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains the PreconditionerLSC class.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 08-11-2010
 */

#include "PreconditionerLSC.hpp"

namespace LifeV
{

PreconditionerLSC::PreconditionerLSC ( boost::shared_ptr<Epetra_Comm> comm ) :
    PreconditionerTeko  (),
    M_precType          ( "" ),
    M_velocityBlockSize ( -1 ),
    M_pressureBlockSize ( -1 ),
    M_comm ( comm )
{

}

PreconditionerLSC::~PreconditionerLSC()
{

}

void
PreconditionerLSC::setDataFromGetPot ( const GetPot& dataFile,
                                       const std::string& section )
{
    this->createParametersList ( M_list, dataFile, section, "LSC" );
    this->setParameters ( M_list );
}

void
PreconditionerLSC::setParameters ( Teuchos::ParameterList& list )
{
    M_precType          = list.get ( "prectype", "LSC" );
}

void
PreconditionerLSC::setFESpace ( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace )
{
    // We setup the size of the blocks
    M_velocityBlockSize = uFESpace->fieldDim() * uFESpace->dof().numTotalDof();
    M_pressureBlockSize = pFESpace->dof().numTotalDof();
}

int
PreconditionerLSC::buildPreconditioner ( matrixPtr_Type& oper )
{
    if ( ( M_velocityBlockSize < 0 ) || ( M_pressureBlockSize < 0 ) )
    {
        std::cout << "You must specify manually the pointers to the FESpaces" << std::endl;
        exit ( -1 );
    }

    // Creating the InverseLibrary from Stratimikos
    RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

    // build the inverse factory needed by the example preconditioner
    RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory ( "Ifpack" );

    // Building the LSC strategy
    RCP<Teko::NS::LSCStrategy> strategy = rcp ( new Teko::NS::InvLSCStrategy ( inverse, true ) );

    // Building the LSC preconditioner factory
    RCP<Teko::BlockPreconditionerFactory> precFact
        = rcp (new Teko::NS::LSCPreconditionerFactory ( strategy ) );

    // Building Block sizes
    std::vector<int> blockSizes;
    blockSizes.push_back ( M_velocityBlockSize );
    blockSizes.push_back ( M_pressureBlockSize );

    // Building the LSC preconditioner
    buildPreconditionerTeko ( precFact, oper, blockSizes );

    return ( EXIT_SUCCESS );
}

void
PreconditionerLSC::createParametersList ( list_Type&         list,
                                          const GetPot&      dataFile,
                                          const std::string& section,
                                          const std::string& /* subSection */ )
{
    bool verbose ( M_comm->MyPID() == 0 );

    //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    //! for more informations on the parameters

    bool displayList = dataFile ( ( section + "/displayList" ).data(), false );

    std::string precType = dataFile ( ( section + "/prectype" ).data(), "LSC" );
    list.set ( "prectype", precType );

    if ( displayList && verbose )
    {
        std::cout << "LSC parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        list.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

Real
PreconditionerLSC::condest()
{
    return 0.0;
}

int
PreconditionerLSC::numBlocksRows() const
{
    return 2;
}

int
PreconditionerLSC::numBlocksCols() const
{
    return 2;
}

} // namespace LifeV
