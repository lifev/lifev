/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-10-08

  Copyright (C) 2010 EPFL

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
   \file LSCPreconditioner.cpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-08
 */

#include "LSCPreconditioner.hpp"

namespace LifeV {

LSCPreconditioner::LSCPreconditioner():
        TekoPreconditioner (),
        M_precType(""),
        M_velocityBlockSize(-1),
        M_pressureBlockSize(-1)
{}

LSCPreconditioner::~LSCPreconditioner()
{}

void LSCPreconditioner::setDataFromGetPot( const GetPot& dataFile,
                                           const std::string& section )
{
    createLSCList(dataFile, section, this->M_List);

    M_velocityBlockSize = this->M_List.get("blocks: velocity block size", -1);
    M_pressureBlockSize = this->M_List.get("blocks: pressure block size", -1);
    M_precType          = this->M_List.get("prectype", "LSC");
}

int LSCPreconditioner::buildPreconditioner(operator_type& oper)
{
    // Creating the InverseLibrary from Stratimikos
    RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

    // build the inverse factory needed by the example preconditioner
    RCP<Teko::InverseFactory> inverse
        = invLib->getInverseFactory("Amesos");

    // Building the LSC strategy
    RCP<Teko::NS::LSCStrategy> strategy = rcp(new Teko::NS::InvLSCStrategy(inverse,true));

    // Building the LSC preconditioner factory
    RCP<Teko::BlockPreconditionerFactory> precFact
        = rcp(new Teko::NS::LSCPreconditionerFactory(strategy));

    // Building Block sizes
    M_velocityBlockSize = this->M_List.get("blocks: velocity block size", -1);
    M_pressureBlockSize = this->M_List.get("blocks: pressure block size", -1);
    std::vector<int> blockSizes;
    blockSizes.push_back(M_velocityBlockSize);
    blockSizes.push_back(M_pressureBlockSize);

    // Building the LSC preconditioner
    buildTekoPreconditioner(precFact,oper,blockSizes);

    return ( EXIT_SUCCESS );
}

Real LSCPreconditioner::Condest()
{
    return 0.0;
}

int LSCPreconditioner::numBlockRow() const
{
    return 2;
}

int LSCPreconditioner::numBlockCol() const
{
    return 2;
}

void createLSCList( const GetPot&              dataFile,
                    const std::string&         section,
                    Teuchos::ParameterList&    list)
{
    //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    //! for more informations on the parameters

    bool displayList = dataFile((section + "/displayList").data(),     false);

    std::string precType = dataFile((section + "/prectype").data(),"LSC");
    list.set("prectype", precType);

    //std::string value1 = dataFile((section + "/LSC/subsection1/value1").data(), "string example");
    //int         value2 = dataFile((section + "/LSC/subsection2/value2").data(), 1);
    //double      value3 = dataFile((section + "/LSC/subsection3/value3").data(), 1.0);
    //list.set("subsection1: value1", value1);
    //list.set("subsection2: value2", value2);
    //list.set("subsection3: value3", value3);

    int velocityBlockSize = dataFile((section + "/LSC/blocks/velocity_block_size").data(), -1);
    int pressureBlockSize = dataFile((section + "/LSC/blocks/pressure_block_size").data(), -1);
    std::cout << "section: " << section << std::endl;

    list.set("blocks: velocity block size", velocityBlockSize);
    list.set("blocks: pressure block size", pressureBlockSize);

    if (displayList) list.print(std::cout);
}

} // namespace LifeV
