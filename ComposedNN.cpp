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
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 29 Jun 2010

    A more detailed description of the file (if necessary)
 */

#include <ComposedNN.hpp>

#include <life/lifefem/bcManage.hpp>

namespace LifeV
{

int ComposedNN::solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver )
{

    boost::shared_ptr<ComposedPreconditioner<Epetra_Operator> > firstCompPrec(new ComposedPreconditioner<Epetra_Operator>(&M_blocks[0]->getMatrixPtr()->Comm()));
    boost::shared_ptr<ComposedPreconditioner<Epetra_Operator> > secondCompPrec(new ComposedPreconditioner<Epetra_Operator>(&M_blocks[0]->getMatrixPtr()->Comm()));

    Chrono chrono;
    int overlapLevel     = M_list.get("overlap level", 2);
    std::string precType = M_list.get("prectype", "Amesos");
    Ifpack factory;

    if(!(M_blockPrecs->getNumber()))
    {
        for(ID k(0); k < M_blocks.size(); ++k)
        {
            M_blockPrecs->displayer().leaderPrint("  M-  Computing double prec. factorization ...        ");
            chrono.start();
            M_prec[k].reset(factory.Create(precType, M_blocks[k]->getMatrixPtr().get(), overlapLevel));
            if ( !M_prec[k].get() )
            {
                ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
            }
            M_prec[k]->SetParameters(M_list);
            M_prec[k]->Initialize();
            M_prec[k]->Compute();
            chrono.stop();
            M_blockPrecs->displayer().leaderPrintMax("done in ", chrono.diff());
        }

        firstCompPrec->push_back(M_prec[0], false, false);
        firstCompPrec->push_back(M_prec[1], false, false);

        secondCompPrec->push_back(M_prec[2], false, false);
        secondCompPrec->push_back(M_prec[3], false, false);

        M_blockPrecs->resetP();
        M_blockPrecs->push_back(firstCompPrec, false, false, false);
        M_blockPrecs->push_back(secondCompPrec, false, false, true /*sum*/);

    }
//     else
//     {
//         for(ID k(0); k < M_blocks.getP().size(); ++k)
//         {
//             if(M_recompute[k])
//             {
//                 M_prec[k].reset(factory.Create(precType, dynamic_cast<Epetra_FECrsMatrix*>(M_blocks.getP()[k].get()), overlapLevel));
//                 if ( !M_prec[k].get() )
//                 {
//                     ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
//                 }
//                 M_prec[k]->SetParameters(M_list);
//                 M_prec[k]->Initialize();
//                 M_prec[k]->Compute();
//                 chrono.stop();
//             }
//         }
//     }
        return linearSolver->solveSystem(rhs, step, boost::static_pointer_cast<Epetra_Operator>(M_blockPrecs));
}


void ComposedNN::coupler(map_shared_ptrtype map,
                         const std::map<ID, ID>& locDofMap,
                         const vector_ptrtype numerationInterface,
                         const Real& timeStep)
{

    UInt totalDofs=map->getMap(Unique)->NumGlobalElements()+1;
    UInt fluidSolid=M_offset[0]+1+M_FESpace[0]->map().getMap(Unique)->NumGlobalElements();

    BlockInterface::swap(M_blocks[1], M_blocks[0]);
    BlockInterface::swap(M_bch[1], M_bch[0]);
    std::vector< fespace_ptrtype > tmpFESpace( M_FESpace );
    BlockInterface::swap(M_FESpace[1], M_FESpace[0]);

    for(ID k=0; k<2; ++k)
    {
        M_blocks[k]->GlobalAssemble();
        matrix_ptrtype block(new matrix_type(*M_blocks[k]));
        M_blocks.push_back(block);
        M_bch.push_back(M_bch[k]);
        M_FESpace.push_back(M_FESpace[k]);
        M_offset.push_back(M_offset[k]);
        M_recompute.push_back(M_recompute[k]);
    }

    matrix_ptrtype coupling(new matrix_type(*map));
    UInt one(1.);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal( one, M_offset[0]+1, fluidSolid );
    coupling->insertValueDiagonal( one, fluidSolid , totalDofs);
    couplingMatrix(coupling, 4, tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 4.);
    coupling->GlobalAssemble();
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal(one, M_offset[1]+1, M_offset[0]+1 );
    coupling->insertValueDiagonal(one,  fluidSolid, totalDofs);
    couplingMatrix(coupling, 8,  tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 4.);
    coupling->GlobalAssemble();
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal( one, M_offset[0]+1, fluidSolid );
    coupling->insertValueDiagonal( one, fluidSolid, totalDofs );
    couplingMatrix(coupling, 2, tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 4.);
    coupling->GlobalAssemble();
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal(one, M_offset[1]+1, M_offset[0]+1 );
    coupling->insertValueDiagonal(-1,  fluidSolid, totalDofs);
    couplingMatrix(coupling, 1,  tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 4.);
    coupling->GlobalAssemble();
    M_coupling.push_back(coupling);

    M_blockPrecs.reset(new ComposedPreconditioner<Epetra_Operator>(&M_blocks[0]->getMatrixPtr()->Comm()));
    M_prec.resize(M_blocks.size());
}

void ComposedNN::setDataFromGetPot(const GetPot& data, const std::string& section)
{
    createIfpackList(data, section, M_list);
}

void ComposedNN::push_back_matrix(const matrix_ptrtype& Mat, const  bool recompute)
{
    Mat->GlobalAssemble();
    *Mat *= 2.;
    super::push_back_matrix(Mat, recompute);
}



void ComposedNN::replace_matrix( matrix_ptrtype& oper, UInt position )
{
    oper->GlobalAssemble();
    *oper *= 2.;
    super::replace_matrix(oper, 1-position );
    super::replace_matrix(oper, 3-position );
}

} // Namespace LifeV
