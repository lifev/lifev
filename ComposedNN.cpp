/* -*- mode: c++ -*- */
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

#include <lifeconfig.h>

#include <life/lifefem/bcManage.hpp>

#include <lifemc/lifesolver/ComposedNN.hpp>

namespace LifeV
{

// ===================================================
//! Public Methods
// ===================================================

int ComposedNN::solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver )
{
    M_firstCompPrec .reset(new composed_prec(M_comm));
    M_secondCompPrec.reset(new composed_prec(M_comm));

    Chrono chrono;
    int overlapLevel     = M_list.get("overlap level", 2);
    std::string precType = M_list.get("prectype", "Amesos");
    Ifpack factory;

    if ( !M_blockPrecs.get() )
        M_blockPrecs.reset(new ComposedOperator< composed_prec >(M_comm));
    if ( true || !set())
    {
        for (ID k(0); k < M_blocks.size(); ++k)
        {
            M_blockPrecs->displayer().leaderPrint("  M-  Computing double prec. factorization ...        ");
            chrono.start();
            M_prec[k].reset(factory.Create(precType, M_blocks[(*M_blockReordering)[k]]->getMatrixPtr().get(), overlapLevel));
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
    }
    else
    {
        for (ID k(0); k < M_blocks.size(); ++k)
        {
            if (M_recompute[(*M_blockReordering)[k]])
            {
                M_blockPrecs->displayer().leaderPrint("  M-  Computing double prec. factorization ...        ");
                chrono.start();
                M_prec[k].reset(factory.Create(precType, M_blocks[(*M_blockReordering)[k]]->getMatrixPtr().get(), overlapLevel));
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
            else
                M_blockPrecs->displayer().leaderPrint("  M-  Reusing double prec. factorization ...        \n");
        }
    }
    M_firstCompPrec->push_back(M_prec[0], false, false);
    M_firstCompPrec->push_back(M_prec[1], false, false);

    M_secondCompPrec->push_back(M_prec[2], false, false);
    M_secondCompPrec->push_back(M_prec[3], false, false);

    //M_blockPrecs->resetP();
    if (!(M_blockPrecs->getNumber()))
    {
        M_blockPrecs->push_back(M_firstCompPrec, false, false, false);
        M_blockPrecs->push_back(M_secondCompPrec, false, false, true /*sum*/);
    }
    else
    {
        M_blockPrecs->replace(M_firstCompPrec, (UInt)0, false, false);
        M_blockPrecs->replace(M_secondCompPrec, (UInt)1, false, false);
    }

    return linearSolver->solveSystem(rhs, step, boost::static_pointer_cast<Epetra_Operator>(M_blockPrecs));
}



void ComposedNN::setDataFromGetPot(const GetPot& data, const std::string& section)
{
    IfpackPreconditioner::createIfpackList( M_list, data, section, "");
}

void ComposedNN::coupler(map_shared_ptrtype& map,
                         const std::map<ID, ID>& locDofMap,
                         const vector_ptrtype& numerationInterface,
                         const Real& timeStep)
{
    UInt totalDofs=map->getMap(Unique)->NumGlobalElements()+1;
    UInt fluidSolid=M_offset[0]+1+M_FESpace[0]->map().getMap(Unique)->NumGlobalElements();

    for (ID k=0; k<2; ++k)
    {
        M_blocks[k]->GlobalAssemble();
        matrixPtr_Type block(new matrix_Type(*M_blocks[k]));
        M_blocks.push_back(block);
        M_bch.push_back(M_bch[k]);
        M_FESpace.push_back(M_FESpace[k]);
        M_offset.push_back(M_offset[k]);
        M_recompute[2+k]=(M_recompute[k]);
    }

    matrixPtr_Type coupling(new matrix_Type(*map));
    UInt one(1.);

    coupling.reset(new matrix_Type(*map, 0));
    coupling->insertValueDiagonal(1., M_offset[fluid]+1, M_offset[solid]+1 );
    coupling->insertValueDiagonal(1.,  fluidSolid, totalDofs);
    couplingMatrix(coupling, (*M_couplingFlags)[0]/*8*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_Type(*map, 0));
    coupling->insertValueDiagonal( 1., M_offset[0]+1, fluidSolid );
    coupling->insertValueDiagonal( 1., fluidSolid , totalDofs);
    couplingMatrix(coupling, (*M_couplingFlags)[1]/*4*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_Type(*map, 0));
    coupling->insertValueDiagonal(1., M_offset[fluid]+1, M_offset[solid]+1 );
    coupling->insertValueDiagonal(-1,  fluidSolid, totalDofs);
    couplingMatrix(coupling, (*M_couplingFlags)[2]/*1*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_Type(*map, 0));
    coupling->insertValueDiagonal( 1., M_offset[0]+1, fluidSolid );
    coupling->insertValueDiagonal( 1., fluidSolid, totalDofs );
    couplingMatrix(coupling, (*M_couplingFlags)[3]/*2*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);
    M_coupling.push_back(coupling);

    M_prec.resize(M_blocks.size());
}

void ComposedNN::applyBoundaryConditions(const Real& time, const UInt i)
{
    M_blocks[i]->openCrsMatrix();
    if ( !M_bch[i]->bdUpdateDone() )
    {
        M_bch[i]->bdUpdate( *M_FESpace[i]->mesh(), M_FESpace[i]->feBd(), M_FESpace[i]->dof() );
        M_bch[i]->setOffset(M_offset[i]);
    }
    bcManageMatrix( *M_blocks[i] , *M_FESpace[i]->mesh(), M_FESpace[i]->dof(), *M_bch[i], M_FESpace[i]->feBd(), 2., time);
}


void ComposedNN::push_back_matrix(const matrixPtr_Type& Mat, const  bool recompute)
{
    Mat->GlobalAssemble();
    *Mat *= 2.;
    super::push_back_matrix(Mat, recompute);
}



void ComposedNN::replace_matrix( const matrixPtr_Type& oper, UInt position )
{
    oper->GlobalAssemble();
    *oper *= 2.;
    M_blocks[position]=oper;
    M_blocks[position]=oper;
}

} // Namespace LifeV
