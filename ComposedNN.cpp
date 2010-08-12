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
    @file ComposedNN.cpp

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 29 Jun 2010
 */

#include <ComposedNN.hpp>

#include <life/lifefem/bcManage.hpp>

namespace LifeV
{

int ComposedNN::solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver )
{
    typedef  ComposedPreconditioner<Ifpack_Preconditioner> composed_prec;
    boost::shared_ptr< composed_prec > firstCompPrec(new composed_prec(M_comm));
    boost::shared_ptr< composed_prec > secondCompPrec(new composed_prec(M_comm));

    Chrono chrono;
    int overlapLevel     = M_list.get("overlap level", 2);
    std::string precType = M_list.get("prectype", "Amesos");
    Ifpack factory;

    if( !M_blockPrecs.get() )
        M_blockPrecs.reset(new ComposedPreconditioner< composed_prec >(M_comm));
    if( !(M_blockPrecs->getNumber()))
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
    }
    else
    {
        for(ID k(0); k < M_blocks.size(); ++k)
        {
            if(M_recompute[k])
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
            else
                M_blockPrecs->displayer().leaderPrint("  M-  Reusing double prec. factorization ...        \n");
        }
    }
        firstCompPrec->push_back(M_prec[0], false, false);
        firstCompPrec->push_back(M_prec[1], false, false);

        secondCompPrec->push_back(M_prec[2], false, false);
        secondCompPrec->push_back(M_prec[3], false, false);

        //M_blockPrecs->resetP();
    if(!(M_blockPrecs->getNumber()))
    {
        M_blockPrecs->push_back(firstCompPrec, false, false, false);
        M_blockPrecs->push_back(secondCompPrec, false, false, true /*sum*/);
    }
    else
    {
        M_blockPrecs->replace(firstCompPrec, (UInt)0, false, false);
        M_blockPrecs->replace(secondCompPrec, (UInt)1, false, false);
    }

    return linearSolver->solveSystem(rhs, step, boost::static_pointer_cast<Epetra_Operator>(M_blockPrecs));
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
    bool tmpRecompute = M_recompute[1];
    M_recompute[1] = M_recompute[0];
    M_recompute[0] = tmpRecompute;

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
    couplingMatrix(coupling, 4, tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);

    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal(one, M_offset[1]+1, M_offset[0]+1 );
    coupling->insertValueDiagonal(one,  fluidSolid, totalDofs);
    couplingMatrix(coupling, 8,  tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);

    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal( one, M_offset[0]+1, fluidSolid );
    coupling->insertValueDiagonal( one, fluidSolid, totalDofs );
    couplingMatrix(coupling, 2, tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);

    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal(one, M_offset[1]+1, M_offset[0]+1 );
    coupling->insertValueDiagonal(-1,  fluidSolid, totalDofs);
    couplingMatrix(coupling, 1,  tmpFESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);

    M_coupling.push_back(coupling);

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



void ComposedNN::replace_matrix( const matrix_ptrtype& oper, UInt position )
{
    oper->GlobalAssemble();
    *oper *= 2.;
    M_blocks[1-position]=oper;
    M_blocks[3-position]=oper;
}

} // Namespace LifeV
