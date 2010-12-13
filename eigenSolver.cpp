/* -*- mode: c++ -*-
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

#include <lifemc/lifealg/eigenSolver.hpp>

#ifdef HAVE_TRILINOS_ANASAZI

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>

#include "Teuchos_RefCountPtrDecl.hpp"
#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// ===================================================
//! Constructors & Destructors
// ===================================================

namespace LifeV
{
EigenSolver::EigenSolver(boost::shared_ptr<Solver> const matrix, Epetra_BlockMap const& block_map, long unsigned int numvec):
        MyProblem( new Anasazi::BasicEigenproblem<DataType, Vector, Solver>(Teuchos::rcp(matrix), M_eigenVectors) ),
        M_eigenVectors(new Epetra_MultiVector(block_map, numvec)),
        MyPL(),
        MySolver()
{
}


// ===================================================
//! Public Methods
// ===================================================

void
EigenSolver::setDataFromGetPot(GetPot const& data , const std::string& section)
{
    int block_size = data((section + "block_size").c_str(), 1);
    int max_blocks = data((section + "max_blocks").c_str(), 20);
    int max_restarts = data((section + "max_restarts").c_str(), 50);
    double tol = data((section + "tol").c_str(), 1e-8);
    std::string which = data((section + "which").c_str(), "ML");
    int nev = data((section + "neval").c_str(), 10);//number of eigenvalues
    bool verbose = data((section + "verbose").c_str(), true);
    bool debug = data((section + "debug").c_str(), true);


    int verbosity = Anasazi::Errors + Anasazi::Warnings;
    if (verbose)
    {
        verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
    }
#ifdef DEBUG
    if (debug)
        verbosity += Anasazi::Debug;
#endif
    //
    // Create parameter list to pass into solver manager
    //

    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Block Size", block_size );
    MyPL.set( "Max Blocks", max_blocks );
    MyPL.set( "Max Restarts", max_restarts );
    MyPL.set( "Tol", tol );
    //Teuchos::RCP<Anasazi::SortManager<Real> > MySort = Teuchos::rcp( new Anasazi::BasicSort<Real>( which.c_str() ) );
    //MyPL.set( "Sort Manager", MySort);

    MyProblem->setNEV(nev);
}


void
EigenSolver
::eigenvalues(std::vector<DataType>& realPart, std::vector<DataType>& imgPart)
{
    Anasazi::Eigensolution<DataType,Vector> sol = MyProblem->getSolution();
    std::vector<Anasazi::Value<DataType> > evals = sol.Evals;
    for (UInt i=0; i<evals.size(); ++i)
    {
        realPart.push_back(evals[i].realpart);
        imgPart.push_back(evals[i].imagpart);
    }
}

int
EigenSolver
::solve()
{
    bool out=MyProblem->setProblem();
    std::cout<<"problem set: "<<out<<std::endl;
    MySolver.reset(new eigensolver_raw_type(MyProblem, MyPL));
    return MySolver->solve();
}
}
#endif
