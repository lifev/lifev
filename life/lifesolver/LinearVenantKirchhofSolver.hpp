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
    @date 01 Sep 2010

    A more detailed description of the file (if necessary)
 */

#ifndef LINEARVENANTKIRCHHOFSOLVER_H
#define LINEARVENANTKIRCHHOFSOLVER_H 1

#include <life/lifecore/life.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>

namespace LifeV
{

template <typename Mesh, typename SolverType = LifeV::SolverTrilinos >
class LinearVenantKirchhofSolver : public VenantKirchhofSolver<Mesh, SolverType>
{
public:

    typedef VenantKirchhofSolver<Mesh, SolverType> super;
    typedef typename super::vector_type vector_type;
    typedef typename super::matrix_ptrtype matrix_ptrtype;
    typedef typename super::bchandler_type bchandler_type;

    LinearVenantKirchhofSolver():
            super()
    {}

    void updateJacobian( vector_type& /*sol*/, matrix_ptrtype& /*jac*/ )
    {
        this->M_Displayer->leaderPrint("  Linear S-  Doing nothing (updating jacobian of a linear system) ...                    ");
    }


    //! solves the tangent problem for newton iterations
    void solveJac( vector_type&       /*step*/,
                   const vector_type& /*res*/,
                   Real&              /*linear_rel_tol*/) {assert(false);}

    void solveJacobian( vector_type&       /*step*/,
                        const vector_type& /*res*/,
                        Real&            /*linear_rel_tol*/,
                        bchandler_type&    /*BCd*/ ) {assert(false);}

};


} // Namespace LifeV

#endif /* LINEARVENANTKIRCHHOFSOLVER_H */
