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
 *  @file
 *  @brief This class implements the Linear St. Venant-kirchhoff solver.
 *
 *  @date 01-09-2010
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
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

    //! @name Type definitions
    //@{

    typedef VenantKirchhofSolver<Mesh, SolverType>     super;
    typedef typename super::vector_Type                vector_Type;
    typedef typename super::matrixPtr_Type             matrixPtr_Type;
    typedef typename super::bchandler_Type             bchandler_Type;

    //@}

    //! @name Constructor
    //@{
    
    LinearVenantKirchhofSolver():
            super()
    {}

    //@}

    //! @name Methods
    //@{
    /**
       These methods are derived from the base class VenantKirchhoff 
     */
    void updateJacobian( vector_Type& /*sol*/, matrixPtr_Type& /*jac*/ )
    {
        this->M_Displayer->leaderPrint("  Linear S-  Doing nothing (updating jacobian of a linear system) ...                    ");
    }


    //! solves the tangent problem for newton iterations
    void solveJac( vector_Type&       /*step*/,
                   const vector_Type& /*res*/,
                   Real&              /*linear_rel_tol*/) {assert(false);}

  void solveJacobian( vector_Type&       /*step*/,
                        const vector_Type& /*res*/,
                        Real&            /*linear_rel_tol*/,
                        bchandler_Type&    /*BCd*/ ) {assert(false);}

    //@}
};


} // Namespace LifeV

#endif /* LINEARVENANTKIRCHHOFSOLVER_H */
