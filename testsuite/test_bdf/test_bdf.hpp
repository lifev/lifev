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
    @brief Simple Fourier test with Dirichlet Boundary condition

    @author Umberto Villa <uvilla@emory.edu>
    @date 14-04-2010

     Solve the problem

           \partial_t u - \nu(t) \Delta u + sigma(t) u = f(t)

                    u = g on the boundary
                        u(t=0) = u0 initial condition

            on a cube
            \nu, \sigma and \source can be function of time
            (which implies that the matrix needs to be reassembled each time)

    Purpose: Test BDF of different order

    The analytical solution used in the test is quadratic both in space and time.
    Therefore a second order approximation in space + BDF2 should provide a solution exact up to
    the tolerance of the linear system.

 */


#ifndef __test_bdf_H
#define __test_bdf_H 1





// ===================================================
//! Includes
// ===================================================
#include <life/lifecore/life.hpp>





/*!
 * \class laplacian
 * \brief LifeV Laplacian test case
 *
 *  @author U. Villa
 *  @see
 */
class test_bdf
//     :
//     public LifeV::Application
{
public:

    /** @name Constructors, destructor
     */
    //@{

    test_bdf( int argc,
              char** argv );

    ~test_bdf()
    {}

    //@}

    /** @name  Methods
     */
    //@{

    void run();
    bool check();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> Members;
};

#endif /* __test_bdf_H */
