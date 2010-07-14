/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
             Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-05-18

  Copyright (C) 2010 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file ethiersteiman.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-05-18
 */

#ifndef __Ethiersteinman_H
#define __Ethiersteinman_H 1

#include <life/lifecore/application.hpp>
#include <life/lifesolver/Oseen.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifemesh/basisElSh.hpp>

#include "ESSteady_function.hpp"
#include "ESUnsteady_function.hpp"

enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class Ethiersteinman
 * \brief 2D/3D Ethiersteinman Simulation class
 *
 *  @author Christophe Prud'homme
 *  @author Gwenol Grandperrin
 *  @see
 */

class Ethiersteinman
//     :
//     public LifeV::Application
{
public:
    typedef LifeV::RegionMesh3D<LifeV::LinearTetra>       mesh_type;
    typedef LifeV::FESpace< mesh_type, LifeV::EpetraMap > fespace_type;
    typedef LifeV::Oseen< mesh_type >                     fluid_type;
    typedef fluid_type::vector_type                       vector_type;
    typedef boost::shared_ptr<vector_type>                vector_ptrtype;
    typedef fluid_type::matrix_type                       matrix_type;

    // Problem definition
    typedef LifeV::EthierSteinmanUnsteady Problem;



    /** @name Typedefs
     */
    //@{

//    typedef LifeV::Application super;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Ethiersteinman( int argc,
              char** argv,
              LifeV::AboutData const& ad,
              LifeV::po::options_description const& od );

    ~Ethiersteinman()
        {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{
    //! Computes L2 errors

    void check();
    void run();

    //@}


private:

    struct RESULT_CHANGED_EXCEPTION
    {
    public:
        RESULT_CHANGED_EXCEPTION()
        {
            std::cout << "Some modifications led to changes in the l2 norm of the solution" << std::endl;
        }
    };


    struct Private;
    boost::shared_ptr<Private> d;
    std::ofstream              out_norm;

    // Tolerance of the test (should be <1)
    // Actually for convTol=1, the test failed
    // if the improvement of accuracy is less
    // than predicted by the theory.
    // convTol lower down the theoretical bounds
    LifeV::Real convTol;

};

#endif /* __Ethiersteinman_H */
