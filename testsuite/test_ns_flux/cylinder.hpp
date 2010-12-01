/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-19

  Copyright (C) 2005 EPFL

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
   \file cylinder.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2008-06-13
 */

#ifndef __Cylinder_H
#define __Cylinder_H 1


enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class Cylinder
 * \brief 3D Cylinder Simulation with fluxes at inflow
 *
 *  @author Simone Deparis
 *  @see
 */
class Cylinder
//     :
//     public LifeV::Application
{
public:


    /** @name Typedefs
     */
    //@{

//    typedef LifeV::Application super;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Cylinder( int argc,
              char** argv );

    ~Cylinder()
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

    void run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> d;
};

#endif /* __Cylinder_H */
