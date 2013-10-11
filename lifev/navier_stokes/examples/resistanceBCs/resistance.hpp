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
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-19
 */

#ifndef __RESISTANCE_H
#define __RESISTANCE_H 1


/*!
 * \class ResistanceTest
 * \brief 2D/3D Cylinder Simulation class
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class ResistanceTest
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

    ResistanceTest ( int argc,
               char** argv );

    ~ResistanceTest()
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
    boost::shared_ptr<Private> parameters;
};
#endif /* __RESISTANCE_H */
