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
    @brief 2D/3D Cylinder Simulation class

    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @date 19-04-2005

 */



#ifndef __Cylinder_H
#define __Cylinder_H 1

enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class Cylinder
 * \brief 2D/3D Cylinder Simulation class
 *
 *  @author Christophe Prud'homme
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
