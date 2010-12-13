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
