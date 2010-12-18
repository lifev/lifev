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

#ifndef __Structure_H
#define __Structure_H 1

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif

enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class Structure
 * \brief 2D/3D Structure Simulation class
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class Structure
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

    Structure( int                                   argc,
               char**                                argv,
               boost::shared_ptr<Epetra_Comm>        structComm );

    ~Structure()
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

    void run()
    {
//             if ( vm().count( "help" ) )
//             {
//                 std::cout << optionsDescription() << "\n";
//                 return;
//             }
        run3d();
        void CheckResults(const LifeV::Real& dispNorm, const LifeV::Real& time);
    }

    //@}



protected:

private:

    /**
     * run the 2D driven cylinder simulation
     */
    void run2d();

    /**
     * run the 3D driven cylinder simulation
     */
    void run3d();
    void CheckResults(const LifeV::Real& dispNorm, const LifeV::Real& time);
    LifeV::UInt RESULT_CHANGED_EXCEPTION(const LifeV::Real time);

private:
    struct Private;
    boost::shared_ptr<Private> parameters;
};

#endif /* __Structure_H */
