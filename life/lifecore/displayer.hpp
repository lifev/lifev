/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Paolo Crosetto <crosetto@iacspc70.epfl.ch>
       Date: 2009-03-02

  Copyright (C) 2009 EPFL, INRIA, Politecnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/**
   \file MonolithicSolid.cpp
   \author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   \date 2009-03-02
 */

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#ifndef _DISPLAYER_H_
#define _DISPLAYER_H_

namespace LifeV
{

/*!
  \class Displayer
  \brief
  This class is used to display messages in parallel simulations: if a communicator is passed to the constructor
  only one processor (the leader) will print out the message. If no communicator is passed to the constructor every processor prints the messages.
*/
class Displayer
{
public:
    Displayer( Epetra_Comm& comm);
    Displayer();
    virtual ~Displayer(){}

    void    leaderPrint(string const message, double const number) const;
    /*! to print one message and one real value
      \param message message to print out
      \param number Real value that we want to print
    */


    void    leaderPrint(string const message) const;

    void    leaderPrintMax(string const message, double const number) const;
    /*!
      Take a Real input value from all processors in the communicator, computes the max, returns the max to all processors of the communicator. Then processor 0 of the communicator prints it.
      \param message message to print out
      \param number Real value that we want to print
    */
    bool    isLeader() const
    {
        if( M_comm  != 0)
            return M_comm->MyPID() == 0;
        else
            return true;
    }
    /*!
      returns the processor 0 of the communicator
    */
    const Epetra_Comm& comm() const {return *M_comm;}

protected:

    Epetra_Comm*                    M_comm;
    int                             M_me;

};

}
#endif
