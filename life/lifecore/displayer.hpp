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
   \file displayer.hpp

   \version 1.0
   \date 2009-03-02
   \author Paolo Crosetto <crosetto@iacspc70.epfl.ch>

   \version 1.4
   \date 2009-03-02
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>

   \brief Template input variables for more general output messages:
   now it is possible to display not only strings but also Int, Real, etc..
 */

#ifndef _DISPLAYER_H_
#define _DISPLAYER_H_

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/life.hpp>
#include <boost/shared_ptr.hpp>

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
    Displayer( const boost::shared_ptr<Epetra_Comm>& comm=boost::shared_ptr<Epetra_Comm>() );
    Displayer( const Displayer& displayer );
	Displayer( const bool verbose);
    virtual ~Displayer() {}

    /*! To print one message.
      \param message1 message to print out
    */
    template <typename T1>
    void leaderPrint( const T1& message1 ) const;

    /*! To print two messages.
      \param message1 first message to print out
      \param message2 second message to print out
    */
    template <typename T1, typename T2>
    void leaderPrint( const T1& message1, const T2& message2 ) const;

    /*! To print three messages.
      \param message1 first message to print out
      \param message2 second message to print out
      \param message3 third message to print out
    */
    template <typename T1, typename T2, typename T3>
    void leaderPrint( const T1& message1, const T2& message2, const T3& message3 ) const;

    /*!
      Take a Real input value from all processors in the communicator, computes the max, returns the max to all processors of the communicator.
      Then processor 0 of the communicator prints it.
      \param message1 message to print out
      \param localMax Int or Real local maximum value that we want to print
    */
    template <typename T1>
    void leaderPrintMax( const T1& message1, const Real& localMax ) const;

    /*!
      Return true if it is process 0 of the communicator
    */
    inline bool isLeader() const
  {
      return M_verbose;
  }

    //Set the communicator
    void SetCommunicator( const boost::shared_ptr<Epetra_Comm>& comm );

    /*!
      Return the communicator
    */
    const boost::shared_ptr<Epetra_Comm>& comm() const { return M_comm; }

protected:

    boost::shared_ptr<Epetra_Comm>					    M_comm;
    bool	           									M_verbose;

};



// IMPLEMENTATION



template <typename T1>
void Displayer::
leaderPrint( const T1& message1 ) const
{
	if ( M_verbose )
		std::cout << message1 << std::flush;
}

template <typename T1, typename T2>
void Displayer::
leaderPrint( const T1& message1, const T2& message2 ) const
{
	if ( M_verbose )
		std::cout << message1 << message2 << std::flush;
}

template <typename T1, typename T2, typename T3>
void Displayer::
leaderPrint( const T1& message1, const T2& message2, const T3& message3 ) const
{
	if ( M_verbose )
		std::cout << message1 << message2 << message3 << std::flush;
}

template <typename T1>
void Displayer::
leaderPrintMax( const T1& message1, const Real& localMax ) const
{
	if ( M_comm.get() )
    {
        Real num( localMax );
        Real globalMax;

		M_comm->MaxAll( &num, &globalMax, 1 );
        if ( M_verbose )
            std::cout << message1 << globalMax << std::endl;
    }
    else
        std::cout << message1 << localMax << std::endl;
}

} // Namespace LifeV

#endif
