//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2009-2010 EPFL, Politecnico di Milano

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
/*!
 * @file
 *
 * @brief Template input variables for more general output messages:
 * now it is possible to display not only strings but also Int, Real, etc..
 *
 * @date 2009-03-02
 * @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 * @maintainer Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef DISPLAYER_H
#define DISPLAYER_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/life.hpp>

namespace LifeV
{

//! Displayer - This class is used to display messages in parallel simulations.
/*!
 * @author Paolo Crosetto, Cristiano Malossi
 *
 * If a communicator is passed to the constructor only one processor (the leader) will print out the message.
 * If no communicator is passed to the constructor every processor prints the messages.
 */
class Displayer
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_Comm                              comm_Type;
    typedef boost::shared_ptr< comm_Type >           commPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    Displayer () {}

    Displayer( const commPtr_Type& comm = commPtr_Type() );

    //! Copy constructor
    /*!
     * @param displayer Displayer
     */
    Displayer( const Displayer& displayer );

    //! Destructor
    virtual ~Displayer() {}

    //@}


    //! @name Methods
    //@{

    //! Print one message.
    /*!
     * @param message1 message to print out
     */
    template <typename T1>
    void leaderPrint( const T1& message1 ) const;

    //! Print two messages.
    /*!
     * @param message1 first message to print out
     * @param message2 second message to print out
     */
    template <typename T1, typename T2>
    void leaderPrint( const T1& message1, const T2& message2 ) const;

    //! Print three messages.
    /*!
     * @param message1 first message to print out
     * @param message2 second message to print out
     * @param message3 third message to print out
     */
    template <typename T1, typename T2, typename T3>
    void leaderPrint( const T1& message1, const T2& message2, const T3& message3 ) const;

    //! Print the maximum value among the processors
    /*!
     * Take a Real input value from all processors in the communicator, computes the max,
     * returns the max to all processors of the communicator.
     * Then processor 0 of the communicator prints it.
     * @param message1 message to print out
     * @param localMax Int or Real local maximum value that we want to print
     */
    template <typename T1>
    void leaderPrintMax( const T1& message1, const Real& localMax ) const;

    //! Determine if it is the leader
    /*!
     * @return true if it is process 0 of the communicator
     */
    const bool& isLeader() const;

    //! @name Set Methods
    //@{

    //! Set the communicator
    /*!
     * @param comm the communicator
     */
    void setCommunicator( const commPtr_Type& comm );

    //! @name Get Methods
    //@{

    //! Get the communicator
    /*!
     * @return the communicator
     */
    const commPtr_Type& comm() const
    {
        return M_comm;
    }

    //@}

protected:

    commPtr_Type					    M_comm;
    bool	           					M_verbose;

};


// ===================================================
// Template implementation
// ===================================================
template <typename T1>
void
Displayer::leaderPrint( const T1& message1 ) const
{
    if ( M_verbose )
        std::cout << message1 << std::flush;
}

template <typename T1, typename T2>
void
Displayer::leaderPrint( const T1& message1, const T2& message2 ) const
{
    if ( M_verbose )
        std::cout << message1 << message2 << std::flush;
}

template <typename T1, typename T2, typename T3>
void
Displayer::leaderPrint( const T1& message1, const T2& message2, const T3& message3 ) const
{
    if ( M_verbose )
        std::cout << message1 << message2 << message3 << std::flush;
}

template <typename T1>
void
Displayer::leaderPrintMax( const T1& message1, const Real& localMax ) const
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

#endif // DISPLAYER_H
