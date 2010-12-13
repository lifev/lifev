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
 *  @file
 *  @brief File containing the CommunicatorEpetra class
 *
 *  @date 2010-12-13
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Tiziano Passerini <tiziano.passerini@gmail.com>
 *  @contributor Umberto Villa <uvilla@emory.edu>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef COMMUNICATOREPETRA_H
#define COMMUNICATOREPETRA_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
    #include <mpi.h>
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV classes
#include <life/lifecore/life.hpp>

namespace LifeV
{

//! CommunicatorEpetra - This class is a wrapper to Epetra_Comm class.
/*!
 * @author Paolo Crosetto, Cristiano Malossi
 *
 * In addition it also has some functionality to
 * If a communicator is passed to the constructor only one processor (the leader) will print out the message.
 * If no communicator is passed to the constructor every processor prints the messages.
 */
class CommunicatorEpetra
{
public:

    //! @name Public Types
    //@{

#ifdef EPETRA_MPI
    typedef Epetra_MpiComm                           communicator_Type;
#else
    typedef Epetra_SerialComm                        communicator_Type;
#endif
    typedef boost::shared_ptr< communicator_Type >   communicatorPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param communicator shared_ptr to communicator
     */
    explicit CommunicatorEpetra( const communicatorPtr_Type& communicator = communicatorPtr_Type() );

    //! Copy constructor
    /*!
     * @param CommunicatorEpetra CommunicatorEpetra
     */
    CommunicatorEpetra( const CommunicatorEpetra& communicatorEpetra );

    //! Destructor
    virtual ~CommunicatorEpetra() {}

    //@}


    //! @name Methods
    //@{

    //! Print one message.
    /*!
     * @param message1 message to print out
     */
    template < typename Message1Type >
    void leaderPrint( const Message1Type& message1 ) const;

    //! Print two messages.
    /*!
     * @param message1 first message to print out
     * @param message2 second message to print out
     */
    template < typename Message1Type, typename Message2Type >
    void leaderPrint( const Message1Type& message1, const Message2Type& message2 ) const;

    //! Print three messages.
    /*!
     * @param message1 first message to print out
     * @param message2 second message to print out
     * @param message3 third message to print out
     */
    template < typename Message1Type, typename Message2Type, typename Message3Type >
    void leaderPrint( const Message1Type& message1, const Message2Type& message2, const Message3Type& message3 ) const;

    //! Print the maximum value among the processors
    /*!
     * Take a Real input value from all processors in the communicator, computes the max,
     * returns the max to all processors of the communicator.
     * Then processor 0 of the communicator prints it.
     * @param message1 message to print out
     * @param localMax Int or Real local maximum value that we want to print
     */
    template < typename Message1Type >
    void leaderPrintMax( const Message1Type& message1, const Real& localMax ) const;

    //! Determine if it is the leader
    /*!
     * @return true if it is process 0 of the communicator
     */
    const bool& isLeader() const { return M_verbose; }

    //@}


    //! @name Set Methods
    //@{

    //! Set the communicator
    /*!
     * @param communicator the communicator
     */
    void setCommunicator( const communicatorPtr_Type& communicator );

    //@}


    //! @name Get Methods
    //@{

    //! Get the const communicator
    /*!
     * @return the const communicator
     */
    const communicatorPtr_Type& communicator() const { return M_communicator; }

    //! Get the communicator
    /*!
     * @return the communicator
     */
    communicatorPtr_Type& communicator()  { return M_communicator; }

    //@}

private :

    communicatorPtr_Type				M_communicator;
    bool	           					M_verbose;
};

// ===================================================
// Template implementation
// ===================================================
template < typename Message1Type >
void
CommunicatorEpetra::leaderPrint( const Message1Type& message1 ) const
{
	if ( M_verbose )
		std::cout << message1 << std::flush;
}

template < typename Message1Type, typename Message2Type >
void
CommunicatorEpetra::leaderPrint( const Message1Type& message1, const Message2Type& message2 ) const
{
	if ( M_verbose )
		std::cout << message1 << message2 << std::flush;
}

template < typename Message1Type, typename Message2Type, typename Message3Type >
void
CommunicatorEpetra::leaderPrint( const Message1Type& message1, const Message2Type& message2, const Message3Type& message3 ) const
{
	if ( M_verbose )
		std::cout << message1 << message2 << message3 << std::flush;
}

template < typename Message1Type >
void
CommunicatorEpetra::leaderPrintMax( const Message1Type& message1, const Real& localMax ) const
{
	if ( M_communicator.get() )
    {
        Real num( localMax );
        Real globalMax;

		M_communicator->MaxAll( &num, &globalMax, 1 );
        if ( M_verbose )
            std::cout << message1 << globalMax << std::endl;
    }
    else
        std::cout << message1 << localMax << std::endl;
}

} // Namespace LifeV

#endif // COMMUNICATOREPETRA_H
