/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-23

 Copyright (C) 2009 EPFL

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
 \file MS_Algorithm_Aitken.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-23
 */

#ifndef __MS_Algorithm_Aitken_H
#define __MS_Algorithm_Aitken_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <life/lifealg/generalizedAitken.hpp>

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV {

//! MS_Algorithm_Aitken - The MultiScale Algorithm implementation of Aitken
/*!
 *  The MS_Algorithm_Aitken is an implementation of MS_Algorithm
 *  which implements the Aitken method.
 *
 *  @author Cristiano Malossi
 */
class MS_Algorithm_Aitken : public virtual MS_Algorithm
{
public:

    typedef MS_Algorithm                  super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm_Aitken();

    //! Copy constructor
    /*!
     * \param algorithm - MS_Algorithm_Aitken
     */
    MS_Algorithm_Aitken( const MS_Algorithm_Aitken& algorithm );

    //! Destructor
    ~MS_Algorithm_Aitken() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param algorithm - MS_Algorithm_Aitken
     */
    MS_Algorithm_Aitken& operator=( const MS_Algorithm_Aitken& algorithm );

    //@}


    //! @name MultiScale PhysicalModel Virtual Functions
    //@{

    //! Setup the data of the algorithm
    void SetupData( const GetPot& DataFile );

    //! Perform sub-iteration on the couplings
    void SubIterate( void );

    //! Display some information about the algorithm
    void ShowMe( void );

    //@}

protected:

    //! @name Protected Methods
    //@{

    //@}

    enum methodType
    {
        Scalar, Vectorial, VectorialBlock
    };

    std::map< std::string, methodType >            M_methodMap;
    methodType                                     M_method;
    bool                                           M_inverseOmega;
    generalizedAitken< VectorType >                M_generalizedAitken;

};

//! Factory create function
inline MS_Algorithm* createAitken()
{
    return new MS_Algorithm_Aitken();
}

} // Namespace LifeV

#endif /* __MS_Algorithm_Aitken_H */
