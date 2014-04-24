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
 @brief Base class for applying cardiac stimulus

 @date 11-2013
 @author Toni Lassila <toni.lassila@epfl.ch>

 @last update 11-2013
 */


#ifndef ELECTROSTIMULUS_HPP_
#define ELECTROSTIMULUS_HPP_

#include <lifev/core/array/VectorEpetra.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
namespace LifeV
{

class ElectroStimulus
{

public:

    //! @name Type definitions
    //@{
    typedef VectorEpetra                    vector_Type;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
    typedef Teuchos::ParameterList          list_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    ElectroStimulus();

    //! Destructor
    virtual ~ElectroStimulus() {};

    //@}

    //! @name Get Methods
    //@{

    //@}

    //! @name Set Methods
    //@{


    //@}

    //! @name Copy Methods
    //@{

    //@}

    //! @name Methods
    //@{
    inline virtual Real appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
    {
        return 0.0;
    }

    virtual void setParameters (list_Type&  list)
    {
	
    }

    virtual void showMe ()
    {

    }

    //@}

private:

};


} // namespace LifeV

#endif /* ELECTROSTIMULUS_HPP_ */
