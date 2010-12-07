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
 *  @brief File containing the BCInterface_OperatorFunctionFile class
 *
 *  @date 26-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface_OperatorFunctionFile_H
#define BCInterface_OperatorFunctionFile_H 1

#include <lifemc/lifesolver/BCInterface_FunctionFile.hpp>
#include <lifemc/lifesolver/BCInterface_OperatorFunction.hpp>

namespace LifeV
{

//! BCInterface_OperatorFunctionFile - LifeV bcFunction wrapper for BCInterface (with Operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface, SpiritParser and and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *  The function string can contain Operator parameters.
 */
template< class physicalSolver_Type >
class BCInterface_OperatorFunctionFile: public virtual BCInterface_FunctionFile< physicalSolver_Type > ,
        public virtual BCInterface_OperatorFunction< physicalSolver_Type >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface_Function< physicalSolver_Type >             super0;
    typedef BCInterface_FunctionFile< physicalSolver_Type >         super1;
    typedef BCInterface_OperatorFunction< physicalSolver_Type >     super2;

    typedef BCInterface_Data                                        data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface_OperatorFunctionFile();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_OperatorFunctionFile( const data_Type& data );

    //! Destructor
    virtual ~BCInterface_OperatorFunctionFile() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const data_Type& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface_OperatorFunctionFile( const BCInterface_OperatorFunctionFile& function );

    virtual BCInterface_OperatorFunctionFile& operator=( const BCInterface_OperatorFunctionFile& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename physicalSolver_Type >
inline BCInterface_Function< physicalSolver_Type >* createBCInterface_OperatorFunctionFile()
{
    return new BCInterface_OperatorFunctionFile< physicalSolver_Type > ();
}

// ===================================================
// Constructors
// ===================================================
template< class physicalSolver_Type >
BCInterface_OperatorFunctionFile< physicalSolver_Type >::BCInterface_OperatorFunctionFile() :
        super0      (),
        super1      (),
        super2      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface_OperatorFunctionFile::BCInterface_OperatorFunctionFile()" << "\n";
#endif

}

template< class physicalSolver_Type >
BCInterface_OperatorFunctionFile< physicalSolver_Type >::BCInterface_OperatorFunctionFile( const data_Type& data ) :
        super0      (),
        super1      (),
        super2      ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface_OperatorFunctionFile::BCInterface_OperatorFunctionFile( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< class physicalSolver_Type >
void
BCInterface_OperatorFunctionFile< physicalSolver_Type >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5024 ) << "BCInterface_OperatorFunctionFile::setData" << "\n";
#endif
    super1::setData( data );

    //super2::setData( data ); Cannot call directly, because it call again BCInterface_Function::setup( data )
    super2::M_flag = data.flag();

    super2::createAccessList( data );
}

} // Namespace LifeV

#endif /* BCInterface_OperatorFunctionFile_H */
