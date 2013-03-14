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
 *  @brief File containing the BCInterfaceFunctionParserFile class
 *
 *  @date 09-07-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */
#ifndef BCInterfaceFunctionParserFile_H
#define BCInterfaceFunctionParserFile_H 1

#include <lifev/bc_interface/core/function/BCInterfaceFunctionParser.hpp>

namespace LifeV
{

//! BCInterfaceFunctionParserFile - LifeV boundary condition function file wrapper for \c BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between the \c BCInterface and the \c Parser. It allows to construct LifeV
 *  functions type for boundary conditions, using a \c GetPot file containing a function string and a
 *  table of discrete data (for example a discrete flow rate or pressure as a function of the time).
 *
 *  See \c BCInterfaceFunctionParser class for more details.
 *
 *  <b>DETAILS</b> <BR>
 *  The constructor of the class takes a string contains the \c GetPot file name.
 *  The \c GetPot file has the following structure:
 *
 *  <ul>
 *      <li> <b>function:</b> contains the expression of the function (as described in the \c BCInterfaceFunctionParser class).
 *      <li> <b>variables:</b> contains the list of variables and coefficients present in the function.
 *                             The first one is the variable and should be sorted in a growing order,
 *                             while all the others are coefficients.
 *      <li> <b>data:</b> contains the discrete list of variable and related coefficients. It must
 *                        have n columns, where n is the number of variables (and coefficients).
 *      <li> <b>scale:</b> Contains n optional coefficients (default = 1), which multiply each column in the data table.
 *      <li> <b>loop:</b> Useful for periodic simulation: at the end of the list, it restarts using the first value instead of extrapolating the last two.
 *  </ul>
 *
 *  <b>NOTE</b> <BR>
 *  During the execution, if the value of the variable (usually the time) is not present in the 'data' table,
 *  the class linearly interpolates the value between the two closest values. Moreover, if the value of the variable is higher
 *  than anyone present in the 'data' table, the class linearly extrapolates the value using the last two values in the table.
 *
 *  <b>EXAMPLE OF DATA FILE</b> <BR>
 *  <CODE>
 *  function  = '(0,0,q)'                    <BR>
 *  loop      = false                        <BR>
 *  variables = 't                q'         <BR>
 *  scale     = '1                1'         <BR>
 *  data      = '0.000000000      1.00       <BR>
 *               0.333333333      2.00       <BR>
 *               0.666666666      3.00       <BR>
 *               1.000000000      4.00'      <BR>
 *  </CODE>
 */
template< typename BcHandlerType, typename PhysicalSolverType >
class BCInterfaceFunctionParserFile: public virtual BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef BcHandlerType                                                          bcHandler_Type;
    typedef PhysicalSolverType                                                     physicalSolver_Type;

    typedef BCInterfaceFunction< bcHandler_Type, physicalSolver_Type >             function_Type;
    typedef BCInterfaceFunctionParser< bcHandler_Type, physicalSolver_Type >       functionParser_Type;

    typedef typename function_Type::data_Type                                      data_Type;
    typedef typename function_Type::dataPtr_Type                                   dataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    explicit BCInterfaceFunctionParserFile();

    //! Destructor
    virtual ~BCInterfaceFunctionParserFile() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data for boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    virtual void setData ( const dataPtr_Type& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionParserFile ( const BCInterfaceFunctionParserFile& function );

    BCInterfaceFunctionParserFile& operator= ( const BCInterfaceFunctionParserFile& function );

    //@}


    //! @name Private methods
    //@{

    //! Linear interpolation (extrapolation) between two values of the data.
    void dataInterpolation();

    //@}

    std::vector< std::string >                   M_variables;
    bool                                         M_loop;
    std::map< std::string, std::vector< Real > > M_data;
    std::vector< Real >::iterator                M_dataIterator;
};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename BcHandlerType, typename PhysicalSolverType >
inline BCInterfaceFunctionParser< BcHandlerType, PhysicalSolverType >* createBCInterfaceFunctionParserFile()
{
    return new BCInterfaceFunctionParserFile< BcHandlerType, PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
BCInterfaceFunctionParserFile< BcHandlerType, PhysicalSolverType >::BCInterfaceFunctionParserFile() :
    function_Type                    (),
    functionParser_Type              (),
    M_variables                      (),
    M_loop                           (),
    M_data                           (),
    M_dataIterator                   ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5022 ) << "BCInterfaceFunctionFile::BCInterfaceFunctionFile()" << "\n";
#endif

}



// ===================================================
// Set Methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
inline void
BCInterfaceFunctionParserFile< BcHandlerType, PhysicalSolverType >::setData ( const dataPtr_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5022 ) << "BCInterfaceFunctionFile::loadData            fileName: " << data->baseString() << "\n";
#endif

    // Create a true copy
    dataPtr_Type dataCopy ( new data_Type ( *data ) );

    std::vector< std::string > stringsVector;
    boost::split ( stringsVector, dataCopy->baseString(), boost::is_any_of ( "[" ) );

    //Load data from file
    GetPot dataFile ( stringsVector[0] );

    //Set variables
    UInt variablesNumber = dataFile.vector_variable_size ( "variables" );

    M_variables.clear();
    M_variables.reserve ( variablesNumber );

    std::vector< Real > scale;
    scale.reserve ( variablesNumber );

    for ( UInt j ( 0 ); j < variablesNumber; ++j )
    {
        M_variables.push_back ( dataFile ( "variables", "unknown", j ) );
        scale.push_back ( dataFile ( "scale", 1.0, j ) );
    }

#ifdef HAVE_LIFEV_DEBUG
    std::stringstream output;
    output << "BCInterfaceFunctionFile::loadData           variables: ";
    for ( UInt j (0); j < variablesNumber; ++j )
    {
        output << M_variables[j] << "  ";
    }

    output << "\n                                                           scale: ";
    for ( UInt j (0); j < variablesNumber; ++j )
    {
        output << scale[j] << "  ";
    }

    debugStream ( 5022 ) << output.str() << "\n";
#endif

    //Load loop flag
    M_loop = dataFile ( "loop", false );

    //Load data
    UInt dataLines = dataFile.vector_variable_size ( "data" ) / variablesNumber;

    M_data.clear();
    for ( UInt j ( 0 ); j < variablesNumber; ++j )
    {
        M_data[M_variables[j]].reserve ( dataLines );
    }

    for ( UInt i ( 0 ); i < dataLines; ++i )
        for ( UInt j ( 0 ); j < variablesNumber; ++j )
        {
            M_data[M_variables[j]].push_back ( scale[j] * dataFile ( "data", 0.0, i * variablesNumber + j ) );
        }

#ifdef HAVE_LIFEV_DEBUG
    output.str ("");
    output << "                                                 loop: " << M_loop << "\n";
    output << "                                                 data:";
    for ( UInt i (0); i < dataLines; ++i )
    {
        if (i > 0)
        {
            output << "                                                                 ";
        }

        for ( UInt j (0); j < variablesNumber; ++j )
        {
            output << " " << M_data[ M_variables[j] ][i];
        }
        output << "\n";
    }
    debugStream ( 5022 ) << output.str();
#endif

    //Initialize iterator
    M_dataIterator = M_data[M_variables[0]].begin();

    //Update the data container (IT IS A COPY!) with the correct base string for the BCInterfaceFunctionParser
    if ( stringsVector.size() < 2 )
    {
        dataCopy->setBaseString ( dataFile ( "function", "Undefined" ) );
    }
    else
    {
        boost::replace_all ( stringsVector[1], "]", "" );
        dataCopy->setBaseString ( dataFile ( ( "function" + stringsVector[1] ).c_str(), "Undefined" ) );
    }

    // Now data contains the real base string
    functionParser_Type::setData ( dataCopy );

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5022 ) << "                                             function: " << dataCopy->baseString() << "\n";
#endif

}


// ===================================================
// Private Methods
// ===================================================
template< typename BcHandlerType, typename PhysicalSolverType >
inline void
BCInterfaceFunctionParserFile< BcHandlerType, PhysicalSolverType >::dataInterpolation()
{
    //Get variable
    Real X = functionParser_Type::M_parser->variable ( M_variables[0] );

    //If it is a loop scale the variable: X = X - (ceil( X / Xmax ) -1) * Xmax
    if ( M_loop )
    {
        X -= ( std::ceil ( X / M_data[M_variables[0]].back() ) - 1 ) * M_data[M_variables[0]].back();
    }

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5022 ) << "                                                    variable: " << X << "\n";
#endif

    //Move Iterator
    for ( ;; )
    {

#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 5022 ) << "                                       iterator  position   : " << static_cast<Real> ( M_dataIterator - M_data[ M_variables[0] ].begin() ) << "\n";
        debugStream ( 5022 ) << "                                       variable (position)  : " << *M_dataIterator << "\n";
        debugStream ( 5022 ) << "                                       variable (position+1): " << * (M_dataIterator + 1) << "\n";
#endif

        if ( X >= *M_dataIterator && X <= * ( M_dataIterator + 1 ) )
        {
            break;
        }

        if ( X > *M_dataIterator )
        {
            if ( M_dataIterator + 1 == M_data[M_variables[0]].end() )
            {
                break;
            }
            else
            {
                ++M_dataIterator;
            }
        }
        else
        {
            if ( M_dataIterator == M_data[M_variables[0]].begin() )
            {
                break;
            }
            else
            {
                --M_dataIterator;
            }
        }
    }

    //Linear interpolation (extrapolation if X > xB)
    Real xA, xB, A, B;
    for ( UInt j ( 1 ), position = static_cast< UInt > ( M_dataIterator - M_data[M_variables[0]].begin() ) ;
            j < static_cast< UInt > ( M_variables.size() ); ++j )
    {
        xA = M_data[M_variables[0]][position];
        xB = M_data[M_variables[0]][position + 1];
        A  = M_data[M_variables[j]][position];
        B  = M_data[M_variables[j]][position + 1];

        functionParser_Type::M_parser->setVariable ( M_variables[j], A + ( B - A ) / ( xB - xA ) * ( X - xA ) );

#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 5022 ) << "                                                          " << M_variables[j] << " = " << A + (B - A) / (xB - xA) * (X - xA) << "\n";
#endif
    }
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionParserFile_H */
