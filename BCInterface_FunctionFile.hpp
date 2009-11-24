//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
//@HEADER

/*!
 *  @file
 *  @brief BCInterface_Function
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 09-07-2009
 */
#ifndef BCInterface_FunctionFile_H
#define BCInterface_FunctionFile_H 1

#include <lifemc/lifefem/BCInterface_Function.hpp>

namespace LifeV {

//! BCInterface_FunctionFile - LifeV bcFunction wrapper for BCInterface
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface and SpiritParser. It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
 *
 *
 *
 *  <b>DETAILS:</b>
 *
 *  The constructor of the class takes a string contains the GetPot file name.
 *  The GetPot file has the following structure:
 *
 *  function:  contains the expression of the function to use (as described in the BCInterfaceFunction class).
 *
 *  variables: contains the list of variables and coefficients present in the function.
 *             The first one is the variable and should be sorted in a growing order,
 *             while all the others are coefficients.
 *
 *  data:      contains the discrete list of variable and related coefficients. It must
 *             have n columns, where n is the number of variables (and coefficients).
 *
 *	scale:     Contains n optional coefficients (Default = 1), which multiply each column
 *	           in the data table.
 *
 *  loop:      Useful for time dependent value: at the end of the list, restart using the first value
 *             instead of extrapolating.
 *
 *	NOTE:
 *	During the execution, if the value of the variable (usually the time) is not present in the 'data' table,
 *	the class linearly interpolates the value between the two closest one. Moreover, if the value of the variable is higher
 *	than anyone present in the 'data' table, the class linearly extrapolates the value using the last two values in the table.
 *
 *   <b>EXAMPLE OF DATAFILE</b>
 *
 *  function	= '(0,0,q)'
 *  variables	=   't				q'
 *  scale		= 	'1				1'
 *  data		=  '0.000000000		1.00
 *  				0.333333333		2.00
 *  				0.666666666		3.00
 *  				1.000000000		4.00'
 */
template< typename Operator >
class BCInterface_FunctionFile: public virtual BCInterface_Function< Operator >
{
public:

    typedef BCInterface_Function< Operator > super;

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    BCInterface_FunctionFile();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_FunctionFile( const BCInterface_Data< Operator >& data );

    //! Copy constructor
    /*!
     * @param function BCInterface_FunctionFile
     */
    BCInterface_FunctionFile( const BCInterface_FunctionFile& function );

    //! Destructor
    virtual ~BCInterface_FunctionFile() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param function BCInterface_FunctionFile
     * @return reference to a copy of the class
     */
    virtual BCInterface_FunctionFile& operator=( const BCInterface_FunctionFile& function );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void SetData( const BCInterface_Data< Operator >& data );

    //! Compare function
    /*!
     * @param data BC data loaded from GetPot file
     * @return true if the functions are equal, false if they aren't
     */
    virtual bool Compare( const BCInterface_Data< Operator >& data );

    //@}

protected:

    std::string M_fileName;

private:

    std::vector< std::string >                   M_variables;
    std::vector< Real >                          M_scale;
    bool                                         M_loop;
    std::map< std::string, std::vector< Real > > M_data;
    std::vector< Real >::iterator                M_dataIterator;

    //! @name Private functions
    //@{

    //! loadData
    inline void LoadData( BCInterface_Data< Operator > data );

    //! Linear interpolation (extrapolation) between two values of the data.
    inline void DataInterpolation();

    //@}
};

//! Factory create function
template< typename Operator >
inline BCInterface_Function< Operator >* createFunctionFile()
{
    return new BCInterface_FunctionFile< Operator > ();
}

// ===================================================
// Constructors
// ===================================================
template< typename Operator >
BCInterface_FunctionFile< Operator >::BCInterface_FunctionFile() :
    BCInterface_Function< Operator > (),
    M_fileName                       (),
    M_variables                      (),
    M_scale                          (),
    M_loop                           (),
    M_data                           (),
    M_dataIterator                   ()
{

#ifdef DEBUG
    Debug( 5022 ) << "BCInterface_FunctionFile::BCInterface_FunctionFile()" << "\n";
#endif

}

template< typename Operator >
BCInterface_FunctionFile< Operator >::BCInterface_FunctionFile( const BCInterface_Data< Operator >& data ) :
    BCInterface_Function< Operator > (),
    M_fileName                       (),
    M_variables                      (),
    M_scale                          (),
    M_loop                           (),
    M_data                           (),
    M_dataIterator                   ()
{

#ifdef DEBUG
    Debug( 5022 ) << "BCInterface_FunctionFile::BCInterface_FunctionFile( data )" << "\n";
#endif

    this->SetData( data );
}

template< typename Operator >
BCInterface_FunctionFile< Operator >::BCInterface_FunctionFile( const BCInterface_FunctionFile& function ) :
    BCInterface_Function< Operator > ( function ),
    M_fileName                       ( function.M_fileName ),
    M_variables                      ( function.M_variables ),
    M_scale                          ( function.M_scale ),
    M_loop                           ( function.M_loop ),
    M_data                           ( function.M_data ),
    M_dataIterator                   ( function.M_dataIterator )
{
}

// ===================================================
// Methods
// ===================================================
template< typename Operator >
BCInterface_FunctionFile< Operator >&
BCInterface_FunctionFile< Operator >::operator=( const BCInterface_FunctionFile& function )
{
    if ( this != &function )
    {
        super::operator=( function );
        M_fileName     = function.M_fileName;
        M_variables    = function.M_variables;
        M_scale        = function.M_scale;
        M_loop         = function.M_loop;
        M_data         = function.M_data;
        M_dataIterator = function.M_dataIterator;
    }

    return *this;
}

template< typename Operator >
void
BCInterface_FunctionFile< Operator >::SetData( const BCInterface_Data< Operator >& data )
{
    M_fileName = data.GetBaseString(); //The base string contains the file name

#ifdef DEBUG
    Debug( 5022 ) << "BCInterface_FunctionFile::setData             fileName: " << M_fileName << "\n";
#endif

    LoadData( data );
}

template< typename Operator >
bool
BCInterface_FunctionFile< Operator >::Compare( const BCInterface_Data< Operator >& data )
{
    return M_fileName.compare( data.GetBaseString() ) == 0 && super::M_comV == data.GetComV();
}

// ===================================================
// Private functions
// ===================================================
template< typename Operator >
inline void
BCInterface_FunctionFile< Operator >::LoadData( BCInterface_Data< Operator > data )
{
    std::vector< std::string > stringsVector;
    boost::split( stringsVector, M_fileName, boost::is_any_of( "[" ) );

    //Load data from file
    GetPot dataFile( stringsVector[0] );

    //Set variables
    UInt variablesNumber = dataFile.vector_variable_size( "variables" );

    M_variables.clear();
    M_variables.reserve( variablesNumber );

    M_scale.clear();
    M_scale.reserve( variablesNumber );

    for ( UInt j( 0 ); j < variablesNumber; ++j )
    {
        M_variables.push_back( dataFile( "variables", "unknown", j ) );
        M_scale.push_back( dataFile( "scale", 1.0, j ) );
    }

#ifdef DEBUG
    std::stringstream output;
    output << "BCInterface_FunctionFile::loadData           variables: ";
    for ( UInt j(0); j < variablesNumber; ++j )
        output << M_variables[j] << "  ";

    output << "\n                                                           scale: ";
    for ( UInt j(0); j < variablesNumber; ++j )
        output << M_scale[j] << "  ";

    Debug( 5022 ) << output.str() << "\n";
#endif

    //Load loop flag
    M_loop = dataFile( "loop", false );

    //Load data
    UInt dataLines = dataFile.vector_variable_size( "data" ) / variablesNumber;

    M_data.clear();
    for ( UInt j( 0 ); j < variablesNumber; ++j )
        M_data[M_variables[j]].reserve( dataLines );

    for ( UInt i( 0 ); i < dataLines; ++i )
        for ( UInt j( 0 ); j < variablesNumber; ++j )
            M_data[M_variables[j]].push_back( M_scale[j] * dataFile( "data", 0.0, i
                    * variablesNumber + j ) );

#ifdef DEBUG
    output.str("");
    output << "                                                 loop: " << M_loop << "\n";
    output << "                                                 data:";
    for ( UInt i(0); i < dataLines; ++i )
    {
        if (i > 0)
            output << "                                                                 ";

        for ( UInt j(0); j < variablesNumber; ++j )
            output << " " << M_data[ M_variables[j] ][i];
        output << "\n";
    }
    Debug( 5022 ) << output.str();
#endif

    //Initialize iterator
    M_dataIterator = M_data[M_variables[0]].begin();

    //Update the data container (IT IS A COPY!) with the correct base string for the BCInterface_Function
    if ( stringsVector.size() < 2 )
        data.SetBaseString( dataFile( "function", "Undefined" ) );
    else
    {
        boost::replace_all( stringsVector[1], "]", "" );
        data.SetBaseString( dataFile( ( "function" + stringsVector[1] ).c_str(), "Undefined" ) );
    }

    // Now data contains the real base string
    super::SetData( data );

#ifdef DEBUG
    Debug( 5022 ) << "                                             function: " << data.GetBaseString() << "\n";
#endif
}

template< typename Operator >
inline void
BCInterface_FunctionFile< Operator >::DataInterpolation()
{
    //Get variable
    Real X = super::M_parser->getVariable( M_variables[0] );

    //If it is a loop scale the variable: X = X - (ceil( X / Xmax ) -1) * Xmax
    if ( M_loop )
        X -= ( std::ceil( X / M_data[M_variables[0]].back() ) - 1 ) * M_data[M_variables[0]].back();

#ifdef DEBUG
    Debug( 5022 ) << "                                                    variable: " << X << "\n";
#endif

    //Move Iterator
    for ( ;; )
    {

#ifdef DEBUG
        Debug( 5022 ) << "                                       iterator  position   : " << static_cast<Real> ( M_dataIterator - M_data[ M_variables[0] ].begin() ) << "\n";
        Debug( 5022 ) << "                                       variable (position)  : " << *M_dataIterator << "\n";
        Debug( 5022 ) << "                                       variable (position+1): " << *(M_dataIterator+1) << "\n";
#endif

        if ( X >= *M_dataIterator && X <= *( M_dataIterator + 1 ) )
            break;

        if ( X > *M_dataIterator )
        {
            if ( M_dataIterator + 1 == M_data[M_variables[0]].end() )
                break;
            else
                ++M_dataIterator;
        }
        else
        {
            if ( M_dataIterator == M_data[M_variables[0]].begin() )
                break;
            else
                --M_dataIterator;
        }
    }

    //Linear interpolation (extrapolation if X > xB)
    Real xA, xB, A, B;
    for ( UInt j( 1 ), position = static_cast< UInt > ( M_dataIterator - M_data[M_variables[0]].begin() ) ;
          j < static_cast< UInt > ( M_variables.size() ); ++j )
    {
        xA = M_data[M_variables[0]][position];
        xB = M_data[M_variables[0]][position + 1];
        A  = M_data[M_variables[j]][position];
        B  = M_data[M_variables[j]][position + 1];

        super::M_parser->setVariable( M_variables[j], A + ( B - A ) / ( xB - xA ) * ( X - xA ) );

#ifdef DEBUG
        Debug( 5022 ) << "                                                          " << M_variables[j] << " = " << A+(B-A)/(xB-xA)*(X-xA) << "\n";
#endif
    }
}

} // Namespace LifeV

#endif /* BCInterface_FunctionFile_H */
