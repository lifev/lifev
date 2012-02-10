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
 *  @brief File containing a class for 0D model data handling.
 *
 *  @version 1.0
 *  @date 16-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalData_H
#define ZeroDimensionalData_H

// LIFEV
#include <lifemc/lifesolver/ZeroDimensionalCircuitData.hpp>

namespace LifeV
{

//! Data container for 0D model
class ZeroDimensionalData
{
public:

    typedef TimeData                                                  time_Type;

    typedef boost::shared_ptr < time_Type >                           timePtr_Type;

    //! Constructor
    explicit ZeroDimensionalData();

    //! Destructor
    virtual ~ZeroDimensionalData();

    //! setup model
    void setup( const GetPot& dataFile, bcInterfacePtr_Type bc, const std::string& section = "0D_Model" );

    //! show some information
    void showMe( ) const;

    //! show variables
    void showMeVariables() ;

    //! set time
    void setTimeData( const timePtr_Type TimeData ) { M_time = TimeData; }

    //! initialize Solution
    void initializeSolution() ;

    //! save solution 
    void saveSolution() ;

    //! update source elements
    void updateBC();

    timePtr_Type dataTime() const { return M_time; }
    
    //! get circuit data container
    zeroDimensionalCircuitDataPtr_Type circuitData() const {return M_circuitData;}

    //!total number of unknowns
    Int unknownCounter() const {return  M_unknownCounter;}

    //! Rhytmos solver data container
    struct SolverData {
        std::string        method;
        int                numberTimeStep;
        double             maxError;
        double             reltol;
        double             abstol;
        int                maxOrder;
        bool               verbose;
        int                verbosLevel;
        bool               useNOX;
        bool               fixTimeStep;
        std::string        extraLSParamsFile;
        std::string        linearSolverParamsFile;
    };
    
    typedef struct SolverData       solverData_Type;

    solverData_Type solverData() const {return M_solverData;}

private:

    void writeHeaders();

    void assignVaribleIndex();

    timePtr_Type                         M_time;
    zeroDimensionalCircuitDataPtr_Type   M_circuitData      ;
    std::ofstream                        M_voltageFileStream;
    std::ofstream                        M_currentFileStream;
    std::ofstream                        M_balanceFileStream;
    Int                                  M_unknownCounter   ;
    OutPutFormat                         M_outPutFormat;
    solverData_Type                      M_solverData;
};

} // LifeV namespace

#endif //ZeroDimensionalData_H
