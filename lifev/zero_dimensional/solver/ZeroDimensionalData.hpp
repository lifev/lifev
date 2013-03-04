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
 *  @version alpha (experimental)
 *
 *  @date 16-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalData_H
#define ZeroDimensionalData_H

// LIFEV
#include <lifev/zero_dimensional/solver/ZeroDimensionalCircuitData.hpp>

namespace LifeV
{

//! Data container for 0D model
class ZeroDimensionalData
{
public:

    // TODO This should be a separate class and not a struct here
    //! Rhytmos solver data container
    struct SolverData
    {
        std::string        method;
        Int                numberTimeStep;
        Real               maxError;
        Real               reltol;
        Real               abstol;
        Int                maxOrder;
        bool               verbose;
        Int                verboseLevel;
        bool               useNOX;
        bool               fixTimeStep;
        std::string        extraLSParamsFile;
        std::string        linearSolverParamsFile;
    };

    typedef struct SolverData                                         solverData_Type;
    typedef TimeData                                                  time_Type;
    typedef boost::shared_ptr < time_Type >                           timePtr_Type;

    //! Constructor
    explicit ZeroDimensionalData();

    //! Destructor
    virtual ~ZeroDimensionalData();

    //! setup model
    void setup ( const GetPot& dataFile, bcPtr_Type bc, const std::string& section = "0D_Model" );

    //! initialize Solution
    void initializeSolution() ;

    //! update source elements
    void updateBC();

    //! save solution
    void saveSolution() ;

    //! show some information
    void showMe() const
    {
        M_circuitData->showMe();
    }

    //! show variables
    void showMeVariables() ;

    //! set time
    void setTimeData ( const timePtr_Type timeData )
    {
        M_time = timeData;
    }

    const timePtr_Type& dataTime() const
    {
        return M_time;
    }

    //! get circuit data container
    zeroDimensionalCircuitDataPtr_Type circuitData() const
    {
        return M_circuitData;
    }

    //!total number of unknowns
    const Int& unknownCounter() const
    {
        return M_unknownCounter;
    }

    const solverData_Type& solverData() const
    {
        return M_solverData;
    }

private:

    // TODO: The output part should be rewritten following the example in the OneDFSI solver
    void writeHeaders();

    void assignVaribleIndex();

    timePtr_Type                         M_time;
    OutPutFormat                         M_outPutFormat;
    zeroDimensionalCircuitDataPtr_Type   M_circuitData;
    std::ofstream                        M_voltageFileStream;
    std::ofstream                        M_currentFileStream;
    std::ofstream                        M_balanceFileStream;
    Int                                  M_unknownCounter;
    solverData_Type                      M_solverData;
};

} // LifeV namespace

#endif //ZeroDimensionalData_H
