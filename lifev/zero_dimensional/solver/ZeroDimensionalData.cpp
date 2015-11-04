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
 *  @date 16-11-2011
 *  @author Mahmoud Jafargholi
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/zero_dimensional/solver/ZeroDimensionalData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
ZeroDimensionalData::ZeroDimensionalData() :
    M_time(), M_outPutFormat ( OutPutFormat ( "20",
                                              "12",
                                              "   ",
                                              5000 ) ), M_circuitData ( new ZeroDimensionalCircuitData ), M_venousPressure ( 0. )
{
}

ZeroDimensionalData::~ZeroDimensionalData()
{
    M_voltageFileStream.close();
    M_currentFileStream.close();
    M_balanceFileStream.close();
}

// ===================================================
// Methods
// ===================================================
void ZeroDimensionalData::setup ( const GetPot& dataFile,
                                  bcPtr_Type bc,
                                  const std::string& section )
{
    if ( !M_time.get() )
        M_time.reset ( new time_Type ( dataFile,
                                       section + "/time_discretization" ) );

    std::string circuirDataFile = dataFile ( ( section + "/CircuitDataFile" ).data(),
                                             "./inputFilexx.dat" );
    GetPot dataFileCircuit ( circuirDataFile );
    std::string circuitFile = dataFileCircuit ( "Files/InputFiles/CircuitFile",
                                                "./circuitFilexx.dat" );
    std::string folder = dataFileCircuit ( "Files/OutputFiles/Folder",
                                           "outputxx" );//TODO create folder output

    std::string voltageFile = folder + "/" + dataFileCircuit ( "Files/OutputFiles/VoltageFile",
                                                               "voltagexx.txt" );
    std::string currentFile = folder + "/" + dataFileCircuit ( "Files/OutputFiles/CurrentFile",
                                                               "currentxx.txt" );
    std::string balanceFile = folder + "/" + dataFileCircuit ( "Files/OutputFiles/BalanceCurrentinNode",
                                                               "balancexx.txt" );
    M_voltageFileStream.open ( voltageFile.c_str() );
    M_currentFileStream.open ( currentFile.c_str() );
    M_balanceFileStream.open ( balanceFile.c_str() );

    //Build the Circuit
    M_circuitData->buildCircuit ( circuitFile.c_str(),
                                  bc );

    //assign varible index
    assignVaribleIndex();

    M_solverData.method = dataFile ( ( section + "/Solver/method" ).data(), "IRK" );
    M_solverData.numberTimeStep = dataFile ( ( section + "/Solver/numberTimeStep" ).data(), 1 );
    M_solverData.maxError = dataFile ( ( section + "/Solver/maxError" ).data(), 1.0 );
    M_solverData.reltol = dataFile ( ( section + "/Solver/reltol" ).data(), 0.1 );
    M_solverData.abstol = dataFile ( ( section + "/Solver/abstol" ).data(), 1.0 );
    M_solverData.maxOrder = dataFile ( ( section + "/Solver/maxOrder" ).data(), 1 );
    M_solverData.verbose = dataFile ( ( section + "/Solver/verbose" ).data(), false );
    M_solverData.verboseLevel = dataFile ( ( section + "/Solver/verboseLevel" ).data(), 0 );
    M_solverData.useNOX = dataFile ( ( section + "/Solver/useNOX" ).data(), false );
    M_solverData.fixTimeStep = dataFile ( ( section + "/Solver/fixTimeStep" ).data(), true );
    M_solverData.extraLSParamsFile = dataFile ( ( section + "/Solver/extraLinearSolverParamsFile" ).data(), "./Extra_AztecOO_Params.xml" );
    M_solverData.linearSolverParamsFile = dataFile ( ( section + "/Solver/linearSolverParamsUsedFile" ).data(), "./lowsf.aztecoo.used.xml" );

    //Set zero initial condition
    // TODO: change to general initial condition
    initializeSolution();

    //Wrire Header OutputFiles
    writeHeaders();
}

void ZeroDimensionalData::initializeSolution()
{
    ptrVecZeroDimensionalElementPtr_Type elementList = M_circuitData->Elements() ->elementList();
    for ( iterZeroDimensionalElement_Type theElement = elementList->begin(); theElement != elementList->end(); theElement++ )
    {
        ( *theElement )->setCurrent ( 0.0 );
        ( *theElement )->setDeltaCurrent ( 0.0 );
    }

    ptrVecZeroDimensionalNodePtr_Type nodeList = M_circuitData->Nodes() ->nodeList();
    for ( iterZeroDimensionalNode_Type theNode = nodeList->begin(); theNode != nodeList->end(); theNode++ )
    {
        ( *theNode )->setVoltage ( 0.0 );
        ( *theNode )->setDeltaVoltage ( 0.0 );
    }
}

//! update source elements
void ZeroDimensionalData::updateBC()
{
    Real time = M_time->time();
    ptrVecZeroDimensionalElementCurrentSourcePtr_Type currentElementList = M_circuitData->Elements()-> currentSourceList();

    for ( iterZeroDimensionalElementCurrentSource_Type theElement = currentElementList->begin(); theElement != currentElementList->end(); theElement++ )
    {
        ( *theElement )->setCurrentByTime ( time );
        std::cout << ( *theElement )->current() << std::endl;
    }

    ptrVecZeroDimensionalElementVoltageSourcePtr_Type voltageElementList = M_circuitData->Elements()-> voltageSourceList();
    for ( iterZeroDimensionalElementVoltageSourcePtr_Type theElement = voltageElementList->begin(); theElement != voltageElementList->end(); theElement++ )
    {
        ( *theElement )->setVoltageByTime ( time );
        std::cout << ( *theElement )->voltage() << std::endl;
    }
}

void ZeroDimensionalData::saveSolution()
{
    ptrVecZeroDimensionalElementPtr_Type elementList = M_circuitData->Elements() ->elementList();
    M_outPutFormat.writeDataFormat ( M_time->time(),
                                     M_currentFileStream,
                                     M_outPutFormat.space );

    //write current
    for ( iterZeroDimensionalElement_Type theElement = elementList->begin(); theElement != elementList->end(); theElement++ )
    {
        M_outPutFormat.writeDataFormat ( ( *theElement )->current(),
                                         M_currentFileStream,
                                         M_outPutFormat.space );
    }

    //write voltage and current balance at each node
    M_outPutFormat.writeNewLine ( M_currentFileStream );
    ptrVecZeroDimensionalNodePtr_Type nodeList = M_circuitData->Nodes() ->nodeList();
    M_outPutFormat.writeDataFormat ( M_time->time(),
                                     M_voltageFileStream,
                                     M_outPutFormat.space );
    M_outPutFormat.writeDataFormat ( M_time->time(),
                                     M_balanceFileStream,
                                     M_outPutFormat.space );

    for ( iterZeroDimensionalNode_Type theNode = nodeList->begin(); theNode != nodeList->end(); theNode++ )
    {
        M_outPutFormat.writeDataFormat ( ( *theNode )->voltage(),
                                         M_voltageFileStream,
                                         M_outPutFormat.space );
        M_outPutFormat.writeDataFormat ( ( *theNode )->currentBalance(),
                                         M_balanceFileStream,
                                         M_outPutFormat.space );
    }

    M_outPutFormat.writeNewLine ( M_voltageFileStream );
    M_outPutFormat.writeNewLine ( M_balanceFileStream );
}


void ZeroDimensionalData::showMeVariables()
{
    std::cout << "This MultiscaleModel0D::ShowVariables----------------------" << std::endl;
    zeroDimensionalNodeSPtr_Type Nodes = M_circuitData->Nodes();
    zeroDimensionalElementSPtr_Type Elements = M_circuitData->Elements();

    for ( iterZeroDimensionalNodeUnknown_Type theUnknownNode = Nodes->unknownNodeList()->begin(); theUnknownNode != Nodes->unknownNodeList()->end(); theUnknownNode++ )
    {
        ( *theUnknownNode )->showMe ( 1 );
    }

    for ( iterZeroDimensionalElementPassiveInductor_Type theInductor = Elements->inductorList() ->begin(); theInductor
            != Elements->inductorList() ->end(); theInductor++ )
    {
        ( *theInductor )->showMe ( 1 );
    }
}

void ZeroDimensionalData::assignVaribleIndex()
{
    zeroDimensionalNodeSPtr_Type Nodes = M_circuitData->Nodes();
    ptrVecZeroDimensionalNodeUnknownPtr_Type unKnownList = Nodes->unknownNodeList();
    zeroDimensionalElementSPtr_Type Elements = M_circuitData->Elements();
    ptrVecZeroDimensionalElementPassiveInductorPtr_Type inductorList = Elements->inductorList();
    M_unknownCounter = 0;

    //set variable index for unknown nodes
    for ( iterZeroDimensionalNodeUnknown_Type theUnknownNode = unKnownList->begin(); theUnknownNode != unKnownList->end(); theUnknownNode++ )
    {
        ( *theUnknownNode )->assignVariableIndex ( M_unknownCounter );
        M_unknownCounter++;
    }

    //set variable index for inductor current
    for ( iterZeroDimensionalElementPassiveInductor_Type theInductor = inductorList->begin(); theInductor != inductorList->end(); theInductor++ )
    {
        ( *theInductor )->assignVariableIndex ( M_unknownCounter );
        M_unknownCounter++;
    }
}

void ZeroDimensionalData::writeHeaders()
{
    //write header for current file
    ptrVecZeroDimensionalElementPtr_Type elementList = M_circuitData->Elements() ->elementList();
    M_outPutFormat.writeDataFormat ( "% time", M_currentFileStream, M_outPutFormat.space );
    for ( iterZeroDimensionalElement_Type theElement = elementList->begin(); theElement != elementList->end(); theElement++ )
    {
        M_outPutFormat.writeDataFormat ( ( *theElement )->id(), M_currentFileStream, M_outPutFormat.space );
    }
    M_outPutFormat.writeNewLine ( M_currentFileStream );

    //write header for voltage and current balance at each node (voltage file and balance file)
    ptrVecZeroDimensionalNodePtr_Type nodeList = M_circuitData->Nodes() ->nodeList();
    M_outPutFormat.writeDataFormat ( "% time", M_balanceFileStream, M_outPutFormat.space );
    M_outPutFormat.writeDataFormat ( "% time", M_voltageFileStream, M_outPutFormat.space );
    for ( iterZeroDimensionalNode_Type theNode = nodeList->begin(); theNode != nodeList->end(); theNode++ )
    {
        M_outPutFormat.writeDataFormat ( ( *theNode )->id(), M_voltageFileStream, M_outPutFormat.space );
        M_outPutFormat.writeDataFormat ( ( *theNode )->id(), M_balanceFileStream, M_outPutFormat.space );
    }
    M_outPutFormat.writeNewLine ( M_voltageFileStream );
    M_outPutFormat.writeNewLine ( M_balanceFileStream );
}

} // LifeV namespace
