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
 *  @brief File containing a class for 1D model data handling.
 *
 *  @version 1.0
 *  @date 01-07-2004
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 12-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <life/lifesolver/ZeroDimensionalCircuitData.hpp>

namespace LifeV
{


// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElement::ZeroDimensionalElement()
{
}
void
ZeroDimensionalElement::showMe(__attribute__((unused)) Int flag)
{
    cout<<"Id = "<<id();
    cout<<"\t type = ";
    cout<<enum2string(type());
}

const std::string
ZeroDimensionalElement::enum2string(const ZeroDimentionalElementType & type)
{
    std::string output;
    switch (type)
        {
            case resistor:
            {
                output = "resistor";
                break;
            }
            case capacitor:
            {
                output = "capacitor";
                break;
            }
            case inductor:
            {
                output = "inductor";
                break;
            }
            case diode:
            {
                output = "diode";
                break;
            }
            case voltageSource:
            {
                output =  "voltageSource";
                break;
            }
            case currentSource:
            {
                output = "currentSource";
                break;
            }
        }
    return (output);
}


// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementPassive::ZeroDimensionalElementPassive():
    M_parameter                 (),
    M_nodeIndex                 ()
{
    M_nodeIndex.reserve(2);
}

void
ZeroDimensionalElementPassive::showMe(Int flag)
{
    ZeroDimensionalElement::showMe(flag);
    cout<<"\t Node1= "<<nodeIndex(0)<<"\t Node2= "<<nodeIndex(1);
}

void
ZeroDimensionalElementPassive::connectElement(zeroDimensionalNodeSPtr_Type & Nodes)
{
    Int Index=M_nodeIndex.at(0);
    zeroDimensionalNodePtr_Type    theNode=Nodes->nodeListAt(Index);

    theNode->setelementListIndex(M_id);
    Int otherEndIndex=M_nodeIndex.at(1);
    theNode->setnodeListIndex(otherEndIndex);

    Index=M_nodeIndex.at(1);
    theNode=Nodes->nodeListAt(Index);

    theNode->setelementListIndex(M_id);
    otherEndIndex=M_nodeIndex.at(0);
    theNode->setnodeListIndex(otherEndIndex);
}

Real
ZeroDimensionalElementPassive::direction(const Int & nodeId)const{
    if (M_nodeIndex.at(0) == nodeId){
        return -1.0;
    } else {
        return  1.0;
    }
}
// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementPassiveResistor::ZeroDimensionalElementPassiveResistor()
{
    M_type      =  resistor                           ;
}

void
ZeroDimensionalElementPassiveResistor::showMe(Int flag)
{
    ZeroDimensionalElementPassive::showMe(flag);
    cout<<"\t parameter = "<<parameter()<<std::endl;
}

void
ZeroDimensionalElementPassiveResistor::buildABC(__attribute__((unused)) matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes)
{
    for (int i = 0; i < 2; i++) {
        const int& theNodeCounter =  (i)%2;
        const int& theOtherNodeCounter =  (i+1)%2;
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< theNodeCounter << "\t"<< theOtherNodeCounter << "\n";
#endif
        const ZeroDimensionalNode& theNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theNodeCounter]));
        const ZeroDimensionalNode& theOtherNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theOtherNodeCounter]));
        if (theNodeTest.type() == unknownNode ){
            const ZeroDimensionalNodeUnknown& theNode = *(Nodes->unknownNodeMapAt(M_nodeIndex[theNodeCounter]));
            const int& equationRow = theNode.equationRow();
            B.addToCoefficient(equationRow,theNode.variableIndex(), -M_parameter);
#ifdef HAVE_LIFEV_DEBUG
          Debug( 8151 )<< equationRow << "\t"<< theNode.variableIndex() << "\t"<< M_parameter << "\t""\n";
#endif
            if (theOtherNodeTest.type() == unknownNode){
                const ZeroDimensionalNodeUnknown& theOtherNode = *(Nodes->unknownNodeMapAt(M_nodeIndex[theOtherNodeCounter]));
                B.addToCoefficient(equationRow, theOtherNode.variableIndex(), M_parameter);
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< equationRow << "\t"<< theOtherNode.variableIndex() << "\n";
#endif
            } else {
                const ZeroDimensionalNodeKnown & theOtherNode = *(Nodes->knownNodeMapAt(M_nodeIndex[theOtherNodeCounter]));
                C[equationRow] += (M_parameter * theOtherNode.voltage());

            }
        }
#ifdef HAVE_LIFEV_DEBUG
        A.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        B.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        C.epetraVector().Print(cout);
        cout<< "------------------\n";
#endif
    }
}

void
ZeroDimensionalElementPassiveResistor::deepUpdate(const ZeroDimensionalNodeS& Nodes){
    M_current       = ( Nodes.nodeListAt(M_nodeIndex.at(0))->voltage()      - Nodes.nodeListAt(M_nodeIndex.at(1))->voltage()       )     * M_parameter;
    M_deltaCurrent  = ( Nodes.nodeListAt(M_nodeIndex.at(0))->deltaVoltage() - Nodes.nodeListAt(M_nodeIndex.at(1))->deltaVoltage()  )     * M_parameter;
}
// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementPassiveDiode::ZeroDimensionalElementPassiveDiode():
    M_alpha                       (),
    M_beta                        (),
    M_forwardBias                 ()
{
    M_type      =  diode;
}

void
ZeroDimensionalElementPassiveDiode::showMe(Int flag)
{
    ZeroDimensionalElementPassiveResistor::showMe(flag);
    cout<<"\t FB = "<<forwardBias()<<"\t alpha = "<<alpha()<<"\t beta = "<<beta()<<std::endl;
}
void
ZeroDimensionalElementPassiveDiode::calculateEffectiveResistance(const double& voltage){

    double e0     =(M_beta*exp(M_alpha*(0-M_forwardBias)));

    if (voltage == 0){
        M_parameter=e0;
    } else {
        double current    = (M_beta*exp(M_alpha*(voltage-M_forwardBias)))-e0;
        M_parameter = abs(current / voltage) ;
    }
}
void
ZeroDimensionalElementPassiveDiode::deepUpdate(const ZeroDimensionalNodeS& Nodes){
    double deltaVoltage = ( Nodes.nodeListAt(M_nodeIndex.at(0))->voltage()      - Nodes.nodeListAt(M_nodeIndex.at(1))->voltage()       );
    calculateEffectiveResistance(deltaVoltage);
    M_current           = deltaVoltage * M_parameter;
    M_deltaCurrent      = ( Nodes.nodeListAt(M_nodeIndex.at(0))->deltaVoltage() - Nodes.nodeListAt(M_nodeIndex.at(1))->deltaVoltage() ) * M_parameter;
}

void
ZeroDimensionalElementPassiveDiode::buildABC(matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes)
{

    const int& theNodeCounter =  0;
    const int& theOtherNodeCounter =  1;
    const ZeroDimensionalNode& theNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theNodeCounter]));
    const ZeroDimensionalNode& theOtherNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theOtherNodeCounter]));
    Real deltaVoltage = theNodeTest.voltage() - theOtherNodeTest.voltage();
    calculateEffectiveResistance(deltaVoltage);
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8151 )<< "------------------\n"<<deltaVoltage<<"------------------\n"<<M_parameter;
#endif
    ZeroDimensionalElementPassiveResistor::buildABC(A,B,C,Nodes);
}

// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementPassiveCapacitor::ZeroDimensionalElementPassiveCapacitor()
{
    M_type      =  capacitor;
}

void
ZeroDimensionalElementPassiveCapacitor::showMe(Int flag)
{
    ZeroDimensionalElementPassive::showMe(flag);
    cout<<"\t parameter = "<<parameter()<<std::endl;
}

void
ZeroDimensionalElementPassiveCapacitor::deepUpdate(const ZeroDimensionalNodeS& Nodes){

    M_current       = ( Nodes.nodeListAt(M_nodeIndex.at(0))->deltaVoltage()      - Nodes.nodeListAt(M_nodeIndex.at(1))->deltaVoltage()       )     * M_parameter;
}

void
ZeroDimensionalElementPassiveCapacitor::buildABC(matrix_Type& A,__attribute__((unused)) matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes)
{
    for (int i = 0; i < 2; i++) {
        const int& theNodeCounter =  (i)%2;
        const int& theOtherNodeCounter =  (i+1)%2;
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< theNodeCounter << "\t"<< theOtherNodeCounter << "\n";
#endif
        const ZeroDimensionalNode& theNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theNodeCounter]));
        const ZeroDimensionalNode& theOtherNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theOtherNodeCounter]));
        if (theNodeTest.type() == unknownNode ){
            const ZeroDimensionalNodeUnknown& theNode = *(Nodes->unknownNodeMapAt(M_nodeIndex[theNodeCounter]));
            const int& equationRow = theNode.equationRow();
            A.addToCoefficient(equationRow,theNode.variableIndex(), -M_parameter);
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< equationRow << "\t"<< theNode.variableIndex() << "\t"<< M_parameter << "\t""\n";
#endif
            if (theOtherNodeTest.type() == unknownNode){
                const ZeroDimensionalNodeUnknown& theOtherNode = *(Nodes->unknownNodeMapAt(M_nodeIndex[theOtherNodeCounter]));
                A.addToCoefficient(equationRow, theOtherNode.variableIndex(), M_parameter);
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< equationRow << "\t"<< theOtherNode.variableIndex() << "\n";
#endif
            } else {
                const ZeroDimensionalNodeKnown & theOtherNode = *(Nodes->knownNodeMapAt(M_nodeIndex[theOtherNodeCounter]));
                C[equationRow] += (M_parameter * theOtherNode.deltaVoltage());
            }
        }
#ifdef HAVE_LIFEV_DEBUG
        A.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        B.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        C.epetraVector().Print(cout);
        cout<< "------------------\n";
#endif

    }
}
// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementPassiveInductor::ZeroDimensionalElementPassiveInductor():
   M_equationRow             (),
   M_variableIndex             ()
{
    M_type      =  inductor;
}

void
ZeroDimensionalElementPassiveInductor::showMe(Int flag)
{
    ZeroDimensionalElementPassive::showMe(flag);
    if (flag==0)
         cout<<"\t parameter = "<<parameter()<<std::endl;
    else
         cout<<"\t EquationRow= "<<M_equationRow<<"\t VariableIndex= "<<M_variableIndex;
}

void
ZeroDimensionalElementPassiveInductor::assignVariableIndex(const Int & index)
{
    M_equationRow       =       index;
    M_variableIndex     =       index;

}
void
ZeroDimensionalElementPassiveInductor::buildABC(matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes)
{
    for (int i = 0; i < 2; i++) {
        const int& theNodeCounter =  (i)%2;
        const ZeroDimensionalNode& theNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theNodeCounter]));
        if (theNodeTest.type() == unknownNode ){
            const ZeroDimensionalNodeUnknown& theNode = *(Nodes->unknownNodeMapAt(M_nodeIndex[theNodeCounter]));
            const int& equationRow = theNode.equationRow();
            B.addToCoefficient(equationRow,M_variableIndex, direction(theNode.id()));
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< equationRow << "\t"<< M_variableIndex << "\t"<< direction(theNode.id()) << "\t""\n";
#endif
        }
#ifdef HAVE_LIFEV_DEBUG
        A.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        B.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        C.epetraVector().Print(cout);
        cout<< "------------------\n";
#endif
    }
    //----------------------------write the Inductor Equations
    A.addToCoefficient(M_equationRow, M_variableIndex, 1.0);
    {
        int i = 0;
        const int& theNodeCounter =  (i)%2;
        const int& theOtherNodeCounter =  (i+1)%2;
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< theNodeCounter << "\t"<< theOtherNodeCounter << "\n";
#endif
        const ZeroDimensionalNode& theNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theNodeCounter]));
        const ZeroDimensionalNode& theOtherNodeTest = *(Nodes->nodeListAt(M_nodeIndex[theOtherNodeCounter]));
        if (theNodeTest.type() == unknownNode ){
            const ZeroDimensionalNodeUnknown& theNode = *(Nodes->unknownNodeMapAt(M_nodeIndex[theNodeCounter]));
            B.addToCoefficient(M_equationRow,theNode.variableIndex(), -M_parameter);
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< theNode.voltage() << "\t"<< theNode.variableIndex() << "\n";
#endif
        }else{
            const ZeroDimensionalNodeKnown & theNode = *(Nodes->knownNodeMapAt(M_nodeIndex[theNodeCounter]));
            C[M_equationRow] += (-M_parameter * theNode.voltage() );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< theNode.voltage() << "\n";
#endif
        }

        if (theOtherNodeTest.type() == unknownNode ){
            const ZeroDimensionalNodeUnknown& theOtherNode = *(Nodes->unknownNodeMapAt(M_nodeIndex[theOtherNodeCounter]));
            B.addToCoefficient(M_equationRow,theOtherNode.variableIndex(), M_parameter);
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< theOtherNode.voltage() << "\t"<< theOtherNode.variableIndex() << "\n";
#endif
        }else{
            const ZeroDimensionalNodeKnown& theOtherNode = *(Nodes->knownNodeMapAt(M_nodeIndex[theOtherNodeCounter]));
            C[M_equationRow] += (+M_parameter * theOtherNode.voltage() );
#ifdef HAVE_LIFEV_DEBUG
            Debug( 8151 )<< theOtherNode.voltage() <<  "\n";
#endif
        }
#ifdef HAVE_LIFEV_DEBUG
        A.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        B.matrixPtr()->Print(cout);
        cout<< "------------------\n";
        C.epetraVector().Print(cout);
        cout<< "------------------\n";
#endif

    }

}


// ==================================================
// Constructors
// ===================================================
ZeroDimensionalElementSource::ZeroDimensionalElementSource():
    M_nodeIndex                             ()

{
}

void
ZeroDimensionalElementSource::showMe(Int flag)
{
    ZeroDimensionalElement::showMe(flag);
    cout<<"\t Node1= "<<nodeIndex()<<std::endl;
}



// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementVoltageSource::ZeroDimensionalElementVoltageSource():
    M_voltage                             (),
    M_deltaVoltage                        ()
{
    M_type      =  voltageSource            ;

}


void
ZeroDimensionalElementVoltageSource::connectElement(zeroDimensionalNodeSPtr_Type & Nodes)
{
    zeroDimensionalNodeKnownPtr_Type    theNode=Nodes->knownNodeMapAt(M_nodeIndex);
    theNode->setelementListIndex(M_id);
    Int otherEndIndex=-1;// in fact this element is not coneccting this node to a new node // Only filling the vector
    theNode->setnodeListIndex(otherEndIndex);
}
void
ZeroDimensionalElementVoltageSource::calculateCurrent(const ZeroDimensionalNodeS& Nodes,const ZeroDimensionalElementS& Elements){
    ZeroDimensionalNodeKnown& theNode =  *Nodes.knownNodeMapAt(M_nodeIndex);
    const vecInt_Type& indexList = theNode.elementListIndex();
    Real current = 0.0 ;
    Int length = indexList.size();
    for (int i = 0; i < length; ++i) {
        const ZeroDimensionalElement& theElement = *Elements.elementListAt(indexList.at(i));
        if (theElement.id()  != M_id){
            Real tmp = theElement.current() * theElement.direction(theNode.id());
#ifdef HAVE_LIFEV_DEBUG
           Debug( 8151 ) << "\n ZeroDimensionalElementVoltageSource::calculateCurrent"<< i <<"\t"<<tmp;
#endif
           current += tmp;
        }
    }
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8151 ) << "\n current is" <<"\t"<<current;
#endif
    M_current = current;
}
// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementCurrentSource::ZeroDimensionalElementCurrentSource()
{
    M_type      =  currentSource;
}


void
ZeroDimensionalElementCurrentSource::connectElement(zeroDimensionalNodeSPtr_Type & Nodes)
{
    Int Index=M_nodeIndex;
    zeroDimensionalNodePtr_Type    theNode=Nodes->nodeListAt(Index);
    theNode->setelementListIndex(M_id);
    Int otherEndIndex=-1;// in fact this element is not coneccting this node to a new node // Only filling the vector
    theNode->setnodeListIndex(otherEndIndex);
}

void
ZeroDimensionalElementCurrentSource::buildABC(__attribute__((unused)) matrix_Type& A,__attribute__((unused)) matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes)
{
    const ZeroDimensionalNode& theNodeTest = *(Nodes->nodeListAt(M_nodeIndex));
    if (theNodeTest.type() == unknownNode ){
        const ZeroDimensionalNodeUnknown& theNode = *(Nodes->unknownNodeMapAt(M_nodeIndex));
        const int& equationRow = theNode.equationRow();
        C[equationRow] += (-M_current) ;
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< "\n"<<equationRow << "\t"<< M_current <<"\n";
        Debug( 8151 )<< "------------------\n"<<M_nodeIndex;
        A.matrixPtr()->Print(cout);
        Debug( 8151 )<< "------------------\n";
        B.matrixPtr()->Print(cout);
        Debug( 8151 )<< "------------------\n";
        C.epetraVector().Print(cout);
        Debug( 8151 )<< "------------------\n";
#endif
    } else {
        std::cerr<<"Error at ZeroDimensionalElementCurrentSource::buildABC, source connected to voltage";
        std::exit(-1);
    }
}


// ===================================================
// Constructors
// ===================================================
ZeroDimensionalNode::ZeroDimensionalNode():
    M_id                          (),
    M_type                        (),
    M_currentBalance              (0),
    M_elementListIndex            (),
    M_nodeListIndex               (),
    M_voltage                     (),
    M_deltaVoltage                ()
{
}

void
ZeroDimensionalNode::showMe(Int flag)
{
    cout<<"Id = "<<id();
    cout<<"\t type = ";
    cout<<enum2string(type())<<"\t";
    if (flag==0)
       {
    for (UInt i=0; i < M_elementListIndex.size();i++)
          cout<<"{"<<M_elementListIndex.at(i)<<","<<M_nodeListIndex.at(i)<<"}  ";
        cout<<std::endl;
       }
}

const std::string
ZeroDimensionalNode::enum2string(const ZeroDimentionalNodeType & type) const
{
    std::string output;
    switch (type)
        {
            case knownNode:
            {
                output = "knownNode";
                break;
            }
            case unknownNode:
            {
                output = "unknownNode";
                break;
            }

        }
    return (output);
}

void
ZeroDimensionalNode::calculateCurrentBalance(const ZeroDimensionalElementS& Elements){
    Real currentBalance = 0.0 ;
    Int length = M_elementListIndex.size();
    for (int i = 0; i < length ; ++i) {
        const ZeroDimensionalElement& theElement = *Elements.elementListAt(M_elementListIndex.at(i));
        currentBalance += theElement.current() * theElement.direction(M_id);
    }
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8151 ) << "\n ZeroDimensionalNode::calculateCurrentBalance"<< currentBalance ;
#endif
    M_currentBalance = currentBalance;
}

// ===================================================
// Constructors
// ===================================================
ZeroDimensionalNodeUnknown::ZeroDimensionalNodeUnknown():
    M_variableIndex                          (),
    M_equationRow                            ()
{
    M_type = unknownNode;
}
void
ZeroDimensionalNodeUnknown::assignVariableIndex(const Int & index)
{
    M_equationRow    =    index;
    M_variableIndex    =    index;

}

void
ZeroDimensionalNodeUnknown::showMe(Int flag)
{
    ZeroDimensionalNode::showMe(flag);
    if (flag==1)
    cout<<"\t EquationRow= "<<M_equationRow<<"\t VariableIndex= "<<M_variableIndex;
}

// ===================================================
// Constructors
// ===================================================
ZeroDimensionalNodeKnown::ZeroDimensionalNodeKnown()
{M_type = knownNode;}
ZeroDimensionalNodeKnown::ZeroDimensionalNodeKnown( const zeroDimensionalElementVoltageSourcePtr_Type& theElement )
{
    M_element= theElement;
    M_type = knownNode;
}

// ===================================================
// Constructors
// ===================================================
ZeroDimensionalElementS::ZeroDimensionalElementS():
    M_elementList                         (new vecZeroDimensionalElementPtr_Type),
    M_resistorList                        (new vecZeroDimensionalElementPassiveResistorPtr_Type),
    M_capacitorList                       (new vecZeroDimensionalElementPassiveCapacitorPtr_Type),
    M_inductorList                        (new vecZeroDimensionalElementPassiveInductorPtr_Type),
    M_diodeList                           (new vecZeroDimensionalElementPassiveDiodePtr_Type),
    M_currentSourceList                   (new vecZeroDimensionalElementCurrentSourcePtr_Type),
    M_voltageSourceList                   (new vecZeroDimensionalElementVoltageSourcePtr_Type),
    M_voltageSourceMap                    (new mapVoltageSource_Type)
{
}

void
ZeroDimensionalElementS::showMe(Int flag)
{
    cout<<"=============== Show all ZeroDimentional Elements ==========="<<std::endl;
    for (iterZeroDimensionalElement_Type theElement = M_elementList->begin();
                    theElement != M_elementList->end(); theElement++) {
        (*theElement)->showMe(flag);
    }
}


// ===================================================
// Constructors
// ===================================================
ZeroDimensionalNodeS::ZeroDimensionalNodeS():
    M_nodeList                            (new vecZeroDimensionalNodePtr_Type),
    M_unknownNodeList                     (new vecZeroDimensionalNodeUnknownPtr_Type),
    M_knownNodeList                       (new vecZeroDimensionalNodeKnownPtr_Type),
    M_knownNodeMap                        (new mapNodeKnown_Type),
    M_unknownNodeMap                      (new mapNodeUnknown_Type)
{
}

void
ZeroDimensionalNodeS::showMe(Int flag)
{
    cout<<"=============== Show all ZeroDimentional Nodes    ==========="<<std::endl;
    for (iterZeroDimensionalNode_Type theNode = M_nodeList->begin();
                    theNode != M_nodeList->end(); theNode++) {
        (*theNode)->showMe(flag);
    }
}

// ===================================================
// Constructors
// ===================================================
ZeroDimensionalCircuitData::ZeroDimensionalCircuitData():
   M_Elements               (new ZeroDimensionalElementS),
   M_Nodes                  (new ZeroDimensionalNodeS),
   M_bc                     ()

{
}

void
ZeroDimensionalCircuitData::showMe(Int flag)
{
    M_Elements->showMe(flag);
    M_Nodes->showMe(flag);
}

// ===================================================
// Methods
// ===================================================
void
ZeroDimensionalCircuitData::buildCircuit(const  char *fileName, bcPtr_Type bc)
{
using namespace std;
    ifstream infile;
    string stringline1;
    string stringline2;
    string stringline3;
    string stringtmp;
    stringstream stringStreamLine1 (stringstream::in | stringstream::out);
    stringstream stringStreamLine2 (stringstream::in | stringstream::out);

    bool boolElementType[ZERO_DIMENTIONAL_DEFINED_ELEMENTS];
    Int numberOfNodes = 0;
    Int numberOfElements = 0;
    Int numberOfTerminalNodes = 0;
    Int ID = 0;
    Int Node1 = 0;
    Int Node2 = 0;
    Real parameter1 = 0.0;
    Real forwardBias = 0.0;
    Real alpha = 0.0;
    Real beta = 0.0;;


    std::vector<ZeroDimentionalNodeType> nodesType;
    vecInt_Type                     nodesConnectingSource;
    vecInt_Type                    terminalNodes;

    infile.open (fileName, ifstream::in);
    if (infile.is_open())
    {
      // ----------Read the first line----------------------------
        getline(infile,stringline1);
        stringStreamLine1 << stringline1;
        stringStreamLine1 >> stringtmp;
        numberOfElements = std::atoi(stringtmp.c_str());
        stringStreamLine1 >> stringtmp;
        numberOfNodes = std::atoi(stringtmp.c_str());
        stringStreamLine1 >> stringtmp;
        numberOfTerminalNodes = std::atoi(stringtmp.c_str());
        getline(infile,stringline2);
        int nodeId = -1;
      // ----------Read the second line---------------------------
        stringStreamLine2 << stringline2;
        for (int i=0;i<numberOfTerminalNodes;i++){
            stringStreamLine2 >> stringtmp;
            nodeId = std::atoi(stringtmp.c_str());
            terminalNodes.push_back(nodeId);
        }
    //---------------------------------------------------------
        for (int i=0;i<numberOfNodes;i++){
            nodesType.push_back(unknownNode);
            nodesConnectingSource.push_back(-1);
        }
      // ----------Read Elements-----------------------------------
        for(int i=0;i<numberOfElements;i++)
        {
            stringstream stringStreamLine3 (stringstream::in | stringstream::out);
            getline(infile,stringline3);
            stringStreamLine3 << stringline3;
            stringStreamLine3 >> stringtmp;
            ID = std::atoi(stringtmp.c_str());
            stringStreamLine3 >> stringtmp;
            boolElementType[0] = stringtmp.compare("resistor");
            boolElementType[1] = stringtmp.compare("capacitor");
            boolElementType[2] = stringtmp.compare("inductor");
            boolElementType[3] = stringtmp.compare("voltageSource");
            boolElementType[4] = stringtmp.compare("currentSource");
            boolElementType[5] = stringtmp.compare("diode");

            stringStreamLine3 >> stringtmp;
            Node1 = std::atoi(stringtmp.c_str());

            if (!boolElementType[0])//resistor
            {
                stringStreamLine3 >> stringtmp;
                Node2 = std::atoi(stringtmp.c_str());
                stringStreamLine3 >> stringtmp;
                parameter1 = std::atof(stringtmp.c_str());
                createElementResistor(ID,Node1,Node2,parameter1);
            }
            if (!boolElementType[1])//capacitor
            {
                stringStreamLine3 >> stringtmp;
                Node2 = std::atoi(stringtmp.c_str());
                stringStreamLine3 >> stringtmp;
                parameter1 = std::atof(stringtmp.c_str());
                createElementCapacitor(ID,Node1,Node2,parameter1);
            }
            if (!boolElementType[2])//inductor
            {
                stringStreamLine3 >> stringtmp;
                Node2 = std::atoi(stringtmp.c_str());
                stringStreamLine3 >> stringtmp;
                parameter1 = std::atof(stringtmp.c_str());
                createElementInductor(ID,Node1,Node2,parameter1);
            }
            if (!boolElementType[5])//diode
            {
                stringStreamLine3 >> stringtmp;
                Node2 = std::atoi(stringtmp.c_str());
                stringStreamLine3 >> stringtmp;
                forwardBias = std::atof(stringtmp.c_str());
                stringStreamLine3 >> stringtmp;
                alpha = std::atof(stringtmp.c_str());
                stringStreamLine3 >> stringtmp;
                beta = std::atof(stringtmp.c_str());
                createElementDiode(ID,Node1,Node2,forwardBias,alpha,beta);
            }
        }
        infile.close();
    }
    else
    {
        cerr << "Error opening circuit file";
        exit(-1);
    }

    for(iterVecInt_Type theNodeId = terminalNodes.begin();
                    theNodeId != terminalNodes.end(); theNodeId++ ){
      // create Source Elements -----------------------------------------------
      switch (bc->bc(*theNodeId).bcType()) {
          case OneDFSI::S://create voltage source
              nodesConnectingSource.at(*theNodeId) = createElementVoltageSource(*theNodeId);
              nodesType.at(*theNodeId)=knownNode;
              break;
          case OneDFSI::Q:// current source
              createElementCurrentSource(*theNodeId);
              break;
          default:
              break;
        }
    }
    //create Nodes --------------------------------------------------------
    for (int i=0;i<numberOfNodes;i++)
    {
        if (nodesType.at(i)==knownNode)
        {
          //Create known Node and connect the voltage source to the node
      createKnownNode(i,M_Elements->voltageSourceMap(nodesConnectingSource.at(i)));
        }
        else
        {
      createUnknownNode(i);
        }
    }
    //connect elements to the nodes
    for (Int i=0;i<M_Elements->elementCounter();i++){
        M_Elements->elementListAt(i)->connectElement(M_Nodes);
    }
    //set BC to source elements
    fixBC(bc);
}

void
ZeroDimensionalCircuitData::createElementResistor(Int ID,Int node1,Int node2,Real parameter)
{
    zeroDimensionalElementPassiveResistorPtr_Type theElement(new ZeroDimensionalElementPassiveResistor());
    theElement->setid(ID);
    if (M_Elements->elementCounter() != ID){
        cerr << "Error: Element Id error at  "<< ID ;
        exit(-1);
    }
    theElement->setnodeIndex(node1);
    theElement->setnodeIndex(node2);
    if (parameter <= 0){
        cerr << "Error: Resistance value <=0, ID =  "<< ID ;
        exit(-1);
    }
    theElement->setparameter(1.0/parameter);
    M_Elements->setelementList(theElement);
    M_Elements->setresistorList(theElement);
}
void
ZeroDimensionalCircuitData::createElementCapacitor(Int ID,Int node1,Int node2,Real parameter)
{
    zeroDimensionalElementPassiveCapacitorPtr_Type theElement(new ZeroDimensionalElementPassiveCapacitor());
    theElement->setid(ID);
    if (M_Elements->elementCounter() != ID){
        cerr << "Error: Element Id error at  "<< ID ;
        exit(-1);
    }
    theElement->setnodeIndex(node1);
    theElement->setnodeIndex(node2);
    theElement->setparameter(parameter);
    M_Elements->setelementList(theElement);
    M_Elements->setcapacitorList(theElement);
}
void
ZeroDimensionalCircuitData::createElementInductor(Int ID,Int node1,Int node2,Real parameter)
{
    zeroDimensionalElementPassiveInductorPtr_Type theElement(new ZeroDimensionalElementPassiveInductor());
    theElement->setid(ID);

    if (M_Elements->elementCounter() != ID){
        cerr << "Error: Element Id error at  "<< ID ;
        exit(-1);
    }
    theElement->setnodeIndex(node1);
    theElement->setnodeIndex(node2);
    theElement->setparameter(1.0/parameter);
    M_Elements->setelementList(theElement);
    M_Elements->setinductorList(theElement);
}
void
ZeroDimensionalCircuitData::createElementDiode(Int ID, Int node1, Int node2, Real forwardBias, Real alpha, Real beta)
{
    zeroDimensionalElementPassiveDiodePtr_Type theElement(new ZeroDimensionalElementPassiveDiode());
    theElement->setid(ID);
    if (M_Elements->elementCounter() != ID){
        cerr << "Error: Element Id error at  "<< ID ;
        exit(-1);
    }
    theElement->setnodeIndex(node1);
    theElement->setnodeIndex(node2);
    theElement->setparameter(0);
    theElement->setforwardBias(forwardBias);
    theElement->setalpha(alpha);
    theElement->setbeta(beta);
    M_Elements->setelementList(theElement);
    M_Elements->setdiodeList(theElement);
}

Int
ZeroDimensionalCircuitData::createElementVoltageSource(Int node1)
{
    zeroDimensionalElementVoltageSourcePtr_Type theElement(new ZeroDimensionalElementVoltageSource());
    theElement->setid(M_Elements->elementCounter());
    theElement->setnodeIndex(node1);
    M_Elements->setelementList(theElement);
    M_Elements->setvoltageSourceList(theElement);
    M_Elements->setvoltageSourceMap(theElement->id(),theElement);
    return theElement->id();
}
void
ZeroDimensionalCircuitData::createElementCurrentSource(Int node1)
{
    zeroDimensionalElementCurrentSourcePtr_Type theElement(new ZeroDimensionalElementCurrentSource());
    theElement->setid(M_Elements->elementCounter());
    theElement->setnodeIndex(node1);
    M_Elements->setelementList(theElement);
    M_Elements->setcurrentSourceList(theElement);
}


void
ZeroDimensionalCircuitData::createUnknownNode(const Int & id)
{
    zeroDimensionalNodeUnknownPtr_Type     theNode(new ZeroDimensionalNodeUnknown());
    theNode->setid(id);
    M_Nodes->setnodeList(theNode);
    M_Nodes->setunknownNodeList(theNode);
    M_Nodes->setunknownNodeMap(theNode->id() , theNode);
}

void
ZeroDimensionalCircuitData::createKnownNode(const Int & id, const zeroDimensionalElementVoltageSourcePtr_Type & theElement)
{
    zeroDimensionalNodeKnownPtr_Type    theNode(new ZeroDimensionalNodeKnown(theElement));
    theNode->setid(id);
    M_Nodes->setnodeList(theNode);
    M_Nodes->setknownNodeList(theNode);
    M_Nodes->setknownNodeMap(theNode->id() , theNode);
}

void
ZeroDimensionalCircuitData::createKnownNode(const Int & id)
{
    zeroDimensionalNodeKnownPtr_Type     theNode(new ZeroDimensionalNodeKnown());
    theNode->setid(id);
    M_Nodes->setnodeList(theNode);
    M_Nodes->setknownNodeList(theNode);
    M_Nodes->setknownNodeMap(theNode->id() , theNode);
}
void
ZeroDimensionalCircuitData::fixBC(bcPtr_Type bc){
    M_bc = bc;
    //set BC hander in source elements

    ptrVecZeroDimensionalElementCurrentSourcePtr_Type currentElementList   = M_Elements->currentSourceList();
    for (iterZeroDimensionalElementCurrentSource_Type theElement = currentElementList->begin(); theElement != currentElementList->end(); theElement++) {
        (*theElement)->setbc(bc);
    }
    ptrVecZeroDimensionalElementVoltageSourcePtr_Type voltageElementList   = M_Elements->voltageSourceList();
    for (iterZeroDimensionalElementVoltageSourcePtr_Type theElement = voltageElementList->begin(); theElement != voltageElementList->end(); theElement++) {
        (*theElement)->setbc(bc);
    }
}

OutPutFormat::OutPutFormat(
                std::string  width,
                std::string  precision,
                std::string  whiteSpace,
                int bufferSize){
    M_width = width;
    M_precision = precision;
    M_whiteSpace = whiteSpace;
    M_buffer = new char[bufferSize];
    M_formatDouble = "% " + M_width + "." + M_precision + "f";
    M_formatInteger = "% " + M_width + "d";
}
void
OutPutFormat::writeDataFormat(const double& number,std::ofstream & stream, const EndLine& flag){
    sprintf (M_buffer, M_formatDouble.data(), number);
    stream << M_buffer;
    switch (flag) {
        case newLine:
            stream << std::endl;
            break;
        case space:
            stream << M_whiteSpace;
            break;
        case nothing:
            break;
        default:
            std::cerr<< "no flag at OutPutFormat::writeDataFormat";
            break;
    }
}
void
OutPutFormat::writeDataFormat(const int& number,std::ofstream & stream, const EndLine& flag){
    sprintf (M_buffer, M_formatInteger.data(), number);
    stream << M_buffer;
    switch (flag) {
        case newLine:
            stream << std::endl;
            break;
        case space:
            stream << M_whiteSpace;
            break;
        case nothing:
            break;
        default:
            std::cerr<< "no flag at OutPutFormat::writeDataFormat";
            break;
    }
}
void
OutPutFormat::writeNewLine(std::ofstream & stream){
    stream << std::endl;
}
OutPutFormat::~OutPutFormat(){
    delete [] M_buffer;
}

void
ZeroDimensionalCircuitData::updateCircuitDataFromY(const double& t, const Epetra_Vector* x,const Epetra_Vector* x_dot){
    const ptrVecZeroDimensionalNodeUnknownPtr_Type& unknownNodeList = M_Nodes->unknownNodeList();

    //update unknown Nodes from solution.
    for (iterZeroDimensionalNodeUnknown_Type theNode = unknownNodeList ->begin();
                    theNode != unknownNodeList->end(); theNode++) {
        const int& variableIndex = (*theNode)->variableIndex();
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< (*theNode)->id()<<"   ---------\n";
        Debug( 8151 )<< variableIndex<<"   ---------\n";
#endif
        (*theNode)->setvoltage((*x)[variableIndex]);
        (*theNode)->setdeltaVoltage((*x_dot)[variableIndex]);
    }

    //unpdate inductor unknown current
    const ptrVecZeroDimensionalElementPassiveInductorPtr_Type& inductorList = M_Elements->inductorList();
    for (iterZeroDimensionalElementPassiveInductor_Type theInductor = inductorList ->begin();
                    theInductor != inductorList->end(); theInductor++) {
        const int& variableIndex = (*theInductor)->variableIndex();
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< (*theInductor)->id()<<"   ---------\n";
        Debug( 8151 )<< variableIndex<<"   ---------\n";
#endif
        (*theInductor)->setcurrent((*x)[variableIndex]);
        (*theInductor)->setdeltaCurrent((*x_dot)[variableIndex]);
    }

    //update BCs by time
    const ptrVecZeroDimensionalNodeKnownPtr_Type& knownNodeList = M_Nodes->knownNodeList();
    for (iterZeroDimensionalNodeKnown_Type theNode = knownNodeList ->begin();
                    theNode != knownNodeList->end(); theNode++) {
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< (*theNode)->id()<<"   ---------\n";

#endif
    (*theNode)->setvoltageByTime(t);
        (*theNode)->setdeltaVoltageByTime(t);
    }
    const ptrVecZeroDimensionalElementCurrentSourcePtr_Type& currentElementList = M_Elements->currentSourceList();
    for (iterZeroDimensionalElementCurrentSource_Type theElement = currentElementList ->begin();
                    theElement != currentElementList->end(); theElement++) {
#ifdef HAVE_LIFEV_DEBUG
        Debug( 8151 )<< (*theElement)->id()<<"   ---------\n";
#endif
        (*theElement)->setcurrentByTime(t);
    }
}
void
ZeroDimensionalCircuitData::deepUpdateFromY(const double& t, const Epetra_Vector& y,const Epetra_Vector& yp){
  //First invoke the updateCircuitDataFromY method (light update)
    updateCircuitDataFromY(t,&y,&yp);
#ifdef HAVE_LIFEV_DEBUG
    y.Print(cout);
    yp.Print(cout);
#endif

    //deep update (mainly currents), iterate over all elements
    const ptrVecZeroDimensionalElementPtr_Type& elementList = M_Elements->elementList();
    for (iterZeroDimensionalElement_Type theElement = elementList ->begin();
                    theElement != elementList->end(); theElement++) {
        (*theElement)->deepUpdate(*M_Nodes);
    }

    // Now we can compute voltage source current
    const ptrVecZeroDimensionalElementVoltageSourcePtr_Type& voltageList = M_Elements->voltageSourceList();
    for (iterZeroDimensionalElementVoltageSourcePtr_Type theElement = voltageList ->begin();
                    theElement != voltageList->end(); theElement++) {
        (*theElement)->calculateCurrent(*M_Nodes,*M_Elements);
    }

    // check the conservation of current at each node.
    const ptrVecZeroDimensionalNodePtr_Type& nodeList = M_Nodes->nodeList();
    for (iterZeroDimensionalNode_Type theNode = nodeList ->begin();
                    theNode != nodeList->end(); theNode++) {
        (*theNode)->calculateCurrentBalance(*M_Elements);
    }

}
void
ZeroDimensionalCircuitData::updateABC(matrix_Type& A,matrix_Type& B,vector_Type& C){
    A.openCrsMatrix();
    B.openCrsMatrix();
    A.matrixPtr()->PutScalar(0.0);
    B.matrixPtr()->PutScalar(0.0);
    C.epetraVector().PutScalar(0.0);

    const ptrVecZeroDimensionalElementPtr_Type& elementList = M_Elements->elementList();

    //iterate over all elements and compute each element's contribution on A,B,C
    for (iterZeroDimensionalElement_Type theElement = elementList ->begin();
                    theElement != elementList->end(); theElement++) {

        (*theElement)->buildABC(A,B,C,M_Nodes);
    }
    A.globalAssemble();
    B.globalAssemble();
}

} // LifeV namespace
