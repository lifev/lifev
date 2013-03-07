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
 *  @brief File containing a class for 0D model circuit data handling.
 *  @version alpha (experimental)
 *
 *  @date 26-09-2011
 *  @author Mahmoud Jafargholi <mahmoud.jafargholi@epfl.ch>
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ZeroDimensionalCircuitData_H
#define ZeroDimensionalCircuitData_H 1

// LIFEV
#include <lifev/core/array/MatrixEpetra.hpp>

// MATHCARD
#include <lifev/zero_dimensional/fem/ZeroDimensionalBCHandler.hpp>
#include <lifev/zero_dimensional/solver/ZeroDimensionalDefinitions.hpp>
#include <lifev/zero_dimensional/fem/ZeroDimensionalBC.hpp>

namespace LifeV
{

// TODO Move forward declarations in ZeroDimensionalDefinitions
// TODO Move type definitions inside classes

//! A container class for all node objects
class ZeroDimensionalNodeS;

//! A container class for all element obkects
class ZeroDimensionalElementS;

typedef boost::shared_ptr< ZeroDimensionalElementS >                    zeroDimensionalElementSPtr_Type;
typedef boost::shared_ptr< ZeroDimensionalNodeS >                       zeroDimensionalNodeSPtr_Type;
typedef std::vector<Int>                                                vecInt_Type;
typedef vecInt_Type::iterator                                           iterVecInt_Type;
typedef ZeroDimensionalBCHandler                                        bc_Type;
typedef boost::shared_ptr< bc_Type >                                    bcPtr_Type;
typedef MatrixEpetra<Real>                                              matrix_Type;
typedef VectorEpetra                                                    vector_Type;
typedef Epetra_Vector                                                   vectorEpetra_Type;
typedef boost::shared_ptr< matrix_Type >                                matrixPtr_Type;
typedef boost::shared_ptr< vector_Type >                                vectorPtr_Type;
typedef boost::shared_ptr<vectorEpetra_Type >                           vectorEpetraPtr_Type;

//! ZeroDimensionalElement - The base element class.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElement
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalElement() {}

    //! Destructor
    virtual ~ZeroDimensionalElement() {}

    //@}


    //! @name Methods
    //@{

    //! Display some information.
    virtual void showMe ( const Int& flag = 0 );

    const std::string enum2string ( const ZeroDimensionalElementType& type );

    //! Connect elements to the nodes.
    /*!
     * After all emenets and nodes are created, each element will call this
     * method to connect itse;f to the nodes.
     */
    virtual void connectElement ( zeroDimensionalNodeSPtr_Type& nodes ) = 0;

    //! Contribution of the element of matrix \bf{A} and \bf{B} and vector \bf{C}.
    /*!
     * After updating the BCs ( or Terminal nodes ) this each element will invoke
     * this method to compute it's contribution on matrices.
     */
    virtual void buildABC ( matrix_Type& /*A*/, matrix_Type& /*B*/, vector_Type& /*C*/, const zeroDimensionalNodeSPtr_Type& /*Nodes*/ ) {}

    //! Compute outputs (currents and voltages) from the solution vector after each succesful iteration.
    /*!
     * After each time step, when Rythmos solver is succesfully finishes, this method will compute
     * finial outputs ( for exmple currents) from the finial solution vector.
     */
    virtual void extractSolution ( const ZeroDimensionalNodeS& /*nodes*/ ) {}

    //! This method specifies the convention of current direction in an element.
    /*!
     * @param A node index connected to the element.
     * @return +1 if the current convention is toward the iniput node and -1 otherwise.
     */
    virtual Real direction ( const Int& nodeId ) const = 0;

    //@}


    //! @name Set Methods
    //@{

    void setId (const Int& id )
    {
        M_id = id;
    }

    void setCurrent (const Real& current )
    {
        M_current = current;
    }

    //! Set derivative of current respect to time.
    void setDeltaCurrent (const Real& deltaCurrent )
    {
        M_deltaCurrent = deltaCurrent;
    }

    //@}


    //! @name Get Methods
    //@{

    const Int& id() const
    {
        return M_id;
    }

    const ZeroDimensionalElementType& type() const
    {
        return M_type;
    }

    const Real& current() const
    {
        return M_current;
    }

    //! Get derivative of current respect to time.
    const Real& deltaCurrent() const
    {
        return M_deltaCurrent;
    }

    //@}

protected:

    Int                             M_id;
    ZeroDimensionalElementType      M_type; //= 'Resistor';%'Capacitor' ,'Inductor','Voltage Source','Current Source' 'Diode'
    Real                            M_current;
    Real                            M_deltaCurrent;
};


// TODO Move type definitions inside classes
typedef boost::shared_ptr<ZeroDimensionalElement>                              zeroDimensionalElementPtr_Type;
typedef std::vector<zeroDimensionalElementPtr_Type>                            vecZeroDimensionalElementPtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalElementPtr_Type>                   ptrVecZeroDimensionalElementPtr_Type;
typedef vecZeroDimensionalElementPtr_Type::iterator                            iterZeroDimensionalElement_Type;

//! ZeroDimensionalElementPassive - A class for passive elements.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementPassive: public ZeroDimensionalElement
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalElementPassive();

    //! Destructor
    virtual ~ZeroDimensionalElementPassive() {}

    //@}


    //! @name Methods
    //@{

    //! Show some information.
    void showMe ( const Int& flag = 0 );

    //! Impleaments the abstarct class for passive elements.
    void connectElement ( zeroDimensionalNodeSPtr_Type& nodes );

    //@}


    //! @name Set Methods
    //@{

    //! set parameter (1/R, 1/L, C, 1/R_{eff})
    void setParameter ( const Real& parameter )
    {
        M_parameter = parameter;
    }

    //! add the node to the list.
    /*!
     * @param node index.
     */
    void setNodeIndex ( const Int& index )
    {
        M_nodeIndex.push_back (index);
    }

    //@}


    //! @name Get Methods
    //@{

    //! get the parameter (1/R, 1/L, C, 1/R_{eff})
    const Real& parameter() const
    {
        return M_parameter;
    }

    //! get the node index connected to the node.
    /*!
     * @param \it{i}th node connected to the elelemt.
     * @return  Index of \it{i}th node connected to the element.
     */
    const Int& nodeIndex ( const Int& position ) const
    {
        return M_nodeIndex.at (position);
    }

    Real direction ( const Int& nodeId ) const;

    //@}

protected:

    //parameter= 'Resistor';%'Capacitor' ,'Inductor','Diode'
    //parameter=  1/R      ;%C            ,1/L      , 1/R_{eff}

    Real        M_parameter;
    vecInt_Type M_nodeIndex; //Index of connected nodes
};




//! ZeroDimensionalElementPassiveResistor - Resistor.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementPassiveResistor: public ZeroDimensionalElementPassive
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Contructor.
    explicit ZeroDimensionalElementPassiveResistor();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveResistor() {}

    //@}


    //! @name Methods
    //@{

    void showMe ( const Int& flag = 0 );

    void buildABC ( matrix_Type& A, matrix_Type& B, vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes );

    void extractSolution ( const ZeroDimensionalNodeS& nodes );

    //@}
};



//! ZerodimentionalElementPassiveDiode - Diode.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementPassiveDiode: public ZeroDimensionalElementPassiveResistor
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalElementPassiveDiode();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveDiode() {}

    //@}


    //! @name Methods
    //@{

    void showMe ( const Int& flag = 0 );

    void extractSolution ( const ZeroDimensionalNodeS& nodes );

    void buildABC ( matrix_Type& A, matrix_Type& B, vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes );

    //@}


    //! @name Set Methods
    //@{

    //!current = beta * exp(alpha * (voltage - forwardBias )) - (beta * exp(alpha * ( - forwardBias )))
    void setalpha (const Real& alpha )
    {
        M_alpha = alpha;
    }

    void setbeta (const Real& beta )
    {
        M_beta = beta;
    }

    void setforwardBias (const Real& forwardBias )
    {
        M_forwardBias = forwardBias;
    }

    //@}


    //! @name Get Methods
    //@{

    const Real& alpha() const
    {
        return M_alpha;
    }

    const Real& beta() const
    {
        return M_beta;
    }

    const Real& forwardBias() const
    {
        return M_forwardBias;
    }

    //@}

protected:

    //! calculate the effective resistance.
    /*!
     * @param voltage difference
     * @return effective ressitance
     */
    void calculateEffectiveResistance (const Real& voltage);

    Real            M_alpha; //current = beta * exp(alpha * (voltage - forwardBias )) - (beta * exp(alpha * ( - forwardBias )))
    Real            M_beta;
    Real            M_forwardBias;
};



//! ZerodimentionalElementPassiveCapacitor - Capacitor.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementPassiveCapacitor: public ZeroDimensionalElementPassive
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalElementPassiveCapacitor();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveCapacitor() {}

    //@}


    //! @name Methods
    //@{

    void showMe ( const Int& flag = 0 );

    void extractSolution ( const ZeroDimensionalNodeS& nodes );

    void buildABC ( matrix_Type& A, matrix_Type& B, vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes );

    //@}
};


//! ZeroDimensionalElementPassiveInductor - Inductor.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementPassiveInductor: public ZeroDimensionalElementPassive
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalElementPassiveInductor();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveInductor() {}

    //@}


    //! @name Methods
    //@{

    //! Set the variable index and equation row index for Inductor.
    /*!
     *Current in Inductor is an unknown.
     */
    virtual void assignVariableIndex ( const Int& index );

    void showMe ( const Int& flag = 0 );

    void buildABC ( matrix_Type& A, matrix_Type& B, vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes );

    //@}


    //! @name Get Methods
    //@{

    //! get equation row for in matrix A,B and C.
    const Int& equationRow() const
    {
        return M_equationRow;
    }

    //! get variable index in solution vector  x  and \dot{x}
    const Int& variableIndex() const
    {
        return M_variableIndex;
    }

    //@}

protected:
    Int            M_equationRow;
    Int            M_variableIndex;
};



//! ZeroDimensionalElementSource - Base class for source elements.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementSource: public ZeroDimensionalElement
{

public:

    //! @name Constructors & Destructor
    //@{

    //!Constructor
    explicit ZeroDimensionalElementSource();

    //! Destructor
    virtual ~ZeroDimensionalElementSource() {}

    //@}


    //! @name Methods
    //@{

    void showMe ( const Int& flag = 0 );

    //@}


    //! @name Set Methods
    //@{

    void setNodeIndex (const Real& index )
    {
        M_nodeIndex = index;
    }

    //! Set BC handler.
    void setBC ( const bcPtr_Type& bc)
    {
        M_bc = bc;
    }

    //@}


    //! @name Get Methods
    //@{

    Int nodeIndex() const
    {
        return M_nodeIndex;
    }

    Real direction ( const Int& /*nodeId*/ ) const
    {
        return -1.0;
    }

    //@}

protected:

    Int            M_nodeIndex; //Index of connected node
    bcPtr_Type     M_bc;
};



//! ZeroDimensionalElementVoltageSource - Voltage Source.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementVoltageSource: public ZeroDimensionalElementSource
{

public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalElementVoltageSource();

    //! Destructor
    virtual ~ZeroDimensionalElementVoltageSource() {}

    //@}


    //! @name Get Methods
    //@{

    void connectElement (zeroDimensionalNodeSPtr_Type& Nodes);

    //! calculate current passing outward in voltage source.
    /*!
     *  This method can be called after all elements invoked extractSolution method.
     */
    void calculateCurrent ( const ZeroDimensionalNodeS& Nodes, const ZeroDimensionalElementS& Elements );

    //@}


    //! @name Set Methods
    //@{

    //! Update voltage source by time.
    void setVoltageByTime ( const Real& time )
    {
        M_voltage = M_bc->bc ( M_nodeIndex ).evaluate (time);
    }

    //! Update \frac{\partial voltage}{\partial t} by time.
    void setDeltaVoltageByTime ( const Real& time )
    {
        M_deltaVoltage = M_bc->bc ( M_nodeIndex + BC_CONSTANT ).evaluate (time);
    }

    //@}


    //! @name Get Methods
    //@{

    const Real& voltage() const
    {
        return M_voltage;
    }

    Real deltaVoltage() const
    {
        return M_deltaVoltage;
    }

    Real voltageByTime (const Real& time) const
    {
        return M_bc->bc ( M_nodeIndex ).evaluate (time);
    }

    Real deltaVoltageByTime (const Real& time) const
    {
        return M_bc->bc ( M_nodeIndex + BC_CONSTANT ).evaluate (time);
    }

    //@}

protected:

    Real    M_voltage                     ; //voltage at time t_{n}
    Real    M_deltaVoltage                ; //\frac{\mathrm{d \text{ voltage}} }{\mathrm{d} t}
};



//! ZeroDimensionalElementCurrentSource - Current Source.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementCurrentSource: public ZeroDimensionalElementSource
{

public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor.
    explicit ZeroDimensionalElementCurrentSource();

    //! Destructor
    virtual ~ZeroDimensionalElementCurrentSource() {}

    //@}


    //! @name Methods
    //@{

    void connectElement ( zeroDimensionalNodeSPtr_Type& Nodes );

    void buildABC ( matrix_Type& A, matrix_Type& B, vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes );

    //@}


    //! @name set Methods
    //@{

    void setCurrentByTime (const Real& time )
    {
        M_current = M_bc->bc ( M_nodeIndex ).evaluate (time);
    }

    //@}


    //! @name Get Methods
    //@{

    Real currentByTime (const Real& time ) const
    {
        return M_bc->bc ( M_nodeIndex ).evaluate (time);
    }

    Real current() const
    {
        return M_current;
    }

    //@}

};













// TODO Move type definitions inside classes
typedef boost::shared_ptr<ZeroDimensionalElementPassiveResistor>        zeroDimensionalElementPassiveResistorPtr_Type;
typedef boost::shared_ptr<ZeroDimensionalElementPassiveCapacitor>       zeroDimensionalElementPassiveCapacitorPtr_Type;
typedef boost::shared_ptr<ZeroDimensionalElementPassiveInductor>        zeroDimensionalElementPassiveInductorPtr_Type;
typedef boost::shared_ptr<ZeroDimensionalElementPassiveDiode>           zeroDimensionalElementPassiveDiodePtr_Type;
typedef boost::shared_ptr<ZeroDimensionalElementCurrentSource>          zeroDimensionalElementCurrentSourcePtr_Type;
typedef boost::shared_ptr<ZeroDimensionalElementVoltageSource>          zeroDimensionalElementVoltageSourcePtr_Type;


typedef std::vector<zeroDimensionalElementPassiveResistorPtr_Type>      vecZeroDimensionalElementPassiveResistorPtr_Type;
typedef std::vector<zeroDimensionalElementPassiveCapacitorPtr_Type>     vecZeroDimensionalElementPassiveCapacitorPtr_Type;
typedef std::vector<zeroDimensionalElementPassiveInductorPtr_Type>      vecZeroDimensionalElementPassiveInductorPtr_Type;
typedef std::vector<zeroDimensionalElementPassiveDiodePtr_Type>         vecZeroDimensionalElementPassiveDiodePtr_Type;
typedef std::vector<zeroDimensionalElementCurrentSourcePtr_Type>        vecZeroDimensionalElementCurrentSourcePtr_Type;
typedef std::vector<zeroDimensionalElementVoltageSourcePtr_Type>        vecZeroDimensionalElementVoltageSourcePtr_Type;


typedef boost::shared_ptr<vecZeroDimensionalElementPassiveResistorPtr_Type>     ptrVecZeroDimensionalElementPassiveResistorPtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalElementPassiveCapacitorPtr_Type>    ptrVecZeroDimensionalElementPassiveCapacitorPtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalElementPassiveInductorPtr_Type>     ptrVecZeroDimensionalElementPassiveInductorPtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalElementPassiveDiodePtr_Type>        ptrVecZeroDimensionalElementPassiveDiodePtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalElementCurrentSourcePtr_Type>       ptrVecZeroDimensionalElementCurrentSourcePtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalElementVoltageSourcePtr_Type>       ptrVecZeroDimensionalElementVoltageSourcePtr_Type;

typedef vecZeroDimensionalElementPassiveResistorPtr_Type::iterator              iterZeroDimensionalElementPassiveResistor_Type;
typedef vecZeroDimensionalElementPassiveCapacitorPtr_Type::iterator             iterZeroDimensionalElementPassiveCapacitor_Type;
typedef vecZeroDimensionalElementPassiveInductorPtr_Type::iterator              iterZeroDimensionalElementPassiveInductor_Type;
typedef vecZeroDimensionalElementPassiveDiodePtr_Type::iterator                 iterZeroDimensionalElementPassiveDiode_Type;
typedef vecZeroDimensionalElementCurrentSourcePtr_Type::iterator                iterZeroDimensionalElementCurrentSource_Type;
typedef vecZeroDimensionalElementVoltageSourcePtr_Type::iterator                iterZeroDimensionalElementVoltageSourcePtr_Type;


//! ZeroDimensionalNode - The base node class.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalNode
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalNode();

    //! Destructor
    virtual ~ZeroDimensionalNode() {}

    //@}


    //! @name Methods
    //@{

    //! Calculate current balance at node.
    /*!
     * After updating current in all elements, we can verify the balance of current flow at each node.
     */
    void calculateCurrentBalance ( const ZeroDimensionalElementS& Elements );

    virtual void showMe ( const Int& flag = 0 );

    //@}


    //! @name Set Methods
    //@{

    const std::string enum2string ( const ZeroDimensionalNodeType& type ) const;

    void setId ( const Int& id )
    {
        M_id = id;
    }

    //! add an element index to the elelemt list.
    void setElementListIndex ( const Int& index )
    {
        M_elementListIndex.push_back (index);
    }

    //! add an node index which is connected by an element in element list.
    /*!
     * Each elelemnt in element list, coonects this node to another node ( except source elementt). nodeList is a container for conecting nodes.
     * If the element connected to this node has only one terminal ( like voltage source and current source), the connecting index would be -1.
     */
    void setNodeListIndex (const Int& index )
    {
        M_nodeListIndex.push_back (index);
    }

    virtual void setVoltage (const Real& voltage )
    {
        M_voltage = voltage;
    }

    virtual void setDeltaVoltage (const Real& deltaVoltage )
    {
        M_deltaVoltage = deltaVoltage;
    }

    //@}


    //! @name Get Methods
    //@{

    const Int& id() const
    {
        return M_id;
    }

    const ZeroDimensionalNodeType& type() const
    {
        return M_type;
    }

    const Int& elementListIndexAt ( const Int& position ) const
    {
        return M_elementListIndex.at (position);
    }

    const vecInt_Type& elementListIndex() const
    {
        return M_elementListIndex;
    }

    const Int& nodeListIndexAt ( const Int& position ) const
    {
        return M_nodeListIndex.at (position);
    }

    virtual const Real& voltage() const
    {
        return M_voltage;
    }

    virtual Real deltaVoltage() const
    {
        return M_deltaVoltage;
    }

    const Real& currentBalance() const
    {
        return M_currentBalance;
    }

    //@}

protected:

    Int                             M_id;
    ZeroDimensionalNodeType         M_type;             //= 'Known';%'Unknown'
    Real                            M_currentBalance;   //sum of currents over all branches
    vecInt_Type                     M_elementListIndex; // List of id(s) of connected Elements to this Node
    vecInt_Type                     M_nodeListIndex;    // List of id(s) of connected Nodes to this Node
    Real                            M_voltage;
    Real                            M_deltaVoltage;
};


//! ZeroDimensionalNodeUnknown - This class defines the unknown node class.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalNodeUnknown: public ZeroDimensionalNode
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalNodeUnknown();

    //! Destructor
    virtual ~ZeroDimensionalNodeUnknown() {}

    //@}


    //! @name Methods
    //@{

    void showMe ( const Int& flag = 0 );

    //@}


    //! @name Set Methods
    //@{

    //! assign the index of the unknown voltage.
    void assignVariableIndex (const Int& index) ;

    //@}


    //! @name Get Methods
    //@{

    const Int& variableIndex() const
    {
        return M_variableIndex;
    }

    const Int& equationRow() const
    {
        return M_equationRow;
    }

    //@}

protected:
    Int                 M_variableIndex; // Index of the variable
    Int                 M_equationRow;   // #Row(equation) in the Matrix
};



//! ZeroDimensionalNodeKnown - This class defines the known node class. A Voltage Source element is connected to this class.
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalNodeKnown: public ZeroDimensionalNode
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Contructor
    explicit ZeroDimensionalNodeKnown();

    //! Contructor.
    /*!
     *@param Voltage Source connected to the knwn node.
     */
    ZeroDimensionalNodeKnown ( const zeroDimensionalElementVoltageSourcePtr_Type& theElement );


    //! Destructor
    virtual ~ZeroDimensionalNodeKnown() {}

    //@}


    //! @name Set Methods
    //@{

    //!Set the VoltageSource Element which is connected to the Node
    void setElement ( const zeroDimensionalElementVoltageSourcePtr_Type& element )
    {
        M_element = element;
    }

    void setVoltageByTime (const Real& time )
    {
        M_voltage = M_element->voltageByTime (time);
        M_element->setVoltageByTime (time);
    }

    void setDeltaVoltageByTime (const Real& time )
    {
        M_deltaVoltage = M_element->deltaVoltageByTime (time);
        M_element->setDeltaVoltageByTime (time);
    }

    const Real& voltage() const
    {
        return M_element->voltage();
    }

    Real  voltageByTime (Real& time) const
    {
        return M_element->voltageByTime (time);
    }

    Real  deltaVoltageByTime (Real& time) const
    {
        return M_element->deltaVoltageByTime (time);
    }

    //@}

protected:

    zeroDimensionalElementVoltageSourcePtr_Type M_element;
};




// TODO Move type definitions inside classes
typedef boost::shared_ptr<ZeroDimensionalNode>              zeroDimensionalNodePtr_Type;
typedef std::vector<zeroDimensionalNodePtr_Type>            vecZeroDimensionalNodePtr_Type;
typedef boost::shared_ptr< vecZeroDimensionalNodePtr_Type > ptrVecZeroDimensionalNodePtr_Type;
typedef vecZeroDimensionalNodePtr_Type::iterator            iterZeroDimensionalNode_Type;

typedef boost::shared_ptr<ZeroDimensionalNodeUnknown>             zeroDimensionalNodeUnknownPtr_Type;
typedef std::vector< zeroDimensionalNodeUnknownPtr_Type >         vecZeroDimensionalNodeUnknownPtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalNodeUnknownPtr_Type>  ptrVecZeroDimensionalNodeUnknownPtr_Type;
typedef vecZeroDimensionalNodeUnknownPtr_Type::iterator           iterZeroDimensionalNodeUnknown_Type;

typedef boost::shared_ptr<ZeroDimensionalNodeKnown>               zeroDimensionalNodeKnownPtr_Type;
typedef std::vector< zeroDimensionalNodeKnownPtr_Type >           vecZeroDimensionalNodeKnownPtr_Type;
typedef boost::shared_ptr< vecZeroDimensionalNodeKnownPtr_Type >  ptrVecZeroDimensionalNodeKnownPtr_Type;
typedef vecZeroDimensionalNodeKnownPtr_Type::iterator             iterZeroDimensionalNodeKnown_Type;

typedef std::map <Int, zeroDimensionalElementVoltageSourcePtr_Type>  mapVoltageSource_Type;
typedef boost::shared_ptr < mapVoltageSource_Type>                   mapVoltageSourcePtr_Type;


//! ZeroDimensionalElementS - Container of elements
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalElementS
{
public:

    //! constructor
    explicit ZeroDimensionalElementS();

    //! Destructor
    virtual ~ZeroDimensionalElementS() {}

    void showMe ( const Int& flag = 0 );

    //! add element to the list.
    void setelementList (const zeroDimensionalElementPtr_Type& theElement )
    {
        M_elementList->push_back (theElement);
    }

    //! get element.
    /*!
     *@param element index
     *@return element
     */
    const zeroDimensionalElementPtr_Type& elementListAt ( const Int& index ) const
    {
        return M_elementList->at (index);
    }

    const ptrVecZeroDimensionalElementPtr_Type& elementList() const
    {
        return M_elementList;
    }

    const ptrVecZeroDimensionalElementPassiveResistorPtr_Type& resistorList() const
    {
        return M_resistorList;
    }

    const ptrVecZeroDimensionalElementPassiveCapacitorPtr_Type& capacitorList() const
    {
        return M_capacitorList;
    }

    const ptrVecZeroDimensionalElementPassiveInductorPtr_Type& inductorList() const
    {
        return M_inductorList;
    }

    const ptrVecZeroDimensionalElementPassiveDiodePtr_Type& diodeList() const
    {
        return M_diodeList;
    }

    const ptrVecZeroDimensionalElementVoltageSourcePtr_Type& voltageSourceList() const
    {
        return M_voltageSourceList;
    }

    const ptrVecZeroDimensionalElementCurrentSourcePtr_Type& currentSourceList() const
    {
        return M_currentSourceList;
    }

    //! total number of elements including sources.
    Int elementCounter() const
    {
        return M_elementList->size();    //TODO Why when I use CONST I get a warning??
    }

    Int resistorCounter() const
    {
        return M_resistorList->size();
    }

    Int capacitorCounter() const
    {
        return M_capacitorList->size();
    }

    Int inductorCounter() const
    {
        return M_inductorList->size();
    }

    Int diodeCounter() const
    {
        return M_diodeList->size();
    }

    Int voltageSourceCounter() const
    {
        return M_voltageSourceList->size();
    }

    Int currentSourceCounter() const
    {
        return M_currentSourceList->size();
    }

    //! add resistor to the resistor list.
    void  setResistorList (const zeroDimensionalElementPassiveResistorPtr_Type& resistorPtr )
    {
        M_resistorList->push_back (resistorPtr);
    }

    //! add capacitor to the capacitor list.
    void  setCapacitorList (const zeroDimensionalElementPassiveCapacitorPtr_Type& capacitorPtr )
    {
        M_capacitorList->push_back (capacitorPtr);
    }

    //! add inductor to the inductor list.
    void  setInductorList (const zeroDimensionalElementPassiveInductorPtr_Type& inductorPtr )
    {
        M_inductorList->push_back (inductorPtr);
    }

    //! add diode to the diode list.
    void  setDiodeList (const zeroDimensionalElementPassiveDiodePtr_Type& diodePtr )
    {
        M_diodeList->push_back (diodePtr);
    }

    //! add currentSource to the current Source list.
    void  setCurrentSourceList (const zeroDimensionalElementCurrentSourcePtr_Type& currentSourcePtr )
    {
        M_currentSourceList->push_back (currentSourcePtr);
    }

    //! add voltgeSource to the voltage source list.
    void  setVoltageSourceList (const zeroDimensionalElementVoltageSourcePtr_Type& voltageSourcePtr )
    {
        M_voltageSourceList->push_back (voltageSourcePtr);
    }

    //! add object to the map from voltage source index to the voltage source object.
    void setVoltageSourceMap (const Int& id, const zeroDimensionalElementVoltageSourcePtr_Type& voltageSource)
    {
        (*M_voltageSourceMap) [id] = voltageSource;
    }

    const zeroDimensionalElementVoltageSourcePtr_Type voltageSourceMap (Int& id) const
    {
        return (*M_voltageSourceMap) [id] ;
    }

protected:

    //!List of Elements Ptr
    ptrVecZeroDimensionalElementPtr_Type                    M_elementList;
    ptrVecZeroDimensionalElementPassiveResistorPtr_Type     M_resistorList;
    ptrVecZeroDimensionalElementPassiveCapacitorPtr_Type    M_capacitorList;
    ptrVecZeroDimensionalElementPassiveInductorPtr_Type     M_inductorList;
    ptrVecZeroDimensionalElementPassiveDiodePtr_Type        M_diodeList;
    ptrVecZeroDimensionalElementCurrentSourcePtr_Type       M_currentSourceList;
    ptrVecZeroDimensionalElementVoltageSourcePtr_Type       M_voltageSourceList;
    mapVoltageSourcePtr_Type                                M_voltageSourceMap;
};

// TODO Move type definitions inside classes
typedef std::map <Int, zeroDimensionalNodeUnknownPtr_Type>                          mapNodeUnknown_Type;
typedef std::map <Int, zeroDimensionalNodeKnownPtr_Type>                            mapNodeKnown_Type;
typedef boost::shared_ptr < mapNodeKnown_Type>                                      mapNodeKnownPtr_Type;
typedef boost::shared_ptr < mapNodeUnknown_Type  >                                  mapNodeUnknownPtr_Type;



//! ZeroDimensionalNodeS - Container of nodes
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalNodeS
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalNodeS();

    //! Destructor
    virtual ~ZeroDimensionalNodeS() {}

    virtual void showMe ( const Int& flag = 0 );

    const zeroDimensionalNodePtr_Type& nodeListAt ( const Int& index ) const
    {
        return M_nodeList->at (index);
    }

    const zeroDimensionalNodeUnknownPtr_Type unknownNodeListAt ( const Int& Index) const
    {
        return M_unknownNodeList->at (Index);
    }

    const zeroDimensionalNodeKnownPtr_Type knownNodeListAt ( const Int& Index) const
    {
        return M_knownNodeList->at (Index);
    }

    const ptrVecZeroDimensionalNodePtr_Type nodeList() const
    {
        return M_nodeList;
    }

    const ptrVecZeroDimensionalNodeUnknownPtr_Type& unknownNodeList() const
    {
        return M_unknownNodeList;
    }

    const ptrVecZeroDimensionalNodeKnownPtr_Type knownNodeList() const
    {
        return M_knownNodeList;
    }

    const zeroDimensionalNodeKnownPtr_Type knownNodeMapAt (Int& id) const
    {
        return (*M_knownNodeMap) [id];
    }

    const zeroDimensionalNodeUnknownPtr_Type unknownNodeMapAt (Int& id) const
    {
        return (*M_unknownNodeMap) [id];
    }

    //! add node to the list
    void setnodeList (const zeroDimensionalNodePtr_Type& theNode )
    {
        M_nodeList->push_back (theNode);
    }

    //! add unknownNode to the unknwnNode List
    void setunknownNodeList ( const zeroDimensionalNodeUnknownPtr_Type& unknownNode)
    {
        M_unknownNodeList->push_back (unknownNode);
    }

    //! add knownNode to the knwnNode List
    void setknownNodeList ( const zeroDimensionalNodeKnownPtr_Type& knownNode)
    {
        M_knownNodeList->push_back (knownNode);
    }

    //! add knownNode to the map. A map from the index (id) to the object.
    void setknownNodeMap ( const Int& id, const zeroDimensionalNodeKnownPtr_Type& knownNode)
    {
        (*M_knownNodeMap) [id] = knownNode;
    }

    //! add unknownNode to the map. A map from the index (id) to the object.
    void setunknownNodeMap ( const Int& id, const zeroDimensionalNodeUnknownPtr_Type& unknownNode)
    {
        (*M_unknownNodeMap) [id] = unknownNode;
    }

    Int unknownNodeCounter() const
    {
        return M_unknownNodeList->size();
    }

    Int knownNodeCounter() const
    {
        return M_knownNodeList->size();
    }

    Int nodeCounter() const
    {
        return M_nodeList->size();
    }

protected:

    //List of Nodes
    ptrVecZeroDimensionalNodePtr_Type               M_nodeList;
    ptrVecZeroDimensionalNodeUnknownPtr_Type        M_unknownNodeList;
    ptrVecZeroDimensionalNodeKnownPtr_Type          M_knownNodeList;
    mapNodeKnownPtr_Type                            M_knownNodeMap;
    mapNodeUnknownPtr_Type                          M_unknownNodeMap;
};



//! ZeroDimensionalCircuitData - Container of circuit data
/*!
 *  @authors Mahmoud Jafargholi
 */
class ZeroDimensionalCircuitData
{
public:

    //! @name Constructors & Destructor
    //@{

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalCircuitData();

    //! Destructor
    virtual ~ZeroDimensionalCircuitData() {}

    //@}

    void showMe ( const Int& flag = 0 );

    //! create the circuit.
    /*!
     * @param circuit file
     * @param BC handler
     */
    void buildCircuit (const char* fileName, bcPtr_Type bc);

    //! get element container object.
    const zeroDimensionalElementSPtr_Type Elements() const
    {
        return M_Elements;
    }

    //! get node container object.
    const zeroDimensionalNodeSPtr_Type Nodes() const
    {
        return M_Nodes;
    }

    //! (shallow) update the circuit data from the solution.
    /*!
     * This method is invoked every iteration before calling the updateABC method.
     * This method updates the circuit data which is dependent on time or solution vector. For example
     * source elements are function of time and diode R_{eff} are function of voltage difference.
     */
    void updateCircuitDataFromY (const Real& t, const Epetra_Vector* y, const Epetra_Vector* yp);

    //! create matrix A,B and C.
    /*!
     * before calling this method, updateCircuitDataFromY method should be invoked.
     */
    void updateABC (matrix_Type& A, matrix_Type& B, vector_Type& C);


    //! (deep) update the circuit data from solution.
    /*!
     * This methed is invoked after Rythoms step is finished. This method computes currents.
     */
    void extractSolutionFromY (const Real& t, const Epetra_Vector& y, const Epetra_Vector& yp);

protected:

    // set BCs to source elements
    void fixBC (bcPtr_Type bc);

    void createElementResistor (Int ID, Int node1, Int node2, Real parameter);
    void createElementCapacitor (Int ID, Int node1, Int node2, Real parameter );
    void createElementInductor (Int ID, Int node1, Int node2, Real parameter);
    void createElementDiode (Int ID, Int node1, Int node2, Real forwardBias, Real alpha, Real beta);
    Int createElementVoltageSource (Int node1);
    void createElementCurrentSource (Int node1);
    void createUnknownNode (const Int& id);
    void createKnownNode (const Int& id);
    void createKnownNode (const Int& id, const zeroDimensionalElementVoltageSourcePtr_Type& theElement);

    zeroDimensionalElementSPtr_Type  M_Elements;
    zeroDimensionalNodeSPtr_Type     M_Nodes;
    bcPtr_Type                       M_bc;
};






typedef boost::shared_ptr< ZeroDimensionalCircuitData > zeroDimensionalCircuitDataPtr_Type;

//! OutPutFormat - Write to output
/*!
 *  @authors Mahmoud Jafargholi
 */
class OutPutFormat
{
public:

    //! Constructor
    explicit OutPutFormat ( std::string  width, std::string  precision, std::string  whiteSpace, Int bufferSize);

    //! Destructor
    virtual ~OutPutFormat();

    enum EndLine { newLine, space, nothing };

    void writeDataFormat (const Real& number, std::ofstream& stream, const EndLine& flag);
    void writeDataFormat (const Int& number, std::ofstream& stream, const EndLine& flag);
    void writeDataFormat (const string& text, std::ofstream& stream, const EndLine& flag);
    void writeNewLine (std::ofstream& stream);

private:

    std::string  M_width;
    UInt         M_integerWidth;
    std::string  M_precision ;
    std::string  M_whiteSpace ;
    std::string  M_formatDouble;
    std::string  M_formatInteger;
    std::string  M_formatString;
    char*         M_buffer;
};

} // LifeV namespace

#endif //ZeroDimensionalCircuitData_H
