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
 *  @date 26-09-2011
 *  @author Mahmoud Jafargholi <mahmoud.jafargholi@epfl.ch>
 *
 *  @contributors Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

//! ZeroDimensionalCircuitData - Class which read and holds all the data for the Zero Dimensional Model Solver.
/*!
 *  @author Mahmoud Jafargholi
 */

#ifndef ZeroDimensionalCircuitData_H
#define ZeroDimensionalCircuitData_H 1

#define BC_CONSTANT 1000

// LIFEV
#include <life/lifearray/MatrixEpetra.hpp>

// MATHCARD
#include <lifemc/lifesolver/BCInterface0D.hpp>
#include <lifemc/lifesolver/MultiscaleModel.hpp>
#include <lifemc/lifesolver/ZeroDimensionalDefinitions.hpp>
#include <lifemc/lifefem/ZeroDimensionalBC.hpp>

#define  ZERO_DIMENTIONAL_DEFINED_ELEMENTS          6
#define  ZERO_DIMENTIONAL_DEFINED_NODES             2




namespace LifeV
{

/**
* ZeroDimentional element type. 
*/
enum ZeroDimentionalElementType
    {
        resistor,
        capacitor,
        inductor,
        diode,
        voltageSource,
        currentSource
    };

/**
* ZeroDimentional node type. 
*/
enum ZeroDimentionalNodeType
    {
        knownNode,
        unknownNode
    };
  

//! A container class for all node objects
class ZeroDimensionalNodeS;

  //! A container class for all element obkects
class ZeroDimensionalElementS;

    typedef boost::shared_ptr< ZeroDimensionalElementS>                         zeroDimensionalElementSPtr_Type;
    typedef boost::shared_ptr< ZeroDimensionalNodeS>                            zeroDimensionalNodeSPtr_Type;
    typedef std::vector<Int>                                                vecInt_Type;
    typedef vecInt_Type::iterator                                           iterVecInt_Type;
    typedef ZeroDimensionalBCHandler                                        bc_Type;
    typedef boost::shared_ptr< bc_Type >                                    bcPtr_Type;
    typedef BCInterface0D< bc_Type, Multiscale::MultiscaleData >            bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >                           bcInterfacePtr_Type;
    typedef MatrixEpetra<double>                                            matrix_Type;
    typedef VectorEpetra                                                    vector_Type;
    typedef Epetra_Vector                                                   vectorEpetra_Type;
    typedef boost::shared_ptr< matrix_Type >                                matrixPtr_Type;
    typedef boost::shared_ptr< vector_Type >                                vectorPtr_Type;
    typedef boost::shared_ptr<vectorEpetra_Type >                           vectorEpetraPtr_Type;


//-----------------------------------------------------------------------
//! ZeroDimensionalElement - The base element class .
/*!
 *  This class is the base class for all elements. 
 */
class ZeroDimensionalElement
{
public:
    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalElement();
    
    //! Destructor
    virtual ~ZeroDimensionalElement() {}

    //@}
    
    //! Display some information.
    virtual void showMe(__attribute__((unused)) Int flag=0) ;

    //! Connect elements to the nodes.
    /*!
     * After all emenets and nodes are created, each element will call this
     * method to connect itse;f to the nodes.
     */  
    virtual void connectElement(zeroDimensionalNodeSPtr_Type & Nodes)=0;
   
    //! Contribution of the element of matrix \bf{A} and \bf{B} and vector \bf{C}.
    /*!
     * After updating the BCs ( or Terminal nodes ) this each element will invoke 
     * this method to compute it's contribution on matrices.
     */
    virtual void buildABC( __attribute__((unused)) matrix_Type& A,__attribute__((unused)) matrix_Type& B, __attribute__((unused)) vector_Type& C, __attribute__((unused)) const zeroDimensionalNodeSPtr_Type& Nodes){};

  //! Compute outputs (currents and voltages) from the solution vector after each succesful iteration.
  /*!
   * After each time step, when Rythmos solver is succesfully finishes, this method will compute 
   * finial outputs ( for exmple currents) from the finial solution vector.
   */
    virtual void deepUpdate(__attribute__((unused)) const ZeroDimensionalNodeS& Nodes){}

  //! This method specifies the convention of current direction in an element.
  /*!
   * @param A node index connected to the element.
   * @return +1 if the current convention is toward the iniput node and -1 otherwise.
   */ 
    virtual Real direction(const Int & nodeId) const=0;
    
    const std::string  enum2string(const ZeroDimentionalElementType & type);

    void setid              (const Int                          & id                ) { M_id          = id                  ; }

    void setcurrent         (const Real                         & current           ) { M_current     = current           ; }

  //! set derivative of current respect to time.
    void setdeltaCurrent    (const Real                         & deltaCurrent      ) { M_deltaCurrent= deltaCurrent           ; }

    const Int                       & id            ()        const { return M_id             ; }

    const ZeroDimentionalElementType& type          ()        const { return M_type           ; }

    const Real                      & current       ()        const { return M_current      ; }

  //! get derivative of current respect to time.
    const Real                      & deltaCurrent  ()        const { return M_deltaCurrent      ; }

protected:

    Int                             M_id                    ;
    ZeroDimentionalElementType      M_type                  ; //= 'Resistor';%'Capacitor' ,'Inductor','Voltage Source','Current Source' 'Diode'
    Real                            M_current               ; 
    Real                            M_deltaCurrent          ; 
};

typedef boost::shared_ptr<ZeroDimensionalElement>                              zeroDimensionalElementPtr_Type;
typedef std::vector<zeroDimensionalElementPtr_Type>                            vecZeroDimensionalElementPtr_Type;
typedef boost::shared_ptr<vecZeroDimensionalElementPtr_Type>                   ptrVecZeroDimensionalElementPtr_Type;
typedef vecZeroDimensionalElementPtr_Type::iterator                            iterZeroDimensionalElement_Type;
//-----------------------------------------------------------------------
//! ZeroDimensionalElementPassive. A class for passive elements.
class ZeroDimensionalElementPassive: public ZeroDimensionalElement
{
public:

  //! Constructor
    explicit ZeroDimensionalElementPassive();

    //! Destructor
    virtual ~ZeroDimensionalElementPassive() {}

  //! Show some information.
    void showMe(Int flag=0);

  //! Impleaments the abstarct class for passive elements.
    void connectElement              (zeroDimensionalNodeSPtr_Type & Nodes);

  //! set parameter (1/R, 1/L, C, 1/R_{eff})
    void setparameter              (const Real                          & parameter                ) { M_parameter          = parameter                  ; }
  //! add the node to the list.
  /*!
   * @param node index.
   */
    void setnodeIndex              (const Int & index1) { M_nodeIndex.push_back(index1); }

  //! get the parameter (1/R, 1/L, C, 1/R_{eff})
    const Real                      & parameter     ()      const { return M_parameter      ; }

  //! get the node index connected to the node.
  /*!
   * @param \it{i}th node connected to the elelemt.
   * @return  Index of \it{i}th node connected to the element.
   */
    const Int                       & nodeIndex     (const Int & position)      const { return M_nodeIndex.at(position)      ; }

    Real direction(const Int & nodeId)const;

protected:

    Real            M_parameter     ;
    //parameter= 'Resistor';%'Capacitor' ,'Inductor','Diode'
    //parameter=  1/R      ;%C            ,1/L      , 1/R_{eff}
    vecInt_Type M_nodeIndex      ; //Index of connected nodes
};
//-----------------------------------------------------------------------
//! ZeroDimentionalElement - Resistor.
class ZeroDimensionalElementPassiveResistor: public ZeroDimensionalElementPassive
{
public:

  //! Contructor.
    explicit ZeroDimensionalElementPassiveResistor();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveResistor() {}

  //    void insertElement(zeroDimensionalElementSPtr_Type & Elements);

    void showMe(Int flag=0);

    void buildABC(matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes);

    void deepUpdate(const ZeroDimensionalNodeS& Nodes);
};

//-----------------------------------------------------------------------
//! ZerodimentionalElement - Diode.
class ZeroDimensionalElementPassiveDiode: public ZeroDimensionalElementPassiveResistor
{
public:
   
    //! Constructor
    explicit ZeroDimensionalElementPassiveDiode();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveDiode() {}

    void showMe(Int flag=0);

    //!current = beta * exp(alpha * (voltage - forwardBias )) - (beta * exp(alpha * ( - forwardBias )))
    void setalpha                  (const Real                          & alpha               ) { M_alpha         = alpha                 ; }

    void setbeta                   (const Real                          & beta                ) { M_beta          = beta                  ; }

    void setforwardBias            (const Real                          & forwardBias         ) { M_forwardBias   = forwardBias           ; }

    const Real                      & alpha       ()      const { return M_alpha      ; }

    const Real                      & beta        ()      const { return M_beta       ; }

    const Real                      & forwardBias ()      const { return M_forwardBias; }

    void deepUpdate(const ZeroDimensionalNodeS& Nodes);

    void buildABC(matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes);

protected:

    Real            M_alpha       ; //current = beta * exp(alpha * (voltage - forwardBias )) - (beta * exp(alpha * ( - forwardBias )))
    Real            M_beta        ;
    Real            M_forwardBias ;
    
  //! calculate the effective resistance.
  /*!
   * @param voltage difference
   * @return effective ressitance
   */
    void calculateEffectiveResistance(const double& voltage);
};

//-----------------------------------------------------------------------
//! Zerodimentional Element - Capacitor.
class ZeroDimensionalElementPassiveCapacitor: public ZeroDimensionalElementPassive
{
public:

  //! Constructor
    explicit ZeroDimensionalElementPassiveCapacitor();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveCapacitor() {}

    void showMe(Int flag=0);

    void deepUpdate(const ZeroDimensionalNodeS& Nodes);

    void buildABC(matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes);

protected:
};


//-----------------------------------------------------------------------
//! ZeroDimentional Element - Inductor.
class ZeroDimensionalElementPassiveInductor: public ZeroDimensionalElementPassive
{
public:

  //! Constructor
    explicit ZeroDimensionalElementPassiveInductor();

    //! Destructor
    virtual ~ZeroDimensionalElementPassiveInductor() {}

  //! Set the variable index and equation row index for Inductor.
  /*! 
   *Current in Inductor is an unknown. 
   */
    virtual void     assignVariableIndex(const Int & index);

    void showMe(Int flag=0);

  //! get equation row for in matrix A,B and C.
    Int   equationRow          ()const                          { return M_equationRow                   ; }

  //! get variable index in solution vector  x  and \dot{x}
    Int   variableIndex        ()const                          { return M_variableIndex                 ; }

    void buildABC(matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes);

protected:
    Int            M_equationRow     ;
    Int            M_variableIndex   ;
};

//-----------------------------------------------------------------------
//! Base class for source elements.
class ZeroDimensionalElementSource: public ZeroDimensionalElement
{

public:

  //!Constructor
    explicit ZeroDimensionalElementSource();

    //! Destructor
    virtual ~ZeroDimensionalElementSource() {}

    void showMe(Int flag=0);

    void setnodeIndex      (const Real & index           ) {        M_nodeIndex = index ; }

    Int  nodeIndex()const                          { return M_nodeIndex         ; }

  //! set BC handler.
    void setbc             (const bcInterfacePtr_Type& bc) {M_bc = bc;}

    Real direction(__attribute__((unused)) const Int & nodeId)const {return -1.0;}

protected:

    Int                        M_nodeIndex     ; //Index of connected node
    bcInterfacePtr_Type        M_bc;
};

//-----------------------------------------------------------------------
//! ZerodimentionalElement - Voltage Source.
class ZeroDimensionalElementVoltageSource: public ZeroDimensionalElementSource
{

public:

  //! Constructor
    explicit ZeroDimensionalElementVoltageSource();

    //! Destructor
    virtual ~ZeroDimensionalElementVoltageSource() {}

    void connectElement    (zeroDimensionalNodeSPtr_Type & Nodes);

  //! Update voltage source by time. 
    void setvoltageByTime                 (const Real& time              ) { M_voltage                   = M_bc->handler()->bc( M_nodeIndex ).evaluate(time);                   ; }

  //! Update \frac{\partial voltage}{\partial t} by time. 
    void setdeltaVoltageByTime            (const Real& time              ) { M_deltaVoltage              = M_bc->handler()->bc( M_nodeIndex + BC_CONSTANT ).evaluate(time);              ; }

    Real voltage()                 const { return M_voltage                  ; }

    Real deltaVoltage()            const { return M_deltaVoltage               ; }

    Real voltageByTime(const Real & time)        const { return M_bc->handler()->bc( M_nodeIndex ).evaluate(time); }

    Real deltaVoltageByTime(const Real & time)   const { return M_bc->handler()->bc( M_nodeIndex + BC_CONSTANT ).evaluate(time); }

  //! calculate current passing outward in voltage source.
  /*!
   *  This method can be called after all elements invoked deepUpdate method.
   */
    void calculateCurrent(const ZeroDimensionalNodeS& Nodes,const ZeroDimensionalElementS& Elements);
protected:

    Real    M_voltage                     ; //voltage at time t_{n}
    Real    M_deltaVoltage                ; //\frac{\mathrm{d \text{ voltage}} }{\mathrm{d} t}
};


//-----------------------------------------------------------------------
//! ZerodimentionalElement - Current Source.
class ZeroDimensionalElementCurrentSource: public ZeroDimensionalElementSource
{

public:

  //! Constructor.
    explicit ZeroDimensionalElementCurrentSource();

    //! Destructor
    virtual ~ZeroDimensionalElementCurrentSource() {}

    void connectElement    (zeroDimensionalNodeSPtr_Type & Nodes);

    void setcurrentByTime               (const Real& time                  ) { M_current                 = M_bc->handler()->bc( M_nodeIndex ).evaluate(time); }

    Real currentByTime           (const Real& time                  ) const { return M_bc->handler()->bc( M_nodeIndex ).evaluate(time);}

    Real current()                 const { return M_current  ; }

    void buildABC(matrix_Type& A,matrix_Type& B,vector_Type& C, const zeroDimensionalNodeSPtr_Type& Nodes);

protected:
};

//-----------------------------------------------------------------------
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
//----------------------------------------------------------------------------

//! ZeroDimensionalNode - The base node class .
/*!
 *  This class is the base class for all nodes. 
 */
class ZeroDimensionalNode
{
public:

  //! Constructor
    explicit ZeroDimensionalNode();

    //! Destructor
    virtual ~ZeroDimensionalNode() {}

    virtual void showMe(Int flag=0);

    const std::string  enum2string(const ZeroDimentionalNodeType & type) const;

    void            setid              (const Int                          & id                ) { M_id          = id                  ;}

  //! add an element index to the elelemt list. 
    void            setelementListIndex(const Int                          & index             ) { M_elementListIndex.push_back(index) ;}

  //! add an node index which is connected by an element in element list.
  /*!
   * Each elelemnt in element list, coonects this node to another node ( except source elementt). nodeList is a container for conecting nodes.
   * If the element connected to this node has only one terminal ( like voltage source and current source), the connecting index would be -1.
   */
    void            setnodeListIndex   (const Int                          & index             ) { M_nodeListIndex.push_back(index)    ;}

    virtual void    setvoltage         (const Real                         & voltage           ) { M_voltage          = voltage        ;}

    virtual void    setdeltaVoltage    (const Real                         & deltaVoltage      ) { M_deltaVoltage     = deltaVoltage   ;}

  //! Calculate current balance at node.
  /*!
   * After updating current in all elements, we can verify the balance of current flow at each node.
   */
    void            calculateCurrentBalance(const ZeroDimensionalElementS& Elements);
    
    const Int                      & id            ()      const { return M_id             ; }

    const ZeroDimentionalNodeType  & type          ()      const { return M_type           ; }

    const Int                      & elementListIndexAt  (const Int & position)      const { return M_elementListIndex.at(position)  ; }

    const vecInt_Type              & elementListIndex   ()      const { return M_elementListIndex  ; }

    const Int                      & nodeListIndexAt    (const Int & position)      const { return M_nodeListIndex.at(position)     ; }

    virtual  const Real              voltage            ()  const { return M_voltage            ; }

    virtual  Real                    deltaVoltage  ()      const { return M_deltaVoltage    ;}

    const Real                     & currentBalance()      const{ return M_currentBalance  ; }

protected:

    Int                             M_id                    ;
    ZeroDimentionalNodeType         M_type                  ; //= 'Known';%'Unknown'
    vecInt_Type                     M_elementListIndex      ; // List of id(s) of connected Elements to this Node
    vecInt_Type                     M_nodeListIndex         ; // List of id(s) of connected Nodes to this Node
    Real                            M_currentBalance        ; //sum of currents over all branches
    Real                            M_voltage               ;
    Real                            M_deltaVoltage          ;
};

//-----------------------------------------------------------------------
//! ZeroDimensionalNodeUnknown.
/*!
 *  This class defines the unknown node class. 
 */
class ZeroDimensionalNodeUnknown: public ZeroDimensionalNode
{
public:

  //! Constructor
    explicit ZeroDimensionalNodeUnknown();

    //! Destructor
    virtual ~ZeroDimensionalNodeUnknown() {}

  //! assign the index of the unknown voltage.
    void assignVariableIndex(const Int & index) ;

    void showMe(Int flag=0);

  //    void setvariableIndex          (const Int       & variableIndex             ) { M_variableIndex         = variableIndex              ;}
  // void setequationRow            (const Int       & equationRow               ) { M_equationRow           = equationRow                ;}

    const Int                      & variableIndex      ()    const { return M_variableIndex      ; }

    const Int                      & equationRow        ()    const { return M_equationRow        ; }

protected:
    Int                 M_variableIndex                ; // Index of the variable
    Int                 M_equationRow                  ; // #Row(equation) in the Matrix
};

//-----------------------------------------------------------------------
//! ZeroDimensionalNodeknown.
/*!
 *  This class defines the known node class.
 *  A Voltage Source element is connected to this class.
 */
class ZeroDimensionalNodeKnown: public ZeroDimensionalNode
{
public:

  //! Contructor
    explicit ZeroDimensionalNodeKnown();

  //! Contructor.
  /*!
   *@param Voltage Source connected to the knwn node.
   */
    ZeroDimensionalNodeKnown( const zeroDimensionalElementVoltageSourcePtr_Type & theElement );


    //! Destructor
    virtual ~ZeroDimensionalNodeKnown() {}

//!Set the VoltageSource Element which is connected to the Node
    void setelement                  (const zeroDimensionalElementVoltageSourcePtr_Type &element){M_element=element                     ;}

    void    setvoltageByTime         (const Real                         & time      ) { M_voltage          = M_element->voltageByTime(time)        ;M_element->setvoltageByTime(time);}

    void    setdeltaVoltageByTime    (const Real                         & time      ) { M_deltaVoltage     = M_element->deltaVoltageByTime(time)   ;M_element->setdeltaVoltageByTime(time);}

     const Real  voltage()     const {return M_element->voltage()          ;}

     Real  voltageByTime(Real& time)                   const { return M_element->voltageByTime(time)                  ; }

     Real  deltaVoltageByTime(Real& time)              const { return M_element->deltaVoltageByTime(time)             ; }

protected:

    zeroDimensionalElementVoltageSourcePtr_Type M_element;
};

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

typedef std::map <int, zeroDimensionalElementVoltageSourcePtr_Type>                 mapVoltageSource_Type;
typedef boost::shared_ptr < mapVoltageSource_Type>                                  mapVoltageSourcePtr_Type;

//! container class for elements
class ZeroDimensionalElementS
{
public:

  //! constructor
    explicit ZeroDimensionalElementS();

    //! Destructor
    virtual ~ZeroDimensionalElementS() {}

    void showMe(Int flag=0);

  //! add element to the list.
    void setelementList(const zeroDimensionalElementPtr_Type& theElement   ) { M_elementList->push_back(theElement);}

  //! get element.
  /*!
   *@param element index
   *@return element
   */
    const zeroDimensionalElementPtr_Type                &elementListAt(const Int & index)      const { return M_elementList->at(index); }

    const ptrVecZeroDimensionalElementPtr_Type          elementList()                         const { return M_elementList; }

    const ptrVecZeroDimensionalElementPassiveResistorPtr_Type      resistorList ()      const { return M_resistorList       ; }

    const ptrVecZeroDimensionalElementPassiveCapacitorPtr_Type     capacitorList ()     const { return M_capacitorList      ; }

    const ptrVecZeroDimensionalElementPassiveInductorPtr_Type      inductorList ()      const { return M_inductorList       ; }

    const ptrVecZeroDimensionalElementPassiveDiodePtr_Type         diodeList ()         const { return M_diodeList          ; }

    const ptrVecZeroDimensionalElementVoltageSourcePtr_Type        voltageSourceList ()const { return M_voltageSourceList  ; }

    const ptrVecZeroDimensionalElementCurrentSourcePtr_Type        currentSourceList ()const { return M_currentSourceList  ; }

  //! total number of elements including sources.
          Int elementCounter        ()                          const { return M_elementList->size()            ; }//TODO Why when I use CONST I get a warning??

          Int resistorCounter       (                    )      const { return M_resistorList->size()        ; }

          Int capacitorCounter      (                    )      const { return M_capacitorList->size()       ; }

          Int inductorCounter       (                    )      const { return M_inductorList->size()        ; }

          Int diodeCounter          (                    )      const { return M_diodeList->size()           ; }

          Int voltageSourceCounter  (                    )      const { return M_voltageSourceList->size()   ; }

          Int currentSourceCounter  (                    )      const { return M_currentSourceList->size()   ; }

  //! add resistor to the resistor list. 
    void  setresistorList        (const zeroDimensionalElementPassiveResistorPtr_Type & resistorPtr      ) { M_resistorList->push_back(resistorPtr); }

  //! add capacitor to the capacitor list. 
    void  setcapacitorList       (const zeroDimensionalElementPassiveCapacitorPtr_Type & capacitorPtr    ) { M_capacitorList->push_back(capacitorPtr); }

  //! add inductor to the inductor list. 
    void  setinductorList        (const zeroDimensionalElementPassiveInductorPtr_Type& inductorPtr       ) { M_inductorList->push_back(inductorPtr); }

  //! add diode to the diode list. 
    void  setdiodeList           (const zeroDimensionalElementPassiveDiodePtr_Type    & diodePtr         ) { M_diodeList->push_back(diodePtr); }

  //! add currentSource to the current Source list. 
    void  setcurrentSourceList   (const zeroDimensionalElementCurrentSourcePtr_Type   & currentSourcePtr ) { M_currentSourceList->push_back(currentSourcePtr); }

  //! add voltgeSource to the voltage source list. 
    void  setvoltageSourceList   (const zeroDimensionalElementVoltageSourcePtr_Type   & voltageSourcePtr ) { M_voltageSourceList->push_back(voltageSourcePtr); }

  //! add object to the map from voltage source index to the voltage source object. 
    void setvoltageSourceMap           (const Int & id, const zeroDimensionalElementVoltageSourcePtr_Type & voltageSource) {(*M_voltageSourceMap)[id]= voltageSource;}

    const zeroDimensionalElementVoltageSourcePtr_Type        voltageSourceMap(Int& id) const {return (*M_voltageSourceMap)[id] ;}

protected:

//!List of Elements Ptr
    ptrVecZeroDimensionalElementPtr_Type                    M_elementList             ; 
    ptrVecZeroDimensionalElementPassiveResistorPtr_Type     M_resistorList            ;
    ptrVecZeroDimensionalElementPassiveCapacitorPtr_Type    M_capacitorList           ;
    ptrVecZeroDimensionalElementPassiveInductorPtr_Type     M_inductorList            ;
    ptrVecZeroDimensionalElementPassiveDiodePtr_Type        M_diodeList               ;
    ptrVecZeroDimensionalElementCurrentSourcePtr_Type       M_currentSourceList       ;
    ptrVecZeroDimensionalElementVoltageSourcePtr_Type       M_voltageSourceList       ;
    mapVoltageSourcePtr_Type                                M_voltageSourceMap        ;
};

//-----------------------------------------------------------------------
typedef std::map <int, zeroDimensionalNodeUnknownPtr_Type>                          mapNodeUnknown_Type;
typedef std::map <int, zeroDimensionalNodeKnownPtr_Type>                            mapNodeKnown_Type;
typedef boost::shared_ptr < mapNodeKnown_Type>                                      mapNodeKnownPtr_Type;
typedef boost::shared_ptr < mapNodeUnknown_Type  >                                  mapNodeUnknownPtr_Type;

//! container class for nodes.
class ZeroDimensionalNodeS
{
public:

  //! Constructor
    explicit ZeroDimensionalNodeS();

    //! Destructor
    virtual ~ZeroDimensionalNodeS() {}

    virtual void showMe(Int flag=0);

    const zeroDimensionalNodePtr_Type&            nodeListAt        (const Int & index)      const {return M_nodeList->at(index)       ;}

    const zeroDimensionalNodeUnknownPtr_Type     unknownNodeListAt (const Int & Index)      const {return M_unknownNodeList->at(Index);}

    const zeroDimensionalNodeKnownPtr_Type       knownNodeListAt   (const Int & Index)      const {return M_knownNodeList->at(Index)  ;}

    const ptrVecZeroDimensionalNodePtr_Type         nodeList        ()        const {return M_nodeList           ;}

    const ptrVecZeroDimensionalNodeUnknownPtr_Type&  unknownNodeList ()        const {return M_unknownNodeList    ;}

    const ptrVecZeroDimensionalNodeKnownPtr_Type    knownNodeList   ()        const {return M_knownNodeList      ;}

    const zeroDimensionalNodeKnownPtr_Type          knownNodeMapAt(Int& id)   const {return (*M_knownNodeMap)[id]   ;}

    const zeroDimensionalNodeUnknownPtr_Type        unknownNodeMapAt(Int& id) const {return (*M_unknownNodeMap)[id] ;}

  //! add node to the list
    void setnodeList                 (const zeroDimensionalNodePtr_Type        & theNode    )    {M_nodeList->push_back(theNode);}

  //! add unknownNode to the unknwnNode List
    void setunknownNodeList          (const zeroDimensionalNodeUnknownPtr_Type & unknownNode)    {M_unknownNodeList->push_back(unknownNode);}

  //! add knownNode to the knwnNode List
    void setknownNodeList            (const zeroDimensionalNodeKnownPtr_Type   & knownNode  )    {M_knownNodeList->push_back(knownNode);}

  //! add knownNode to the map. A map from the index (id) to the object.
    void setknownNodeMap             (const Int & id, const zeroDimensionalNodeKnownPtr_Type   & knownNode)   {(*M_knownNodeMap)[id]  = knownNode;}

  //! add unknownNode to the map. A map from the index (id) to the object.
    void setunknownNodeMap           (const Int & id, const zeroDimensionalNodeUnknownPtr_Type & unknownNode) {(*M_unknownNodeMap)[id]= unknownNode;}

    Int unknownNodeCounter    ()                       const {return M_unknownNodeList->size()       ; }

    Int knownNodeCounter      ()                       const {return M_knownNodeList->size()         ; }

    Int nodeCounter           ()                       const {return M_nodeList->size()              ; }

protected:

//List of Nodes
    ptrVecZeroDimensionalNodePtr_Type               M_nodeList            ; 
    ptrVecZeroDimensionalNodeUnknownPtr_Type        M_unknownNodeList     ;
    ptrVecZeroDimensionalNodeKnownPtr_Type          M_knownNodeList       ;
    mapNodeKnownPtr_Type                            M_knownNodeMap;
    mapNodeUnknownPtr_Type                          M_unknownNodeMap;
};
//-----------------------------------------------------------------------

//! Data container for circuit data.
class ZeroDimensionalCircuitData
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ZeroDimensionalCircuitData();

    //! Destructor
    virtual ~ZeroDimensionalCircuitData() {}

    //@}
  
  void showMe(Int flag=0);

  //! create the circuit.
  /*!
   * @param circuit file
   * @param BC handler
   */
  void buildCircuit (const char *fileName, bcInterfacePtr_Type bc);

  //! get element container object. 
    const zeroDimensionalElementSPtr_Type          Elements()  const {return M_Elements;}

  //! get node container object.
    const zeroDimensionalNodeSPtr_Type             Nodes   ()  const {return M_Nodes   ;}

  //! (shallow) update the circuit data from the solution.
  /*!
   * This method is invoked every iteration before calling the updateABC method.
   * This method updates the circuit data which is dependent on time or solution vector. For example
   * source elements are function of time and diode R_{eff} are function of voltage difference.
   */
    void updateCircuitDataFromY(const double& t, const Epetra_Vector* y,const Epetra_Vector* yp);

  //! create matrix A,B and C.
  /*!
   * before calling this method, updateCircuitDataFromY method should be invoked.
   */
    void updateABC(matrix_Type& A,matrix_Type& B,vector_Type& C);


  //! (deep) update the circuit data from solution.
  /*!
   * This methed is invoked after Rythoms step is finished. This method computes currents.  
   */
    void deepUpdateFromY(const double& t, const Epetra_Vector& y,const Epetra_Vector& yp);

protected:

  // set BCs to source elements
    void fixBC(bcInterfacePtr_Type bc);

    void  createElementResistor(Int ID,Int node1,Int node2,Real parameter);

    void  createElementCapacitor(Int ID,Int node1,Int node2,Real parameter );

    void  createElementInductor(Int ID,Int node1,Int node2,Real parameter);

    void  createElementDiode(Int ID, Int node1, Int node2, Real forwardBias, Real alpha, Real beta);

    Int   createElementVoltageSource(Int node1);

    void  createElementCurrentSource(Int node1);

    void  createUnknownNode(const Int& id);

    void  createKnownNode  (const Int& id);

    void  createKnownNode  (const Int& id, const zeroDimensionalElementVoltageSourcePtr_Type & theElement);

    zeroDimensionalElementSPtr_Type                           M_Elements;
    zeroDimensionalNodeSPtr_Type                              M_Nodes;
    bcInterfacePtr_Type                                       M_bc;
};
typedef boost::shared_ptr< ZeroDimensionalCircuitData > zeroDimensionalCircuitDataPtr_Type;

//! write data to a file.
class OutPutFormat {
public:

  //! Constructor
        explicit OutPutFormat(
                    std::string  width,
                    std::string  precision,
                    std::string  whiteSpace,
                    int bufferSize);

  //! Destructor
        virtual ~OutPutFormat();

        enum EndLine {
            newLine,
            space,
            nothing
        };

        void writeDataFormat(const double& number, std::ofstream & stream, const EndLine& flag);

        void writeDataFormat(const int& number,std::ofstream & stream, const EndLine& flag);

        void writeNewLine(std::ofstream & stream);

private:

        std::string  M_width;
        std::string  M_precision ;
        std::string  M_whiteSpace ;
        std::string  M_formatDouble;
        std::string  M_formatInteger;
        char         *M_buffer;
};

} // LifeV namespace

#endif //ZeroDimensionalCircuitData_H
