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
    @brief File contains classes to holds the FE vectors used for prescribing boundary conditions

	@author Miguel Fernandez <miguel.fernandez@inria.fr>
    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @author Vincent Martin <vincent.martin@inria.fr>
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @maintainer Mauro Perego <perego.mauro@gmail.com>

    */

#ifndef BCVECTOR_H
#define BCVECTOR_H 1

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifefem/dofInterfaceBase.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifecore/factory.hpp>



namespace LifeV
{

//! BCVectorBase - class that holds the FE vectors used for prescribing boundary conditions.
/*!
   @author Miguel Fernandez
   @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   @author Vincent Martin <vincent.martin@inria.fr>

  This class holds the FE vectors used for prescribing  boundary conditions.
  The FE vectors given by the user must have the dimension of the total DOFs, although only DOFs on the boundary are considered.

  In the case of Essential boundary condition, we want to prescribe u = v on part of the boundary ( u is the solution, v the given FE vector).
  In the case of Natural boundary condition, depending on the type, we want to add to the right hand side of the equation one of the following terms:

  type 0: 	v, 				in this case v is a quantity integrated on the boundary (i.e. a residual) <br>
  type 1: ( v, n phi)_bd  	with v scalar and phi vector <br>
  type 2: ( v n, phi)_bd  	v vector, phi scalar <br>
  type 3: ( v, phi)_bd   	v and phi can be both vectors or scalars  (not yet implemented)

  here ( . , . )_bd  denote the L2 inner product on the boundary, n is the normal, v the given FE vector

  This class holds data structure also for Resistance and Flux boundary conditions

  This is the base class for other BCVectorXXX classes.
  Inheritance is used to hold specific boundary condition data.
*/

class BCVectorBase
{
public:

    //! @name Constructors and Destructor
    //@{


    //! Empty Constructor
    /*!
     * The user must call setBCVector(..).
     */
    BCVectorBase();


    //! Constructor
    /*!
      @param rightHandSideVector The given Finite Element vector holding data to prescribe on boundary
      @param numberOfTotalDof number of total dof in the vector of data
      @param type The type can assume the following values (0, 1, 2); see BCVector class description for their meaning
    */

    BCVectorBase( const EpetraVector& rightHandSideVector, const UInt numberOfTotalDof, UInt type = 0 );


    //! Copy Constructor
    BCVectorBase( const BCVectorBase& bcVectorBase);


    //! Destructor
    virtual ~BCVectorBase()
    {
        // nothing to be done here
    }


    //@}

    //! @name Operators
    //@{

    //! Assignment operator
    virtual BCVectorBase& operator=( BCVectorBase const& );

    //! Return the value of the selected component of rightHandSideVector at position globalDofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    virtual Real operator() ( const ID& globalDofId, const ID& component ) const;

    //@}


    //! @name Methods

    //!  Return the value of the selected component of the boundary mass coefficient vector at position dofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    virtual Real MixteVec( const ID& globalDofId, const ID& component ) const;


    //!  Return the value of the selected component of the beta coefficient vector at position dofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    virtual Real BetaVec( const ID& globalDofId, const ID& component ) const;


    //!  Return the value of the selected component of the gamma coefficient vector at position dofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    virtual Real GammaVec( const ID& globalDofId, const ID& component ) const;

    //! showMe
    /*!
     * @param verbose The verbosity
     * @param out The output stream (default: cout)
     */
    virtual std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const = 0;


    //@}

    //! @name Get Methods
    //@{




    //! Return the number of total DOF
    inline UInt nbTotalDOF() const { return M_numberOfTotalDof; }


    //! Return the type of conditions (see BCVector class description)
    inline UInt type() const { return M_type; }


    //! determine whether the BCVector is updated
    inline bool isFinalized() const { return M_finalized;}


    //! determine whether the boundary mass coefficient for Robin bc is a Vector
    inline bool ismixteVec() const  {return M_isRobinBdMassCoeffAVector;}


    //true if beta coefficient is a Vector
    inline bool isbetaVec() const  {return M_isBetaCoeffAVector;}


    //true if gamma coefficient is a Vector
    inline bool isgammaVec() const {return M_isGammaCoeffAVector;}


    //! Return the value of the boundary mass coefficient of Robin conditions
    inline Real mixteCoef() const { return M_robinBoundaryMassCoeff; }


    //! Return the value of the resistance coefficient
    inline Real resistanceCoef() const { return M_resistanceCoeff; }


    //! Return the value of the beta coefficient
    inline Real betaCoef() const { return M_betaCoeff; }


    //! Return the value of the gamma coefficient
    inline Real gammaCoef() const { return M_gammaCoeff; }

    //@}


    //! @name Set Methods
    //@{

    //! set the finalize status
    /*!
      @param isFinalized Boolean value saying whether BCVectorBase have been updated
    */
    inline void setFinalized( bool isFinalized ) { M_finalized = isFinalized; }


    //! set M_isRobinBdMassCoeffAVector
    /*!
      @param isRobinBdMassCoeffAVector Boolean value saying whether the boundary mass coefficient is a vector
     */
    inline void setmixte( bool isRobinBdMassCoeffAVector ) { M_isRobinBdMassCoeffAVector = isRobinBdMassCoeffAVector; }


    //! set M_isBetaCoeffAVector
    /*!
      @param isRobinBdMassCoeffAVector Boolean value saying whether the beta coefficient is a vector
     */
    inline void setisbeta( bool isBetaCoeffAVector) { M_isBetaCoeffAVector = isBetaCoeffAVector; }


    //! set M_isGammaCoeffAVector
    /*!
      @param isGammaCoeffAVector Boolean value saying whether the gamma coefficient is a vector
     */
    inline void setisgamma( bool isGammaCoeffAVector ) { M_isGammaCoeffAVector = isGammaCoeffAVector; }


    //! set the boundary mass coefficient of Robin bc
    /*!
      @param robinBoundaryMassCoeff The boundary mass coefficient of robin conditions
     */
    inline void setMixteCoef( const Real& robinBoundaryMassCoeff ) { M_robinBoundaryMassCoeff = robinBoundaryMassCoeff;}


    //! set the Resistance coefficient
    /*!
      @param robinBoundaryMassCoeff The boundary mass coefficient of robin conditions
     */

    inline void setResistanceCoef( const Real& resistanceCoeff ) { M_resistanceCoeff = resistanceCoeff; }


    //! set the Beta coefficient FE vector
    /*!
      @param betaCoeff The beta coefficient
     */

    inline void setBetaCoef( const Real& betaCoeff ) { M_betaCoeff = betaCoeff;}


    //! set the Gamma coefficient FE vector
    /*!
      @param gammaCoeff The gamma coefficient
     */

    void setGammaCoef( const Real& gammaCoeff ) { M_gammaCoeff = gammaCoeff; }


    //! set the boundary mass coefficient FE vector for Robin boundary conditions
    /*!
      @param robinBoundaryMassCoeff The boundary mass coefficient of robin conditions
     */

    void setMixteVec( const EpetraVector& robinBoundaryMassCoeffVector );



    //! set the beta coefficient FE vector
    /*!
      @param betaCoeffVector The beta coefficient FE vector
     */
    void setBetaVec( const EpetraVector& betaCoeffVector );


    //! set the gamma coefficient FE vector
    /*!
      @param gammaCoeffVector The gamaa coefficient FE vector
     */
    void setGammaVec( const EpetraVector& gammaCoeffVector );


    //! set the right hand side FE vector
    /*!
      @param righHandSideVector
      @param numberOfTotalDOF
      @param type
     */
    void setVector( const EpetraVector& righHandSideVector, UInt numberOfTotalDOF, UInt type=0 );

    //@}

protected:

    //! The pointer to  FE vector for the right hand side part of the equation
    const EpetraVector* M_rightHandSideVectorPtr;

    //! The pointer to FE Vector holding the robin boundary Mass coefficients
    const EpetraVector* M_robinBoundaryMassCoeffVectorPtr;
    const EpetraVector* M_betaCoeffVectorPtr;
    const EpetraVector* M_gammaCoeffVectorPtr;

    //! Number of total dof in the vector of data
    UInt M_numberOfTotalDof;

    //! Coefficient for boundary mass term in Robin conditions
    Real M_robinBoundaryMassCoeff;

    //! Coefficient for Resistance coefficient
    Real M_resistanceCoeff;


    //! Coefficient for the beta coefficient
    Real M_betaCoeff;

    //! Coefficient for gamma coefficient
    Real M_gammaCoeff;

    //! boolean determining whether the boundary mass coefficient is a FE Vector
    bool M_isRobinBdMassCoeffAVector;

    //! boolean determining whether the boundary mass coefficient is a FE Vector
    bool M_isBetaCoeffAVector;

    //! boolean determining whether the boundary mass coefficient is a FE Vector
    bool M_isGammaCoeffAVector;


    //!  Type of boundary condition; see the BCBase class description
    UInt M_type;

private:

    //! true when the BCVector is updated
    bool M_finalized;

};



// ============ BCVector ================


//! BCVector - class that holds the FE vectors used for prescribing boundary conditions.
/*!
   @author Miguel Fernandez
   @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   @author Vincent Martin <vincent.martin@inria.fr>

  This class holds the FE vectors used for prescribing  boundary conditions. It is derived from the class BCVectorBase
  The FE vectors given by the user must have the dimension of the total DOFs, although only DOFs on the boundary are considered.

  In the case of Essential boundary condition, we want to prescribe u = v on part of the boundary ( u is the solution, v the given FE vector).
  In the case of Natural boundary condition, depending on the type, we want to add to the right hand side of the equation one of the following terms:

  type 0: @c	v, 				in this case v is a quantity integrated on the boundary (i.e. a residual) <br>
  type 1: @c ( v, n phi)_bd  	with v scalar and phi vector <br>
  type 2: @c ( v n, phi)_bd  	v vector, phi scalar  <br>
  type 3: @c ( v, phi)_bd   	v and phi can be both vectors or scalars  (not yet implemented)

  here ( . , . )_bd  denote the the L2 inner product on the boundary, n is the normal, v the given FE vector and phi the FE test function

  This class holds data structure also for Robin, Resistance and Flux boundary conditions.

  Since the FE vector v is used in the right hand side of the equation, the associate code variable is @c M_rightHandSideVector  (EpetraVector)

  In Robin boundary conditions we add to the matrix the term: <br>
  @c (coeff u, phi)_bd   <br>
  and to the right hand side of the equation the term: <br>
  @ ( v, phi)_bd     <br>
  The code variables associated with coeff are M_boundaryMassCoeff and M_boundaryMassCoeffVector, the former being a scalar (Real), the latter an FE Vector (EpetraVector)

   This is the base class for other BCVectorInterface class.
  Inheritance is used to hold specific boundary condition data.
*/

class BCVector:
        public BCVectorBase
{
public:

    //! @name Public Types
    //@{
    //! super class
    typedef BCVectorBase super; //deprecated
    typedef BCVectorBase bcVectorBase_Type;

    //@}

    //! @name Constructors and Destructor
    //@{

    //! Default Constructor
    /*!
     * The user must call setVector(..)
     */
    BCVector() {}

    //! Constructor
    /*!
    	@param rightHandSideVector The given Finite Element vector holding data to prescribe on boundary
    	@param numberOfTotalDof number of total dof in the vector of data
    	@param type The type can assume the following values (0, 1, 2); see BCVector class description for their meaning
    */
    BCVector( EpetraVector& rightHandSideVector, UInt const numberOfTotalDof, UInt type=0 );


    //Copy Constructor
    BCVector( const BCVector& bcVector );

    //! Destructor
    virtual ~BCVector( ) {}

    //@}


    //! @name Operators
    //@{

    //! Assignment operator for BCVector
    BCVector & operator=( const BCVector & bcVector );

    //@}

    //! @name Methods
    //@{

    //! showMe
    /*!
    * @param verbose The verbosity
    * @param out The output stream (default: cout)
    */
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;

    //@}
};

// ============ BCVectorInterface ================

//! BCVectorInterface - class that holds the FE vectors used for prescribing boundary conditions on Interfaces.
/*!
   @author Miguel Fernandez
   @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   @author Vincent Martin <vincent.martin@inria.fr>

  This class holds the FE vectors used for prescribing  boundary conditions. It is derived from the class BCVectorBase
  The FE vectors given by the user must have the dimension of the total DOFs, although only DOFs on the boundary are considered.

  In the case of Essential boundary condition, we want to prescribe u = v on part of the boundary ( u is the solution, v the given FE vector).
  In the case of Natural boundary condition, depending on the type, we want to add to the right hand side of the equation one of the following terms:

  type 0: 	v, 				in this case v is a quantity integrated on the boundary (i.e. a residual) <br>
  type 1: ( v, n phi)_bd  	with v scalar and phi vector  <br>
  type 2: ( v n, phi)_bd  	v vector, phi scalar  <br>
  type 3: ( v, phi)_bd   	v and phi can be both vectors or scalars  (not yet implemented)

  here ( . , . )_bd  denote the L2 inner product on the boundary, n is the normal, v the given FE vector

  This class holds data structure also for Robin, Resistance and Flux boundary conditions.

  Since the FE vector v is used in the right hand side of the equation, the associate code variable is @c M_rightHandSideVector  (EpetraVector)

  In Robin boundary conditions we add to the matrix the term: <br>
  @c (coeff u, phi)_bd   <br>
  and to the right hand side of the equation the term: <br>
  @ ( v, phi)_bd     <br>
  The code variables associated with coeff are M_boundaryMassCoeff and M_boundaryMassCoeffVector, the former being a scalar (Real), the latter an FE Vector (EpetraVector)

  This is the base class for other BCVectorInterface class.
  Inheritance is used to hold specific boundary condition data.
*/

class BCVectorInterface
        :
        public BCVectorBase
{
public:

    //! @name Public Types
    //@{


    typedef BCVectorBase super;  //deprecated
    typedef boost::shared_ptr<DofInterfaceBase> dof_interface_type; //deprecated

    typedef BCVectorBase bcVectorBase_Type;
    typedef boost::shared_ptr<DofInterfaceBase> dofInterfacePtr_Type;

    //@}


    //! @name Constructors and Destructor
    //@{


    //! Default Constructor
    /*!
     * The user must call setVector(..)
     */
    BCVectorInterface () {}

    //! Constructor
    /*!
      @param rightHandSideVector The given Finite Element vector holding data to prescribe on boundary
      @param numberOfTotalDof Number of total dof in the vector of data
      @param interfaceDofPtr The pointer to the container of connections between the DOFs on two matching meshes
      @param type The type can assume the following values (0, 1, 2); see BCVector class description for their meaning
    */
    BCVectorInterface( const EpetraVector& rightHandSideVector, UInt numberOfTotalDof, const dofInterfacePtr_Type& interfaceDofPtr, UInt type=0 );


    //! Copy Constructor
    BCVectorInterface( const BCVectorInterface & bcVectorInterface );


    //!Destructor
    virtual ~BCVectorInterface() {}


    //@}



    //! @name Operators
    //@{

    //! Assignment operator for BCVectorInterface
    BCVectorInterface & operator=( const BCVectorInterface & bcVectorInterface );


    //! Return the value of the selected component of rightHandSideVector at position globalDofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    Real operator() ( const ID& globalDofId, const ID& component ) const;



    //@}



    //! setup after default constructor
    /*!
          @param rightHandSideVector The given Finite Element vector holding data to prescribe on boundary
          @param numberOfTotalDof Number of total dof in the vector of data
          @param interfaceDofPtr The pointer to the container of connections between the DOFs on two matching meshes
          @param type The type can assume the following values (0, 1, 2); see BCVectorInterface class description for their meaning
        */
    void setup ( const EpetraVector& rightHandSideVector, UInt numberOfTotalDof, const dofInterfacePtr_Type& interfaceDofPtr, UInt type=0 );


    //! set the BC vector (after default construction)
    /*!
      @param rightHandSideVector The given Finite Element vector holding data to prescribe on boundary
      @param numberOfTotalDof Number of total dof in the vector of data
      @param interfaceDofPtr The pointer to the container of connections between the DOFs on two matching meshes
      @param type The type can assume the following values (0, 1, 2); see BCVectorInterface class description for their meaning
    */
    void setVector( const EpetraVector& rightHandSideVector, UInt numberOfTotalDof, const dofInterfacePtr_Type& interfaceDofPtr, UInt type=0);



    //! @name Methods
    //@{

    //!  Return the value of the selected component of the boundary mass coefficient vector at position dofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    Real MixteVec( const ID& globalDofId, const ID& component ) const;


    //!  Return the value of the selected component of the beta coefficient vector at position dofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    Real BetaVec( const ID& globalDofId, const ID& component ) const;


    //!  Return the value of the selected component of the gamma coefficient vector at position dofID
    /*!
      @param globalDofId The global DOF id
      @param component The vector component
    */
    Real GammaVec( const ID& globalDofId, const ID& component ) const;



    //! showMe
    /*!
    * @param verbose The verbosity
    * @param out The output stream (default: cout)
    */
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;


    //! Return reference to  DofInterfaceBase object, the container of connection of DOFs
    inline DofInterfaceBase const & dofInterface() const { return *M_interfaceDofPtr; }

    //@}

protected:

    //! DofInterfaceBase object holding the connections between the interface dofs
    dofInterfacePtr_Type M_interfaceDofPtr;

};
typedef LifeV::singleton< LifeV::factoryClone< BCVectorBase > > FactoryCloneBCVector;
}
#endif
