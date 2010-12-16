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
    @brief classes to handle boundary conditions

    @author M.A. Fernandez
    @author M.Prosi
    @contributor Lucia Mirabella <lucia.mirabell@gmail.com>
    @maintainer Lucia Mirabella <lucia.mirabell@gmail.com>

    @date 06-2002

    @date 11-2002

  This file contains the classes which may be used to store boundary
  conditions. A boundary condition object will have the following
  elements:
<ol>
  <li> a name identifying a specific BC,

  <li> a flag identifying a specific part of the mesh boundary,

  <li> a type (Essential, Natural, Mixte, Flux, Resistance),

  <li> a mode of implementation (Scalar, Full, Component, Normal,
     Tangential, Resistance, Directional),

  <li> a functor holding the data function,

  <li> a bool vector  describing the components involved in this boundary condition

  <li> a list of pointers to identifiers allowing the user to know to
     which DOF the boundary condition applies.
</ol>

 */

#ifndef BCCOND_H
#define BCCOND_H

#include <set>
#include <map>

#include <boost/shared_ptr.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifemesh/identifier.hpp>
#include <life/lifemesh/markers.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/currentFE.hpp>
#include <life/lifefem/currentBdFE.hpp>

#include <life/lifefem/bcVector.hpp>
#include <life/lifefem/bcFunction.hpp>
#include <life/lifearray/EpetraVector.hpp>


namespace LifeV
{

/*! @enum BCType
	Boundary condition basic types: Natural, Mixte, Flux, Resistance, Periodic, Essential, EssentialEdges, EssentialVertices
 */
enum BCType
{
    Natural, 			/*!< Neumann boundary conditions */
    Mixte, 				/*!< Robin boundary conditions */
    Flux, 				/*!< Flux boundary conditions */
    Resistance,			/*!< Resistance boundary conditions */
    Periodic,			/*!< Periodic boundary conditions */
    Essential, 			/*!< Dirichlet boundary conditions */
    EssentialEdges, 	/*!< Dirichlet boundary conditions on edges */
    EssentialVertices, 	/*!< Dirichlet boundary conditions on vertices */
};

/*! @enum BCMode
  	Type for boundary conditions application modes
 */
enum BCMode
{
    Scalar, 	/*!< To be used for scalar problems */
    Full, 		/*!< To be used for vector problems, when the boundary condition involves all components*/
    Component, 	/*!< To be used for vector problems, when the boundary condition DOESN'T involve all components*/
    Normal, 	/*!< To be used for vector problems, when the boundary condition involve the normal component*/
    Tangential,	/*!< To be used for vector problems, when the boundary condition involve the tangential component*/
    Directional /*!< To be used for vector problems, when the boundary condition involve a specific direction*/
};



/*! Type of the name of the Boundary conditions
 */
typedef std::string BCName;
const BCName nullBCName; 		//!< Empty string



//! BCBase - Base class which holds the boundary condition information
/*!
    @author M.A. Fernandez
    @author M.Prosi
    @see

  For each boundary condition the user must give
<ol>
  <li> a name,

  <li> a mesh flag,

  <li> a type,

  <li> a mode,

  <li> a data BCFunction,

  <li> three (or two in 2D) bools describing the components involved in
  this boundary condition.
</ol>
  Finally the list of pointers to identifiers will be updated in the
  Dof class (\c BCHandler::bdUpdate method).

  \warning The idea is to not use inheritance from this class

 */
class BCBase
{
public:

    //! BCHandle is a friend class of BCBase
    friend class BCHandler;

    //! @name Public Types
    //@{

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    BCBase();

    //! Constructor for BCBase
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Scalar, Full,
       Component, Normal, Tangential, Directional
       @param bcFunction the function holding the user defined function defining the boundary condition
       @param components vector of IDs storing the list of components involved in this boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionBase& bcFunction,
            const std::vector<ID>& components );

    //! Constructor for BCBase without specifying components for Scalar, Tangential or Normal mode problems
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Scalar, Normal, Tangential
       @param bcFunction the BCFunctionBase holding the function defining the boundary condition
       involved in this boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionBase& bcFunction );

    //! Constructor for BCBase without specifying components for without list of components for Full mode problems
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Full
       @param bcFunction BCFunctionBase holding the function defining the boundary condition
       @param numberOfComponents number of components involved
       in this boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionBase& bcFunction,
            const UInt& numberOfComponents );

    //! Constructor for BCBase to prescribe a boundary condition from a vector of dof values
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Scalar, Full,
       Component, Normal, Tangential, Directional
       @param vector the vector containing the dof values to be prescribed as boundary data
       @param components vector of IDs storing the list of components involved in this boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCVectorBase& vector,
            const std::vector<ID>& components );

    //! Constructor for BCBase to prescribe a boundary condition from a vector of dof values  without specifying components for Scalar, Tangential or Normal mode problems
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Scalar, Full,
       Component, Normal, Tangential, Directional
       @param bcVector the vector containing the dof values to be prescribed as boundary data
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCVectorBase& bcVector );

    //! Constructor for BCBase to prescribe a boundary condition from a vector of dof values  without specifying components for Full mode problems
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Scalar, Full,
       Component, Normal, Tangential, Directional
       @param bcVector the vector containing the dof values to be prescribed as boundary data
       @param numberOfComponents number of components involved in this boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCVectorBase& bcVector,
            const UInt& numberOfComponents );

    //! Constructor for BCBase. The BC function depends on a generic FE vector (e.g. the solution at the previous time step)
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Scalar, Full,
       Component, Normal, Tangential, Directional
       @param bcFunctionFEVectorDependent the BCFunctionUDepBase holding the function (depending on a generic finite element vector ) defining the boundary condition
       @param components vector of IDs storing the list of components involved in this boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionUDepBase& bcFunctionFEVectorDependent,
            const std::vector<ID>& components );

    //! Constructor for BCBase without specifying components for Scalar, Tangential or Normal mode problems. The BC function depends on a generic FE vector (e.g. the solution at the previous time step)
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode: Scalar, Normal, Tangential
       @param bcFunctionFEVectorDependent the BCFunctionUDepBase holding the function (depending on a generic finite element vector ) defining the boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionUDepBase& bcFunctionFEVectorDependent);

    //! Constructor for BCBase without specifying components for Full mode problems. The BC function depends on a generic FE vector (e.g. the solution at the previous time step)
    /*!
       @param name the name of the boundary condition
       @param flag the mesh flag identifying the part of the mesh where
       the boundary condition applies
       @param type the boundary condition type: Natural, Essential, Mixte, Flux, Resistance
       @param mode the boundary condition mode:  Full
       @param bcFunctionFEVectorDependent the BCFunctionUDepBase holding the function (depending on a generic finite element vector ) defining the boundary condition
       @param numberOfComponents number of components involved in this boundary condition
     */
    BCBase( const std::string& name,
            const entityFlag_Type& flag,
            const BCType& type,
            const BCMode& mode,
            BCFunctionUDepBase& bcFunctionFEVectorDependent,
            const UInt& numberOfComponents );

    //! Copy constructor for BCBase
    /*!
     @param bcBase a BCBase object
     @warning This is not a copy constructor since the lists are built empty
     */
    BCBase( const BCBase& bcBase );

    //! Destructor
    ~BCBase();
    //@}


    //! @name Methods
    //@{
    //! Returns the index of the component of the solution associated to the iComponent-th component prescribed in the boundary condition at hand
    /*!
    	Example: the solution has 3 components and we prescribe a boundary condition on component 1 and 3.
    	Then, component(2) returns 3, since 3 is the index of the 2nd BC prescribed.

       @param iComponent the "local" component (from 1 to numberOfComponents)
       @return the index of the component of the solution associated to the iComponent-th component prescribed in the boundary condition at hand
     */
    ID component( const ID i ) const;

    //! Returns true if mixte (in BC Vector ) is a EpetraVector (mixteVec), false if it is scalar (default alphaCoef=1)
    /*!
       @return true if mixte (in BC Vector ) is a EpetraVector (mixteVec), false if it is scalar
     */
    bool ismixteVec() const;

    //! Returns true if beta (in BC Vector ) is a EpetraVector (betaVec) (default betaCoef=1)
    /*!
       @return true if beta (in BC Vector ) is a EpetraVector (betaVec)
     */
    bool isbetaVec() const;

    //! Returns true if gamma (in BC Vector ) is a EpetraVector (gammaVec) (default gammaCoef=1)
    /*!
       @return true if gamma (in BC Vector ) is a EpetraVector (gammaVec)
     */
    bool isgammaVec() const;

    //! Returns the value of the mixte coefficient vector (in BC Vector)
    /*!
       corresponding to DOF iDof and component iComponent
       @param iDof DOF we are looking for in MixteVec
       @param iComponent component we are looking for in MixteVec
       @return value of the mixte coefficient vector (in BC Vector) corresponding to iDof and iComponent
     */
    Real MixteVec( const ID& iDof, const ID& iComponent ) const;

    //! Returns the value of the beta coefficient vector (in BC Vector)
    /*!
       corresponding to DOF iDof and component iComponent
       @param iDof DOF we are looking for in BetaVec
       @param iComponent component we are looking for in BetaVec
       @return value of the Beta coefficient vector (in BC Vector) corresponding to iDof and iComponent
     */
    Real BetaVec( const ID& iDof, const ID& iComponent ) const;

    //! Returns the value of the gamma coefficient vector (in BC Vector)
    /*!
       corresponding to DOF iDof and component iComponent
       @param iDof DOF we are looking for in GammaVec
       @param iComponent component we are looking for in GammaVec
       @return value of the Gamma coefficient vector (in BC Vector) corresponding to iDof and iComponent
     */
    Real GammaVec( const ID& iDof, const ID& iComponent ) const;


    //! Returns a pointer to the BCFunctionBase object
    /*!
       @return pointer to the BCFunctionBase object
     */
    const BCFunctionBase* pointerToFunctor() const;

    //! Returns a pointer to the BCFunctionUDepBase object
    /*!
       @return pointer to the BCFunctionUDepBase object
     */
    const BCFunctionUDepBase* pointerToFunctorUDep() const;

    //! Returns a pointer to the BCVector object
    /*!
       @return pointer to the BCVector object
     */
    const BCVectorBase* pointerToBCVector() const;

    //! Adds a new identifier to the list
    /*!
       @param identifierToAddPtr pointer to the IdentifierBase object to be added
     */
    void addIdentifier( IdentifierBase* identifierToAddPtr );

    //! Adds a new identifier to the list of IdGlobal
    /*!
       @param identifierToAddPtr pointer to the IdentifierBase object to be added
     */
    void addIdentifierIdGlobal( IdentifierBase* identifierToAddPtr);

    //! Returns the size of the identifiers list
    /*!
       @return the size of the identifiers list
     */
    UInt list_size() const;

    //! Returns the size of global identifiers list
    /*!
       @return the size of global identifiers list
     */
    UInt list_size_IdGlobal() const;

    //! Returns element of M_IdGlobal vector
    int IdGlobal( int id) const;

    //! Method that writes info in output
    /*!
       @param verbose to specify the level of verbosity (false by default)
       @param outStream to specify the output stream (std::cout by default)
     */
    std::ostream & showMe( bool verbose = false, std::ostream & outStream = std::cout ) const;
    //@}


    //! @name Operators
    //@{

    //! The assignment operator for BCBase
    /*!
        @param bcBase a BCBase object
        @return Reference to a new BCBase with the same
                content of bcBase
     */
    BCBase & operator=( const BCBase& bcBase);

    //! Returns a pointer to the (i)-th element of the list of identifiers
    /*!
       The list of identifiers has to be finalized before calling this operator.
       @param i index of the element in the list of identifier that we want to be returned (starting from 0)
     */
    const IdentifierBase* operator[] ( const ID& i ) const;

    //! Returns a pointer to the (i-1)-th element of the list of identifiers
    /*!
       The list of identifiers has to be finalized before calling this operator.
       @param i index of the element in the list of identifier that we want to be returned (starting from 1)
     */
    const IdentifierBase* operator() ( const ID& i ) const;


    //! Overloading function operator by calling the BCFunctionBase user specified function
    /*!
       @param t time
       @param x coordinate
       @param y coordinate
       @param z coordinate
       @param iComponent component of the vector function
       @return iComponent of the user defined function evaluated in (t,x,y,z)
     */
    Real operator() ( const Real& t, const Real& x, const Real& y,
                      const Real& z, const ID& iComponent ) const;

    //! Overloading function operator by calling the BCFunctionUDepBase user specified function
    /*!
       @param t time
       @param x coordinate
       @param y coordinate
       @param z coordinate
       @param iComponent component of the vector function
       @param u value of the FE vector in t, x, y, z, component iComp
       @return i-component of the user defined function evaluated in (t,x,y,z,u)
     */
    Real operator() ( const Real& t, const Real& x, const Real& y,
                      const Real& z, const ID& iComponent, const Real& u ) const;


    //! Overloading function operator by querying the BCVector in DOF iDof and component iComponent
    /*!
       @param iDof global dof index
       @param iComponent component index
     */
    Real operator() ( const ID& iDof, const ID& iComponent ) const;

    //! Overloading "less-than" operator between BCBase objects
    /*!
       The "smaller" (or weaker) boundary conditions is the one to be applied first
       @param bcBase1 first BCBase to compare
       @param bcBase2 second BCBase to compare
       @return True if bcBase1 is smaller, False if bcBase1 is bigger
     */
    friend bool operator<( const BCBase& bcBase1, const BCBase& bcBase2 )
    {
        if (bcBase1.type() == bcBase2.type())
            return (bcBase1.flag() < bcBase2.flag());
        else
            return ( bcBase1.type() < bcBase2.type() );
    }


    //! Overloading "is-equal" operator for BCBase objects
    /*!
       Check if the flag of bcBase is equal to flag argument
       @param bcBase BCBase to check
       @param flag EntityFlag to be compared with bcBase flag
       @return True if bcBase's flag is equal to flag
     */
    friend bool operator==( const BCBase& bcBase, const entityFlag_Type flag )
    {
        return bcBase.flag() == flag;
    }



    //@}

    //! @name Set Methods
    //@{
    //! set BCVectorBase boundary condition
    /*!
       @param bcVector to be set in BCBase class
     */
    void setBCVector( const BCVectorBase& bcVector );

    //! set BCFunctionBase boundary condition
    /*!
       @param bcFunction to be set in BCBase class
     */
    void setBCFunction( const BCFunctionBase& bcFunction );

    //! set BCFunctionUDepBase boundary condition
    /*!
       @param bcFunctionFEVectorDependent to be set in BCBase class
     */
    void setBCFunction( const BCFunctionUDepBase& bcFunctionFEVectorDependent );

    //! Set the BC offset
    /*!
       @param bcOffset to be set in BCBase class
     */
    void setOffset(int bcOffset) {M_offset = bcOffset;}
    //@}

    //! @name Get Methods
    //@{

    //! Returns the boundary condition name
    /*!
       @return boundary condition name
     */
    std::string name() const;

    //! Returns the flag associated to the boundary condition name
    /*!
       @return boundary condition flag
     */
    entityFlag_Type flag() const;

    //! Returns the boundary condition type
    /*!
       @return boundary condition type
     */
    BCType type() const;

    //! Returns the boundary condition mode
    /*!
       @return boundary condition mode
     */
    BCMode mode() const;

    //! Returns the number of components involved in this boundary condition
    /*!
       @return number of components prescribed by this boundary condition
     */
    UInt numberOfComponents() const;


    //! Returns the offset associated to this boundary condition
    /*!
       @return offset associated to this boundary condition
     */
    const int& offset() const {return M_offset;}

    //! Returns the value of the mixte coefficient (in BC Vector)
    /*!
       @return value of the mixte coefficient (in BC Vector)
     */
    Real mixteCoef() const;

    //! Returns the value of the resistance coefficient (in BC Vector)
    /*!
       @return value of the resistance coefficient (in BC Vector)
     */
    Real resistanceCoef() const;

    //! Returns the value of the beta coefficient (in BC Vector)
    /*!
       @return value of the beta coefficient (in BC Vector)
     */
    Real betaCoef() const;

    //! Returns the value of the gamma coefficient (in BC Vector)
    /*!
       @return value of the gamma coefficient (in BC Vector)
     */
    Real gammaCoef() const;

    //! Returns True if a FE BCVector has been provided to the class, False otherwise
    /*!
       @return True if FE BCVector has been provided to the class, False otherwise
     */
    bool dataVector() const;

    //! Returns whether the list is finalized and the vector of ID's is then accessible.
    /*!
      @return M_finalized private member
     */
    bool finalized() const;

    //! Returns whether the list is finalized and the vector of IDGlobal's is then accessible.
    /*!
      @return M_finalizedIdGlobal private member
     */
    bool finalizedIdGlobal() const;


    //! Returns True if the BCBase is based on a BCFunctionUDepBase function, False otherwise
    /*!
       @return True if the BCBase is based on a BCFunctionUDepBase function, False otherwise
     */
    bool isUDep() const;


    //@}
private:

    std::string                           M_name; 	//!< name of the boundary condition

    entityFlag_Type                       M_flag;  //!< flag identifying a specific part of the mesh boundary

    BCType                                M_type;  //!< the boundary condition type

    BCMode                                M_mode;  //!< the boundary condition mode of application

    std::vector<ID>                       M_components;		//! the list of components involved in this BC

    boost::shared_ptr<BCFunctionBase>     M_bcFunction;   //!< Pointer to a user defined BC function

    boost::shared_ptr<BCFunctionUDepBase> M_bcFunctionFEVectorDependent; //!< Pointer to a user defined BC function (depending on a generic FE vector)

    boost::shared_ptr<BCVectorBase >      M_bcVector;   	//!< Pointer to a user given BC vector

    bool                                  M_isStored_BcVector; //! True if a FE BCVector has been provided

    bool                                  M_isStored_BcFunctionVectorDependent; //!< True if the BCBase is based on a BCFunctionUDepBase function, False otherwise

    std::set<boost::shared_ptr<IdentifierBase>, identifierComp> M_idSet; //!< set of pointers to identifiers allowing the user to get hold the DOF to which the BC applies

    std::set<boost::shared_ptr<IdentifierBase>, identifierComp> M_idGlobalSet; //!< set of pointers to identifiers allowing the user to get hold the global DOF to which the BC applies

    std::vector<int> M_IdGlobalVector; //!< vector of IdGlobal

    std::vector<boost::shared_ptr<IdentifierBase> > M_idList; //!< container for id's when the list is finalized

    std::vector<boost::shared_ptr<IdentifierBase> > M_idGlobalList; //!< container for global id's when the list is finalized

    int M_offset; //!< boundary condition offset

    bool M_finalized; //!< True, when M_idList is finalized

    bool M_finalizedIdGlobal; //!< True, when M_idGlobalList is finalized

    void finalize(); //!< Transfer information from M_idSet to M_idList

    void  finalizeIdGlobal(); //!< Transfer information from M_idGlobalSet to M_idGlobalList
};





}

#endif
