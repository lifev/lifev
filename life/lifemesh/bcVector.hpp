/*!
  \file bcVector.h
  \brief classes to handle data vectors for boundary conditions.  
  \version 1.0
  \author M.A. Fernandez 
  \date 11/2002 
  
  \version 1.1
  \author V. Martin 
  \date 02/2003 

  This file contains the classes which may be used to store boundary
  conditions. 
*/

#ifndef __BCVECTOR__
#define __BCVECTOR__

#include "lifeV.hpp"
#include "dofInterfaceBase.hpp"
#include "vecUnknown.hpp"

// ============ BCVector_Interface ================ 

/*! 

 \class BCVector_Interface

 Class that holds a user data vector for boundary conditions

 The data funcitions given by the user must have the following declaration
 Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
*/

class BCVector_Interface {
 public:

  //! Constructor
  /*!
    \param vec data vector holding data
    \param nbTotalDof number of total dof in the vector of data
    \param dofIn dofInterfaceBase object holding the connections between the interface dofs of the
           data vector and those of the associated to the boundary conditions
  */
  BCVector_Interface( Vector& vec, UInt nbTotalDof, DofInterfaceBase& dofIn);
   
  //! Default Constructor (the user must call setBCVector(..))
  BCVector_Interface ();

  //! set the BC vector (after default construction)
  void setBCVector_Interface( Vector& vec, UInt nbTotalDof, DofInterfaceBase& dofIn );

  //! set the Mixte coefficient
  void setBCVec_MixteCoef( const Real& coef );

  //! This method returns the value to be imposed in the component iComp of the dof iDof
  /*!
    \param iDof the number of the Dof
    \param iComp the number of the component
  */
   Real operator()(const ID& iDof, const ID& iComp) const;

   //! Return the value of the Mixte coefficient
   Real MixteCoef() const;

   //! Assignment operator for BCVector_Interface
   BCVector_Interface & operator=(const BCVector_Interface & BCv);


   //! Output
   ostream &  showMe(bool verbose=false, ostream & out=cout) const;

 protected:

   //! The data vector 
   Vector* _vec;
   
   //! DofInterfaceBase object holding the connections between the interface dofs
   DofInterfaceBase* _dofIn;

   //! Number of total dof in the vector of data
   UInt _nbTotalDof;

   //! Coefficient for mixte boundary conditions (Robin) 
   /*! For the moment, it is the same for all the entries of the data vector.
    */
   Real _MixteCoef;

   //! true when the BCVector is updated (and can be used)
   bool _finalized;

};

#endif
