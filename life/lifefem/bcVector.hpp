/*
 This file is part of the LifeV library

 Authors: Miguel Fernandez
          Vincent Martin
          Christophe Prud'homme <christophe.prudhomme@epfl.ch>

 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file bcVector.hpp
  \brief classes to handle data vectors for boundary conditions.
  \author M.A. Fernandez
  \author V. Martin
  \author C. Prud'homme

  This file contains the classes which may be used to store boundary
  conditions.
*/

#ifndef __BCVECTOR__
#define __BCVECTOR__

#include "lifeV.hpp"
#include "dofInterfaceBase.hpp"
#include "vecUnknown.hpp"
#include <boost/function.hpp>

#include <singleton.hpp>
#include <factory.hpp>



namespace LifeV
{
/*!

 \class BCVectorBase

 Base class that holding data vector for boundary conditions

*/
class BCVectorBase
{
public:

    //! Constructor
    /*!
      \param vec data vector holding data
      \param nbTotalDof number of total dof in the vector of data
    */
    BCVectorBase( Vector& vec, const UInt nbTotalDof );


    //! Constructor
    /*!
      \param vec data vector holding data
      \param nbTotalDof number of total dof in the vector of data
      \param type must be
      -# 0:  boundary integration done (ex. residual of a variational problem)
      -# 1:  needs boundary integration of \f$\lambda n \cdot  \mathbf{\phi}_i\f$
      \warning (implemented only for the Natural BC AM 10/2004)
      -# 2:  needs boundary integration of \f$\mathbf{\lambda} \cdot n \phi_i\f$
      \warning (not yet implemented AM 10/2004)
    */
    BCVectorBase( Vector& vec, const UInt nbTotalDof, UInt type );


    //! Default Constructor (the user must call setBCVector(..))
    BCVectorBase();

    //! Do nothing destructor
    virtual ~BCVectorBase()
        {
            // nothing to be done here
        }

    //@{

    //! assignement operator
    virtual BCVectorBase& operator=( BCVectorBase const& );

    //! This method returns the value to be imposed in the component iComp of the dof iDof
    /*!
      \param iDof the number of the Dof
      \param iComp the number of the component
    */
    virtual Real operator() ( const ID& iDof, const ID& iComp ) const;

    //! Return the value of the Mixte coefficient vector to be imposed in the component iComp of the dof iDof
    virtual Real MixteVec( const ID& iDof, const ID& iComp ) const;

    //@}



    //@{

    //! \return true if finalized, false otherwise
    bool isFinalized() const
        {
            return _M_finalized;
        }

    //! return the total number of DOF
    UInt nbTotalDOF() const
        {
            return _M_nbTotalDof;
        }

    /*!
      Return the type of the kind of information in BCVector

      -# 0:  boundary integration done (ex. residual of a variational problem)

      -# 1: needs boundary integration of \f$\lambda n \cdot \mathbf{\phi}_i\f$
      \warning (implemented only for the Natural BC - AM 10/2004)

      -# 2:  needs boundary integration of \f$\mathbf{\lambda} \cdot n \phi_i\f$
      \warning (not yet implemented - AM 10/2004)
    */
    UInt type() const
        {
            return _M_type;
        }


    //! Return the value of the Mixte coefficient
    Real mixteCoef() const
        {
            return _M_mixteCoef;
        }


    //@}


    //@{

    //! finalize the BC
    void setFinalized( bool __v )
        {
            _M_finalized = __v;
        }

    //! set the Mixte coefficient
    void setMixteCoef( const Real& coef )
        {
            _M_mixteCoef = coef;
        }

    //! set the Mixte coefficient data vector
    void setMixteVec( Vector& vec_mixte )
        {
	    _M_vec_mixte= &vec_mixte;
	}


    //! set the vector
    void setVector( Vector& __vec, UInt nbDOF );

    //@}

    //@{

    //! Output
    virtual std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const = 0;

    //@}

protected:

    //! The data vector
    Vector* _M_vec;

    //! The data vector of the mixte coefficient
    Vector* _M_vec_mixte;

    //! Number of total dof in the vector of data
    UInt _M_nbTotalDof;

    //! Coefficient for mixte boundary conditions (Robin)
    /*! For the moment, it is the same for all the entries of the data vector.
     */
    Real _M_mixteCoef;


    /*!
      Type of boundary condition

      -# boundary integration done (ex. residual of a variational problem)

      -# needs boundary integration of \f$\lambda n \cdot  \mathbf{\phi}_i\f$

      -# needs boundary integration of \f$\mathbf{\lambda} \cdot n \phi_i\f$
    */
    UInt _M_type;

private:

    //! true when the BCVector is updated (and can be used)
    bool _M_finalized;

};



// ============ BCVector ================

/*!

 \class BCVector

 Class that holds a user data vector for boundary conditions

*/

class BCVector:
            public BCVectorBase
{
public:

    //! super class
    typedef BCVectorBase super;

    //! Constructor
    /*!
      \param vec data vector holding data
      \param nbTotalDof number of total dof in the vector of data
             data vector and those of the associated to the boundary conditions
    */
    BCVector( Vector& vec, UInt nbTotalDof );

    //! Constructor
    /*!
      \param vec data vector holding data
      \param nbTotalDof number of total dof in the vector of data
      \param type must be
       -# boundary integration done (ex. residual of a variational problem)
       -# needs boundary integration of \f$\lambda n \cdot  \mathbf{\phi}_i\f$ \warning (not yet implemented - AM 10/2004)
       -# needs boundary integration of \f$\mathbf{\lambda} \cdot n \phi_i\f$ \warning (implemented only for the Natural BC - AM 10/2004)
    */
    BCVector( Vector& vec, UInt const nbTotalDof, UInt type );

    //! Default Constructor (the user must call setVector(..))
    BCVector();

    //! Assignment operator for BCVectorInterface
    BCVector & operator=( const BCVector & BCv );

    //! Output
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;
};

// ============ BCVectorInterface ================

/*!

 \class BCVectorInterface

 Class that holds a user data vector for boundary conditions on
 interfaces

 The data functions given by the user must have the following
 declaration

 \verbatim
 Real g( const Real& time, const Real& x, const Real& y,
         const Real& z, const ID& icomp)
 \endverbatim
*/

class BCVectorInterface
    :
    public BCVectorBase
{
public:

    //! super class
    typedef BCVectorBase super;

    //! Constructor
    /*!
      \param vec data vector holding data
      \param nbTotalDof number of total dof in the vector of data
      \param dofIn dofInterfaceBase object holding the connections between the interface dofs of the
             data vector and those of the associated to the boundary conditions
    */
    BCVectorInterface( Vector& vec, UInt nbTotalDof, DofInterfaceBase& dofIn );

    //! Default Constructor (the user must call setBCVector(..))
    BCVectorInterface ();

    //! set the BC vector (after default construction)
    void setVector( Vector& vec, UInt nbTotalDof, DofInterfaceBase& dofIn );

    /*!
      This method returns the value to be imposed in the component iComp of the dof iDof.

      \param iDof the number of the Dofin
      \param iComp the number of the component
    */
    Real operator() ( const ID& iDof, const ID& iComp ) const;

    //! This method returns the value of the mixte coefficient to be imposed in the component iComp of the dof iDof
    Real MixteVec( const ID& iDof, const ID& iComp ) const;

    //! Assignment operator for BCVectorInterface
    BCVectorInterface & operator=( const BCVectorInterface & BCv );

    //! getter
    DofInterfaceBase const & dofInterface() const
        {
            return *_M_dofIn;
        }

    //! Output
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;

protected:

    //! DofInterfaceBase object holding the connections between the interface dofs
    DofInterfaceBase* _M_dofIn;

};
typedef LifeV::singleton< LifeV::factoryClone< BCVectorBase > > FactoryCloneBCVector;
}
#endif
