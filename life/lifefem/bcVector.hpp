/*
 This file is part of the LifeV library

 Authors: Miguel Fernandez
          Vincent Martin
          Christophe Prud'homme <christophe.prudhomme@epfl.ch>

 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifefem/dofInterfaceBase.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifecore/factory.hpp>



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
      \param type must be
      -# 0:  boundary integration done (ex. residual of a variational problem)
      -# 1:  needs boundary integration of \f$\lambda n \cdot  \mathbf{\phi}_i\f$
      \warning (implemented only for the Natural BC AM 10/2004)
      -# 2:  needs boundary integration of \f$\mathbf{\lambda} \cdot n \phi_i\f$
    */

    BCVectorBase( const EpetraVector& vec, const UInt nbTotalDof, UInt type = 0 );


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

    //! Return the value of the Beta coefficient vector to be imposed in the component iComp of the dof iDof
    virtual Real BetaVec( const ID& iDof, const ID& iComp ) const;

    //! Return the value of the Gamma coefficient vector to be imposed in the component iComp of the dof iDof
    virtual Real GammaVec( const ID& iDof, const ID& iComp ) const;

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

    //! Return the value of the Mixte coefficient
    Real resistanceCoef() const
    {
        return _M_resistanceCoef;
    }

    //! Return the value of the beta coefficient
    Real betaCoef() const
    {
        return _M_betaCoef;
    }
    //! Return the value of the gamma coefficient
    Real gammaCoef() const
    {
        return _M_gammaCoef;
    }
    //@}


    //@{

    //! finalize the BC
    void setFinalized( bool __v )
    {
        _M_finalized = __v;
    }

    //true if Mixte coefficient is a Vector
    bool ismixteVec()const  {return _M_ismixteVec;}

    //true if beta coefficient is a Vector
    bool isbetaVec() const  {return _M_isbetaVec;}

    //true if gamma coefficient is a Vector
    bool isgammaVec() const {return _M_isgammaVec;}

    //setting modes of coefficent

    void setmixte(bool _s) {_M_ismixteVec = _s;}

    void setisbeta(bool _s) {_M_isbetaVec = _s;}

    void setisgamma(bool _s) {_M_isgammaVec =_s;}

    //! set the Mixte coefficient
    void setMixteCoef( const Real& coef )
    {
        _M_mixteCoef = coef;
    }

    //! set the Resistance coefficient
    void setResistanceCoef( const Real& coef )
    {
        _M_resistanceCoef = coef;
    }

    //! set the Mixte coefficient data vector
    void setMixteVec( EpetraVector& vec_mixte )
    {
        _M_ismixteVec = true;
        _M_vec_mixte= &vec_mixte;
    }

    //! set the Beta coefficient data vector
    void setBetaCoef( const Real& coef )
    {
        _M_betaCoef = coef;
    }

    //! set the Gamma coefficient data vector
    void setGammaCoef( const Real& coef )
    {
        _M_gammaCoef = coef;
    }


    //! set the beta coefficient data vector
    void setBetaVec( EpetraVector& vec_beta )
    {
        _M_isbetaVec = true;
        _M_vec_beta= &vec_beta;
    }

    //! set the gamma coefficient data vector
    void setGammaVec( EpetraVector& vec_gamma )
    {
        _M_isgammaVec = true;
        _M_vec_gamma= &vec_gamma;
    }

    //! set the vector
    void setVector( EpetraVector& __vec, UInt nbDOF, UInt type=0 );

    //@}

    //@{

    //! Output
    virtual std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const = 0;

    //@}

protected:

    //! The data vector
    const EpetraVector* _M_vec;

    //! The data vector of the mixte coefficient
    EpetraVector* _M_vec_mixte;
    EpetraVector* _M_vec_beta;
    EpetraVector* _M_vec_gamma;

    //! Number of total dof in the vector of data
    UInt _M_nbTotalDof;

    //! Coefficient for mixte boundary conditions (Robin)
    /*! For the moment, it is the same for all the entries of the data vector.
     */
    Real _M_mixteCoef;

    //! Coefficient for mixte boundary conditions (Resistance)
    /*! For the moment, it is the same for all the entries of the data vector.
     */
    Real _M_resistanceCoef;


    //! Coefficient for mixte boundary conditions (Robin)
    /*! For the moment, it is the same for all the entries of the data vector.
     */
    Real _M_betaCoef;

    //! Coefficient for  boundary conditions (Natural type 0)
    /*! For the moment, it is the same for all the entries of the data vector.
     */
    Real _M_gammaCoef;

    bool _M_ismixteVec;
    bool _M_isbetaVec;
    bool _M_isgammaVec;


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
      \param type must be
      -# boundary integration done (ex. residual of a variational problem)
      -# needs boundary integration of \f$\lambda n \cdot  \mathbf{\phi}_i\f$ \warning (not yet implemented - AM 10/2004)
      -# needs boundary integration of \f$\mathbf{\lambda} \cdot n \phi_i\f$ \warning (implemented only for the Natural BC - AM 10/2004)
    */
    BCVector( EpetraVector& vec, UInt const nbTotalDof, UInt type=0 );

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

//    typedef FSIOperator::vector_type vector_type;
    typedef BCVectorBase super;
    typedef boost::shared_ptr<DofInterfaceBase> dof_interface_type;

    //! Constructor
    /*!
      \param vec data vector holding data
      \param nbTotalDof number of total dof in the vector of data
      \param dofIn dofInterfaceBase object holding the connections between the interface dofs of the
      data vector and those of the associated to the boundary conditions
      \param type must be
      -# 0:  boundary integration done (ex. residual of a variational problem)
      -# 1:  needs boundary integration of \f$\lambda n \cdot  \mathbf{\phi}_i\f$
      \warning: implemented only for the Natural BC
      -# 2:  needs boundary integration of \f$\mathbf{\lambda} \cdot n \phi_i\f$
    */
    BCVectorInterface( const EpetraVector& vec, UInt nbTotalDof, dof_interface_type dofIn, UInt type=0 );

    //! Default Constructor (the user must call setBCVector(..))
    BCVectorInterface ();

    //! setup after default constructor

    void setup ( const EpetraVector& vec, UInt nbTotalDof, dof_interface_type dofIn, UInt type=0 );

    //! set the BC vector (after default construction)
    void setVector( EpetraVector& vec, UInt nbTotalDof, dof_interface_type dofIn, UInt type=0);

    /*!
      This method returns the value to be imposed in the component iComp of the dof iDof.

      \param iDof the number of the Dofin
      \param iComp the number of the component
    */
    Real operator() ( const ID& iDof, const ID& iComp ) const;

    //! This method returns the value of the mixte coefficient to be imposed in the component iComp of the dof iDof
    Real MixteVec( const ID& iDof, const ID& iComp ) const;

    //! This method returns the value of the beta coefficient to be imposed in the component iComp of the dof iDof
    Real BetaVec( const ID& iDof, const ID& iComp ) const;

    //! This method returns the value of the gamma coefficient to be imposed in the component iComp of the dof iDof
    Real GammaVec( const ID& iDof, const ID& iComp ) const;

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
    dof_interface_type _M_dofIn;

};
typedef LifeV::singleton< LifeV::factoryClone< BCVectorBase > > FactoryCloneBCVector;
}
#endif
