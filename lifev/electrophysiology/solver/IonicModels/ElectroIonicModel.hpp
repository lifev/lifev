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
  @file ElectroIonicModel
  @brief Base class for ionic models

  This class is a  abstract base class for all ionic models.

  Assume that the ionic model is written in the form

  \f[ \dfrac{\partial \mathbf{v} }{\partial t} = f(\mathbf{v), t\f]

  If you wish to implement a new ionic model you should create a new
  class which inherits from this class and implement the abstract methods:

  Given \f[ \mathbf{v} \f] compute \f[ f(\mathbf{v}, t) \f]
  void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);
  In principle this method should call the following two methods:

  Given \f[ \mathbf{v} \f] compute \f[ f(\mathbf{v}, t) \f]  only for the first
  equation, that is only for the  transmembrane potential
  Real computeLocalPotentialRhs ( const std::vector<Real>& v );

  Given \f[ \mathbf{v} \f] compute \f[ f(\mathbf{v}, t) \f] for all the variables
  excpet the first one, that is all variables except the transmembrane potential
  void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

  Print out informations about the ionic model
  void showMe();

  Some ionic models may be solved with the Rush Larsen scheme.
  The implementation of this method is not mandatory.
  For biophysically detailed models, it may be convenient,
  to have it.
  void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt ) {}

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 02-2012
 */


#ifndef _ELECTROIONICMODEL_H_
#define _ELECTROIONICMODEL_H_

#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/electrophysiology/stimulus/ElectroStimulus.hpp>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

namespace LifeV
{
class ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{

    typedef VectorEpetra                            vector_Type;

    typedef boost::shared_ptr<VectorEpetra>         vectorPtr_Type;

    typedef boost::shared_ptr<VectorElemental>      elvecPtr_Type;

    typedef RegionMesh<LinearTetra>                 mesh_Type;

    typedef MatrixEpetra<Real>                      matrix_Type;

    typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;

    typedef FESpace<mesh_Type, MapEpetra>           feSpace_Type;

    typedef boost::shared_ptr<feSpace_Type>         feSpacePtr_Type;

    typedef boost::function < Real (const Real& t,
                                    const Real& x,
                                    const Real& y,
                                    const Real& z,
                                    const ID& i) >   function_Type;
    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
     */
    ElectroIonicModel();

    //! Constructor
    /*!
     * @param n number of equations in the ionic model
     */
    ElectroIonicModel ( int n );

    //! Constructor
    /*!
     *  If the number of gating variables is unknown use the method ElectroIonicModel ( int n );
     */
    /*!
     * @param n number of equations in the ionic model
     * @param g number of gating variables in the ionic model
     */
    ElectroIonicModel ( int n, int g );

    //! Copy Constructor
    /*!
     * @param Ionic an ionic model
     */
    ElectroIonicModel ( const ElectroIonicModel& Ionic );

    //! Destructor
    virtual ~ElectroIonicModel() {};

    //@}

    //! @name Methods
    //@{

    //! returns the number of equations of the ionic model
    /*!
     * @param
     */
    inline const short int Size() const
    {
        return M_numberOfEquations;
    }

    //! returns the number of gating variables in the ionic model
    /*!
     * @param
     */
    inline const short int numberOfGatingVariables() const
    {
        return M_numberOfGatingVariables;
    }

    //! returns the value of the membrane capacitance in the model
    /*!
     * @param
     */
    inline const Real membraneCapacitance() const
    {
        return M_membraneCapacitance;
    }

    //! returns the value of the applied current in the model/point
    /*!
     * @param
     */
    inline const Real appliedCurrent() const
    {
        return M_appliedCurrent;
    }

    //! returns the pointer to the applied current FE vector in the 3D case
    /*!
     * @param
     */
    inline vectorPtr_Type appliedCurrentPtr()
    {
        return M_appliedCurrentPtr;
    }

    //! returns the vector with the resting values of the variables in the ionic model
    /*!
     * @param
     */
    inline const std::vector<Real> restingConditions() const
    {
        return M_restingConditions;
    }

    //! returns the function describing the pacing protocol for the ionic model
    /*!
     * @param
     */
    inline const function_Type pacaingProtocol() const
    {
        return M_pacingProtocol;
    }

    //! set the membrane capacitance in the ionic model
    /*!
     * @param p membrane capacitance
     */
    inline void setMembraneCapacitance ( const Real p )
    {
        M_membraneCapacitance = p;
    }

    //! set the applied current in the ionic model/point
    /*!
     * @param p applied current magnitude
     */
    inline void setAppliedCurrent (const Real p)
    {
        M_appliedCurrent = p;
    }

    //! set the pointer to the applied current in the 3D ionic model
    /*!
     * @param p applied current magnitude
     */
    inline void setAppliedCurrentPtr (const vectorPtr_Type p)
    {
        this->M_appliedCurrentPtr = p;
    }

    //! set the pointer to the applied current in the 3D ionic model
    /*!
     * @param p applied current magnitude
     */
    inline void setAppliedCurrent (const vector_Type& p)
    {
        M_appliedCurrentPtr.reset ( new vector_Type ( p ) );
    }

    //! Interpolate the function f on the FE space feSpacePtr at time time
    /*!
     *  This method is a wrapper of the interpolation method from the FESpace class
     */
    /*!
     * @param f boost function f(t,x,y,x,ID)
     * @param feSpacePtr pointer to the finite element space on which we interpolate
     * @param time time at which we evaluate the boost function f
     */
    inline void setAppliedCurrentFromFunction (function_Type& f, feSpacePtr_Type feSpacePtr,
                                               Real time = 0.0)
    {

        feSpacePtr -> interpolate (
            static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (f),
            *M_appliedCurrentPtr, time);
    }

    //! Interpolate the function of the electro stimulus
    /*!
     *  This method is a wrapper of the interpolation method from the FESpace class
     */
    /*!
     * @param stimulus pacing protocol defined as ElectroStimulus boost function f(t,x,y,x,ID)
     * @param feSpacePtr pointer to the finite element space on which we interpolate
     * @param time time at which we evaluate the boost function f
     */
    inline void setAppliedCurrentFromElectroStimulus ( ElectroStimulus& stimulus, feSpacePtr_Type feSpacePtr, Real time = 0.0)
    {

        // boost::ref() is needed here because otherwise a copy of the base object is reinstantiated
        function_Type f = boost::bind (&ElectroStimulus::appliedCurrent, boost::ref (stimulus), _1, _2, _3, _4, _5 );

        feSpacePtr -> interpolate (
            static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (f),
            *M_appliedCurrentPtr, time);
    }

    //! Set the pacing protocol as boost function
    /*!
     * @param pacingProtocol  boost function defining the pacing protocol
     */
    inline void setPacingProtocol ( function_Type pacingProtocol )
    {
        M_pacingProtocol = pacingProtocol;
    }

    //! Set component of the resting conditions
    /*!
     * @param value value of the resting condition
     * @param j number of the component to be changed
     */
    inline void setRestingCondtions (Real value, int j)
    {
        M_restingConditions[j] = value;
    }

    //! Set resting conditions
    /*!
     * @param restingCondition values of the resting condition
     */
    inline void setRestingCondtions (std::vector<Real>& restingConditions)
    {
        M_restingConditions = restingConditions;
    }

    //! Simple wrapper to add the applied current
    /*!
     * @param rhs right hand side of the voltage equation
     */
    inline void addAppliedCurrent (Real& rhs)
    {
        rhs += M_appliedCurrent;
    }

    //! Simple wrapper to add the applied current
    /*!
     * @param rhs right hand side of the ionic model (in a point)
     */
    inline void addAppliedCurrent (std::vector<Real>& rhs)
    {
        rhs[0] += M_appliedCurrent;
    }

    //@}

    //! This methods computes the Jacobian numerically
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param h differentiation step
     */
    virtual matrix_Type getJac (const vector_Type& v, Real h = 1.0e-8);

    //! This methods computes the Jacobian numerically
    /*!
     * @param v vector of  the  state variables
     * @param h differentiation step
     */
    virtual std::vector< std::vector<Real> > getJac (const std::vector<Real>& v, Real h = 1.0e-8);

    //! This methods computes the right hand side of the gating variables in the 3D case
    /*!
     *  In 3D the gating variables are still treated nodewise (that is, as a local 0D system)
     *  It computes the right hand side of all  state variables except for the voltage equation.
     */
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param rhs vector of pointers to the right hand side vectors of each variable
     */
    virtual void computeGatingRhs ( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );

    //! This methods computes the right hand side of the state variables that are not gating variables in the 3D case
    /*!
     *  This method is to be used together with Rush-Larsen time advancing scheme.
     *  In 3D the non gating variables are still treated nodewise (that is, as a local 0D system).
     *  It computes the right hand side of all non gating  state variables except for the voltage equation.
     */
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param rhs vector of pointers to the right hand side vectors of each variable
     */
    virtual void computeNonGatingRhs ( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );

    //! Compute the new value of the gating variables in 3D with the Rush Larsen method specified in the 0D version of the ionic model
    /*!
     *  This method wraps the 0D model to be used in 3D
     */
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param dt time step
     */
    virtual void computeGatingVariablesWithRushLarsen ( std::vector<vectorPtr_Type>& v, const Real dt );

    //! Compute the right hand side of the ionic model in 3D
    /*!
     *  This method wraps the 0D model to be used in 3D
     */
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param rhs vector of right hand side state variables
     */
    virtual void computeRhs ( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );

    //! Compute the right hand side of the voltage equation linearly interpolating the ionic currents
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param rhs vector of right hand side state variables
     * @param massMatrix mass matrix of the monodomain system (may be lumped)
     */
    virtual void computePotentialRhsICI ( const std::vector<vectorPtr_Type>& v,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          matrix_Type&                        massMatrix );

    //! Compute the right hand side of the voltage equation using SVI
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param rhs vector of right hand side state variables
     * @param uFESpace finite element space of the voltage

     */
    virtual void computePotentialRhsSVI ( const std::vector<vectorPtr_Type>& v,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          FESpace<mesh_Type, MapEpetra>&  uFESpace );

    //! Compute the right hand side of the voltage equation using SVI   specifying the quadrature rule
    /*!
     * @param v vector of pointers to the  state variables vectors
     * @param rhs vector of right hand side state variables
     * @param uFESpace finite element space of the voltage
     * @param qr quadrature rule for the integration
     */
    virtual void computePotentialRhsSVI ( const std::vector<vectorPtr_Type>& v,
                                          std::vector<vectorPtr_Type>&        rhs,
                                          FESpace<mesh_Type, MapEpetra>&  uFESpace,
                                          const QuadratureRule& qr );

    //! Initialize the ionic model with a given vector of state variable (0D version)
    /*!
     * @param v vector of state variables initial conditions
     */
    virtual void initialize ( std::vector<Real>& v );

    //initialize with given conditions
    //! Initialize the ionic model in 3D with a given vector of state variable vector pointers
    /*!
     * @param v vector of state variables initial conditions
     */
    virtual void initialize ( std::vector<vectorPtr_Type>& v );


    //! @name Overloads
    //@{

    //! Assignment operator
    /*!
     * @param Ionic ionic model
     */
    ElectroIonicModel& operator= ( const ElectroIonicModel& Ionic );

    //@}



    /////////////////////////////////////////////////////////////
    /// A new ionic model must  have the following methods    ///
    /////////////////////////////////////////////////////////////

    //! This methods contains the actual evaluation of the rhs of all state variables except for the voltage equation (0D version)
    /*!
     * @param v vector of state variables including the voltage (with n elements)
     * @param rhs vector of right hand side state variables excluding the voltage (with n-1 elements)
     */
    virtual void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs ) = 0;

    //! This methods contains the actual evaluation of the rhs of the voltage equation only (0D version)
    /*!
     * @param v vector of state variables including the voltage (with n elements)
     */
    virtual Real computeLocalPotentialRhs ( const std::vector<Real>& v ) = 0;

    //! This methods contains the actual evaluation of the rhs of all state variablesin the model (0D version)
    /*!
     *  Although this method can just wrap the  computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs ) and
     *  the computeLocalPotentialRhs ( const std::vector<Real>& v ) methods, for efficiency it may be better to duplicate the code.
     */
    /*!
     * @param v vector of state variables including the voltage (with n elements)
     * @param rhs vector of right hand side state variables including the voltage (with n elements)
     */
    virtual void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs) = 0;

    //! This methods shows the parameters of the ionic model
    /*!
     * @param
     */
    virtual void showMe() = 0;

    //! This methods contains the actual evaluation of the rhs of the voltage equation only (0D version)
    /*!
     *  Overload this method in order to solve the ionic model with the Rush-Larsen method in the monodomain solver
     */
    /*!
     * @param v vector of state variables
     * @param dt time step of the simulation
     */
    virtual void computeGatingVariablesWithRushLarsen ( std::vector<Real>& /*v*/, const Real /*dt*/ ) {}


    //! In the case this method is improperly used, it should use this default implementation
    /*!
     *  This method should be used together with the Rush-Larsen time advancing scheme
     */
    /*!
     * @param v vector of state variables including the voltage (with n elements)
     * @param rhs vector of right hand side state variables excluding the voltage and the gating variables (with n-g-1 elements)
     */
    virtual void computeNonGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs )
    {
        if (M_numberOfGatingVariables == 0 )
        {
            computeGatingRhs (v, rhs);
        }
    };

    //@}

protected:

    //Number of equations in the model
    short int  M_numberOfEquations;

    //Number of gating variables in the model
    short int  M_numberOfGatingVariables;

    //Resting conditions or default initial conditions of the ionic model
    std::vector<Real> M_restingConditions;

    //Value of the membrane capacitance in the ionic model
    Real M_membraneCapacitance;

    //Applied current in the 0D version and applied current in a point in the 3D version
    Real M_appliedCurrent;

    //Pointer to the applied current FE vector for 3D simulations
    vectorPtr_Type M_appliedCurrentPtr;

    //Function describing the pacing protocol of the model - NEEDS TO BE CONFIRMED
    function_Type M_pacingProtocol;


};


}

#endif
