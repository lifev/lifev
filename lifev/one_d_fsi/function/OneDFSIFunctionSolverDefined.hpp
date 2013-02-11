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
 *  @brief File containing some functions for the boundary conditions of 1D models.
 *
 *  @version 1.0
 *  @date 01-08-2006
 *  @author Lucia Mirabella  <lucia.mirabella@gmail.com>
 *  @author Tiziano Passerini <tiziano.passerini@gmail.com>
 *
 *  @version 2.0
 *  @date 20-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef OneDFSIBCFunctionDefault_H
#define OneDFSIBCFunctionDefault_H

#include <lifev/one_d_fsi/function/OneDFSIFunction.hpp>
#include <lifev/one_d_fsi/solver/OneDFSIData.hpp>
#include <lifev/one_d_fsi/solver/OneDFSIFlux.hpp>
#include <lifev/one_d_fsi/solver/OneDFSISource.hpp>

namespace LifeV
{

//! OneDFSIModelBCFunctionDefault - Base class for deriving specific 1D boundary functions
/*!
 *  @author Lucia Mirabella, Tiziano Passerini, Cristiano Malossi
 *  @see Equations and networks of 1-D models \cite FormaggiaLamponi2003
 *  @see Geometrical multiscale coupling of 1-D models \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D
 *
 *  This class provide a general interface for implementing some specific boundary conditions
 *  for the 1D segment.
 */
class OneDFSIFunctionSolverDefined
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDFSIFunction                         bcFunction_Type;
    typedef boost::shared_ptr<bcFunction_Type>      bcFunctionPtr_Type;

    typedef OneDFSIFlux                             flux_Type;
    typedef boost::shared_ptr< flux_Type >          fluxPtr_Type;

    typedef OneDFSISource                           source_Type;
    typedef boost::shared_ptr< source_Type >        sourcePtr_Type;

    typedef OneDFSIData                             data_Type;
    typedef data_Type::mesh_Type                    mesh_Type;

    typedef data_Type::container2D_Type             container2D_Type;

    typedef SolverAmesos                            linearSolver_Type;
    typedef linearSolver_Type::vector_type          vector_Type;
    typedef boost::shared_ptr< vector_Type >        vectorPtr_Type;
    typedef boost::array< vectorPtr_Type, 2 >       vectorPtrContainer_Type;

    typedef linearSolver_Type::matrix_type          matrix_Type;

    typedef std::map< std::string, vectorPtr_Type > solution_Type;
    typedef boost::shared_ptr< solution_Type >      solutionPtr_Type;

    typedef OneDFSI::bcLine_Type                    bcLine_Type;
    typedef OneDFSI::bcSide_Type                    bcSide_Type;
    typedef OneDFSI::bcType_Type                    bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     */
    explicit OneDFSIFunctionSolverDefined ( const bcSide_Type& bcSide, const bcType_Type& bcType );

    //! Copy constructor
    /*!
     * @param bcFunctionDefault OneDFSIFunctionSolverDefined
     */
    explicit OneDFSIFunctionSolverDefined ( const OneDFSIFunctionSolverDefined& bcFunctionDefault );

    //! Destructor
    virtual ~OneDFSIFunctionSolverDefined() {}

    //@}


    //! @name Methods
    //@{

    //! Operator()
    /*!
     *  Evaluate the function.
     *
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @return the value of the function.
     */
    virtual Real operator() ( const Real& time, const Real& timeStep );

    //@}


    //! @name Set Methods
    //@{

    //! Set the flux and the source classes for the problem
    /*!
     *  @param fluxPtr pointer to the flux term of the problem.
     *  @param sourcePtr pointer to the source term of the problem.
     */
    void setFluxSource ( const fluxPtr_Type& fluxPtr, const sourcePtr_Type& sourcePtr );

    //! Set the solution of the problem
    /*!
     *  @param solutionPtr pointer to the solution of the problem.
     */
    void setSolution ( const solutionPtr_Type& solutionPtr )
    {
        M_solutionPtr = solutionPtr;
    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Automatically identify the boundary node.
    virtual void setupNode();

    //@}

    fluxPtr_Type                             M_fluxPtr;
    sourcePtr_Type                           M_sourcePtr;
    solutionPtr_Type                         M_solutionPtr;

    UInt                                     M_bcNode;
    bcSide_Type                              M_bcSide;
    bcType_Type                              M_bcType;
};



//! OneDFSIFunctionSolverDefinedRiemann - Class which implements Riemann boundary conditions for the 1D segment
/*!
 *  @author Lucia Mirabella, Tiziano Passerini
 *
 *  \cond \TODO Add the equation and some descriptions \endcond
 */
class OneDFSIFunctionSolverDefinedRiemann : public OneDFSIFunctionSolverDefined
{
public:

    //! @name Type definitions and Enumerators
    //@{

    typedef OneDFSIFunctionSolverDefined                    super;
    typedef super::container2D_Type                         container2D_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     */
    explicit OneDFSIFunctionSolverDefinedRiemann ( const bcSide_Type& bcSide, const bcType_Type& bcType );

    //! Copy constructor
    /*!
     * @param bcFunctionRiemann OneDFSIFunctionSolverDefinedRiemann
     */
    explicit OneDFSIFunctionSolverDefinedRiemann ( const OneDFSIFunctionSolverDefinedRiemann& bcFunctionRiemann );

    //! Destructor
    virtual ~OneDFSIFunctionSolverDefinedRiemann() {}

    //@}


    //! @name Methods
    //@{

    //! Operator()
    /*!
     *  Evaluate the function.
     *
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @return the value of the function.
     */
    virtual Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Update the boundary condition variables.
    void updateBCVariables();

    //@}

    //! Value of U at the boundary
    container2D_Type                         M_bcU;

    //! Value of W at the boundary
    container2D_Type                         M_bcW;
};


//! OneDFSIFunctionSolverDefinedCompatibility - Class which implements Compatibility boundary conditions for the 1D segment
/*!
 *  @author Lucia Mirabella, Tiziano Passerini, Cristiano Malossi
 *
 *  The compatibility equations are derived using the pseudo-characteristic teory:
 *
 *  \f[
 *  \mathbf L(\mathbf U^{n+1}-\mathbf U^{n}_\star + \mathbf U^0 - \mathbf U^0_\star) =
 *  \Delta t \left( \mathbf \Lambda \displaystyle\frac{\partial \mathbf L}{\partial z}(\mathbf U^n_\star -
 *  \mathbf U^0_\star) - \mathbf L \mathbf B(\mathbf U^n_\star) + \mathbf L \mathbf B(\mathbf U^0_\star) \right).
 *  \f]
 */
class OneDFSIFunctionSolverDefinedCompatibility : public OneDFSIFunctionSolverDefinedRiemann
{
public:

    //! @name Type definitions
    //@{

    typedef OneDFSIFunctionSolverDefinedRiemann              super;

    typedef super::fluxPtr_Type                              fluxPtr_Type;
    typedef super::sourcePtr_Type                            sourcePtr_Type;
    typedef super::solutionPtr_Type                          solutionPtr_Type;

    typedef super::mesh_Type                                 mesh_Type;
    typedef super::container2D_Type                          container2D_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     */
    explicit OneDFSIFunctionSolverDefinedCompatibility ( const bcSide_Type& bcSide,  const bcType_Type& bcType );

    //! Copy constructor
    /*!
     * @param bcFunctionCompatibility OneDFSIFunctionSolverDefinedCompatibility
     */
    explicit OneDFSIFunctionSolverDefinedCompatibility ( const OneDFSIFunctionSolverDefinedCompatibility& bcFunctionCompatibility );

    //! Destructor
    virtual ~OneDFSIFunctionSolverDefinedCompatibility() {}

    //@}


    //! @name Methods
    //@{

    //! Operator()
    /*!
     *  Evaluate the function.
     *
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @return the value of the function.
     */
    virtual Real operator() ( const Real& /*time*/, const Real& timeStep )
    {
        return computeRHS ( timeStep );
    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Automatically identify the boundary node.
    void setupNode();

    //! Compute the rhs
    /*!
     *  @param timeStep the time step.
     *  @return rhs of the problem.
     */
    Real computeRHS ( const Real& timeStep );

    //! Compute the current eigenvalues and eigenvectors
    void computeEigenValuesVectors();

    //! Compute the rhs
    /*!
     *  @param eigenvalue eigenvalue
     *  @param eigenvector eigenvector
     *  @param deltaEigenvector derivative of the eigenvector
     *  @param timeStep the time step.
     *  @return rhs of the problem
     */
    Real evaluateRHS ( const Real& eigenvalue, const container2D_Type& eigenvector,
                       const container2D_Type& deltaEigenvector, const Real& timeStep );

    //! Compute the current CFL
    /*!
     *  @param eigenvalue eigenvalue
     *  @param timeStep the time step.
     *  @return CFL
     */
    Real computeCFL ( const Real& eigenvalue, const Real& timeStep ) const;

    //! Scalar product between 2 2D vectors
    /*!
     *  @pararm vector1 first vector
     *  @pararm vector2 second vector
     *  @return scalar product
     */
    Real scalarProduct ( const container2D_Type& vector1, const container2D_Type& vector2 )
    {
        return vector1[0] * vector2[0] + vector1[1] * vector2[1];
    }

    //@}

    //! ID of the boundary edge
    UInt                                               M_bcElement;
    //! Dof of the internal node adjacent to the boundary
    UInt                                               M_bcInternalNode;

    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    container2D_Type                                   M_eigenvalues;
    container2D_Type                                   M_deltaEigenvalues;

    //! Left eigen vectors for the two eigen values
    container2D_Type                                   M_leftEigenvector1;
    container2D_Type                                   M_leftEigenvector2;
    container2D_Type                                   M_deltaLeftEigenvector1;
    container2D_Type                                   M_deltaLeftEigenvector2;
};


//! OneDFSIFunctionSolverDefinedAbsorbing - Class which implements absorbing boundary conditions for the 1D segment
/*!
 *  @author Maria Rita de Luca
 *
 *  \cond \TODO Add the equation and some descriptions \endcond
 */
class OneDFSIFunctionSolverDefinedAbsorbing : public OneDFSIFunctionSolverDefinedCompatibility
{
public:

    //! @name Type definitions
    //@{

    typedef OneDFSIFunctionSolverDefinedCompatibility        super;

    typedef super::fluxPtr_Type                              fluxPtr_Type;
    typedef super::sourcePtr_Type                            sourcePtr_Type;
    typedef super::solutionPtr_Type                          solutionPtr_Type;

    typedef super::mesh_Type                                 mesh_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     */
    explicit OneDFSIFunctionSolverDefinedAbsorbing ( const bcSide_Type& bcSide, const bcType_Type& bcType ) : super ( bcSide, bcType ) {}

    //! Copy constructor
    /*!
     * @param bcFunctionAbsorbing OneDFSIFunctionSolverDefinedAbsorbing
     */
    explicit OneDFSIFunctionSolverDefinedAbsorbing ( const OneDFSIFunctionSolverDefinedAbsorbing& bcFunctionAbsorbing ) : super ( bcFunctionAbsorbing ) {}

    //! Destructor
    virtual ~OneDFSIFunctionSolverDefinedAbsorbing() {}

    //@}


    //! @name Methods
    //@{

    //! Operator()
    /*!
     *  Evaluate the function.
     *
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @return the value of the function.
     */
    Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:

    //! Set the value of the resistance
    /*!
     * For absorbing BC do nothing.
     * @param resistance value of the resistance
     */
    virtual void resistance ( Real& /*resistance*/ )
    {
        /*Do nothing => absorbing!*/
    }

    //! Venous pressure
    /*!
     * For absorbing BC the venous pressure is equal to the external pressure.
     * @return venous pressure.
     */
    virtual Real venousPressure()
    {
        return M_fluxPtr->physics()->externalPressure();
    }

};


//! OneDFSIFunctionSolverDefinedResistance - Class which implements resistance boundary conditions for the 1D segment
/*!
 *  @author Lucia Mirabella, Tiziano Passerini
 *
 *  \cond \TODO Add the equation and some descriptions \endcond
 */
class OneDFSIFunctionSolverDefinedResistance : public OneDFSIFunctionSolverDefinedAbsorbing
{
public:

    //! @name Type definitions
    //@{

    typedef OneDFSIFunctionSolverDefinedAbsorbing            super;

    typedef super::fluxPtr_Type                              fluxPtr_Type;
    typedef super::sourcePtr_Type                            sourcePtr_Type;
    typedef super::solutionPtr_Type                          solutionPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     *  @param resistance the terminal resistance.
     */
    explicit OneDFSIFunctionSolverDefinedResistance ( const bcSide_Type& bcSide,  const bcType_Type& bcType, const Real& resistance );

    //! Copy constructor
    /*!
     * @param bcFunctionResistance OneDFSIFunctionSolverDefinedResistance
     */
    explicit OneDFSIFunctionSolverDefinedResistance ( const OneDFSIFunctionSolverDefinedResistance& bcFunctionResistance );

    //! Destructor
    virtual ~OneDFSIFunctionSolverDefinedResistance() {}

    //@}

protected:

    //! Set the value of the resistance
    /*!
     * @param resistance value of the resistance
     */
    void resistance ( Real& resistance )
    {
        resistance = M_resistance;
    }

    //! Venous pressure
    /*!
     * @return venous pressure.
     */
    Real venousPressure()
    {
        return M_fluxPtr->physics()->venousPressure();
    }

    Real M_resistance;
};



//! OneDFSIFunctionSolverDefinedWindkessel3 - Class which implements windkessel RCR boundary conditions for the 1D segment
/*!
 *
 *  \cond \TODO Description should be reordered using latex etc... \endcond
 *  \cond \TODO This method has not been tested yet \endcond
 *
 *   *   Q -> ---R1-------R2---
 *         ^      |       ^
 *         P      C       Pv
 *         ^      |       ^
 *
 *  The holding equation is:
 *
 *  P + C R2 dP/dt = (R1 + R2 ) * Q + R1 R2 C dQ/dt + Pv
 *
 *  where the "venous" pressure Pv is taken constant and equal to 5mmHg (6666 dyn/cm^2).
 *
 *  You can solve this ODE in analytical fashion and obtain:
 *
 *  P(t) = P(0) + [ \int_0^t ( pv + (R1+R2) Q(s) + R1*R2*C dQ(s)/ds ) exp( s / (R2*C) ) ds ] * exp( - t / (R2*C) )
 *
 *  then you just need to exploit numerical integration rules.
 *  Here a simple trapezoidal rule is used, with a first order approximation of derivative
 *
 *  dQ(s)/ds = Q(t(n+1)) - Q(t(n)) * (1/dt)
 *
 *  @author Lucia Mirabella, Tiziano Passerini
 */
class OneDFSIFunctionSolverDefinedWindkessel3 : public OneDFSIFunctionSolverDefinedCompatibility
{
public:

    //! @name Type definitions
    //@{

    typedef OneDFSIFunctionSolverDefinedCompatibility    super;

    typedef super::fluxPtr_Type                          fluxPtr_Type;
    typedef super::sourcePtr_Type                        sourcePtr_Type;
    typedef super::solutionPtr_Type                      solutionPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     *  @param bcLine the line of the boundary condition (first or second).
     *  @param bcType the type of the boundary condition (\f$Q\f$, \f$A\f$, \f$P\f$, \f$S\f$, \f$W_1\f$, \f$W_2\f$).
     *  @param resistance1 the first terminal resistance.
     *  @param resistance2 the second terminal resistance.
     *  @param compliance the compliance.
     *  @param absorbing is an absorbing boundary condition
     *  @param venousPressure the venous pressure
     */
    explicit OneDFSIFunctionSolverDefinedWindkessel3 ( const bcSide_Type& bcSide, const bcType_Type& bcType,
                                                       const Real& resistance1, const Real& resistance2,
                                                       const Real& compliance,
                                                       const bool& absorbing = false,
                                                       const Real& venousPressure = 6666. );

    //! Copy constructor
    /*!
     * @param bcFunctionWindkessel3 OneDFSIFunctionSolverDefinedWindkessel3
     */
    explicit OneDFSIFunctionSolverDefinedWindkessel3 ( const OneDFSIFunctionSolverDefinedWindkessel3& bcFunctionWindkessel3 );

    //! Destructor
    virtual ~OneDFSIFunctionSolverDefinedWindkessel3() {}

    //@}


    //! @name Methods
    //@{

    //! Operator()
    /*!
     *  Evaluate the function.
     *
     *  @param time the current time.
     *  @param timeStep the time step.
     *  @return the value of the function.
     */
    Real operator() ( const Real& time, const Real& timeStep );

    //@}

protected:

    Real M_resistance1;
    Real M_resistance2;
    Real M_compliance;
    bool M_absorbing;
    Real M_venousPressure;

    // Initial value of the pressure
    Real M_P0;

    // Q at previous time step
    Real M_Q_tn;

    // dQdt at the previous time step
    Real M_dQdt_tn;

    // Integral of the pressure at the previous time step
    Real M_integral_tn;
};

}

#endif // OneDFSIBCFunctionDefault_H
