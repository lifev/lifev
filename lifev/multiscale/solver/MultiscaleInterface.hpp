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
 *  @brief File containing the Multiscale Interface for Fluid problems
 *
 *  @date 31-03-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleInterface_H
#define MultiscaleInterface_H 1

#include <lifev/multiscale/solver/MultiscaleDefinitions.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleInterface - The multiscale interface for fluid problems
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleInterface class provides a general and abstract interface for the
 *  coupling of fluid problems.
 */
class MultiscaleInterface
{
public:

    //! @name Type definitions
    //@{

    typedef boost::function< Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! The main constructor.
    explicit MultiscaleInterface() {}

    //! Destructor
    virtual ~MultiscaleInterface() {}

    //@}


    //! @name Multiscale Interface Virtual Methods
    //@{

    //! Impose the flow rate on a specific interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    virtual void imposeBoundaryFlowRate( const multiscaleID_Type& boundaryID, const function_Type& function ) = 0;

    //! Impose the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    virtual void imposeBoundaryMeanNormalStress( const multiscaleID_Type& boundaryID, const function_Type& function ) = 0;

    //! Impose the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    virtual void imposeBoundaryMeanTotalNormalStress( const multiscaleID_Type& boundaryID, const function_Type& function ) = 0;

    //! Impose the area on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    virtual void imposeBoundaryArea( const multiscaleID_Type& boundaryID, const function_Type& function ) = 0;

    //! Get the flow rate on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return flow rate value
     */
    virtual Real boundaryFlowRate( const multiscaleID_Type& boundaryID ) const = 0;

    //! Get the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean normal stress value
     */
    virtual Real boundaryMeanNormalStress( const multiscaleID_Type& boundaryID ) const = 0;

    //! Get the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean total normal stress value
     */
    virtual Real boundaryMeanTotalNormalStress( const multiscaleID_Type& boundaryID ) const = 0;

    //! Get the area on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return area value
     */
    virtual Real boundaryArea( const multiscaleID_Type& boundaryID ) const = 0;

    //! Get the variation of the flow rate (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    virtual Real boundaryDeltaFlowRate( const multiscaleID_Type& boundaryID, bool& solveLinearSystem ) = 0;

    //! Get the variation of the integral of the mean normal stress (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean normal stress
     */
    virtual Real boundaryDeltaMeanNormalStress( const multiscaleID_Type& boundaryID, bool& solveLinearSystem ) = 0;

    //! Get the variation of the integral of the mean total normal stress (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean total normal stress
     */
    virtual Real boundaryDeltaMeanTotalNormalStress( const multiscaleID_Type& boundaryID, bool& solveLinearSystem ) = 0;

    //! Get the variation of the integral of the area (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the area
     */
    virtual Real boundaryDeltaArea( const multiscaleID_Type& boundaryID, bool& solveLinearSystem ) = 0;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleInterface( const MultiscaleInterface& interface );

    MultiscaleInterface& operator=( const MultiscaleInterface& interface );

    //@}
};

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleInterface_H */
