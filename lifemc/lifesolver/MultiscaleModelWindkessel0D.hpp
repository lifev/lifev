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
 *  @brief File containing the Multiscale Windkessel 0D
 *
 *  @date 08-02-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @author Mahmoud Jafargholi <mahmoud.jafargholi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleModelWindkessel0D_H
#define MultiscaleModelWindkessel0D_H 1

// Mathcard includes
#include <lifemc/lifesolver/MultiscaleModel.hpp>
#include <lifemc/lifesolver/MultiscaleInterfaceFluid.hpp>

#include <lifemc/lifesolver/BCInterface0D.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelWindkessel0D - Multiscale model for Windkessel 0D terminals
/*!
 *  @author Cristiano Malossi, Mahmoud Jafargholi
 *
 *  The MultiscaleModelWindkessel0D class is an implementation of the multiscaleModel_Type
 *  for 1D Fluid problem.
 */
class MultiscaleModelWindkessel0D: public virtual multiscaleModel_Type,
                                   public virtual MultiscaleInterfaceFluid
{
public:

    //! @name Type definitions
    //@{

    typedef ZeroDimensionalBCHandler                                        bc_Type;
    typedef boost::shared_ptr< bc_Type >                                    bcPtr_Type;

    typedef BCInterface0D< bc_Type, MultiscaleModelWindkessel0D >           bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >                           bcInterfacePtr_Type;

    typedef bc_Type::bcSide_Type                                            bcSide_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelWindkessel0D();

    //! Destructor
    virtual ~MultiscaleModelWindkessel0D() {}

    //@}


    //! @name MultiscaleModel Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //@}


    //! @name MultiscaleInterfaceFluid Methods
    //@{

    //! Impose the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    void imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function ) const;

    //! Impose the integral of the normal stress on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    void imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function ) const;

    //! Get the flux on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flux value
     */
    Real boundaryFlowRate( const bcFlag_Type& /*flag*/ ) const { return M_flowRate; }

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param stressType Type of approximation for the stress
     * @return stress value
     */
    Real boundaryStress( const bcFlag_Type& flag ) const;

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @param stressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem );

    //@}


    //! @name Get Methods
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface() { return *M_bc; }

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real boundaryDensity( const bcFlag_Type& /*flag*/) const { return M_globalData->fluidDensity(); }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real boundaryViscosity( const bcFlag_Type& /*flag*/) const { return M_globalData->fluidViscosity(); }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real boundaryPressure( const bcFlag_Type& /*flag*/ ) const { return M_pressure; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelWindkessel0D( const MultiscaleModelWindkessel0D& model );

    MultiscaleModelWindkessel0D& operator=( const MultiscaleModelWindkessel0D& model );

    //@}


    //! @name Private Methods
    //@{

    //! Setup the global data of the model.
    /*!
     * In particular, it replaces the default local values with the ones in the global container.
     * If a value is already specified in the data file, do not perform the replacement.
     *
     * @param fileName File name of the specific model.
     */
    void setupGlobalData();

    void setupExporterImporter();

    //! Initialize the solution.
    void initializeSolution();

    //! Solving for the pressure.
    /*!
     * @return the computed pressure
     */
    Real solveForPressure();    // Solve the equations for Pressure
    Real solveForFlowRate();    // Solve the equations for FlowRate
    Real tangentSolveForFlowRate(); //solve dQ/dP
    Real tangentSolveForPressure(); //solve dP/dQ


    //! Convert the flag from a bcFlag type to a bcSide type
    /*!
     * @param flag boundary condition flag
     * @return boundary condition side.
     */
    bcSide_Type flagConverter( const bcFlag_Type& flag ) const { return (flag == 0) ? OneDimensional::left : OneDimensional::right; }

    void solveLinearModel( bool& solveLinearSystem );

    //@}

    std::ofstream          M_outputFile;

    bcInterfacePtr_Type    M_bc;

    Real                   M_pressure_tn;      // Pressure (P1) @ t=t(n)
    Real                   M_flowRate_tn;      // flowRate (Q1) @ t=t(n)

    Real                   M_pressure;         // Pressure (P2) @ t=t(n+1)
    Real                   M_flowRate;         // flowRate (Q2) @ t=t(n+1)

    Real                   M_tangentPressure;  // Tangent pressure
    Real                   M_tangentFlowRate;  // Tangent flowRate

    Real                   M_resistance1;      // Resistance 1 (R1)
    Real                   M_resistance2;      // Resistance 2 (R2)
    Real                   M_capacitance;      // capacitance  (C)

    Real                   M_venousPressure;   // Back Pressure (Pv)
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelWindkessel0D()
{
    return new MultiscaleModelWindkessel0D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelWindkessel0D_H */
