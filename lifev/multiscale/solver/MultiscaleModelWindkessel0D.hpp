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
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleModelWindkessel0D_H
#define MultiscaleModelWindkessel0D_H 1

#include <lifev/bc_interface/fem/BCInterface0D.hpp>

#include <lifev/multiscale/solver/MultiscaleModel.hpp>
#include <lifev/multiscale/solver/MultiscaleInterface.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelWindkessel0D - Multiscale model for Windkessel 0D terminals
/*!
 *  @author Cristiano Malossi, Mahmoud Jafargholi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleModelWindkessel0D class is an implementation of the multiscaleModel_Type
 *  for 0D problems.
 */
class MultiscaleModelWindkessel0D: public virtual multiscaleModel_Type,
    public virtual MultiscaleInterface
{
public:

    //! @name Type definitions
    //@{

    typedef ZeroDimensionalBCHandler                                        bc_Type;
    typedef boost::shared_ptr< bc_Type >                                    bcPtr_Type;

    typedef BCInterface0D< bc_Type, MultiscaleGlobalData >                  bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >                           bcInterfacePtr_Type;

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
    void setupData ( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Update the solution.
    void updateSolution();

    //! Save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //! Return a specific scalar quantity to be used for a comparison with a reference value.
    /*!
     * This method is meant to be used for night checks.
     * @return reference quantity.
     */
    Real checkSolution() const;

    //@}


    //! @name MultiscaleInterface Methods
    //@{

    //! Impose the flow rate on a specific interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryFlowRate ( const multiscaleID_Type& boundaryID, const function_Type& function );

    //! Impose the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryMeanNormalStress ( const multiscaleID_Type& boundaryID, const function_Type& function );

    //! Impose the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * Note: mean total normal stress cannot be imposed at the interfaces of this model.
     *
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryMeanTotalNormalStress ( const multiscaleID_Type& /*boundaryID*/, const function_Type& /*function*/ )
    {
        multiscaleErrorCheck ( ModelInterface, "Invalid interface [MeanTotalNormalStress] for model type [" + enum2String ( M_type, multiscaleModelsMap ) + "]", M_comm->MyPID() == 0 );
    }

    //! Impose the area on a specific boundary interface of the model
    /*!
     * Note: area cannot be imposed at the interfaces of this model.
     *
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryArea ( const multiscaleID_Type& /*boundaryID*/, const function_Type& /*function*/ )
    {
        multiscaleErrorCheck ( ModelInterface, "Invalid interface [Area] for model type [" + enum2String ( M_type, multiscaleModelsMap ) + "]", M_comm->MyPID() == 0 );
    }

    //! Get the flow rate on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return flow rate value
     */
    Real boundaryFlowRate ( const multiscaleID_Type& boundaryID ) const
    {
        return ( boundaryFlag ( boundaryID ) == 0 ) ? M_flowRateLeft : 0;
    }

    //! Get the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean normal stress value
     */
    Real boundaryMeanNormalStress ( const multiscaleID_Type& boundaryID ) const
    {
        return -boundaryPressure ( boundaryID );
    }

    //! Get the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * Note: returns always a NaN since the mean total normal stress is not defined in the windkessel model
     *
     * @param boundaryID ID of the boundary interface
     * @return mean total normal stress value
     */
    Real boundaryMeanTotalNormalStress ( const multiscaleID_Type& /*boundaryID*/ ) const
    {
        return NaN;
    }

    //! Get the area on a specific boundary interface of the model
    /*!
     *  Note: returns always a NaN since the area is not defined in the windkessel model
     *
     * @param boundaryID ID of the boundary interface
     * @return area value
     */
    Real boundaryArea ( const multiscaleID_Type& /*boundaryID*/ ) const
    {
        return NaN;
    }

    //! Get the variation of the flow rate (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the mean normal stress (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean normal stress
     */
    Real boundaryDeltaMeanNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the mean total normal stress (on a specific boundary interface) using the linear model
    /*!
     *  Note: returns always a NaN since the mean total normal stress is not defined in the windkessel model.
     *
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean total normal stress
     */
    Real boundaryDeltaMeanTotalNormalStress ( const multiscaleID_Type& /*boundaryID*/, bool& /*solveLinearSystem*/ )
    {
        return NaN;
    }

    //! Get the variation of the integral of the area (on a specific boundary interface) using the linear model
    /*!
     *  Note: returns always a NaN since the area is not defined in the windkessel model.
     *
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the area
     */
    Real boundaryDeltaArea ( const multiscaleID_Type& /*boundaryID*/, bool& /*solveLinearSystem*/ )
    {
        return NaN;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface()
    {
        return *M_bc;
    }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param boundaryID ID of the boundary interface
     * @return pressure value
     */
    Real boundaryPressure ( const multiscaleID_Type& boundaryID ) const
    {
        return ( boundaryFlag ( boundaryID ) == 0 ) ? M_pressureLeft : M_pressureRight;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelWindkessel0D ( const MultiscaleModelWindkessel0D& model );

    MultiscaleModelWindkessel0D& operator= ( const MultiscaleModelWindkessel0D& model );

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
    void setupGlobalData ( const std::string& fileName );

    //! Initialize the solution.
    void initializeSolution();

    void setupExporterImporter();

    //! Solving for the flow rate.
    /*!
     * TODO ADD THE EQUATIONS
     * @return the computed flow rate
     */
    Real solveForFlowRate();

    //! Solving for the pressure.
    /*!
     * TODO ADD THE EQUATIONS
     * @return the computed pressure
     */
    Real solveForPressure();

    //! Solve the tangent problem
    /*!
     * @param solveLinearSystem if true the system as already been solved.
     */
    void solveLinearModel ( bool& solveLinearSystem );

    //! Solve the tangent problem for the flow rate
    /*!
     * TODO ADD THE EQUATIONS
     * @return \f$ dQ/dP \f$
     */
    Real tangentSolveForFlowRate();

    //! Solve the tangent problem for the flow rate
    /*!
     * TODO ADD THE EQUATIONS
     * @return \f$ dP/dQ \f$
     */
    Real tangentSolveForPressure();

    //@}

    std::ofstream          M_outputFile;

    bcInterfacePtr_Type    M_bc;

    Real                   M_pressureLeft_tn;  // pressure left (P1) @ t=t(n)
    Real                   M_flowRateLeft_tn;  // flowRate left (Q1) @ t=t(n)

    Real                   M_pressureLeft;     // pressure left (P2) @ t=t(n+1)
    Real                   M_flowRateLeft;     // flowRate left (Q2) @ t=t(n+1)

    Real                   M_pressureRight;    // pressure right (usually venous pressure)

    Real                   M_tangentPressureLeft;  // Tangent pressure left
    Real                   M_tangentFlowRateLeft;  // Tangent flowRate left

    Real                   M_resistance1;      // Resistance 1 (R1)
    Real                   M_resistance2;      // Resistance 2 (R2)
    Real                   M_capacitance;      // capacitance  (C)
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelWindkessel0D()
{
    return new MultiscaleModelWindkessel0D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelWindkessel0D_H */
