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
class MultiscaleModelWindkessel0D: public virtual multiscaleModel_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelWindkessel0D();

    //! Destructor
    virtual ~MultiscaleModelWindkessel0D() {}

    //@}


    //! @name Multiscale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial system (matrix and vectors).
    void buildSystem();

    //! Update the system for (matrix and vectors).
    void updateSystem();

    //! Solve the problem.
    void solveSystem();

    //! save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //@}


    //! @name Methods
    //@{

    //! Setup the linear model
    void setupLinearModel();

    //! Update the linear system matrix and vectors
    void updateLinearModel();

    //! Solve the linear problem
    void solveLinearModel( bool& solveLinearSystem );

    //@}


    //! @name Get Methods (couplings)
    //@{

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real boundaryDensity( const bcFlag_Type& /*flag*/) const { return 0; }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real boundaryViscosity( const bcFlag_Type& /*flag*/) const { return 0; }

    //! Get the flux on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flux value
     */
    Real boundaryFlowRate( const bcFlag_Type& flag ) const { return 0; }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real boundaryPressure( const bcFlag_Type& flag ) const { return 0; }

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param stressType Type of approximation for the stress
     * @return stress value
     */
    Real boundaryStress( const bcFlag_Type& flag, const stress_Type& stressType = Pressure ) const;

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem );

    //! Get the variation of the pressure (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the pressure
     */
    Real boundaryDeltaPressure( const bcFlag_Type& flag, bool& solveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @param stressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem, const stress_Type& stressType = Pressure );

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
    void setupGlobalData( const std::string& fileName );

    //! Initialize the solution.
    void initializeSolution();

    //! Solving for the pressure.
    /*!
     * @return the computed pressure
     */
    Real solveForPressure();    // Solve the equations for Pressure
    Real solveForFlowRate();    // Solve the equations for FlowRate

    Real rightHandSideQ  (Real R1, Real R2 , Real RC, Real T, Real P1, Real P2, Real dP, Real Pv1, Real dt) //Dummy function
    {
        return ((exp( ( R1 + R2 ) * RC * (-T)) )*( RC * ( (P2-P1)*(T)/dt+ P1 ) + dP/R1 - (Pv1)*RC ));
    }

    Real rightHandSideP (Real R1,Real R2,Real RC,Real T,Real Q1,Real Q2,Real dQ,Real Pv1,Real dt) //Dummy function
    {
        return (exp( RC * (-T))*( (R1+R2) * RC *( (Q2 -Q1) *(T) /dt + Q1 ) + R1 * dQ + Pv1*RC ));
    }

    //@}

     Real M_pressure_tn;    //Pressure (P1) @ t=t(n)
     Real M_flowRate_tn;    //flowRate (Q1) @ t=t(n)

     Real M_pressure;       //Pressure (P2) @ t=t(n+1)
     Real M_flowRate;       //flowRate (Q2) @ t=t(n+1)

     Real M_backPressure;    //Back Pressure (Pv)

     Real M_dt      ;       //time step
     Real M_nIntegration;   //Number of Integration steps in each time step

     Real M_resistance1;    //Resistance 1 (R1)
     Real M_resistance2;    //Resistance 2 (R2)
     Real M_capacitance;    //capacitance  (C)

     UInt M_solveForPressure;      //1= solve for Pressure, 0=solve for flowRate;


};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelWindkessel0D()
{
    return new MultiscaleModelWindkessel0D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelWindkessel0D_H */
