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
 *  @brief File containing the MultiScale Model MultiScale
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_Model_MultiScale_H
#define MS_Model_MultiScale_H 1

#include <lifemc/lifesolver/MS_PhysicalModel.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>
#include <lifemc/lifesolver/MS_Model_1D.hpp>
#include <lifemc/lifesolver/MS_Model_FSI3D.hpp>

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Coupling_BoundaryCondition.hpp>
#include <lifemc/lifesolver/MS_Coupling_Stress.hpp>
#include <lifemc/lifesolver/MS_Coupling_FluxStress.hpp>

namespace LifeV
{

//! MS_Model_MultiScale - MultiScale model
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Model_MultiScale class is an implementation of the MS_PhysicalModel
 *  for a general multiscale problem.
 */
class MS_Model_MultiScale: public virtual MS_PhysicalModel
{
public:

    //! @name Constructors & Destructor
    //@{

    typedef MS_PhysicalModel                                            super;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Model_MultiScale();

    //! Destructor
    ~MS_Model_MultiScale();

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Load data and create the models and the couplings
    /*!
     * @param FileName Name of data file.
     */
    void SetupData( const std::string& FileName );

    //! Setup the model
    void SetupModel();

    //! Build the system matrix and vectors
    void BuildSystem();

    //! Update the system matrix and vectors
    void UpdateSystem();

    //! Solve the problem
    void SolveSystem();

    //! save the solution
    void SaveSolution();

    //! Display some information about the multiscale problem (call after SetupProblem)
    void ShowMe();

    //@}


    //! @name Methods
    //@{

    //! Build the global map for the coupling vectors
    /*!
     * @param couplingMap Global coupling map
     */
    void CreateCouplingMap( EpetraMap& couplingMap );

    //! Initialize coupling variables for the first time step
    void InitializeCouplingVariables();

    //! Extrapolate coupling variables for the next time step
    void ExtrapolateCouplingVariables();

    //! Import the values of the coupling variables
    void ImportCouplingVariables( const MS_Vector_Type& CouplingVariables );

    //! Export the values of the coupling variables
    void ExportCouplingVariables( MS_Vector_Type& CouplingVariables );

    //! Export the values of the coupling residuals
    void ExportCouplingResiduals( MS_Vector_Type& CouplingResiduals );

    //! Export the Jacobian matrix
    /*!
     * @param Jacobian Matrix
     */
    void ExportJacobian( MS_Matrix_Type& Jacobian );

    //@}


    //! @name Get Methods
    //@{

    //! Get the number of the coupling variables
    /*!
     * @return number of the coupling variables
     */
    UInt GetCouplingVariablesNumber();

    //@}

private:

    // Models & Couplings
    MS_ModelsVector_Type        M_modelsList;
    MS_CouplingsVector_Type     M_couplingsList;
};

//! Factory create function
inline MS_PhysicalModel* MS_createMultiScale()
{
    return new MS_Model_MultiScale();
}

} // Namespace LifeV

#endif /* MS_Model_MultiScale_H */
