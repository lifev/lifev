//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Model MultiScale
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-03-2009
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

namespace LifeV {

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

    typedef MS_PhysicalModel                                            super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Model_MultiScale();

    //! Copy constructor
    /*!
     * @param multiscale MS_Model_MultiScale
     */
    MS_Model_MultiScale( const MS_Model_MultiScale& multiscale );

    //! Destructor
    ~MS_Model_MultiScale();

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param multiscale MS_Model_MultiScale
     * @return reference to a copy of the class
     */
    MS_Model_MultiScale& operator=( const MS_Model_MultiScale& multiscale );

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Load data and create the models and the couplings
    /*!
     * @param FileName Name of data file.
     */
    void SetupData( const std::string& FileName );

    //! Setup the global data of the model.
    /*!
     * @param PhysicalData Global data container.
     */
    void SetupGlobalData( const boost::shared_ptr< MS_PhysicalData >& PhysicalData );

    //! Setup the model
    void SetupModel();

    //! Build the system matrix and vectors
    void BuildSystem();

    //! Update the system matrix and vectors
    void UpdateSystem();

    //! Solve the problem
    void SolveSystem();

    //! Save the solution
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

    //! Setup parameters for the implicit coupling
    void InitializeCouplingVariables();

    //! Import the values of the coupling variables
    void ImportCouplingVariables( const VectorType& CouplingVariables );

    //! Export the values of the coupling variables
    void ExportCouplingVariables( VectorType& CouplingVariables );

    //! Export the values of the coupling residuals
    void ExportCouplingResiduals( VectorType& CouplingResiduals );

    //! Export the Jacobian matrix
    /*!
     * @param Jacobian Matrix
     */
    void ExportJacobian( MatrixType& Jacobian );

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

    //! @name Private Methods
    //@{

    inline void loadModels( const std::string& FileName );
    inline void loadCouplings( const std::string& FileName );
    inline void loadGeometry( const std::string& FileName );

    template< typename number >
    inline std::vector< number > string2numVect( const std::string& string );

    //@}

    // Models & Couplings
    ModelsVector_Type        M_modelsList;
    CouplingsVector_Type     M_couplingsList;
};

//! Factory create function
inline MS_PhysicalModel* createMultiScale()
{
    return new MS_Model_MultiScale();
}

} // Namespace LifeV

#endif /* MS_Model_MultiScale_H */
