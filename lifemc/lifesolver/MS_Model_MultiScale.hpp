/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-03-12

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Model_MultiScale.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-03-12
 */

#ifndef __MS_Model_MultiScale_H
#define __MS_Model_MultiScale_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <life/lifecore/life.hpp>
#include <life/lifearray/tab.hpp>

#include <lifemc/lifesolver/MS_PhysicalData.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Coupling_BoundaryCondition.hpp>
#include <lifemc/lifesolver/MS_Coupling_Stress.hpp>
#include <lifemc/lifesolver/MS_Coupling_FluxStress.hpp>

#include <boost/array.hpp>
#include <boost/algorithm/string.hpp>

namespace LifeV {

//! MS_Model_MultiScale - MultiScale model
/*!
 *  The MS_Model_MultiScale class is an implementation of the MS_PhysicalModel
 *  for a general multiscale problem.
 *
 *  @author Cristiano Malossi
 */
class MS_Model_MultiScale: public virtual MS_PhysicalModel
{
public:

    typedef MS_PhysicalModel                                            super;

    typedef singleton< factory< MS_PhysicalModel, modelsTypes > >       FactoryModels;
    typedef singleton< factory< MS_PhysicalCoupling, couplingsTypes > > FactoryCouplings;

    typedef boost::shared_ptr< MS_PhysicalModel >                       PhysicalModel_ptr;
    typedef boost::shared_ptr< MS_PhysicalCoupling >                    PhysicalCoupling_ptr;

    //! @name Constructors, Destructor
    //@{

    //! Constructor
    MS_Model_MultiScale();

    //! Copy constructor
    /*!
     * \param multiscale - MS_Model_MultiScale
     */
    MS_Model_MultiScale( const MS_Model_MultiScale& multiscale );

    //! Destructor
    ~MS_Model_MultiScale() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param multiscale - MS_Model_MultiScale
     */
    MS_Model_MultiScale& operator=( const MS_Model_MultiScale& multiscale );

    //! Build the global map for the coupling vectors
    /*!
     * \param couplingVariablesGlobalNumber - Global number of coupling variables
     * \param MyGlobalElements - Epetra_Map MyGlobalElements
     */
    void BuildCouplingVectorsMap( UInt&             couplingVariablesGlobalNumber,
                                  std::vector<int>& MyGlobalElements );

    //! Setup parameters for the implicit coupling
    void InitializeCouplingVariables( void );

    //! Import the values of the coupling variables
    void ImportCouplingVariables( const VectorType& CouplingVariables );

    //! Export the values of the coupling variables
    void ExportCouplingVariables( VectorType& CouplingVariables );

    //! Export the values of the coupling residuals
    void ExportCouplingResiduals( VectorType& CouplingResiduals );

    //! Export the values of the Jacobian product
    /*!
     * \param deltaCouplingVariables - variation of the coupling variables
     */
    void ExportJacobianProduct( const VectorType& deltaCouplingVariables, VectorType& JacobianProduct );

    //@}


    //! @name Get Methods
    //@{

    //! Get the number of the coupling variables
    UInt GetCouplingVariablesNumber( void );

    //@}


    //! @name MultiScale Physical Model
    //@{

    //! Load data and create the models and the couplings
    void SetupData( void );

    //! Setup the model
    void SetupModel( void );

    //! Build the system matrix and vectors
    void BuildSystem( void );

    //! Update the system matrix and vectors
    void UpdateSystem( void );

    //! Solve the problem
    void SolveSystem( void );

    //! Save the solution
    void SaveSolution( void );

    //! Display some information about the multiscale problem (call after SetupProblem)
    void ShowMe( void );

    //@}

private:

    //! @name Private Methods
    //@{

    inline void loadModels();
    inline void loadCouplings();
    inline void loadGeometry();

    template< typename number >
    inline std::vector< number > string2numVect( const std::string& string );

    //@}

    // Models & Connections
    std::map< UInt, PhysicalModel_ptr >    M_models;
    std::map< UInt, PhysicalCoupling_ptr > M_couplings;

    // Models & Connections size
    UInt                                   M_modelsNumber;
    UInt                                   M_couplingsNumber;
};

//! Factory create function
inline MS_PhysicalModel* createMultiScale()
{
    return new MS_Model_MultiScale();
}

} // Namespace LifeV

#endif /* __MS_Model_MultiScale_H */
