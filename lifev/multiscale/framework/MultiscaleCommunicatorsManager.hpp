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
 *  @brief File containing the Multiscale Communicators Manager
 *
 *  @date 13-04-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleCommunicatorsManager_H
#define MultiscaleCommunicatorsManager_H 1

#include <lifev/multiscale/framework/MultiscaleDefinitions.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCommunicatorsManager - The Multiscale Communicators Manager
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleCommunicatorsManager class partitions a communicator among different models.
 */
class MultiscaleCommunicatorsManager
{
public:

    //! @name Type definitions
    //@{

    typedef std::map< UInt, multiscaleCommPtr_Type >                modelsCommunicatorContainer_Type;
    typedef modelsCommunicatorContainer_Type::const_iterator        modelsCommunicatorContainerIterator_Type;
    typedef std::vector< UInt >                                     modelsID_Type;
    typedef modelsID_Type::iterator                                 modelsIDIterator_Type;
    typedef std::vector< Real >                                     modelsLoad_Type;
    typedef modelsLoad_Type::iterator                               modelsLoadIterator_Type;
    typedef std::vector< std::vector< Int > >                       modelsProcessesList_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleCommunicatorsManager();

    //! Destructor
    virtual ~MultiscaleCommunicatorsManager() {}

    //@}


    //! @name Methods
    //@{

    //! Split the communicator among the models
    void splitCommunicator();

    //! Determine if the model is owned by the process
    /*!
     * @param modelID ID of the model.
     * @return true if the model is owned by the process, false otherwise
     */
    bool myModel ( const UInt& modelID ) const;

    //! Determine the number of model owned by the process
    /*!
     * @return number of model owned by the process
     */
    UInt myModelsNumber() const
    {
        return M_commContainer.size();
    }

    //! Display some information about the communicators
    void showMe();

    //@}


    //! @name Set Methods
    //@{

    //! Set the main epetra communicator
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator ( const multiscaleCommPtr_Type& comm )
    {
        M_comm = comm;
    }

    //! Add a group of models
    /*!
     * This method add a group of models for the forthcoming partitioning of the communicator.
     *
     * @param load percentage load of the model (-1 means each model on a different processor).
     * @param modelsIDList list of models.
     */
    void addGroup ( const Real& load, const modelsID_Type& modelsID );

    //@}


    //! @name Get Methods
    //@{

    //! Get the communicator of a specific model
    /*!
     * @param modelID ID of the model.
     * @return communicator.
     */
    const multiscaleCommPtr_Type& modelCommunicator ( const UInt& modelID ) const
    {
        return M_commContainer.find ( modelID )->second;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCommunicatorsManager ( const MultiscaleCommunicatorsManager& solver );

    MultiscaleCommunicatorsManager& operator= ( const MultiscaleCommunicatorsManager& solver );

    //@}


    //! @name Private Methods
    //@{

    void parallelProcessesDistribution ( std::vector<Real>& localNumberOfProcesses, const Int& numberOfProcesses );

    void parallelProcessesAssignment ( std::vector< std::vector< Int > >& parallelProcesses, const std::vector<Real>& localNumberOfProcesses, const Int& numberOfProcesses );

    //! Round a real number to the closest integer
    /*!
     * NOTE: x.5 is rounded to x+1;
     * @param value Real value
     * @return rounded integer value
     */
    Int roundToInteger ( const Real& value ) const
    {
        return static_cast<Int> ( std::floor ( value + 0.5 ) );
    }

    //@}

    // Main Communicator
    multiscaleCommPtr_Type              M_comm;

    // Models communicators
    modelsCommunicatorContainer_Type    M_commContainer;

    // Serial models data
    modelsID_Type                       M_serialModelsID;
    modelsProcessesList_Type            M_serialProcesses;

    // Parallel models data
    modelsID_Type                       M_parallelModelsID;
    modelsLoad_Type                     M_parallelModelsLoad;
    modelsProcessesList_Type            M_parallelProcesses;
};

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleCommunicatorsManager_H */
