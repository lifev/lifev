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

#include <lifemc/lifesolver/MultiscaleDefinitions.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleCommunicatorsManager - The Multiscale Communicators Manager
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleCommunicatorsManager class partition a communicator among different models
 *  and couplings.
 */
class MultiscaleCommunicatorsManager
{
public:

    //! @name Type definitions
    //@{

    typedef std::vector< multiscaleCommPtr_Type >                modelsCommunicatorContainer_Type;

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
    void splitCommunicators();

    //! Display some information about the Multiscale problem (should be called after setupProblem)
    void showMe();

    //@}


    //! @name Set Methods
    //@{

    //! Set the epetra communicator for the Multiscale problem
    /*!
     * @param comm Epetra communicator
     */
    void setCommunicator( const multiscaleCommPtr_Type& comm ) { M_comm = comm; }

    //! Add a group of models
    /*!
     * This method add a group of models for the forthcoming partitioning of the communicator.
     *
     * @param load percentage load of the model (-1 means each model on a different processor).
     * @param modelsIDList list of models.
     */
    void addGroup( const Real& load, const std::vector< UInt >& modelsID )
    {
        if ( load < 0 )
            M_serialModelsID.insert( M_serialModelsID.end(), modelsID.begin(), modelsID.end() );
        else
        {
            M_parallelModelsID.insert( M_parallelModelsID.end(), modelsID.begin(), modelsID.end() );

            std::vector< Real > loadVector( modelsID.size(), load );
            M_parallelModelsLoad.insert( M_parallelModelsLoad.end(), loadVector.begin(), loadVector.end() );
        }
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the communicator of a specific model
    /*!
     * @param modelID ID of the model.
     * @return communicator.
     */
    const multiscaleCommPtr_Type& modelCommunicator( const UInt& modelID ) { return M_commContainer[modelID];}

    //! Get the communicator of a specific coupling
    /*!
     * @param listOfModels vector of models ID.
     * @return communicator.
     */
    const multiscaleCommPtr_Type& couplingCommunicator( const std::vector< UInt >& listOfModels ) {}

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleCommunicatorsManager( const MultiscaleCommunicatorsManager& solver );

    MultiscaleCommunicatorsManager& operator=( const MultiscaleCommunicatorsManager& solver );

    //@}

    // Main Communicator
    multiscaleCommPtr_Type              M_comm;

    // Models communicator
    modelsCommunicatorContainer_Type    M_commContainer;

    std::vector< UInt >                 M_serialModelsID;
    std::vector< UInt >                 M_parallelModelsID;
    std::vector< Real >                 M_parallelModelsLoad;
};

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleCommunicatorsManager_H */
