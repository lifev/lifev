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
  @file IonicMitchellSchaeffer
  @brief Ionic model based on Mitchell-Schaeffer model.

  model as in
  Mitchell, Colleen C., and David G. Schaeffer.
  "A two-current model for the dynamics of cardiac membrane."
  Bulletin of mathematical biology 65.5 (2003): 767-793.

  @date 03-2013
  @author Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>

  @contributors
  @mantainer Luis Miguel De Oliveira Vilaca <luismiguel.deoliveiravilaca@epfl.ch>
  @last update 03-2013
 */


#ifndef _IONICMITCHELLSCHAEFFER_H_
#define _IONICMITCHELLSCHAEFFER_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.


class IonicMitchellSchaeffer : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroIonicModel super;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
    typedef RegionMesh<LinearTetra> mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicMitchellSchaeffer();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicMitchellSchaeffer ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicMitchellSchaeffer object
     */
    IonicMitchellSchaeffer ( const IonicMitchellSchaeffer& model );
    //! Destructor
    virtual ~IonicMitchellSchaeffer() {}

    //@}

    //! @name Overloads
    //@{

    IonicMitchellSchaeffer& operator= ( const IonicMitchellSchaeffer& model );

    //@}

    //! @name Setters and getters
    //@{

    inline const Real& vGate() const
    {
        return M_vGate;
    }

    inline void setVGate (const Real& vGate)
    {
        this->M_vGate = vGate;
    }

    inline const Real& tauClose() const
    {
        return M_tauClose;
    }

    inline void setTauClose (const Real& tauClose)
    {
        this->M_tauClose = tauClose;
    }

    inline const Real& tauOpen() const
    {
        return M_tauOpen;
    }

    inline void setTauOpen (const Real& tauOpen)
    {
        this->M_tauOpen = tauOpen;
    }

    inline const Real& tauIn() const
    {
        return M_tauIn;
    }

    inline void setTauIn (const Real& tauIn)
    {
        this->M_tauIn = tauIn;
    }

    inline const Real& tauOut() const
    {
        return M_tauOut;
    }

    inline void setTauOut (const Real& tauOut)
    {
        this->M_tauOut = tauOut;
    }

    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );


    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v);
    Real computeLocalGatingRhs ( const std::vector<Real>& v );

    //! Display information about the model
    void showMe();


    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters

    Real M_vGate;
    Real M_tauClose;
    Real M_tauOpen;
    Real M_tauIn;
    Real M_tauOut;


    //@}

}; // class IonicMitchellSchaeffer

}

#endif
