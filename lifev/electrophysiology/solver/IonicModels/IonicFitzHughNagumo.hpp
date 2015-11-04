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
  @file IonicFitzHughNagumo
  @brief Ionic model of FitzHugh-Nagumo

  Model taken as in
  Franzone, Piero Colli, et al.
  "Adaptivity in space and time for reaction-diffusion systems in electrocardiology."
  SIAM Journal on Scientific Computing 28.3 (2006): 942-962.


  @date 01-2013
  @author Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>

  @contributors
  @mantainer Giacomo Rosilho de Souza <giacomo.rosilhodesouza@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICFITZHUGHNAGUMO_H_
#define _IONICFITZHUGHNAGUMO_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.

//! @name Type definitions
//@{
typedef ElectroIonicModel super;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef std::shared_ptr<vector_Type> vectorPtr_Type;
typedef RegionMesh<LinearTetra> mesh_Type;
//@}

class IonicFitzHughNagumo : public virtual ElectroIonicModel
{

public:
    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicFitzHughNagumo();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicFitzHughNagumo ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicFitzHughNagumo object
     */
    IonicFitzHughNagumo ( const IonicFitzHughNagumo& model );
    //! Destructor
    virtual ~IonicFitzHughNagumo() {}

    //@}

    //! @name Overloads
    //@{

    IonicFitzHughNagumo& operator= ( const IonicFitzHughNagumo& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& G()        const
    {
        return M_G;
    }
    inline const Real& Vth()      const
    {
        return M_Vth;
    }
    inline const Real& Vp()    const
    {
        return M_Vp;
    }
    inline const Real& Eta1()    const
    {
        return M_Eta1;
    }
    inline const Real& Eta2()    const
    {
        return M_Eta2;
    }
    inline const Real& Eta3()    const
    {
        return M_Eta3;
    }
    inline const Real& Eta()    const
    {
        return M_Eta;
    }
    inline const Real& Gamma()    const
    {
        return M_Gamma;
    }

    inline void setG      ( const Real& G )
    {
        this->M_G       = G;
    }
    inline void setVth    ( const Real& Vth )
    {
        this->M_Vth     = Vth;
    }
    inline void setVp     ( const Real& Vp )
    {
        this->M_Vp      = Vp;
        this->M_Eta     = M_Eta2 / M_Vp;
    }
    inline void setEta1   ( const Real& Eta1 )
    {
        this->M_Eta1    = Eta1;
    }
    inline void setEta2   ( const Real& Eta2 )
    {
        this->M_Eta2    = Eta2;
        this->M_Eta     = M_Eta2 / M_Vp;
        this->M_Gamma   = M_Eta2 * M_Eta3;
    }
    inline void setEta3   ( const Real& Eta3 )
    {
        this->M_Eta3    = Eta3;
        this->M_Gamma   = M_Eta2 * M_Eta3;
    }

    void setup ( Teuchos::ParameterList& parameterList );
    //Compute the rhs on a single node or for the 0D case
    //this is to make visible the computeRhs defined in the base class, otherwise it is hidden by the ones defined here
    using ElectroIonicModel::computeRhs;
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);
    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);
    void computeRhs ( const VectorSmall<2>& v, VectorSmall<2>& rhs);

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    //compute the Jacobian
    matrix_Type getJac (const vector_Type& v, Real h = 1.0e-8);
    std::vector< std::vector<Real> > getJac (const std::vector<Real>& v, Real h = 1.0e-8);
    MatrixSmall<2, 2> getJac (const VectorSmall<2>& v, Real h = 1.0e-8);

    //! Display information about the model
    void showMe();


private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_G;
    Real M_Vth;
    Real M_Vp;
    Real M_Eta1;
    Real M_Eta2;
    Real M_Eta3;
    Real M_Eta;
    Real M_Gamma;


    //@}

}; // class IonicFitzHughNagumo


inline ElectroIonicModel* createIonicFitzHughNagumo()
{
    return new IonicFitzHughNagumo();
}

namespace
{
static bool register_IonicFitzHughNagumo = ElectroIonicModel::IonicModelFactory::instance().registerProduct ("FitzHughNagumo", &createIonicFitzHughNagumo );
}

}

#endif
