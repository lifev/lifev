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
  @file IonicGoldbeter
  @brief Intracellular Calcium model from Goldbeter et al. (1990). By "potential" we
  refer to the cytosolic calcium concentration, whereas the gating variable
   represents the sarcoplasmic calcium concentration

    Model as in
    Goldbeter, Albert, GENEVItVE Dupont, and Michael J. Berridge.
    "Minimal model for signal-induced Ca2+ oscillations and for their frequency encoding through protein phosphorylation."
    Proceedings of the National Academy of Sciences 87.4 (1990): 1461-1465.

  @date 09-2013

  @author Ricardo Ruiz <ricardo.ruizbaier@unil.ch>

  @contributors
  @mantainer Ricardo Ruiz <ricardo.ruizbaier@epfl.ch>
  @last update 09-2013
 */


#ifndef _IONICGOLDBETER_H_
#define _IONICGOLDBETER_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>


#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{


class IonicGoldbeter : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroIonicModel super;
    typedef std::shared_ptr<VectorEpetra> vectorPtr_Type;
    typedef RegionMesh<LinearTetra> mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicGoldbeter();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicGoldbeter ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IntracellularCalciumGoldbeter object
     */
    IonicGoldbeter ( const IonicGoldbeter& model );
    //! Destructor
    virtual ~IonicGoldbeter() {}

    //@}

    //! @name Overloads
    //@{

    IonicGoldbeter& operator= ( const IonicGoldbeter& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& Nu1()     const
    {
        return M_nu1;
    }
    inline const Real& Nu2()     const
    {
        return M_nu2;
    }
    inline const Real& Nu3()     const
    {
        return M_nu3;
    }
    inline const Real& Nu4()     const
    {
        return M_nu4;
    }
    inline const Real& Nu5()     const
    {
        return M_nu5;
    }
    inline const Real& K1()      const
    {
        return M_k1;
    }
    inline const Real& K2()      const
    {
        return M_k2;
    }
    inline const Real& K3()      const
    {
        return M_k3;
    }
    inline const Real& K4()      const
    {
        return M_k4;
    }

    inline void setNu1      ( const Real& nu1 )
    {
        this->M_nu1       = nu1;
    }
    inline void setNu2      ( const Real& nu2 )
    {
        this->M_nu2       = nu2;
    }
    inline void setNu3      ( const Real& nu3 )
    {
        this->M_nu3       = nu3;
    }
    inline void setNu4      ( const Real& nu4 )
    {
        this->M_nu4       = nu4;
    }
    inline void setNu5      ( const Real& nu5 )
    {
        this->M_nu5       = nu5;
    }
    inline void setK1        ( const Real& k1 )
    {
        this->M_k1         = k1;
    }
    inline void setK2        ( const Real& k2 )
    {
        this->M_k2         = k2;
    }
    inline void setK3        ( const Real& k3 )
    {
        this->M_k3         = k3;
    }
    inline void setK4        ( const Real& k4 )
    {
        this->M_k4         = k4;
    }

    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);
    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    //Compute the rhs on a mesh/ 3D case
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Iapp, std::vector<vectorPtr_Type>& rhs );
    //

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v);

    //    void computePotentialRhs(     const std::vector<vectorPtr_Type>& v,
    //                      const VectorEpetra& Iapp,
    //                      std::vector<vectorPtr_Type>& rhs,
    //                      FESpace<mesh_Type, MapEpetra>& uFESpace );

    //! Display information about the model
    void showMe();

    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_nu1;
    Real M_nu2;
    Real M_nu3;
    Real M_nu4;
    Real M_nu5;
    Real M_k1;
    Real M_k2;
    Real M_k3;
    Real M_k4;

    //@}

}; // class IntracellularCalciumGoldbeter

}

#endif
