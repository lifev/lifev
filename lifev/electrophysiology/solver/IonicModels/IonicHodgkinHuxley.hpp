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
  @file IonicHodgkinHuxley
  @brief Ionic model of Hodgkin and Huxley

    Model as in
    Keener, James, and James Sneyd.
    Mathematical Physiology: I: Cellular Physiology. Vol. 1.
    Springer, 2010.

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICHODGKINHUXLEY_H_
#define _IONICHODGKINHUXLEY_H_

//Include the base class
#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <Teuchos_ParameterList.hpp>

namespace LifeV
{
//! IonicModel - This class implements the Hodgkin-Huxley model.

class IonicHodgkinHuxley : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroIonicModel                         super;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicHodgkinHuxley();


    /*!
     * @param list of parameters in an xml file
     */
    IonicHodgkinHuxley ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicHodgkinHuxley object
     */
    IonicHodgkinHuxley ( const IonicHodgkinHuxley& model );
    //! Destructor
    virtual ~IonicHodgkinHuxley() {}

    //@}

    //! @name Overloads
    //@{

    IonicHodgkinHuxley& operator= ( const IonicHodgkinHuxley& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& gNa()             const
    {
        return M_gNa;
    }

    inline const Real& gK()             const
    {
        return M_gK;
    }

    inline const Real& gL()             const
    {
        return M_gL;
    }


    inline const Real& vL()             const
    {
        return M_vL;
    }

    inline const Real& vK()             const
    {
        return M_vK;
    }

    inline const Real& vNa()             const
    {
        return M_vNa;
    }




    inline void setGNa           ( const Real& p )
    {
        this->M_gNa        = p;
    }

    inline void setGK           ( const Real& p )
    {
        this->M_gK        = p;
    }

    inline void setGL           ( const Real& p )
    {
        this->M_gL        = p;
    }

    inline void setVL           ( const Real& p )
    {
        this->M_vL        = p;
    }

    inline void setVNa          ( const Real& p )
    {
        this->M_vNa        = p;
    }

    inline void setVK           ( const Real& p )
    {
        this->M_vK         = p;
    }

    //@}


    //! @name Methods
    //@{
    //! Compute the rhs of the gating variables a single node or for the 0D case
    /*!
     * @param v vector with the variables (V, M, N, H)
     * @param rhs vector where we will insert the rhs of the equations of the gating variables
     */
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    //! compute the rhs of the potential equation on a single node or for the 0D case
    /*!
     * @param v vector with the variables (V, M, N, H)
     */
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    //! Compute the rhs on a single node or for the 0D case
    /*!
     * @param v vector with the variables (V, M, N, H)
     * @param rhs vector where we will insert the rhs of the equations
     */
    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    //! compute the rhs of the gating variables with the RushLarsen scheme
    /*!
     * @param v vector with the variables (V, M, N, H)
     * @param dt timestep
     */
    void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt );

    //! Display information about the model
    void showMe();

    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_gNa;
    Real M_gK;
    Real M_gL;
    Real M_vNa;
    Real M_vK;
    Real M_vL;

    //@}

}; // class IonicMinimalModel



}// end namespace LifeV

#endif
