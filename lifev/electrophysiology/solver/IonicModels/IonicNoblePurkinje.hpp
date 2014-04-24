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
  @file IonicNoblePurkinje
  @brief Ionic model of Noble for Purkinje cells

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


#ifndef _IONICNOBLEPURKINJE_H_
#define _IONICNOBLEPURKINJE_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.

class IonicNoblePurkinje : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroIonicModel                         super;
    typedef boost::shared_ptr<VectorEpetra>         vectorPtr_Type;
    typedef boost::shared_ptr<VectorElemental>  elvecPtr_Type;
    typedef RegionMesh<LinearTetra>                 mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicNoblePurkinje();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicNoblePurkinje ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicNoblePurkinje object
     */
    IonicNoblePurkinje ( const IonicNoblePurkinje& model );
    //! Destructor
    virtual ~IonicNoblePurkinje() {}

    //@}

    //! @name Overloads
    //@{

    IonicNoblePurkinje& operator= ( const IonicNoblePurkinje& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& gi()             const
    {
        return M_gi;
    }


    inline const Real& Cm()             const
    {
        return M_Cm;
    }

    inline const Real& vK()             const
    {
        return M_vK;
    }

    inline const Real& vNa()             const
    {
        return M_vNa;
    }

    inline const Real& Itotal()             const
    {
        return M_Itot;
    }


    inline void setgi           ( const Real& p )
    {
        this->M_gi       = p;
    }


    inline void setCm           ( const Real& p )
    {
        this->M_Cm        = p;
    }

    inline void setvNa          ( const Real& p )
    {
        this->M_vNa        = p;
    }

    inline void setvK           ( const Real& p )
    {
        this->M_vK         = p;
    }


    inline static Real GeneralFunctionAlphaAndBeta (const Real& V, const Real& C1, const Real& C2, const Real& C3, const Real& C4, const Real& C5, const Real& V_0)
    {
        if (! (C1 != 0 && C2 == 0) )
        {
            return (C1 * std::exp ( (V - V_0) / C2) + C3 * (V - V_0) ) / (1 + C4 * std::exp ( (V - V_0) / C5) ) ;
        }
        else
        {
            return (C1 + C3 * (V - V_0) ) / (1 + C4 * std::exp ( (V - V_0) / C5) ) ;
        }
    }


    inline static Real mInf (const Real& V)
    {
        Real alpham = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.1, -1, -15, -48);
        Real betam  = GeneralFunctionAlphaAndBeta (V, 0, 1, -0.12, -1, 5, -8);
        Real taum = alpham + betam;
        return alpham / (taum);
    }

    inline static Real hInf (const Real& V)
    {
        Real alphah = GeneralFunctionAlphaAndBeta (V, 0.17, -20, 0, 0, 1, -90);
        Real betah  = GeneralFunctionAlphaAndBeta (V, 1, 0, 0, 1, -10, -42);
        Real tauh = alphah + betah;
        return alphah / (tauh);
    }

    inline static Real nInf (const Real& V)
    {
        Real alphan = GeneralFunctionAlphaAndBeta (V, 0, 1, 0.0001, -1, -10, -50);
        Real betan  = GeneralFunctionAlphaAndBeta (V, 0.002, -80, 0, 0, 1, -90);
        Real taun = alphan + betan;
        return alphan / (taun);
    }
    //inline const short int& Size() const { return M_numberOfEquations; }
    //@}


    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt );


    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_gi;
    Real M_vNa;
    Real M_vK;
    Real M_Cm;
    Real M_Itot;


    //! Xb states == equivalent to the number of equations
    //short int M_numberOfEquations;

    //@}

}; // class IonicMinimalModel



}

#endif
