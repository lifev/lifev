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
  @file IonicLuoRudyI
  @brief Ionic model Luo-Rudy I

	model as in
	Luo, Ching-hsing, and Yoram Rudy.
	"A model of the ventricular cardiac action potential.
	Depolarization, repolarization, and their interaction."
	Circulation research 68.6 (1991): 1501-1526.

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICLUORUDYI_H_
#define _IONICLUORUDYI_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.

class IonicLuoRudyI : public virtual ElectroIonicModel
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
    IonicLuoRudyI();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicLuoRudyI ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicLuoRudyI object
     */
    IonicLuoRudyI ( const IonicLuoRudyI& model );
    //! Destructor
    virtual ~IonicLuoRudyI() {}

    //@}

    //! @name Overloads
    //@{

    IonicLuoRudyI& operator= ( const IonicLuoRudyI& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real ENa()             const
    {
        return M_ENa;
    }
    inline const Real gNa()             const
    {
        return M_gNa;
    }

    //    //Slow Inward Current
    inline const Real gsi()             const
    {
        return M_gsi;
    }

    //    //Time dependent potassium current
    inline const Real K0()             const
    {
        return M_K0;
    }
    inline const Real gK()             const
    {
        return M_gK;
    }
    inline const Real EK()             const
    {
        return M_EK;
    }

    //    //Time independent potassium current
    inline const Real gK1()             const
    {
        return M_gK1;
    }
    inline const Real EK1()             const
    {
        return M_EK1;
    }

    //    //Plateau Current
    inline const Real EKp()             const
    {
        return M_EKp;
    }
    inline const Real gKp()             const
    {
        return M_gKp;
    }

    //    //Background Current
    inline const Real gb()             const
    {
        return M_gb;
    }


    //set methods
    inline void setENa  ( const Real p )
    {
        M_ENa = p;
    }
    inline void setGNa  ( const Real p )
    {
        M_gNa = p;
    }

    //    //Slow Inward Current
    inline void setGsi  ( const Real p )
    {
        M_gsi = p;
    }

    //    //Time dependent potassium current
    inline void setK0  ( const Real p )
    {
        M_K0 = p;
        M_gK = computeGK (M_K0);
        M_gK1 = computeGK1 (M_K0);
    }
    inline void setGK  ( const Real p )
    {
        M_gK = p;
    }
    inline void setEK  ( const Real p )
    {
        M_EK = p;
    }

    //    //Time independent potassium current
    inline void setGK1  ( const Real p )
    {
        M_gK1 = p;
    }
    inline void setEK1  ( const Real p )
    {
        M_EK1 = p;
    }

    //    //Plateau Current
    inline void setEKp  ( const Real p )
    {
        M_EKp = p;
    }
    inline void setGKp  ( const Real p )
    {
        M_gKp = p;
    }

    //    //Background Current
    inline void setGb  ( const Real p )
    {
        M_gb = p;
    }

    //inline const short int& Size() const { return M_numberOfEquations; }
    //@}

    //! @name Methods
    //@{

    //Fast inward current
    //===================
    //gating variable m
    inline Real am ( Real V )
    {
        return ( 0.32 * ( V + 47.13 ) / ( 1.0 - std::exp ( - 0.1 * ( V + 47.13 ) ) ) );
    }
    inline Real bm ( Real V )
    {
        return ( 0.08 * std::exp ( - V / 11.0 ) );
    }
    inline Real tm (Real V)
    {
        return ( 1.0 / ( am (V) + bm (V) ) );
    }
    inline Real minf (Real V)
    {
        return ( am (V) * tm (V) );
    }
    inline Real dm (Real V, Real m)
    {
        return ( am (V) * (1.0 - m ) - bm (V) * m );
    }

    //gating variable m
    inline Real ah ( Real V )
    {
        if ( V < - 40.0 )
        {
            return ( 0.135 * std::exp ( - ( V + 80 ) / 6.8 ) );
        }
        else
        {
            return 0.0;
        }
    }
    inline Real bh ( Real V )
    {
        if ( V < - 40.0 )
        {
            return ( ( 3.56  * std::exp ( 0.079 * V ) + 3.1 * 1e5 * std::exp ( 0.35 * V) ) );
        }
        else
        {
            return ( 1.0 / ( 0.13 * ( 1 + std::exp ( - ( V + 10.66 ) / 11.1 ) ) ) );
        }
    }
    inline Real th (Real V)
    {
        return ( 1.0 / ( ah (V) + bh (V) ) );
    }
    inline Real hinf (Real V)
    {
        return ( ah (V) * th (V) );
    }
    inline Real dh (Real V, Real h)
    {
        return ( ah (V) * (1.0 - h ) - bh (V) * h );
    }

    //gating variable j
    inline Real aj ( Real V )
    {
        if ( V < - 40.0 )
            return ( ( -1.2714 * 1e5 * std::exp ( 0.2444 * V ) - 3.474 * 1e-5 * std::exp ( - 0.04391 * V ) )
                     * ( V + 37.78 ) / ( 1 + std::exp ( 0.311 * ( V + 79.23 ) ) ) );
        else
        {
            return 0.0;
        }
    }
    inline Real bj ( Real V )
    {
        if ( V < - 40.0 )
        {
            return ( 0.1212 * std::exp ( - 0.01052 * V ) / ( 1 + std::exp ( -0.1378 * ( V + 40.14 ) ) ) );
        }
        else
        {
            return ( 0.3 * std::exp ( -2.535 * 1e-7 * V ) / ( 1 + std::exp ( -0.1 * ( V + 32 ) ) ) );
        }
    }
    inline Real tj (Real V)
    {
        return ( 1.0 / ( aj (V) + bj (V) ) );
    }
    inline Real jinf (Real V)
    {
        return ( aj (V) * tj (V) );
    }
    inline Real dj (Real V, Real j)
    {
        return ( aj (V) * (1.0 - j ) - bj (V) * j );
    }

    //Fast Inward Current
    inline Real INa (Real V, Real m, Real h, Real j)
    {
        return ( M_gNa * m * m * m * h * j * ( V - M_ENa ) );
    }

    //    // slow inward current %%
    //    //=========================
    inline Real Esi ( Real Ca )
    {
        return ( 7.7 - 13.0287 * std::log (Ca) );
    }
    //gating variable d
    inline Real ad ( Real V )
    {
        return ( 0.095 * std::exp ( - 0.01 * ( V - 5 ) ) / ( 1 + std::exp ( - 0.072 * ( V - 5 ) ) ) );
    }
    inline Real bd ( Real V )
    {
        return ( 0.07 * std::exp ( - 0.017 * ( V + 44) ) / ( 1 + std::exp ( 0.05 * ( V + 44) ) ) );
    }
    inline Real td (Real V)
    {
        return ( 1.0 / ( ad (V) + bd (V) ) );
    }
    inline Real dinf (Real V)
    {
        return ( ad (V) * td (V) );
    }
    inline Real dd (Real V, Real d)
    {
        return ( ad (V) * (1.0 - d ) - bd (V) * d );
    }

    //gating variable f
    inline Real af ( Real V )
    {
        return ( 0.012 * std::exp ( - 0.008 * ( V + 28 ) ) / ( 1 + std::exp ( 0.15 * ( V + 28 ) ) ) );
    }
    inline Real bf ( Real V )
    {
        return ( 0.0065 * std::exp ( - 0.02 * ( V + 30 ) ) / ( 1 + std::exp ( - 0.2  * ( V + 30 ) ) ) );
    }
    inline Real tf (Real V)
    {
        return ( 1.0 / ( af (V) + bf (V) ) );
    }
    inline Real finf (Real V)
    {
        return ( af (V) * tf (V) );
    }
    inline Real df (Real V, Real f)
    {
        return ( af (V) * (1.0 - f ) - bf (V) * f );
    }

    //Slow inward current
    inline Real Isi (Real V, Real d, Real f, Real Ca)
    {
        return ( M_gsi * d * f * ( V - Esi (Ca) ) );
    }

    //Calcium concentration
    inline Real dCa (Real V, Real d, Real f, Real Ca)
    {
        return ( - 1e-4 * Isi (V, d, f, Ca) + 0.07 * ( 1e-4 - Ca ) );
    }

    //Time dependent potassium current %%
    //===================================
    inline Real computeGK (Real K0)
    {
        return ( 0.282 * std::sqrt ( K0 / 5.4) );
    }

    //gating variable X
    inline Real aX ( Real V )
    {
        return ( 0.0005 * std::exp (  0.083 * ( V + 50 ) ) / ( 1 + std::exp ( 0.057 * ( V + 50 ) ) ) );
    }
    inline Real bX ( Real V )
    {
        return ( 0.0013 * std::exp ( - 0.06 * ( V + 20 ) ) / ( 1 + std::exp ( - 0.04  * ( V + 20 ) ) ) );
    }
    inline Real tX (Real V)
    {
        return ( 1.0 / ( aX (V) + bX (V) ) );
    }
    inline Real Xinf (Real V)
    {
        return ( aX (V) * tX (V) );
    }
    inline Real dX (Real V, Real X)
    {
        return ( aX (V) * (1.0 - X ) - bX (V) * X );
    }

    //gating Xi
    inline Real Xi (Real V)
    {
        if ( V > -100 )
        {
            return ( 2.837 * ( std::exp ( 0.04 * ( V + 77 ) ) - 1 ) / ( ( V + 77 ) * std::exp ( 0.04 * ( V + 35 ) ) ) );
        }
        else
        {
            return 1.0;
        }
    }
    //Time dependent potassium current
    inline Real IK (Real V, Real X)
    {
        return M_gK * X * Xi (V) * ( V - M_EK );
    }

    //Time independent potassium current
    //===================================
    inline Real computeGK1 (Real K0)
    {
        return ( 0.6047 * std::sqrt ( K0 / 5.4) );
    }

    //gating variable K1
    inline Real aK1 ( Real V )
    {
        return ( 1.02 / ( 1 + std::exp ( 0.2385 * ( V - M_EK1 - 59.2915 ) ) ) );
    }
    inline Real bK1 ( Real V )
    {
        return ( ( 0.49124 * std::exp ( 0.08032 * ( V - M_EK1 + 5.4760 ) )
                   + exp ( 0.06175 * ( V - M_EK1 - 594.31 ) ) ) / ( 1 + std::exp ( - 0.5143  * ( V - M_EK1 + 4.753 ) ) ) );
    }
    inline Real tK1 (Real V)
    {
        return ( 1.0 / ( aK1 (V) + bK1 (V) ) );
    }
    inline Real K1inf (Real V)
    {
        return ( aK1 (V) * tK1 (V) );
    }

    //Time independent potassium current
    inline Real IK1 (Real V)
    {
        return M_gK1 * K1inf (V) * ( V - M_EK1 );
    }

    //Plateau Current
    //====================
    inline Real Kp ( Real V )
    {
        return ( 1.0 / ( 1 + std::exp ( ( 7.488 - V ) / 5.98 ) ) );
    }
    inline Real IKp (Real V)
    {
        return M_gKp * Kp (V) * ( V - M_EKp );
    }

    //Background Current
    //====================
    inline Real Ib (Real V)
    {
        return M_gb * ( V + 59.87 );
    }

    //Total current
    //=============
    inline Real Itot (Real V, Real m, Real h, Real j, Real d, Real f, Real X, Real Ca)
    {
        return ( INa (V, m, h, j) + Isi (V, d, f, Ca) + IK (V, X) + IK1 (V) + IKp (V) + Ib (V) );
    }
    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeNonGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );


    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt );

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    //Fast Sodium Current
    Real M_ENa;
    Real M_gNa;

    //Slow Inward Current
    Real M_gsi;

    //Time dependent potassium current
    Real M_K0;
    Real M_gK;
    Real M_EK;

    //Time independent potassium current
    Real M_gK1;
    Real M_EK1;

    //Plateau Current
    Real M_EKp;
    Real M_gKp;
    //Background Current
    Real M_gb;
    //! Xb states == equivalent to the number of equations
    //short int M_numberOfEquations;

    //@}

}; // class IonicLuoRudyI



}

#endif
