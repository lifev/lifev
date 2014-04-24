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
  @file IonicMinimalModel
  @brief Ionic model of Bueno-Orovio et al.

  Model as in
  Bueno-Orovio, Alfonso, Elizabeth M. Cherry, and Flavio H. Fenton.
  "Minimal model for human ventricular action potentials in tissue."
  Journal of theoretical biology 253.3 (2008): 544-560.

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICMINIMALMODEL_H_
#define _IONICMINIMALMODEL_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.

class IonicMinimalModel : public virtual ElectroIonicModel
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
    IonicMinimalModel();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicMinimalModel ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicMinimalModel object
     */
    IonicMinimalModel ( const IonicMinimalModel& model );
    //! Destructor
    virtual ~IonicMinimalModel() {}

    //@}

    //! @name Overloads
    //@{

    IonicMinimalModel& operator= ( const IonicMinimalModel& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& Uo()             const
    {
        return M_uo;
    }
    inline const Real& Uu()             const
    {
        return M_uu;
    }
    inline const Real& Tetav()      const
    {
        return M_tetav;
    }
    inline const Real& Tetaw()      const
    {
        return M_tetaw;
    }
    inline const Real& Tetavm()         const
    {
        return M_tetavm;
    }
    inline const Real& Tetao()      const
    {
        return M_tetao;
    }
    inline const Real& Tauv1()      const
    {
        return M_tauv1;
    }
    inline const Real& Tauv2()      const
    {
        return M_tauv2;
    }
    inline const Real& Tauvp()      const
    {
        return M_tauvp;
    }
    inline const Real& Tauw1()      const
    {
        return M_tauw1;
    }
    inline const Real& Tauw2()      const
    {
        return M_tauw2;
    }
    inline const Real& Kw()             const
    {
        return M_kw;
    }
    inline const Real& Uw()             const
    {
        return M_uw;
    }
    inline const Real& Tauwp()          const
    {
        return M_tauwp;
    }
    inline const Real& Taufi()      const
    {
        return M_taufi;
    }
    inline const Real& Tauo1()      const
    {
        return M_tauo1;
    }
    inline const Real& Tauo2()      const
    {
        return M_tauo2;
    }
    inline const Real& Tauso1()         const
    {
        return M_tauso1;
    }
    inline const Real& Tauso2()         const
    {
        return M_tauso2;
    }
    inline const Real& Kso()            const
    {
        return M_kso;
    }
    inline const Real& Uso()            const
    {
        return M_uso;
    }
    inline const Real& Taus1()      const
    {
        return M_taus1;
    }
    inline const Real& Taus2()      const
    {
        return M_taus2;
    }
    inline const Real& Ks()             const
    {
        return M_ks;
    }
    inline const Real& Us()         const
    {
        return M_us;
    }
    inline const Real& Tausi()      const
    {
        return M_tausi;
    }
    inline const Real& Tauinf()         const
    {
        return M_tauwinf;
    }
    inline const Real& Winfstart()  const
    {
        return M_winfstar;
    }


    inline void setUo           ( const Real& p )
    {
        this->M_uo        = p;
    }
    inline void setUu           ( const Real& p )
    {
        this->M_uu        = p;
    }
    inline void setTetav        ( const Real& p )
    {
        this->M_tetav     = p;
    }
    inline void setTetaw        ( const Real& p )
    {
        this->M_tetaw     = p;
    }
    inline void setTetavm       ( const Real& p )
    {
        this->M_tetavm    = p;
    }
    inline void setTetao        ( const Real& p )
    {
        this->M_tetao     = p;
    }
    inline void setTauv1        ( const Real& p )
    {
        this->M_tauv1     = p;
    }
    inline void setTauv2        ( const Real& p )
    {
        this->M_tauv2     = p;
    }
    inline void setTauvp        ( const Real& p )
    {
        this->M_tauvp     = p;
    }
    inline void setTauw1        ( const Real& p )
    {
        this->M_tauw1     = p;
    }
    inline void setTauw2        ( const Real& p )
    {
        this->M_tauw2     = p;
    }
    inline void setKw           ( const Real& p )
    {
        this->M_kw        = p;
    }
    inline void setUw           ( const Real& p )
    {
        this->M_uw        = p;
    }
    inline void setTauwp        ( const Real& p )
    {
        this->M_tauwp     = p;
    }
    inline void setTaufi        ( const Real& p )
    {
        this->M_taufi     = p;
    }
    inline void setTauo1        ( const Real& p )
    {
        this->M_tauo1     = p;
    }
    inline void setTauo2        ( const Real& p )
    {
        this->M_tauo2     = p;
    }
    inline void setTauso1       ( const Real& p )
    {
        this->M_tauso1    = p;
    }
    inline void setTauso2       ( const Real& p )
    {
        this->M_tauso2    = p;
    }
    inline void setKso      ( const Real& p )
    {
        this->M_kso       = p;
    }
    inline void setUso      ( const Real& p )
    {
        this->M_uso       = p;
    }
    inline void setTaus1        ( const Real& p )
    {
        this->M_taus1     = p;
    }
    inline void setTaus2        ( const Real& p )
    {
        this->M_taus2     = p;
    }
    inline void setKs           ( const Real& p )
    {
        this->M_ks        = p;
    }
    inline void setUs           ( const Real& p )
    {
        this->M_us        = p;
    }
    inline void setTausi        ( const Real& p )
    {
        this->M_tausi     = p;
    }
    inline void setTauinf       ( const Real& p )
    {
        this->M_tauwinf   = p;
    }
    inline void setWinfstart    ( const Real& p )
    {
        this->M_winfstar  = p;
    }



    //inline const short int& Size() const { return M_numberOfEquations; }
    //@}

    //! @name Methods
    //@{

    inline static Real Heaviside ( const Real& x)
    {
        if ( x > 0 )
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    //
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
    Real M_uo;
    Real M_uu;
    Real M_tetav;
    Real M_tetaw;
    Real M_tetavm;
    Real M_tetao;
    Real M_tauv1;
    Real M_tauv2;
    Real M_tauvp;
    Real M_tauw1;
    Real M_tauw2;
    Real M_kw;
    Real M_uw;
    Real M_tauwp;
    Real M_taufi;
    Real M_tauo1;
    Real M_tauo2;
    Real M_tauso1;
    Real M_tauso2;
    Real M_kso;
    Real M_uso;
    Real M_taus1;
    Real M_taus2;
    Real M_ks;
    Real M_us;
    Real M_tausi;
    Real M_tauwinf;
    Real M_winfstar;


    //! Xb states == equivalent to the number of equations
    //short int M_numberOfEquations;

    //@}

}; // class IonicMinimalModel



class MMTanhFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<1>& p)
    {
        return std::tanh ( p[0] );
    }
    return_Type operator() (const Real& p)
    {
        return std::tanh ( p );
    }


    MMTanhFunctor() {}
    MMTanhFunctor (const MMTanhFunctor&) {}
    ~MMTanhFunctor() {}
};

class MMHFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<1>& p)
    {
        if (p[0] > 0.0)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    return_Type operator() (const Real& p)
    {
        if (p > 0.0)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }


    MMHFunctor() {}
    MMHFunctor (const MMHFunctor&) {}
    ~MMHFunctor() {}
};

class MMSV
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<1>& p)
    {
        std::cout.precision (15);
        std::cout << "\nvalue: " << p[0];
    }
    return_Type operator() (const Real& p)
    {
        std::cout.precision (15);
        std::cout << "\nvalue: " << p;
    }


    MMSV() {}
    MMSV (const MMSV&) {}
    ~MMSV() {}
};



}

#endif
