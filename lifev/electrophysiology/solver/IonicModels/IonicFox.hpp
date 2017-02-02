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
  @file IonicFoxModel
  @brief Ionic model based on Fox model.

    Ionic model as in

    Fox, Jeffrey J., Jennifer L. McHarg, and Robert F. Gilmour Jr.
    "Ionic mechanism of electrical alternans."
    American Journal of Physiology-Heart and Circulatory Physiology 51.2 (2002): H516.

  @date 04-2013
  @author Marie Dupraz <dupraz.marie@gmail.com>

  @contributors
  @mantainer Marie Dupraz <dupraz.marie@gmail.com>
  @last update 04-2013
 */


#ifndef FOX_HPP_INCLUDED
#define FOX_HPP_INCLUDED

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <cmath>


namespace LifeV
{
//! IonicModel - This class implements an ionic model.

class IonicFox : public virtual ElectroIonicModel
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
    IonicFox();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicFox ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicFox object
     */
    IonicFox ( const IonicFox& model );

    //! Destructor
    virtual ~IonicFox() {}

    //@}

    //! @name Overloads
    //@{

    IonicFox& operator= ( const IonicFox& model );

    //@}

    //! @name Setters and getters
    //@{

    inline const Real& areaCap() const
    {
        return M_ACap;
    }
    inline void setACap (const Real& areaCap)
    {
        this->M_ACap = areaCap;
    }

    inline const Real& volMyo() const
    {
        return M_VMyo;
    }
    inline void setVMyo (const Real& volMyo)
    {
        this->M_VMyo = volMyo;
    }

    inline const Real& volSR() const
    {
        return M_Vsr;
    }
    inline void setVsr (const Real& volSR)
    {
        this->M_Vsr = volSR;
    }

    inline const Real& volUp() const
    {
        return M_Vup;
    }

    inline void setVup (const Real& volUp)
    {
        this->M_Vup = volUp;
    }

    inline const Real& concNa0() const
    {
        return M_NaO;
    }
    inline void setNa0 (const Real& concNa0)
    {
        this->M_NaO = concNa0;
    }

    inline const Real& concCa0() const
    {
        return M_CaO;
    }
    inline void setCa0 (const Real& concCa0)
    {
        this->M_CaO = concCa0;
    }

    inline const Real& concK0() const
    {
        return M_K0;
    }
    inline void setK0 (const Real& concK0)
    {
        this->M_K0 = concK0;
    }

    inline const Real& concNaIn() const
    {
        return M_NaIn;
    }
    inline void setNaIn (const Real& concNaIn)
    {
        this->M_NaIn = concNaIn;
    }

    inline const Real& concKIn() const
    {
        return M_KIn;
    }
    inline void setKIn (const Real& concKIn)
    {
        this->M_KIn = concKIn;
    }

    inline const Real& cmdnTot() const
    {
        return M_CmdnTot;
    }
    inline void setCmdnTot (const Real& cmdnTot)
    {
        this->M_CmdnTot = cmdnTot;
    }

    inline const Real& csqnTot() const
    {
        return M_CsqnTot;
    }
    inline void setCsqnTot (const Real& csqnTot)
    {
        this->M_CsqnTot = csqnTot;
    }

    inline const Real& constmCmdn() const
    {
        return M_KmCmdn;
    }
    inline void setKmCmdn (const Real& constmCmdn)
    {
        this->M_KmCmdn = constmCmdn;
    }

    inline const Real& constmCsqn() const
    {
        return M_KmCsqn;
    }
    inline void setKmCsqn (const Real& constmCsqn)
    {
        this->M_KmCsqn = constmCsqn;
    }

    inline const Real& capMem() const
    {
        return M_Cm;
    }
    inline void setCapMem (const Real& capMem)
    {
        this->M_Cm = capMem;
    }

    inline const Real& farad() const
    {
        return M_F;
    }
    inline void setFarad (const Real& farad)
    {
        this->M_F = farad;
    }

    inline const Real& temp() const
    {
        return M_T;
    }
    inline void setTemp (const Real& temp)
    {
        this->M_T = temp;
    }

    inline const Real& gasConst() const
    {
        return M_R;
    }
    inline void setR (const Real& gasConst)
    {
        this->M_R = gasConst;
    }

    inline const Real& maxCondNa() const
    {
        return M_GNa;
    }
    inline void setGNa (const Real& maxCondNa)
    {
        this->M_GNa = maxCondNa;
    }

    inline const Real& maxCondKp() const
    {
        return M_GKp;
    }
    inline void setGkp (const Real& maxCondKp)
    {
        this->M_GKp = maxCondKp;
    }

    inline const Real& maxCondK1() const
    {
        return M_GK1;
    }
    inline void setGk1 (const Real& maxCondK1)
    {
        this->M_GK1 = maxCondK1;
    }

    inline const Real& maxCondKr() const
    {
        return M_GKr;
    }
    inline void setGkr (const Real& maxCondKr)
    {
        this->M_GKr = maxCondKr;
    }

    inline const Real& maxCondKs() const
    {
        return M_GKs;
    }
    inline void setGks (const Real& maxCondKs)
    {
        this->M_GKs = maxCondKs;
    }

    inline const Real& maxCondt0() const
    {
        return M_Gt0;
    }
    inline void setGkt0 (const Real& maxCondt0)
    {
        this->M_Gt0 = maxCondt0;
    }

    inline const Real& kNaCa() const
    {
        return M_kNaCa;
    }
    inline void setKNaCa (const Real& kNaCa)
    {
        this->M_kNaCa = kNaCa;
    }

    inline const Real& constmfCa() const
    {
        return M_KmfCa;
    }
    inline void setKmfCa (const Real& constmfCa)
    {
        this->M_KmfCa = constmfCa;
    }

    inline const Real& constmK1() const
    {
        return M_KmK1;
    }
    inline void setKmK1 (const Real& constmK1)
    {
        this->M_KmK1 = constmK1;
    }

    inline const Real& constmNa() const
    {
        return M_KmNa;
    }
    inline void setKmNa (const Real& constmNa)
    {
        this->M_KmNa = constmNa;
    }

    inline const Real& constmCa() const
    {
        return M_KmCa;
    }
    inline void setKmCa (const Real& constmCa)
    {
        this->M_KmCa = constmCa;
    }

    inline const Real& kSat() const
    {
        return M_kSat;
    }
    inline void setKSat (const Real& kSat)
    {
        this->M_kSat = kSat;
    }

    inline const Real& eta() const
    {
        return M_eta;
    }
    inline void setEta (const Real& eta)
    {
        this->M_eta = eta;
    }

    inline const Real& courNaK() const
    {
        return M_INaK;
    }
    inline void setINaK (const Real& courNaK)
    {
        this->M_INaK = courNaK;
    }

    inline const Real& constmNai() const
    {
        return M_KmNai;
    }
    inline void setKmNai (const Real& constmNai)
    {
        this->M_KmNai = constmNai;
    }

    inline const Real& constmK0() const
    {
        return M_KmK0;
    }
    inline void setKmK0 (const Real& constmK0)
    {
        this->M_KmK0 = constmK0;
    }

    inline const Real& permCa() const
    {
        return M_PCa;
    }
    inline void setPnsCa (const Real& permCa)
    {
        this->M_PCa = permCa;
    }

    inline const Real& constmpCa() const
    {
        return M_KmpCa;
    }
    inline void setKmpCa (const Real& constmpCa)
    {
        this->M_KmpCa = constmpCa;
    }

    inline const Real& courpCa() const
    {
        return M_IpCa;
    }
    inline void setIpCa (const Real& courpCa)
    {
        this->M_IpCa = courpCa;
    }

    inline const Real& maxCondCab() const
    {
        return M_GCab;
    }
    inline void setGCab (const Real& maxCondCab)
    {
        this->M_GCab = maxCondCab;
    }

    inline const Real& maxCondNab() const
    {
        return M_GNab;
    }
    inline void setGNab (const Real& maxCondNab)
    {
        this->M_GNab = maxCondNab;
    }

    inline const Real& constmUp() const
    {
        return M_KmUp;
    }
    inline void setKmUp (const Real& constmUp)
    {
        this->M_KmUp = constmUp;
    }

    inline const Real& permCaK() const
    {
        return M_PCaK;
    }
    inline void setPCaK (const Real& permCaK)
    {
        this->M_PCaK = permCaK;
    }

    inline const Real& permrel() const
    {
        return M_Prel;
    }
    inline void setPrel (const Real& permrel)
    {
        this->M_Prel = permrel;
    }

    inline const Real& permleak() const
    {
        return M_Pleak;
    }
    inline void setPleak (const Real& permleak)
    {
        this->M_Pleak = permleak;
    }

    inline const Real& courCaHalf() const
    {
        return M_ICaHalf;
    }
    inline void setICaHalf (const Real& courCaHalf)
    {
        this->M_ICaHalf = courCaHalf;
    }

    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );


    //Compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    std::vector<Real> computeLocalGatingRhs ( const std::vector<Real>& v );

    void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt );

    void computeNonGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

    std::vector<Real> computeLocalSubSysCaRhs ( const std::vector<Real>& v );

    //Compute the ionic currents (Luo and Rudy)
    std::vector<Real> fastINa ( const std::vector<Real>& v );

    std::vector<Real> rapidIK ( const std::vector<Real>& v );

    std::vector<Real> transOutIK ( const std::vector<Real>& v );

    std::vector<Real> slowIK ( const std::vector<Real>& v );

    std::vector<Real> timeDIK ( const std::vector<Real>& v );

    Real timeIIK1 ( const std::vector<Real>& v );

    Real plaIKp ( const std::vector<Real>& v );

    Real exINaCa ( const std::vector<Real>& v );

    Real pumpINaK ( const std::vector<Real>& v );

    Real pumpIpCa ( const std::vector<Real>& v );

    Real backICab ( const std::vector<Real>& v );

    Real backINab ( const std::vector<Real>& v);


    //! Display information about the model
    void showMe();

    //@}

private:
    //! Cell Geometry Parameters (4)

    Real M_ACap;
    Real M_VMyo;
    Real M_Vsr;
    Real M_Vup;

    //! Standard Ionic Concentrations (5)

    Real M_NaO;
    Real M_CaO;
    Real M_K0;
    Real M_NaIn;
    Real M_KIn;

    //! Buffering Parameters (4)

    Real M_CmdnTot;
    Real M_CsqnTot;
    Real M_KmCmdn;
    Real M_KmCsqn;

    //! Membrane Current Parameters (24)

    Real M_Cm;
    Real M_F;
    Real M_T;
    Real M_R;
    Real M_GNa;
    Real M_GKp;
    Real M_kNaCa;
    Real M_KmNa;
    Real M_KmCa;
    Real M_kSat;
    Real M_eta;
    Real M_INaK;
    Real M_KmNai;
    Real M_KmK0;
    Real M_GK1;
    Real M_GKr;
    Real M_GKs;
    Real M_Gt0;
    Real M_IpCa;
    Real M_KmPCa;
    Real M_GCab;
    Real M_GNab;
    Real M_KmK1;
    Real M_KmfCa;

    //! SR Parameters (1)

    Real M_KmUp;

    //! L-type Ca2+ Channel Parameters (6)

    Real M_PCa;
    Real M_PCaK;
    Real M_Prel;
    Real M_Pleak;
    Real M_ICaHalf;
    Real M_KmpCa;


    //@}



}; // class IonicFox


}

#endif // IonicFox_HPP_INCLUDED
