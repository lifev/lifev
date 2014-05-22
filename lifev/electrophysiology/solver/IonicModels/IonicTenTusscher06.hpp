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
  @file IonicTenTusscher06
  @brief Ionic model of ten Tusscher-Panfilov 2006

  model as in
  Ten Tusscher, K. H. W. J., and A. V. Panfilov.
  "Alternans and spiral breakup in a human ventricular."
  Am J Physiol Heart Circ Physiol 291 (2006): H1088-H1100.


  Implementation derived from the original code of TenTusscher freely available at
  http://www-binf.bio.uu.nl/khwjtuss/SourceCodes/

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

    Note that the model is not solved correctly using Rush-Larsen
    with forward Euler. In the original code of Ten Tusscher
    the solution of CaSS and CaSR is computed using a specific
    algorithm (which I do not know).
    Therefore I consider all the variables except for the
    potential as gating variable and use the specific
    method in their original code to make it work.
    Therefore the number of gating variable is augmented.

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICTENTUSSCHER06_H_
#define _IONICTENTUSSCHER06_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.

class IonicTenTusscher06 : public virtual ElectroIonicModel
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
    IonicTenTusscher06();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicTenTusscher06 ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicTenTusscher06 object
     */
    IonicTenTusscher06 ( const IonicTenTusscher06& model );
    //! Destructor
    virtual ~IonicTenTusscher06() {}

    //@}

    //! @name Overloads
    //@{

    IonicTenTusscher06& operator= ( const IonicTenTusscher06& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline  Real  getknak() const
    {
        return knak;
    }
    inline  Real  getKmNa() const
    {
        return KmNa;
    }
    inline  Real  getKmK() const
    {
        return KmK;
    }
    inline  Real  getknaca() const
    {
        return knaca;
    }
    inline  Real  getKmNai() const
    {
        return KmNai;
    }
    inline  Real  getKmCa() const
    {
        return KmCa;
    }
    inline  Real  getksat() const
    {
        return ksat;
    }
    inline  Real  getn() const
    {
        return n;
    }


    inline  Real  getKo() const
    {
        return Ko;
    }
    inline  Real  getCao() const
    {
        return Cao;
    }
    inline  Real  getNao() const
    {
        return Nao;
    }


    inline  Real  getBufc() const
    {
        return Bufc;
    }
    inline  Real  getKbufc() const
    {
        return Kbufc;
    }
    inline  Real  getBufsr() const
    {
        return Bufsr;
    }
    inline  Real  getKbufsr() const
    {
        return Kbufsr;
    }
    inline  Real  getBufss() const
    {
        return Bufss;
    }
    inline  Real  getKbufss() const
    {
        return Kbufss;
    }

    inline  Real  getVmaxup() const
    {
        return Vmaxup;
    }
    inline  Real  getKup() const
    {
        return Kup;
    }
    inline  Real  getVrel() const
    {
        return Vrel;
    }
    inline  Real  getk1_() const
    {
        return k1_;
    }
    inline  Real  getk2_() const
    {
        return k2_;
    }
    inline  Real  getk3() const
    {
        return k3;
    }
    inline  Real  getk4() const
    {
        return k4;
    }
    inline  Real  getEC() const
    {
        return EC;
    }
    inline  Real  getmaxsr() const
    {
        return maxsr;
    }
    inline  Real  getminsr() const
    {
        return minsr;
    }
    inline  Real  getVleak() const
    {
        return Vleak;
    }
    inline  Real  getVxfer() const
    {
        return Vxfer;
    }


    inline  Real  getpKNa() const
    {
        return pKNa;
    }


    inline  Real  getRTONF() const
    {
        return RTONF;
    }
    inline  Real  getF() const
    {
        return F;
    }
    inline  Real  getR() const
    {
        return R;
    }
    inline  Real  getT() const
    {
        return T;
    }


    inline  Real  getGkr() const
    {
        return Gkr;
    }
    inline  Real  getGks() const
    {
        return Gks;
    }
    inline  Real  getGK1() const
    {
        return GK1;
    }
    inline  Real  getGto() const
    {
        return Gto;
    }
    inline  Real  getGNa() const
    {
        return GNa;
    }
    inline  Real  getGbNa() const
    {
        return GbNa;
    }
    inline  Real  getGCaL() const
    {
        return GCaL;
    }
    inline  Real  getGbCa() const
    {
        return GbCa;
    }
    inline  Real  getGpCa() const
    {
        return GpCa;
    }
    inline  Real  getKpCa() const
    {
        return KpCa;
    }
    inline  Real  getGpK() const
    {
        return GpK;
    }


    inline  Real  getVc() const
    {
        return Vc;
    }
    inline  Real  getVsr() const
    {
        return Vsr;
    }
    inline  Real  getVss() const
    {
        return Vss;
    }

    inline void setknak ( const Real p )
    {
        knak = p;
    }
    inline void setKmNa ( const Real p )
    {
        KmNa = p;
    }
    inline void setKmK ( const Real p )
    {
        KmK = p;
    }
    inline void setknaca ( const Real p )
    {
        knaca = p;
    }
    inline void setKmNai ( const Real p )
    {
        KmNai = p;
    }
    inline void setKmCa ( const Real p )
    {
        KmCa = p;
    }
    inline void setksat ( const Real p )
    {
        ksat = p;
    }
    inline void setn ( const Real p )
    {
        n = p;
    }


    inline void setKo ( const Real p )
    {
        Ko = p;
    }
    inline void setCao ( const Real p )
    {
        Cao = p;
    }
    inline void setNao ( const Real p )
    {
        Nao = p;
    }


    inline void setBufc ( const Real p )
    {
        Bufc = p;
    }
    inline void setKbufc ( const Real p )
    {
        Kbufc = p;
    }
    inline void setBufsr ( const Real p )
    {
        Bufsr = p;
    }
    inline void setKbufsr ( const Real p )
    {
        Kbufsr = p;
    }
    inline void setBufss ( const Real p )
    {
        Bufss = p;
    }
    inline void setKbufss ( const Real p )
    {
        Kbufss = p;
    }

    inline void setVmaxup ( const Real p )
    {
        Vmaxup = p;
    }
    inline void setKup ( const Real p )
    {
        Kup = p;
    }
    inline void setVrel ( const Real p )
    {
        Vrel = p;
    }
    inline void setk1_ ( const Real p )
    {
        k1_ = p;
    }
    inline void setk2_ ( const Real p )
    {
        k2_ = p;
    }
    inline void setk3 ( const Real p )
    {
        k3 = p;
    }
    inline void setk4 ( const Real p )
    {
        k4 = p;
    }
    inline void setEC ( const Real p )
    {
        EC = p;
    }
    inline void setmaxsr ( const Real p )
    {
        maxsr = p;
    }
    inline void setminsr ( const Real p )
    {
        minsr = p;
    }
    inline void setVleak ( const Real p )
    {
        Vleak = p;
    }
    inline void setVxfer ( const Real p )
    {
        Vxfer = p;
    }


    inline void setpKNa ( const Real p )
    {
        pKNa = p;
    }


    inline void setF ( const Real p )
    {
        F = p;
        computeRTONF();
    }
    inline void setR ( const Real p )
    {
        R = p;
        computeRTONF();
    }
    inline void setT ( const Real p )
    {
        T = p;
        computeRTONF();
    }

    inline void computeRTONF()
    {
        RTONF = R * T / F;
    }

    inline void setGkr ( const Real p )
    {
        Gkr = p;
    }
    inline void setGks ( const Real p )
    {
        Gks = p;
    }
    inline void setGK1 ( const Real p )
    {
        GK1 = p;
    }
    inline void setGto ( const Real p )
    {
        Gto = p;
    }
    inline void setGNa ( const Real p )
    {
        GNa = p;
    }
    inline void setGbNa ( const Real p )
    {
        GbNa = p;
    }
    inline void setGCaL ( const Real p )
    {
        GCaL = p;
    }
    inline void setGbCa ( const Real p )
    {
        GbCa = p;
    }
    inline void setGpCa ( const Real p )
    {
        GpCa = p;
    }
    inline void setKpCa ( const Real p )
    {
        KpCa = p;
    }
    inline void setGpK ( const Real p )
    {
        GpK = p;
    }


    inline void setVc ( const Real p )
    {
        Vc = p;
    }
    inline void setVsr ( const Real p )
    {
        Vsr = p;
    }
    inline void setVss ( const Real p )
    {
        Vss = p;
    }


    inline void computeinverseVcF2()
    {
        inverseVcF2 = 1. / ( 2. * Vc * F ) ;
    }
    inline void computeinverseVcF()
    {
        inverseVcF = 1. / ( Vc * F );
    }
    inline void computeinversevssF2()
    {
        inversevssF2 = 1. / ( 2. * Vss * F );
    }


    //inline const short int& Size() const { return M_numberOfEquations; }

    inline Real Ek ( Real Ki)
    {
        return RTONF * std::log ( ( Ko / Ki ) );
    }
    inline Real Ena (Real Nai)
    {
        return RTONF * std::log ( ( Nao / Nai ) );
    }
    inline Real Eks (Real Ki, Real Nai)
    {
        return RTONF * std::log ( ( Ko + pKNa * Nao ) / ( Ki + pKNa * Nai ) );
    }
    inline Real Eca (Real Cai)
    {
        return 0.5 * RTONF * std::log ( ( Cao / Cai ) );
    }
    inline Real Ak1 (Real V, Real Ki)
    {
        return 0.1 / ( 1. + std::exp ( 0.06 * ( V - Ek (Ki) - 200. ) ) );
    }
    inline Real Bk1 (Real V, Real Ki)
    {
        return ( 3. * std::exp ( 0.0002 * ( V - Ek (Ki) + 100. ) ) +
                 std::exp ( 0.1 * ( V - Ek (Ki) - 10. ) ) )
               / ( 1. + std::exp ( -0.5 * ( V - Ek (Ki) ) ) );
    }
    inline Real rec_iK1 ( Real V, Real Ki )
    {
        return Ak1 (V, Ki) / ( Ak1 (V, Ki) + Bk1 (V, Ki) );
    }
    inline Real rec_iNaK (Real V)
    {
        return ( 1. / ( 1. + 0.1245 * std::exp ( -0.1 * V * F / ( R * T ) )
                        + 0.0353 * std::exp (- V * F / ( R * T ) ) ) );
    }
    inline Real rec_ipK (Real V)
    {
        return 1. / ( 1. + std::exp ( ( 25. - V ) / 5.98 ) );
    }

    //@}

    //! @name Methods
    //@{



    inline Real INa (Real V, Real m, Real h, Real j, Real Nai)
    {
        return GNa * m * m * m * h * j * ( V - Ena (Nai) );
    }
    inline Real ICaL (Real V, Real d, Real f, Real f2, Real fcass, Real CaSS)
    {
        return GCaL * d * f *  f2 * fcass * 4. * ( V - 15. ) * ( F * F / ( R * T ) ) *
               ( 0.25 * std::exp ( 2. * ( V - 15 ) * F / ( R * T ) ) * CaSS - Cao )
               / ( std::exp ( 2. * ( V - 15. ) * F / ( R * T ) ) - 1.);

    }
    inline Real Ito (Real V, Real r, Real s, Real Ki)
    {
        return Gto * r * s * ( V - Ek (Ki) );
    }
    inline Real IKr (Real V, Real xr1, Real xr2, Real Ki)
    {
        return Gkr * std::sqrt (Ko / 5.4) * xr1 * xr2 * ( V - Ek (Ki) );
    }
    inline Real IKs (Real V, Real xs, Real Ki, Real Nai)
    {
        return Gks * xs * xs * ( V - Eks (Ki, Nai) );

        std::cout << "\n****  Iks   ******* ";
        std::cout << "\nGks: " << Gks;
        std::cout << "\nxs: " << xs;
        std::cout << "\nV: " << V;
        std::cout << "\nKi: " << Ki;
        std::cout << "\nNai: " << Nai;
        std::cout << "\nEks: " << Eks (Ki, Nai);
        std::cout << "\n*********** ";


    }
    inline Real IK1 (Real V, Real Ki)
    {
        return GK1 * rec_iK1 (V, Ki) * ( V - Ek (Ki) );
    }
    inline Real INaCa (Real V, Real Nai, Real Cai)
    {
        return knaca * ( 1. / ( KmNai * KmNai * KmNai + Nao * Nao * Nao ) ) * ( 1. / ( KmCa + Cao ) ) *
               ( 1. / ( 1. + ksat * std::exp ( ( n - 1. ) * V * F / ( R * T ) ) ) ) *
               (std::exp ( n * V * F / ( R * T ) ) * Nai * Nai * Nai * Cao -
                std::exp ( ( n - 1. ) * V * F / ( R * T ) ) * Nao * Nao * Nao * Cai * 2.5);
    }
    inline Real INaK (Real V, Real Nai)
    {
        return knak * ( Ko / ( Ko + KmK ) ) * ( Nai / ( Nai + KmNa ) ) * rec_iNaK (V);
    }
    inline Real IpCa (Real Cai)
    {
        return GpCa * Cai / ( KpCa + Cai );
    }
    inline Real IpK (Real V, Real Ki)
    {
        return GpK * rec_ipK (V) * ( V - Ek (Ki) );
    }
    inline Real IbNa (Real V, Real Nai)
    {
        return GbNa * ( V - Ena (Nai) );
    }
    inline Real IbCa (Real V, Real Cai)
    {
        return GbCa * ( V - Eca (Cai) );
    }
    inline Real Itot (Real V, Real m, Real h, Real j, Real d, Real f, Real f2, Real fcass,
                      Real r, Real s, Real xr1, Real xr2, Real xs, Real Nai, Real Ki,
                      Real Cai, Real CaSS)
    {
        return IKr (V, xr1, xr2, Ki)    +
               IKs (V, xs, Ki, Nai)   +
               IK1 (V, Ki)   +
               Ito (V, r, s, Ki)   +
               INa (V, m, h, j, Nai)   +
               IbNa (V, Nai)  +
               ICaL (V, d, f, f2, fcass, CaSS)  +
               IbCa (V, Cai) +
               INaK (V, Nai)  +
               INaCa (V, Nai, Cai) +
               IpCa (Cai)  +
               IpK (V, Ki);
    }

    //update concentrations

    inline Real kCaSR (Real CaSR)
    {
        return maxsr - ( ( maxsr - minsr ) / ( 1. + ( EC / CaSR ) * ( EC / CaSR ) ) );
    }
    inline Real k1 (Real CaSR)
    {
        return k1_ / kCaSR (CaSR);
    }
    inline Real k2 (Real CaSR)
    {
        return k2_ * kCaSR (CaSR);
    }
    inline Real dRR (Real CaSR, Real CaSS, Real RR)
    {
        return k4 * ( 1 - RR ) - k2 (CaSR) * CaSS * RR;
    }
    inline Real solveRR (Real CaSR, Real CaSS, Real RR, Real HT)
    {
        return  RR + HT * dRR (CaSR, CaSS, RR);
    }
    inline Real OO (Real CaSR, Real CaSS, Real RR)
    {
        return k1 (CaSR) * CaSS * CaSS * RR / ( k3 + k1 (CaSR) * CaSS * CaSS );
    }

    inline Real Irel (Real CaSR, Real CaSS, Real RR)
    {
        return Vrel * OO (CaSR, CaSS, RR) * ( CaSR - CaSS );
    }
    inline Real Ileak (Real CaSR, Real Cai)
    {
        return Vleak * ( CaSR - Cai);
    }
    inline Real Iup (Real Cai)
    {
        return Vmaxup / ( 1. + ( ( Kup * Kup ) / ( Cai * Cai ) ) );
    }
    inline Real Ixfer (Real CaSS, Real Cai)
    {
        return Vxfer * ( CaSS - Cai);
    }
    inline Real CaCSQN (Real CaSR)
    {
        return Bufsr * CaSR / ( CaSR + Kbufsr );
    }
    inline Real dCaSR (Real Cai, Real CaSR, Real CaSS, Real RR, Real HT = 1.0)
    {
        return HT * ( Iup (Cai) - Irel (CaSR, CaSS, RR) - Ileak (CaSR, Cai) );
    }
    inline Real bjsr (Real Cai, Real CaSR, Real CaSS, Real RR, Real HT)
    {
        return Bufsr - CaCSQN (CaSR) - dCaSR (Cai, CaSR, CaSS, RR, HT) - CaSR + Kbufsr;
    }
    inline Real cjsr (Real Cai, Real CaSR, Real CaSS, Real RR, Real HT)
    {
        return Kbufsr * ( CaCSQN (CaSR) + dCaSR (Cai, CaSR, CaSS, RR, HT) + CaSR);
    }
    inline Real solveCaSR (Real Cai, Real CaSR, Real CaSS, Real RR, Real HT)
    {
        return ( std::sqrt ( bjsr (Cai, CaSR, CaSS, RR, HT) * bjsr (Cai, CaSR, CaSS, RR, HT)
                             + 4. * cjsr (Cai, CaSR, CaSS, RR, HT) )
                 - bjsr (Cai, CaSR, CaSS, RR, HT) ) / 2.;
    }

    inline Real CaSSBuf (Real CaSS)
    {
        return Bufss * CaSS / ( CaSS + Kbufss);
    }
    inline Real dCaSS (Real Cai, Real CaSR, Real CaSS, Real RR, Real V, Real d, Real f, Real f2, Real fcass, Real HT = 1.0)
    {
        return  HT * ( - Ixfer (CaSS, Cai) * ( Vc / Vss )
                       + Irel (CaSR, CaSS, RR) * ( Vsr / Vss )
                       + ( - ICaL (V, d, f, f2, fcass, CaSS) * inversevssF2 * M_cellularCapacitance ) );
    };
    inline Real bcss (Real Cai, Real CaSR, Real CaSS, Real RR, Real V, Real d, Real f, Real f2, Real fcass, Real HT)
    {
        return Bufss - CaSSBuf (CaSS) - dCaSS (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, HT) - CaSS + Kbufss;
    }
    inline Real ccss (Real Cai, Real CaSR, Real CaSS, Real RR, Real V, Real d, Real f, Real f2, Real fcass, Real HT)
    {
        return Kbufss * ( CaSSBuf (CaSS) + dCaSS (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, HT) + CaSS);
    }
    inline Real solveCaSS (Real Cai, Real CaSR, Real CaSS, Real RR, Real V, Real d, Real f, Real f2, Real fcass, Real HT)
    {
        return ( std::sqrt ( bcss (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, HT) * bcss (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, HT)
                             + 4. * ccss (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, HT) )
                 - bcss (Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, HT) ) / 2.;
    }

    inline Real CaBuf (Real Cai)
    {
        return Bufc * Cai / ( Cai + Kbufc);
    }
    inline Real dCai (Real V, Real Nai, Real Cai, Real CaSR, Real CaSS, Real HT = 1.0)
    {
        return HT * ( ( - ( IbCa (V, Cai) + IpCa (Cai) - 2. * INaCa (V, Nai, Cai) )
                        * inverseVcF2 * M_cellularCapacitance )
                      - ( Iup (Cai) - Ileak (CaSR, Cai) ) * ( Vsr / Vc )
                      + Ixfer (CaSS, Cai) );
    }
    inline Real bc (Real V, Real Nai, Real Cai, Real CaSR, Real CaSS, Real HT)
    {
        return Bufc - CaBuf (Cai) - dCai (V, Nai, Cai, CaSR, CaSS, HT) - Cai + Kbufc;
    }
    inline Real cc (Real V, Real Nai, Real Cai, Real CaSR, Real CaSS, Real HT)
    {
        return Kbufc * ( CaBuf (Cai) + dCai (V, Nai, Cai, CaSR, CaSS, HT) + Cai);
    }
    inline Real solveCai (Real V, Real Nai, Real Cai, Real CaSR, Real CaSS, Real HT)
    {
        return ( std::sqrt ( bc (V, Nai, Cai, CaSR, CaSS, HT) * bc (V, Nai, Cai, CaSR, CaSS, HT)
                             + 4. * cc (V, Nai, Cai, CaSR, CaSS, HT) )
                 - bc (V, Nai, Cai, CaSR, CaSS, HT) ) / 2.;
    }


    inline Real dNai (Real V, Real m, Real h, Real j, Real Nai, Real Cai)
    {
        return - ( INa (V, m, h, j, Nai)
                   + IbNa (V, Nai)
                   + 3. * INaK (V, Nai)
                   + 3. * INaCa (V, Nai, Cai) ) * inverseVcF * M_cellularCapacitance;
    }
    inline Real solveNai (Real V, Real m, Real h, Real j, Real Nai, Real Cai, Real HT)
    {
        return Nai + HT * dNai (V, m, h, j, Nai, Cai);
    }

    inline Real dKi (Real V, Real r, Real s, Real xr1, Real xr2, Real xs, Real Nai, Real Ki)
    {
        return - (- M_appliedCurrent
                  + IK1 (V, Ki)
                  + Ito (V, r, s, Ki)
                  + IKr (V, xr1, xr2, Ki)
                  + IKs (V, xs, Ki, Nai)
                  - 2. * INaK (V, Nai)
                  + IpK (V, Ki) ) * inverseVcF * M_cellularCapacitance;
    }
    inline Real solveKi (Real V, Real r, Real s, Real xr1, Real xr2, Real xs, Real Nai, Real Ki, Real HT)
    {
        return Ki + HT * dKi (V, r, s, xr1, xr2, xs, Nai, Ki);
    }


    inline Real AM (Real V)
    {
        return 1. / ( 1. + std::exp ( ( -60. -  V ) / 5. ) );
    }
    inline Real BM (Real V)
    {
        return 0.1 / ( 1. + std::exp ( ( V + 35. ) / 5. ) )
               + 0.10 / ( 1. + std::exp ( ( V - 50. ) / 200. ) );
    }
    inline Real TAU_M (Real V)
    {
        return AM (V) * BM (V);
    }
    inline Real M_INF (Real V)
    {
        return 1. / ( ( 1. + std::exp ( ( -56.86 - V ) / 9.03 ) )
                      * ( 1. + std::exp ( ( - 56.86 - V ) / 9.03 ) ) );
    }
    inline Real AH (Real V)
    {
        if ( V >= - 40.)
        {
            return 0.0;
        }
        else
        {
            return   ( 0.057 * std::exp ( - ( V + 80. ) / 6.8 ) );
        }
    }
    inline Real BH (Real V)
    {
        if ( V >= - 40.)
        {
            return ( 0.77 / ( 0.13 * ( 1. + std::exp ( - ( V + 10.66 ) / 11.1 ) ) ) );
        }
        else
        {
            return ( 2.7 * std::exp ( 0.079 * V ) + ( 3.1e5 ) * std::exp ( 0.3485 * V ) );
        }
    }
    inline Real TAU_H (Real V)
    {
        return 1.0 / (AH (V) + BH (V) );
    }

    inline Real H_INF (Real V)
    {
        return 1. / ( ( 1. + std::exp ( ( V + 71.55 ) / 7.43 ) )
                      * ( 1. + std::exp ( ( V + 71.55 ) / 7.43 ) ) );
    }

    inline Real AJ (Real V)
    {
        if ( V >= - 40.)
        {
            return 0.0;
        }
        else
            return ( ( ( -2.5428e4 ) * std::exp ( 0.2444 * V ) - ( 6.948e-6 ) *
                       std::exp ( -0.04391 * V ) ) * ( V + 37.78 ) /
                     ( 1. + std::exp ( 0.311 * ( V + 79.23 ) ) ) );

    }
    inline Real BJ (Real V)
    {
        if ( V >= - 40.)
        {
            return ( 0.6 * std::exp ( ( 0.057 ) * V ) / ( 1. + std::exp ( -0.1 * ( V + 32. ) ) ) );
        }
        else
            return ( 0.02424 * std::exp ( -0.01052 * V )
                     / ( 1. + std::exp ( -0.1378 * ( V + 40.14 ) ) ) );
    }
    inline Real TAU_J (Real V)
    {
        return 1.0 / ( AJ (V) + BJ (V) );
    }
    inline Real J_INF (Real V)
    {
        return H_INF (V);
    }

    inline Real Xr1_INF (Real V)
    {
        return 1. / ( 1. + std::exp ( ( -26. - V ) / 7. ) );
    }
    inline Real axr1 (Real V)
    {
        return 450. / ( 1. + std::exp ( ( -45. - V ) / 10. ) );
    }
    inline Real bxr1 (Real V)
    {
        return 6. / ( 1. + std::exp ( ( V - ( -30. ) ) / 11.5 ) );
    }
    inline Real TAU_Xr1 (Real V)
    {
        return axr1 (V) * bxr1 (V);
    }
    inline Real Xr2_INF (Real V)
    {
        return 1. / ( 1. + std::exp ( ( V - (-88.) ) / 24. ) );
    }
    inline Real axr2 (Real V)
    {
        return 3. / ( 1. + std::exp ( ( -60. - V ) / 20. ) );
    }
    inline Real bxr2 (Real V)
    {
        return 1.12 / ( 1. + std::exp ( ( V - 60. ) / 20. ) );
    }
    inline Real TAU_Xr2 (Real V)
    {
        return axr2 (V) * bxr2 (V);
    }

    inline Real Xs_INF (Real V)
    {
        return 1. / ( 1. + std::exp ( ( -5. - V ) / 14. ) );
    }
    inline Real Axs (Real V)
    {
        return ( 1400. / (std::sqrt ( 1. + std::exp ( ( 5. - V ) / 6. ) ) ) );
    }
    inline Real Bxs (Real V)
    {
        return ( 1. / ( 1. + std::exp ( ( V - 35. ) / 15.) ) );
    }
    inline Real TAU_Xs (Real V)
    {
        return Axs (V) * Bxs (V) + 80.;
    }

    inline Real R_INF (Real V)
    {
        return 1. / ( 1. + std::exp ( ( 20. - V ) / 6. ) );
    }
    inline Real S_INF (Real V)
    {
        if (flag == Epi)
        {
            return 1. / ( 1. + std::exp ( ( V + 20. ) / 5. ) );
        }
        else if (flag == Endo)
        {
            return 1. / ( 1. + std::exp ( ( V + 28. ) / 5. ) );
        }
        else
        {
            return 1. / ( 1. + std::exp ( ( V + 20. ) / 5. ) );
        }
    }
    inline Real TAU_R (Real V)
    {
        return 9.5 * std::exp ( - ( V + 40. ) * ( V + 40. ) / 1800. ) + 0.8;
    }
    inline Real TAU_S (Real V)
    {
        if (flag == Epi)
            return 85. * std::exp ( - ( V + 45. ) * ( V + 45. ) / 320. )
                   + 5. / ( 1. + std::exp ( ( V - 20. ) / 5. ) ) + 3.;
        else if (flag == Endo)
        {
            return 1000. * std::exp ( - ( V + 67. ) * ( V + 67. ) / 1000. ) + 8.;
        }
        else
            return 85. * std::exp ( - ( V + 45. ) * ( V + 45. ) / 320. )
                   + 5. / ( 1. + std::exp ( ( V - 20. ) / 5. ) ) + 3.;
    }


    inline Real D_INF (Real V)
    {
        return 1. / ( 1. + std::exp ( ( - 8. - V ) / 7.5 ) );
    }
    inline Real Ad (Real V)
    {
        return 1.4 / ( 1. + std::exp ( ( - 35. - V ) / 13. ) ) + 0.25;
    }
    inline Real Bd (Real V)
    {
        return 1.4 / ( 1. + std::exp ( ( V + 5. ) / 5. ) );
    }
    inline Real Cd (Real V)
    {
        return 1. / ( 1. + std::exp ( ( 50. - V ) / 20. ) );
    }
    inline Real TAU_D (Real V)
    {
        return Ad (V) * Bd (V) + Cd (V);
    }
    inline Real F_INF (Real V)
    {
        return 1. / ( 1. + std::exp ( ( V + 20. ) / 7. ) );
    }
    inline Real Af (Real V)
    {
        return 1102.5 * std::exp ( - ( V + 27. ) * ( V + 27. ) / 225. );
    }
    inline Real Bf (Real V)
    {
        return 200. / ( 1. + std::exp ( ( 13. - V ) / 10. ) );
    }
    inline Real Cf (Real V)
    {
        return ( 180. / ( 1. + std::exp ( ( V + 30. ) / 10. ) ) ) + 20.;
    }
    inline Real TAU_F (Real V)
    {
        return Af (V) + Bf (V) + Cf (V);
    }
    inline Real F2_INF (Real V)
    {
        return 0.67 / ( 1. + std::exp ( ( V + 35. ) / 7. ) ) + 0.33;
    }
    inline Real Af2 (Real V)
    {
        return 600. * std::exp ( - ( V + 25. ) * ( V + 25. ) / 170. );
    }
    inline Real Bf2 (Real V)
    {
        return 31. / ( 1. + std::exp ( ( 25. - V ) / 10. ) );
    }
    inline Real Cf2 (Real V)
    {
        return 16. / ( 1. + std::exp ( ( V + 30. ) / 10. ) );
    }
    inline Real TAU_F2 (Real V)
    {
        return Af2 (V) + Bf2 (V) + Cf2 (V);
    }
    inline Real FCaSS_INF (Real CaSS)
    {
        return 0.6 / ( 1. + ( CaSS / 0.05 ) * ( CaSS / 0.05 ) ) + 0.4;
    }
    inline Real TAU_FCaSS (Real CaSS)
    {
        return 80. / ( 1. + ( CaSS / 0.05 ) * ( CaSS / 0.05 ) ) + 2.;
    }


    inline Real solveM (Real V, Real m, Real HT)
    {
        return m + HT * (M_INF (V) - m) / TAU_M (V);
    }
    inline Real solveH (Real V, Real h, Real HT)
    {
        return h +  HT * (H_INF (V) - h) / TAU_H (V);
    }
    inline Real solveJ (Real V, Real j, Real HT)
    {
        return j +  HT * (J_INF (V) - j) / TAU_J (V);
    }
    inline Real solveXr1 (Real V, Real xr1, Real HT)
    {
        return xr1 +  HT * ( Xr1_INF (V) - xr1) / TAU_Xr1 (V);
    }
    inline Real solveXr2 (Real V, Real xr2, Real HT)
    {
        return xr2 +  HT * ( Xr2_INF (V) - xr2 ) / TAU_Xr2 (V);
    }
    inline Real solveXs (Real V, Real xs, Real HT)
    {
        return xs + HT * ( Xs_INF (V) - xs ) / TAU_Xs (V);
    }
    inline Real solveS (Real V, Real s, Real HT)
    {
        return s +  HT * ( S_INF (V) - s ) / TAU_S (V) ;
    }
    inline Real solveR (Real V, Real r, Real HT)
    {
        return r +  HT * ( R_INF (V) - r) / TAU_R (V);
    }
    inline Real solveD (Real V, Real d, Real HT)
    {
        return d +  HT * ( D_INF (V) - d ) / TAU_D (V);
    }
    inline Real solveF (Real V, Real f, Real HT)
    {
        return f + HT * (F_INF (V) - f) / TAU_F (V);
    }
    inline Real solveF2 (Real V, Real f2, Real HT)
    {
        return f2 + HT * (F2_INF (V) - f2) / TAU_F2 (V);
    }
    inline Real solveFCaSS (Real CaSS, Real fcass, Real HT)
    {
        return fcass + HT * (FCaSS_INF (CaSS) - fcass) / TAU_FCaSS (CaSS);
    }


    inline Real dM (Real V, Real m)
    {
        return (M_INF (V) - m) / TAU_M (V);
    }
    inline Real dH (Real V, Real h)
    {
        return (H_INF (V) - h) / TAU_H (V);
    }
    inline Real dJ (Real V, Real j)
    {
        return (J_INF (V) - j) / TAU_J (V);
    }
    inline Real dXr1 (Real V, Real xr1)
    {
        return ( Xr1_INF (V) - xr1) / TAU_Xr1 (V);
    }
    inline Real dXr2 (Real V, Real xr2)
    {
        return ( Xr2_INF (V) - xr2 ) / TAU_Xr2 (V);
    }
    inline Real dXs (Real V, Real xs)
    {
        return ( Xs_INF (V) - xs ) / TAU_Xs (V);
    }
    inline Real dS (Real V, Real s)
    {
        return ( S_INF (V) - s ) / TAU_S (V) ;
    }
    inline Real dR (Real V, Real r)
    {
        return ( R_INF (V) - r) / TAU_R (V);
    }
    inline Real dD (Real V, Real d)
    {
        return ( D_INF (V) - d ) / TAU_D (V);
    }
    inline Real dF (Real V, Real f)
    {
        return (F_INF (V) - f) / TAU_F (V);
    }
    inline Real dF2 (Real V, Real f2)
    {
        return (F2_INF (V) - f2) / TAU_F2 (V);
    }
    inline Real dFCaSS (Real CaSS, Real fcass)
    {
        return (FCaSS_INF (CaSS) - fcass) / TAU_FCaSS (CaSS);
    }


    //update voltage
    inline Real solveV (Real V, Real m, Real h, Real j, Real d, Real f, Real f2, Real fcass,
                        Real r, Real s, Real xr1, Real xr2, Real xs, Real Nai, Real Ki,
                        Real Cai, Real CaSS , Real HT)
    {
        return V + HT * ( M_appliedCurrent / M_membraneCapacitance + (- Itot (V, m, h, j, d, f, f2, fcass, r, s, xr1, xr2, xs, Nai, Ki, Cai, CaSS ) ) );
    }

    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeNonGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    //
    void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt );


    //! Display information about the model
    void showMe();

    void showCurrents (std::vector<Real>& v);
    void showCurrents (Real V, Real m, Real h, Real j, Real d, Real f, Real f2, Real fcass,
                       Real r, Real s, Real xr1, Real xr2, Real xs, Real Nai, Real Ki,
                       Real Cai, Real CaSS);

    void solveOneStep (std::vector<Real>& v, Real dt);
    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

    enum WallFlag
    {
        Endo, Epi, MCell
    };

    inline WallFlag getFlag()
    {
        return flag;
    }


private:
    //! Model Parameters

    //Parameters
    Real knak;
    Real KmNa;
    Real KmK;
    Real knaca;
    Real KmNai;
    Real KmCa;
    Real ksat;
    Real n;


    Real Ko;
    Real Cao;
    Real Nao;


    Real Bufc;
    Real Kbufc;
    Real Bufsr;
    Real Kbufsr;
    Real Bufss;
    Real Kbufss;

    Real Vmaxup;
    Real Kup;
    Real Vrel;
    Real k1_;
    Real k2_;
    Real k3;
    Real k4;
    Real EC;
    Real maxsr;
    Real minsr;
    Real Vleak;
    Real Vxfer;


    Real pKNa;

    Real RTONF;
    Real F;
    Real R;
    Real T;


    Real Gkr;
    Real Gks;
    Real GK1;
    Real Gto;
    Real GNa;
    Real GbNa;
    Real GCaL;
    Real GbCa;
    Real GpCa;
    Real KpCa;
    Real GpK;


    Real Vc;
    Real Vsr;
    Real Vss;

    Real inverseVcF2;
    Real inverseVcF;
    Real inversevssF2;

    WallFlag flag;

    Real M_cellularCapacitance;
    //@}

}; // class IonicTenTusscher06

inline ElectroIonicModel* createIonicTenTusscher06()
{
    return new IonicTenTusscher06();
}

namespace
{
static bool register_IonicMinimalModel = ElectroIonicModel::IonicModelFactory::instance().registerProduct ("TenTusscher06", &createIonicTenTusscher06 );
}

}

#endif
