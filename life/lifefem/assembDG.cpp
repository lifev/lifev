//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief A short description of the file content

    @author Simone Deparis <simone.deparis@epfl.ch>
    @date 09 Mar 2010

    A more detailed description of the file (if necessary)
 */

#include <assembDG.hpp>

namespace LifeV {

void
compute_vec_DG_BF(const BCHandler& BCh, ElemVec& bfvec, const CurrentBFDG& bfDG, int iblock)
{

    int i, ig, icoor;

    ElemVec::vector_view vec = bfvec.block(iblock);

    Real s, x, y, z;

    const BCBase& CurrBC = BCh.GetBCWithFlag(bfDG.marker);
    if(bfDG.bcType == Essential){
        //Dirichlet boundary conditions
        for(i = 0; i < bfDG.nbNodeAd; i++){
            s = 0.;
            for(ig = 0; ig < bfDG.nbQuadPt; ig++){
                bfDG.coorQuadPt(x, y, z, ig);

                for(icoor = 0; icoor < bfDG.nbCoorAd; icoor++){
                    s += - (bfDG.phiDerAd(i, icoor, ig) * bfDG.normal(icoor, ig))
                        * CurrBC(0., x, y, z, iblock) * bfDG.weightMeas(ig);
                }
            }
            vec(i) += s;
        }
    }else{
        //Von Neumann boundary conditions
        for(i = 0; i < bfDG.nbNodeAd; i++){
            s = 0.;
            for(ig = 0; ig < bfDG.nbQuadPt; ig++){
                bfDG.coorQuadPt(x, y, z, ig);
                s += bfDG.phiAd(i, ig) * CurrBC(0., x, y, z, iblock) * bfDG.weightMeas(ig);
            }
            vec(i) += s;
        }
    }
}

} // Namespace LifeV
