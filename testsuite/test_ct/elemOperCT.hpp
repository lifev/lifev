#ifndef __ELEMOPER_CT_HH
#define __ELEMOPER_CT_HH

#include <life/lifecore/life.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/currentFE.hpp>

namespace LifeV
{

void source_pdivv(Real alpha, ElemVec& pLoc, ElemVec& elvec, 
                  const CurrentFE& fe_p, const CurrentFE& fe_u, const int iblock)
{
    int i, j, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;

    for (i=0; i < fe_u.nbNode; i++)
    {
        s = 0;
	for (iq = 0; iq < fe_u.nbQuadPt; ++iq)
	    for (j = 0; j < fe_p.nbNode; ++j)
	        s += pLoc[j]*fe_p.phi(j,iq)*fe_u.phiDer(i,iblock,iq)*fe_u.weightDet(iq);
        vec( i ) += s*alpha; 
    }
}

}

#endif
