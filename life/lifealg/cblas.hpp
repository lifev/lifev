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
    @file
    @brief Fortran BLAS function with undersore handling

    @author Keita Teranishi
    @author Jeff Horner
    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @date 2004

    @contributor Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @maintainer Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>

 */

#ifndef CBLAS_F77_H
#define CBLAS_F77_H

#include <life/lifecore/FortranWrapper.hpp>

#ifdef __cplusplus
#include <cstdio>
extern "C"
{
#endif
    int F77NAME(xerbla)(char*, int*);
    /*
     * Level 1 Fortran Prototypes
     */

    /* Single Precision */

    void F77NAME(srot)(const int*, float *, const int*, float *, const int*, const float *, const float *);
    void F77NAME(srotg)(float *,float *,float *,float *);
    void F77NAME(srotm)( const int*, float *, const int*, float *, const int*, const float *);
    void F77NAME(srotmg)(float *,float *,float *,const float *, float *);
    void F77NAME(sswap)( const int*, float *, const int*, float *, const int*);
    void F77NAME(scopy)( const int*, const float *, const int*, float *, const int*);
    void F77NAME(saxpy)( const int*, const float *, const float *, const int*, float *, const int*);
    void F77NAME(sdot_sub)(const int*, const float *, const int*, const float *, const int*, float *);
    void F77NAME(sdsdot_sub)( const int*, const float *, const float *, const int*, const float *, const int*, float *);
    void F77NAME(sscal)( const int*, const float *, float *, const int*);
    void F77NAME(snrm2_sub)( const int*, const float *, const int*, float *);
    void F77NAME(sasum_sub)( const int*, const float *, const int*, float *);
    void F77NAME(isamax_sub)( const int*, const float * , const int*, const int*);

    /* Double Precision */

    void F77NAME(drot)(const int*, double *, const int*, double *, const int*, const double *, const double *);
    void F77NAME(drotg)(double *,double *,double *,double *);
    void F77NAME(drotm)( const int*, double *, const int*, double *, const int*, const double *);
    void F77NAME(drotmg)(double *,double *,double *,const double *, double *);
    void F77NAME(dswap)( const int*, double *, const int*, double *, const int*);
    void F77NAME(dcopy)( const int*, const double *, const int*, double *, const int*);
    void F77NAME(daxpy)( const int*, const double *, const double *, const int*, double *, const int*);
    void F77NAME(dswap)( const int*, double *, const int*, double *, const int*);
    void F77NAME(dsdot_sub)(const int*, const float *, const int*, const float *, const int*, double *);
    void F77NAME(ddot_sub)( const int*, const double *, const int*, const double *, const int*, double *);
    void F77NAME(dscal)( const int*, const double *, double *, const int*);
    void F77NAME(dnrm2_sub)( const int*, const double *, const int*, double *);
    void F77NAME(dasum_sub)( const int*, const double *, const int*, double *);
    void F77NAME(idamax_sub)( const int*, const double * , const int*, const int*);

    /* Single Complex Precision */

    void F77NAME(cswap)( const int*, void *, const int*, void *, const int*);
    void F77NAME(ccopy)( const int*, const void *, const int*, void *, const int*);
    void F77NAME(caxpy)( const int*, const void *, const void *, const int*, void *, const int*);
    void F77NAME(cswap)( const int*, void *, const int*, void *, const int*);
    void F77NAME(cdotc_sub)( const int*, const void *, const int*, const void *, const int*, void *);
    void F77NAME(cdotu_sub)( const int*, const void *, const int*, const void *, const int*, void *);
    void F77NAME(cscal)( const int*, const void *, void *, const int*);
    void F77NAME(icamax_sub)( const int*, const void *, const int*, const int*);
    void F77NAME(csscal)( const int*, const float *, void *, const int*);
    void F77NAME(scnrm2_sub)( const int*, const void *, const int*, float *);
    void F77NAME(scasum_sub)( const int*, const void *, const int*, float *);

    /* Double Complex Precision */

    void F77NAME(zswap)( const int*, void *, const int*, void *, const int*);
    void F77NAME(zcopy)( const int*, const void *, const int*, void *, const int*);
    void F77NAME(zaxpy)( const int*, const void *, const void *, const int*, void *, const int*);
    void F77NAME(zswap)( const int*, void *, const int*, void *, const int*);
    void F77NAME(zdotc_sub)( const int*, const void *, const int*, const void *, const int*, void *);
    void F77NAME(zdotu_sub)( const int*, const void *, const int*, const void *, const int*, void *);
    void F77NAME(zdscal)( const int*, const double *, void *, const int*);
    void F77NAME(zscal)( const int*, const void *, void *, const int*);
    void F77NAME(dznrm2_sub)( const int*, const void *, const int*, double *);
    void F77NAME(dzasum_sub)( const int*, const void *, const int*, double *);
    void F77NAME(izamax_sub)( const int*, const void *, const int*, const int*);

    /*
     * Level 2 Fortran Prototypes
     */

    /* Single Precision */

    void F77NAME(sgemv)(char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(sgbmv)(char*, const int*, const int*, const int*, const int*, const float *,  const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(ssymv)(char*, const int*, const float *, const float *, const int*, const float *,  const int*, const float *, float *, const int*);
    void F77NAME(ssbmv)(char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(sspmv)(char*, const int*, const float *, const float *, const float *, const int*, const float *, float *, const int*);
    void F77NAME(strmv)( char*, char*, char*, const int*, const float *, const int*, float *, const int*);
    void F77NAME(stbmv)( char*, char*, char*, const int*, const int*, const float *, const int*, float *, const int*);
    void F77NAME(strsv)( char*, char*, char*, const int*, const float *, const int*, float *, const int*);
    void F77NAME(stbsv)( char*, char*, char*, const int*, const int*, const float *, const int*, float *, const int*);
    void F77NAME(stpmv)( char*, char*, char*, const int*, const float *, float *, const int*);
    void F77NAME(stpsv)( char*, char*, char*, const int*, const float *, float *, const int*);
    void F77NAME(sger)( const int*, const int*, const float *, const float *, const int*, const float *, const int*, float *, const int*);
    void F77NAME(ssyr)(char*, const int*, const float *, const float *, const int*, float *, const int*);
    void F77NAME(sspr)(char*, const int*, const float *, const float *, const int*, float *);
    void F77NAME(sspr2)(char*, const int*, const float *, const float *, const int*, const float *, const int*,  float *);
    void F77NAME(ssyr2)(char*, const int*, const float *, const float *, const int*, const float *, const int*,  float *, const int*);

    /* Double Precision */

    void F77NAME(dgemv)(char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(dgbmv)(char*, const int*, const int*, const int*, const int*, const double *,  const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(dsymv)(char*, const int*, const double *, const double *, const int*, const double *,  const int*, const double *, double *, const int*);
    void F77NAME(dsbmv)(char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(dspmv)(char*, const int*, const double *, const double *, const double *, const int*, const double *, double *, const int*);
    void F77NAME(dtrmv)( char*, char*, char*, const int*, const double *, const int*, double *, const int*);
    void F77NAME(dtbmv)( char*, char*, char*, const int*, const int*, const double *, const int*, double *, const int*);
    void F77NAME(dtrsv)( char*, char*, char*, const int*, const double *, const int*, double *, const int*);
    void F77NAME(dtbsv)( char*, char*, char*, const int*, const int*, const double *, const int*, double *, const int*);
    void F77NAME(dtpmv)( char*, char*, char*, const int*, const double *, double *, const int*);
    void F77NAME(dtpsv)( char*, char*, char*, const int*, const double *, double *, const int*);
    void F77NAME(dger)( const int*, const int*, const double *, const double *, const int*, const double *, const int*, double *, const int*);
    void F77NAME(dsyr)(char*, const int*, const double *, const double *, const int*, double *, const int*);
    void F77NAME(dspr)(char*, const int*, const double *, const double *, const int*, double *);
    void F77NAME(dspr2)(char*, const int*, const double *, const double *, const int*, const double *, const int*,  double *);
    void F77NAME(dsyr2)(char*, const int*, const double *, const double *, const int*, const double *, const int*,  double *, const int*);

    /* Single Complex Precision */

    void F77NAME(cgemv)(char*, const int*, const int*, const void *, const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(cgbmv)(char*, const int*, const int*, const int*, const int*, const void *,  const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(chemv)(char*, const int*, const void *, const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(chbmv)(char*, const int*, const int*, const void *, const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(chpmv)(char*, const int*, const void *, const void *, const void *, const int*, const void *, void *, const int*);
    void F77NAME(ctrmv)( char*, char*, char*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ctbmv)( char*, char*, char*, const int*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ctpmv)( char*, char*, char*, const int*, const void *, void *, const int*);
    void F77NAME(ctrsv)( char*, char*, char*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ctbsv)( char*, char*, char*, const int*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ctpsv)( char*, char*, char*, const int*, const void *, void *,const int*);
    void F77NAME(cgerc)( const int*, const int*, const void *, const void *, const int*, const void *, const int*, void *, const int*);
    void F77NAME(cgeru)( const int*, const int*, const void *, const void *, const int*, const void *, const int*, void *,  const int*);
    void F77NAME(cher)(char*, const int*, const float *, const void *, const int*, void *, const int*);
    void F77NAME(cher2)(char*, const int*, const void *, const void *, const int*, const void *, const int*, void *, const int*);
    void F77NAME(chpr)(char*, const int*, const float *, const void *, const int*, void *);
    void F77NAME(chpr2)(char*, const int*, const float *, const void *, const int*, const void *, const int*, void *);

    /* Double Complex Precision */

    void F77NAME(zgemv)(char*, const int*, const int*, const void *, const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(zgbmv)(char*, const int*, const int*, const int*, const int*, const void *,  const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(zhemv)(char*, const int*, const void *, const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(zhbmv)(char*, const int*, const int*, const void *, const void *, const int*, const void *, const int*, const void *, void *, const int*);
    void F77NAME(zhpmv)(char*, const int*, const void *, const void *, const void *, const int*, const void *, void *, const int*);
    void F77NAME(ztrmv)( char*, char*, char*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ztbmv)( char*, char*, char*, const int*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ztpmv)( char*, char*, char*, const int*, const void *, void *, const int*);
    void F77NAME(ztrsv)( char*, char*, char*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ztbsv)( char*, char*, char*, const int*, const int*, const void *, const int*, void *, const int*);
    void F77NAME(ztpsv)( char*, char*, char*, const int*, const void *, void *,const int*);
    void F77NAME(zgerc)( const int*, const int*, const void *, const void *, const int*, const void *, const int*, void *, const int*);
    void F77NAME(zgeru)( const int*, const int*, const void *, const void *, const int*, const void *, const int*, void *,  const int*);
    void F77NAME(zher)(char*, const int*, const double *, const void *, const int*, void *, const int*);
    void F77NAME(zher2)(char*, const int*, const void *, const void *, const int*, const void *, const int*, void *, const int*);
    void F77NAME(zhpr)(char*, const int*, const double *, const void *, const int*, void *);
    void F77NAME(zhpr2)(char*, const int*, const double *, const void *, const int*, const void *, const int*, void *);

    /*
     * Level 3 Fortran Prototypes
     */

    /* Single Precision */

    void F77NAME(sgemm)(char*, char*, const int*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(ssymm)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(ssyrk)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, float *, const int*);
    void F77NAME(ssyr2k)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(strmm)(char*, char*, char*, char*, const int*, const int*, const float *, const float *, const int*, float *, const int*);
    void F77NAME(strsm)(char*, char*, char*, char*, const int*, const int*, const float *, const float *, const int*, float *, const int*);

    /* Double Precision */

    void F77NAME(dgemm)(char*, char*, const int*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(dsymm)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);

    /**
      blas routine declaration for the symmetric rank k operations  :
      C := alpha*A*A' + beta*C,
    */
    void F77NAME(dsyrk)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, double *, const int*);

    void F77NAME(dsyr2k)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(dtrmm)(char*, char*, char*, char*, const int*, const int*, const double *, const double *, const int*, double *, const int*);
    void F77NAME(dtrsm)(char*, char*, char*, char*, const int*, const int*, const double *, const double *, const int*, double *, const int*);

    /* Single Complex Precision */

    void F77NAME(cgemm)(char*, char*, const int*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(csymm)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(chemm)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(csyrk)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, float *, const int*);
    void F77NAME(cherk)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, float *, const int*);
    void F77NAME(csyr2k)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(cher2k)(char*, char*, const int*, const int*, const float *, const float *, const int*, const float *, const int*, const float *, float *, const int*);
    void F77NAME(ctrmm)(char*, char*, char*, char*, const int*, const int*, const float *, const float *, const int*, float *, const int*);
    void F77NAME(ctrsm)(char*, char*, char*, char*, const int*, const int*, const float *, const float *, const int*, float *, const int*);

    /* Double Complex Precision */

    void F77NAME(zgemm)(char*, char*, const int*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(zsymm)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(zhemm)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(zsyrk)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, double *, const int*);
    void F77NAME(zherk)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, double *, const int*);
    void F77NAME(zsyr2k)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(zher2k)(char*, char*, const int*, const int*, const double *, const double *, const int*, const double *, const int*, const double *, double *, const int*);
    void F77NAME(ztrmm)(char*, char*, char*, char*, const int*, const int*, const double *, const double *, const int*, double *, const int*);
    void F77NAME(ztrsm)(char*, char*, char*, char*, const int*, const int*, const double *, const double *, const int*, double *, const int*);
}
#endif /*  CBLAS_F77_H */
