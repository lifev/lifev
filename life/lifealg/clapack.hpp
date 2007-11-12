#ifndef __CLAPACK_H
#define __CLAPACK_H

#include <life/lifealg/cblas.hpp>

typedef int integer;
typedef unsigned long uinteger;
typedef char *address;
typedef short int shortint;
typedef float Real;
typedef double doubleReal;
typedef struct { Real r, i; } Complex;
typedef struct { doubleReal r, i; } doubleComplex;
typedef int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;

#define VOID void

#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef Real (*R_fp)(...);
typedef doubleReal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef Real (*R_fp)();
typedef doubleReal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif

#undef VOID

/*
  ====================================================
         LAPACK and BLAS routine declaration
  ====================================================
  Use of fortran_wrap for fortran declarations.

  The 3 following functions are already declared by Aztec
  with a difference  (supplementary argument(s) : strlen("L")...)
  That's why we put their declaration in a comment.

  Please note that the Aztec declarations requires C arrays for all
  the arguments (double and integer are not declared with references).
  In order to be coherent with
  this (wierd?) convention, we work only with arrays. Therefore, there
  is no difference in the declarations between a true array (such as "A")
  and a true integer (such as "LDA").

  (1) Lapack routine declaration for matrix factorization using choleski:

  SUBROUTINE_F77  F77NAME(dpotrf)
             (char * UPLO,
              I_F77 & N, R8_F77 * A, I_F77 & LDA, I_F77 & INFO);

  (2) Blas 3 routine declaration for general matrix * matrix product :
  C := alpha*op(A)*op(B) + beta*C

  SUBROUTINE_F77  F77NAME(dgemm)
             (char * TRANSA, char * TRANSB,
              I_F77 & M, I_F77 & N, I_F77 & K, R8_F77 & ALPHA,
              R8_F77 * A, I_F77 & LDA, R8_F77 * B, I_F77 & LDB,
          R8_F77 & BETA, R8_F77 * C, I_F77 & LDC);


  (3) Blas 2 routine declaration for general matrix * vector product :
  y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y.

  SUBROUTINE_F77  F77NAME(dgemv)
             (char * TRANS,
              I_F77 * M, I_F77 * N, R8_F77 * ALPHA, R8_F77 *  A, I_F77 * LDA,
              R8_F77 * X, I_F77 * INCX, R8_F77 *  BETA, R8_F77 * Y, I_F77 * INCY );

*/
#ifdef __cplusplus
#include <cstdio>
extern "C" {
#endif
/* Subroutine */ int F77NAME(cbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, Real *d__, Real *e, Complex *vt, integer *ldvt,
    Complex *u, integer *ldu, Complex *c__, integer *ldc, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, Complex *ab, integer *ldab, Real *d__,
    Real *e, Complex *q, integer *ldq, Complex *pt, integer *ldpt,
    Complex *c__, integer *ldc, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, integer *ipiv, Real *anorm, Real *rcond,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, Real *r__, Real *c__, Real *rowcnd, Real
    *colcnd, Real *amax, integer *info);

/* Subroutine */ int F77NAME(cgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, Complex *ab, integer *ldab, Complex *afb, integer *
    ldafb, integer *ipiv, Complex *b, integer *ldb, Complex *x, integer *
    ldx, Real *ferr, Real *berr, Complex *work, Real *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, Complex *ab, integer *ldab, integer *ipiv, Complex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(cgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, Complex *ab, integer *ldab, Complex *afb,
     integer *ldafb, integer *ipiv, char *equed, Real *r__, Real *c__,
    Complex *b, integer *ldb, Complex *x, integer *ldx, Real *rcond, Real
    *ferr, Real *berr, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, Complex *ab, integer *ldab, integer *ipiv, Complex
    *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, Real *scale, integer *m, Complex *v, integer *ldv,
    integer *info);

/* Subroutine */ int F77NAME(cgebal)(char *job, integer *n, Complex *a, integer *lda,
    integer *ilo, integer *ihi, Real *scale, integer *info);

/* Subroutine */ int F77NAME(cgebd2)(integer *m, integer *n, Complex *a, integer *lda,
     Real *d__, Real *e, Complex *tauq, Complex *taup, Complex *work,
    integer *info);

/* Subroutine */ int F77NAME(cgebrd)(integer *m, integer *n, Complex *a, integer *lda,
     Real *d__, Real *e, Complex *tauq, Complex *taup, Complex *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgecon)(char *norm, integer *n, Complex *a, integer *lda,
     Real *anorm, Real *rcond, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgeequ)(integer *m, integer *n, Complex *a, integer *lda,
     Real *r__, Real *c__, Real *rowcnd, Real *colcnd, Real *amax,
    integer *info);

/* Subroutine */ int F77NAME(cgees)(char *jobvs, char *sort, L_fp select, integer *n,
    Complex *a, integer *lda, integer *sdim, Complex *w, Complex *vs,
    integer *ldvs, Complex *work, integer *lwork, Real *rwork, logical *
    bwork, integer *info);

/* Subroutine */ int F77NAME(cgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, Complex *a, integer *lda, integer *sdim, Complex *
    w, Complex *vs, integer *ldvs, Real *rconde, Real *rcondv, Complex *
    work, integer *lwork, Real *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cgeev)(char *jobvl, char *jobvr, integer *n, Complex *a,
    integer *lda, Complex *w, Complex *vl, integer *ldvl, Complex *vr,
    integer *ldvr, Complex *work, integer *lwork, Real *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, Complex *a, integer *lda, Complex *w, Complex *vl,
    integer *ldvl, Complex *vr, integer *ldvr, integer *ilo, integer *ihi,
     Real *scale, Real *abnrm, Real *rconde, Real *rcondv, Complex *work,
    integer *lwork, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgegs)(char *jobvsl, char *jobvsr, integer *n, Complex *
    a, integer *lda, Complex *b, integer *ldb, Complex *alpha, Complex *
    beta, Complex *vsl, integer *ldvsl, Complex *vsr, integer *ldvsr,
    Complex *work, integer *lwork, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgegv)(char *jobvl, char *jobvr, integer *n, Complex *a,
    integer *lda, Complex *b, integer *ldb, Complex *alpha, Complex *beta,
     Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, Complex *
    work, integer *lwork, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgehd2)(integer *n, integer *ilo, integer *ihi, Complex *
    a, integer *lda, Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgehrd)(integer *n, integer *ilo, integer *ihi, Complex *
    a, integer *lda, Complex *tau, Complex *work, integer *lwork, integer
    *info);

/* Subroutine */ int F77NAME(cgelq2)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgelqf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgels)(char *trans, integer *m, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *b, integer *ldb, Complex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgelsx)(integer *m, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *jpvt, Real *rcond,
     integer *rank, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgelsy)(integer *m, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *jpvt, Real *rcond,
     integer *rank, Complex *work, integer *lwork, Real *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeql2)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgeqlf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgeqp3)(integer *m, integer *n, Complex *a, integer *lda,
     integer *jpvt, Complex *tau, Complex *work, integer *lwork, Real *
    rwork, integer *info);

/* Subroutine */ int F77NAME(cgeqpf)(integer *m, integer *n, Complex *a, integer *lda,
     integer *jpvt, Complex *tau, Complex *work, Real *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeqr2)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgeqrf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgerfs)(char *trans, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, integer *ipiv, Complex *
    b, integer *ldb, Complex *x, integer *ldx, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgerq2)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgerqf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgesc2)(integer *n, Complex *a, integer *lda, Complex *
    rhs, integer *ipiv, integer *jpiv, Real *scale);

/* Subroutine */ int F77NAME(cgesv)(integer *n, integer *nrhs, Complex *a, integer *
    lda, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, integer *
    ipiv, char *equed, Real *r__, Real *c__, Complex *b, integer *ldb,
    Complex *x, integer *ldx, Real *rcond, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgetc2)(integer *n, Complex *a, integer *lda, integer *
    ipiv, integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(cgetf2)(integer *m, integer *n, Complex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgetrf)(integer *m, integer *n, Complex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgetri)(integer *n, Complex *a, integer *lda, integer *
    ipiv, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgetrs)(char *trans, integer *n, integer *nrhs, Complex *
    a, integer *lda, integer *ipiv, Complex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, Real *lscale, Real *rscale, integer *m, Complex *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(cggbal)(char *job, integer *n, Complex *a, integer *lda,
    Complex *b, integer *ldb, integer *ilo, integer *ihi, Real *lscale,
    Real *rscale, Real *work, integer *info);

/* Subroutine */ int F77NAME(cgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, integer *n, Complex *a, integer *lda, Complex *b, integer *
    ldb, integer *sdim, Complex *alpha, Complex *beta, Complex *vsl,
    integer *ldvsl, Complex *vsr, integer *ldvsr, Complex *work, integer *
    lwork, Real *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, char *sense, integer *n, Complex *a, integer *lda, Complex *b,
     integer *ldb, integer *sdim, Complex *alpha, Complex *beta, Complex *
    vsl, integer *ldvsl, Complex *vsr, integer *ldvsr, Real *rconde, Real
    *rcondv, Complex *work, integer *lwork, Real *rwork, integer *iwork,
    integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cggev)(char *jobvl, char *jobvr, integer *n, Complex *a,
    integer *lda, Complex *b, integer *ldb, Complex *alpha, Complex *beta,
     Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, Complex *
    work, integer *lwork, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
     Complex *alpha, Complex *beta, Complex *vl, integer *ldvl, Complex *
    vr, integer *ldvr, integer *ilo, integer *ihi, Real *lscale, Real *
    rscale, Real *abnrm, Real *bbnrm, Real *rconde, Real *rcondv, Complex
    *work, integer *lwork, Real *rwork, integer *iwork, logical *bwork,
    integer *info);

/* Subroutine */ int F77NAME(cggglm)(integer *n, integer *m, integer *p, Complex *a,
    integer *lda, Complex *b, integer *ldb, Complex *d__, Complex *x,
    Complex *y, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, Complex *a, integer *lda, Complex *b, integer *ldb,
     Complex *q, integer *ldq, Complex *z__, integer *ldz, integer *info);

/* Subroutine */ int F77NAME(cgglse)(integer *m, integer *n, integer *p, Complex *a,
    integer *lda, Complex *b, integer *ldb, Complex *c__, Complex *d__,
    Complex *x, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cggqrf)(integer *n, integer *m, integer *p, Complex *a,
    integer *lda, Complex *taua, Complex *b, integer *ldb, Complex *taub,
    Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cggrqf)(integer *m, integer *p, integer *n, Complex *a,
    integer *lda, Complex *taua, Complex *b, integer *ldb, Complex *taub,
    Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, Complex *a, integer *
    lda, Complex *b, integer *ldb, Real *alpha, Real *beta, Complex *u,
    integer *ldu, Complex *v, integer *ldv, Complex *q, integer *ldq,
    Complex *work, Real *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(cggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, Complex *a, integer *lda, Complex *b, integer
    *ldb, Real *tola, Real *tolb, integer *k, integer *l, Complex *u,
    integer *ldu, Complex *v, integer *ldv, Complex *q, integer *ldq,
    integer *iwork, Real *rwork, Complex *tau, Complex *work, integer *
    info);

/* Subroutine */ int F77NAME(cgtcon)(char *norm, integer *n, Complex *dl, Complex *
    d__, Complex *du, Complex *du2, integer *ipiv, Real *anorm, Real *
    rcond, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgtrfs)(char *trans, integer *n, integer *nrhs, Complex *
    dl, Complex *d__, Complex *du, Complex *dlf, Complex *df, Complex *
    duf, Complex *du2, integer *ipiv, Complex *b, integer *ldb, Complex *
    x, integer *ldx, Real *ferr, Real *berr, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cgtsv)(integer *n, integer *nrhs, Complex *dl, Complex *
    d__, Complex *du, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, Complex *dl, Complex *d__, Complex *du, Complex *dlf, Complex *
    df, Complex *duf, Complex *du2, integer *ipiv, Complex *b, integer *
    ldb, Complex *x, integer *ldx, Real *rcond, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cgttrf)(integer *n, Complex *dl, Complex *d__, Complex *
    du, Complex *du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgttrs)(char *trans, integer *n, integer *nrhs, Complex *
    dl, Complex *d__, Complex *du, Complex *du2, integer *ipiv, Complex *
    b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgtts2)(integer *itrans, integer *n, integer *nrhs,
    Complex *dl, Complex *d__, Complex *du, Complex *du2, integer *ipiv,
    Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(chbev)(char *jobz, char *uplo, integer *n, integer *kd,
    Complex *ab, integer *ldab, Real *w, Complex *z__, integer *ldz,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(chbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    Complex *ab, integer *ldab, Real *w, Complex *z__, integer *ldz,
    Complex *work, integer *lwork, Real *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, Complex *ab, integer *ldab, Complex *q, integer *ldq,
    Real *vl, Real *vu, integer *il, integer *iu, Real *abstol, integer *
    m, Real *w, Complex *z__, integer *ldz, Complex *work, Real *rwork,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, Complex *ab, integer *ldab, Complex *bb, integer *ldbb,
    Complex *x, integer *ldx, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(chbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, Complex *ab, integer *ldab, Complex *bb, integer *ldbb,
    Real *w, Complex *z__, integer *ldz, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, Complex *ab, integer *ldab, Complex *bb,
    integer *ldbb, Complex *q, integer *ldq, Real *vl, Real *vu, integer *
    il, integer *iu, Real *abstol, integer *m, Real *w, Complex *z__,
    integer *ldz, Complex *work, Real *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(chbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    Complex *ab, integer *ldab, Real *d__, Real *e, Complex *q, integer *
    ldq, Complex *work, integer *info);

/* Subroutine */ int F77NAME(checon)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, Real *anorm, Real *rcond, Complex *work, integer *
    info);

/* Subroutine */ int F77NAME(cheev)(char *jobz, char *uplo, integer *n, Complex *a,
    integer *lda, Real *w, Complex *work, integer *lwork, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cheevd)(char *jobz, char *uplo, integer *n, Complex *a,
    integer *lda, Real *w, Complex *work, integer *lwork, Real *rwork,
    integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(cheevr)(char *jobz, char *range, char *uplo, integer *n,
    Complex *a, integer *lda, Real *vl, Real *vu, integer *il, integer *
    iu, Real *abstol, integer *m, Real *w, Complex *z__, integer *ldz,
    integer *isuppz, Complex *work, integer *lwork, Real *rwork, integer *
    lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(cheevx)(char *jobz, char *range, char *uplo, integer *n,
    Complex *a, integer *lda, Real *vl, Real *vu, integer *il, integer *
    iu, Real *abstol, integer *m, Real *w, Complex *z__, integer *ldz,
    Complex *work, integer *lwork, Real *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(chegs2)(integer *itype, char *uplo, integer *n, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chegst)(integer *itype, char *uplo, integer *n, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chegv)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, Real *w,
    Complex *work, integer *lwork, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(chegvd)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, Real *w,
    Complex *work, integer *lwork, Real *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chegvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    Real *vl, Real *vu, integer *il, integer *iu, Real *abstol, integer *
    m, Real *w, Complex *z__, integer *ldz, Complex *work, integer *lwork,
     Real *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(cherfs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, integer *ipiv, Complex *
    b, integer *ldb, Complex *x, integer *ldx, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(chesv)(char *uplo, integer *n, integer *nrhs, Complex *a,
     integer *lda, integer *ipiv, Complex *b, integer *ldb, Complex *work,
     integer *lwork, integer *info);

/* Subroutine */ int F77NAME(chesvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, integer *
    ipiv, Complex *b, integer *ldb, Complex *x, integer *ldx, Real *rcond,
     Real *ferr, Real *berr, Complex *work, integer *lwork, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chetf2)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(chetrd)(char *uplo, integer *n, Complex *a, integer *lda,
     Real *d__, Real *e, Complex *tau, Complex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(chetrf)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(chetri)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, Complex *work, integer *info);

/* Subroutine */ int F77NAME(chetrs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, integer *ipiv, Complex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(chgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, Complex *a, integer *lda, Complex *b,
    integer *ldb, Complex *alpha, Complex *beta, Complex *q, integer *ldq,
     Complex *z__, integer *ldz, Complex *work, integer *lwork, Real *
    rwork, integer *info);

/* Subroutine */ int F77NAME(chpcon)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, Real *anorm, Real *rcond, Complex *work, integer *info);

/* Subroutine */ int F77NAME(chpev)(char *jobz, char *uplo, integer *n, Complex *ap,
    Real *w, Complex *z__, integer *ldz, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chpevd)(char *jobz, char *uplo, integer *n, Complex *ap,
    Real *w, Complex *z__, integer *ldz, Complex *work, integer *lwork,
    Real *rwork, integer *lrwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(chpevx)(char *jobz, char *range, char *uplo, integer *n,
    Complex *ap, Real *vl, Real *vu, integer *il, integer *iu, Real *
    abstol, integer *m, Real *w, Complex *z__, integer *ldz, Complex *
    work, Real *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chpgst)(integer *itype, char *uplo, integer *n, Complex *
    ap, Complex *bp, integer *info);

/* Subroutine */ int F77NAME(chpgv)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *ap, Complex *bp, Real *w, Complex *z__, integer *ldz,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(chpgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *ap, Complex *bp, Real *w, Complex *z__, integer *ldz,
    Complex *work, integer *lwork, Real *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chpgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, Complex *ap, Complex *bp, Real *vl, Real *vu,
    integer *il, integer *iu, Real *abstol, integer *m, Real *w, Complex *
    z__, integer *ldz, Complex *work, Real *rwork, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chprfs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *afp, integer *ipiv, Complex *b, integer *ldb, Complex *x,
     integer *ldx, Real *ferr, Real *berr, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chpsv)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chpsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *ap, Complex *afp, integer *ipiv, Complex *b, integer *
    ldb, Complex *x, integer *ldx, Real *rcond, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(chptrd)(char *uplo, integer *n, Complex *ap, Real *d__,
    Real *e, Complex *tau, integer *info);

/* Subroutine */ int F77NAME(chptrf)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(chptri)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, Complex *work, integer *info);

/* Subroutine */ int F77NAME(chptrs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, Complex *h__, integer *ldh, Complex *w, Complex *
    vl, integer *ldvl, Complex *vr, integer *ldvr, integer *mm, integer *
    m, Complex *work, Real *rwork, integer *ifaill, integer *ifailr,
    integer *info);

/* Subroutine */ int F77NAME(chseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, Complex *h__, integer *ldh, Complex *w, Complex *z__,
    integer *ldz, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(clabrd)(integer *m, integer *n, integer *nb, Complex *a,
    integer *lda, Real *d__, Real *e, Complex *tauq, Complex *taup,
    Complex *x, integer *ldx, Complex *y, integer *ldy);

/* Subroutine */ int F77NAME(clacgv)(integer *n, Complex *x, integer *incx);

/* Subroutine */ int F77NAME(clacon)(integer *n, Complex *v, Complex *x, Real *est,
    integer *kase);

/* Subroutine */ int F77NAME(clacp2)(char *uplo, integer *m, integer *n, Real *a,
    integer *lda, Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(clacpy)(char *uplo, integer *m, integer *n, Complex *a,
    integer *lda, Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(clacrm)(integer *m, integer *n, Complex *a, integer *lda,
     Real *b, integer *ldb, Complex *c__, integer *ldc, Real *rwork);

/* Subroutine */ int F77NAME(clacrt)(integer *n, Complex *cx, integer *incx, Complex *
    cy, integer *incy, Complex *c__, Complex *s);

/* Subroutine */ int F77NAME(claed0)(integer *qsiz, integer *n, Real *d__, Real *e,
    Complex *q, integer *ldq, Complex *qstore, integer *ldqs, Real *rwork,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(claed7)(integer *n, integer *cutpnt, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, Real *d__, Complex *
    q, integer *ldq, Real *rho, integer *indxq, Real *qstore, integer *
    qptr, integer *prmptr, integer *perm, integer *givptr, integer *
    givcol, Real *givnum, Complex *work, Real *rwork, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(claed8)(integer *k, integer *n, integer *qsiz, Complex *
    q, integer *ldq, Real *d__, Real *rho, integer *cutpnt, Real *z__,
    Real *dlamda, Complex *q2, integer *ldq2, Real *w, integer *indxp,
    integer *indx, integer *indxq, integer *perm, integer *givptr,
    integer *givcol, Real *givnum, integer *info);

/* Subroutine */ int F77NAME(claein)(logical *rightv, logical *noinit, integer *n,
    Complex *h__, integer *ldh, Complex *w, Complex *v, Complex *b,
    integer *ldb, Real *rwork, Real *eps3, Real *smlnum, integer *info);

/* Subroutine */ int F77NAME(claesy)(Complex *a, Complex *b, Complex *c__, Complex *
    rt1, Complex *rt2, Complex *evscal, Complex *cs1, Complex *sn1);

/* Subroutine */ int F77NAME(claev2)(Complex *a, Complex *b, Complex *c__, Real *rt1,
    Real *rt2, Real *cs1, Complex *sn1);

/* Subroutine */ int F77NAME(clags2)(logical *upper, Real *a1, Complex *a2, Real *a3,
    Real *b1, Complex *b2, Real *b3, Real *csu, Complex *snu, Real *csv,
    Complex *snv, Real *csq, Complex *snq);

/* Subroutine */ int F77NAME(clagtm)(char *trans, integer *n, integer *nrhs, Real *
    alpha, Complex *dl, Complex *d__, Complex *du, Complex *x, integer *
    ldx, Real *beta, Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(clahef)(char *uplo, integer *n, integer *nb, integer *kb,
     Complex *a, integer *lda, integer *ipiv, Complex *w, integer *ldw,
    integer *info);

/* Subroutine */ int F77NAME(clahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, Complex *h__, integer *ldh, Complex *w,
    integer *iloz, integer *ihiz, Complex *z__, integer *ldz, integer *
    info);

/* Subroutine */ int F77NAME(clahrd)(integer *n, integer *k, integer *nb, Complex *a,
    integer *lda, Complex *tau, Complex *t, integer *ldt, Complex *y,
    integer *ldy);

/* Subroutine */ int F77NAME(claic1)(integer *job, integer *j, Complex *x, Real *sest,
     Complex *w, Complex *gamma, Real *sestpr, Complex *s, Complex *c__);

/* Subroutine */ int F77NAME(clals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, Complex *b, integer *ldb, Complex *bx,
    integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, Real *givnum, integer *ldgnum, Real *poles, Real *
    difl, Real *difr, Real *z__, integer *k, Real *c__, Real *s, Real *
    rwork, integer *info);

/* Subroutine */ int F77NAME(clalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, Complex *b, integer *ldb, Complex *bx, integer *ldbx,
    Real *u, integer *ldu, Real *vt, integer *k, Real *difl, Real *difr,
    Real *z__, Real *poles, integer *givptr, integer *givcol, integer *
    ldgcol, integer *perm, Real *givnum, Real *c__, Real *s, Real *rwork,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(clapll)(integer *n, Complex *x, integer *incx, Complex *
    y, integer *incy, Real *ssmin);

/* Subroutine */ int F77NAME(clapmt)(logical *forwrd, integer *m, integer *n, Complex
    *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(claqgb)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, Real *r__, Real *c__, Real *rowcnd, Real
    *colcnd, Real *amax, char *equed);

/* Subroutine */ int F77NAME(claqge)(integer *m, integer *n, Complex *a, integer *lda,
     Real *r__, Real *c__, Real *rowcnd, Real *colcnd, Real *amax, char *
    equed);

/* Subroutine */ int F77NAME(claqhb)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, Real *s, Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(claqhe)(char *uplo, integer *n, Complex *a, integer *lda,
     Real *s, Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(claqhp)(char *uplo, integer *n, Complex *ap, Real *s,
    Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(claqp2)(integer *m, integer *n, integer *offset, Complex
    *a, integer *lda, integer *jpvt, Complex *tau, Real *vn1, Real *vn2,
    Complex *work);

/* Subroutine */ int F77NAME(claqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, Complex *a, integer *lda, integer *jpvt, Complex *
    tau, Real *vn1, Real *vn2, Complex *auxv, Complex *f, integer *ldf);

/* Subroutine */ int F77NAME(claqsb)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, Real *s, Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(claqsp)(char *uplo, integer *n, Complex *ap, Real *s,
    Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(claqsy)(char *uplo, integer *n, Complex *a, integer *lda,
     Real *s, Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(clar1v)(integer *n, integer *b1, integer *bn, Real *
    sigma, Real *d__, Real *l, Real *ld, Real *lld, Real *gersch, Complex
    *z__, Real *ztz, Real *mingma, integer *r__, integer *isuppz, Real *
    work);

/* Subroutine */ int F77NAME(clar2v)(integer *n, Complex *x, Complex *y, Complex *z__,
     integer *incx, Real *c__, Complex *s, integer *incc);

/* Subroutine */ int F77NAME(clarcm)(integer *m, integer *n, Real *a, integer *lda,
    Complex *b, integer *ldb, Complex *c__, integer *ldc, Real *rwork);

/* Subroutine */ int F77NAME(clarf)(char *side, integer *m, integer *n, Complex *v,
    integer *incv, Complex *tau, Complex *c__, integer *ldc, Complex *
    work);

/* Subroutine */ int F77NAME(clarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, Complex *v, integer *ldv,
    Complex *t, integer *ldt, Complex *c__, integer *ldc, Complex *work,
    integer *ldwork);

/* Subroutine */ int F77NAME(clarfg)(integer *n, Complex *alpha, Complex *x, integer *
    incx, Complex *tau);

/* Subroutine */ int F77NAME(clarft)(char *direct, char *storev, integer *n, integer *
    k, Complex *v, integer *ldv, Complex *tau, Complex *t, integer *ldt);

/* Subroutine */ int F77NAME(clarfx)(char *side, integer *m, integer *n, Complex *v,
    Complex *tau, Complex *c__, integer *ldc, Complex *work);

/* Subroutine */ int F77NAME(clargv)(integer *n, Complex *x, integer *incx, Complex *
    y, integer *incy, Real *c__, integer *incc);

/* Subroutine */ int F77NAME(clarnv)(integer *idist, integer *iseed, integer *n,
    Complex *x);

/* Subroutine */ int F77NAME(clarrv)(integer *n, Real *d__, Real *l, integer *isplit,
    integer *m, Real *w, integer *iblock, Real *gersch, Real *tol,
    Complex *z__, integer *ldz, integer *isuppz, Real *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(clartg)(Complex *f, Complex *g, Real *cs, Complex *sn,
    Complex *r__);

/* Subroutine */ int F77NAME(clartv)(integer *n, Complex *x, integer *incx, Complex *
    y, integer *incy, Real *c__, Complex *s, integer *incc);

/* Subroutine */ int F77NAME(clarz)(char *side, integer *m, integer *n, integer *l,
    Complex *v, integer *incv, Complex *tau, Complex *c__, integer *ldc,
    Complex *work);

/* Subroutine */ int F77NAME(clarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, Complex *v,
    integer *ldv, Complex *t, integer *ldt, Complex *c__, integer *ldc,
    Complex *work, integer *ldwork);

/* Subroutine */ int F77NAME(clarzt)(char *direct, char *storev, integer *n, integer *
    k, Complex *v, integer *ldv, Complex *tau, Complex *t, integer *ldt);

/* Subroutine */ int F77NAME(clascl)(char *type__, integer *kl, integer *ku, Real *
    cfrom, Real *cto, integer *m, integer *n, Complex *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(claset)(char *uplo, integer *m, integer *n, Complex *
    alpha, Complex *beta, Complex *a, integer *lda);

/* Subroutine */ int F77NAME(clasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, Real *c__, Real *s, Complex *a, integer *lda);

/* Subroutine */ int F77NAME(classq)(integer *n, Complex *x, integer *incx, Real *
    scale, Real *sumsq);

/* Subroutine */ int F77NAME(claswp)(integer *n, Complex *a, integer *lda, integer *
    k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(clasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     Complex *a, integer *lda, integer *ipiv, Complex *w, integer *ldw,
    integer *info);

/* Subroutine */ int F77NAME(clatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, Complex *ab, integer *ldab, Complex *
    x, Real *scale, Real *cnorm, integer *info);

/* Subroutine */ int F77NAME(clatdf)(integer *ijob, integer *n, Complex *z__, integer
    *ldz, Complex *rhs, Real *rdsum, Real *rdscal, integer *ipiv, integer
    *jpiv);

/* Subroutine */ int F77NAME(clatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, Complex *ap, Complex *x, Real *scale, Real *cnorm,
     integer *info);

/* Subroutine */ int F77NAME(clatrd)(char *uplo, integer *n, integer *nb, Complex *a,
    integer *lda, Real *e, Complex *tau, Complex *w, integer *ldw);

/* Subroutine */ int F77NAME(clatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, Complex *a, integer *lda, Complex *x, Real *scale,
     Real *cnorm, integer *info);

/* Subroutine */ int F77NAME(clatrz)(integer *m, integer *n, integer *l, Complex *a,
    integer *lda, Complex *tau, Complex *work);

/* Subroutine */ int F77NAME(clatzm)(char *side, integer *m, integer *n, Complex *v,
    integer *incv, Complex *tau, Complex *c1, Complex *c2, integer *ldc,
    Complex *work);

/* Subroutine */ int F77NAME(clauu2)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(clauum)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpbcon)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, Real *anorm, Real *rcond, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cpbequ)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, Real *s, Real *scond, Real *amax, integer *info);

/* Subroutine */ int F77NAME(cpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Complex *ab, integer *ldab, Complex *afb, integer *ldafb,
    Complex *b, integer *ldb, Complex *x, integer *ldx, Real *ferr, Real *
    berr, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cpbstf)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Complex *ab, integer *ldab, Complex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, Complex *ab, integer *ldab, Complex *afb, integer *
    ldafb, char *equed, Real *s, Complex *b, integer *ldb, Complex *x,
    integer *ldx, Real *rcond, Real *ferr, Real *berr, Complex *work,
    Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cpbtf2)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbtrf)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Complex *ab, integer *ldab, Complex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cpocon)(char *uplo, integer *n, Complex *a, integer *lda,
     Real *anorm, Real *rcond, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cpoequ)(integer *n, Complex *a, integer *lda, Real *s,
    Real *scond, Real *amax, integer *info);

/* Subroutine */ int F77NAME(cporfs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, Complex *b, integer *ldb,
     Complex *x, integer *ldx, Real *ferr, Real *berr, Complex *work,
    Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cposv)(char *uplo, integer *n, integer *nrhs, Complex *a,
     integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, char *
    equed, Real *s, Complex *b, integer *ldb, Complex *x, integer *ldx,
    Real *rcond, Real *ferr, Real *berr, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cpotf2)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotrf)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotri)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotrs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cppcon)(char *uplo, integer *n, Complex *ap, Real *anorm,
     Real *rcond, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cppequ)(char *uplo, integer *n, Complex *ap, Real *s,
    Real *scond, Real *amax, integer *info);

/* Subroutine */ int F77NAME(cpprfs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *afp, Complex *b, integer *ldb, Complex *x, integer *ldx,
    Real *ferr, Real *berr, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cppsv)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *ap, Complex *afp, char *equed, Real *s, Complex *b,
    integer *ldb, Complex *x, integer *ldx, Real *rcond, Real *ferr, Real
    *berr, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cpptrf)(char *uplo, integer *n, Complex *ap, integer *
    info);

/* Subroutine */ int F77NAME(cpptri)(char *uplo, integer *n, Complex *ap, integer *
    info);

/* Subroutine */ int F77NAME(cpptrs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptcon)(integer *n, Real *d__, Complex *e, Real *anorm,
    Real *rcond, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cptrfs)(char *uplo, integer *n, integer *nrhs, Real *d__,
     Complex *e, Real *df, Complex *ef, Complex *b, integer *ldb, Complex
    *x, integer *ldx, Real *ferr, Real *berr, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cptsv)(integer *n, integer *nrhs, Real *d__, Complex *e,
    Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptsvx)(char *fact, integer *n, integer *nrhs, Real *d__,
     Complex *e, Real *df, Complex *ef, Complex *b, integer *ldb, Complex
    *x, integer *ldx, Real *rcond, Real *ferr, Real *berr, Complex *work,
    Real *rwork, integer *info);

/* Subroutine */ int F77NAME(cpttrf)(integer *n, Real *d__, Complex *e, integer *info);

/* Subroutine */ int F77NAME(cpttrs)(char *uplo, integer *n, integer *nrhs, Real *d__,
     Complex *e, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptts2)(integer *iuplo, integer *n, integer *nrhs, Real *
    d__, Complex *e, Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(crot)(integer *n, Complex *cx, integer *incx, Complex *
    cy, integer *incy, Real *c__, Complex *s);

/* Subroutine */ int F77NAME(cspcon)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, Real *anorm, Real *rcond, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cspmv)(char *uplo, integer *n, Complex *alpha, Complex *
    ap, Complex *x, integer *incx, Complex *beta, Complex *y, integer *
    incy);

/* Subroutine */ int F77NAME(cspr)(char *uplo, integer *n, Complex *alpha, Complex *x,
     integer *incx, Complex *ap);

/* Subroutine */ int F77NAME(csprfs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *afp, integer *ipiv, Complex *b, integer *ldb, Complex *x,
     integer *ldx, Real *ferr, Real *berr, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cspsv)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *ap, Complex *afp, integer *ipiv, Complex *b, integer *
    ldb, Complex *x, integer *ldx, Real *rcond, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(csptrf)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(csptri)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, Complex *work, integer *info);

/* Subroutine */ int F77NAME(csptrs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(csrot)(integer *n, Complex *cx, integer *incx, Complex *
    cy, integer *incy, Real *c__, Real *s);

/* Subroutine */ int F77NAME(csrscl)(integer *n, Real *sa, Complex *sx, integer *incx);

/* Subroutine */ int F77NAME(cstedc)(char *compz, integer *n, Real *d__, Real *e,
    Complex *z__, integer *ldz, Complex *work, integer *lwork, Real *
    rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(cstein)(integer *n, Real *d__, Real *e, integer *m, Real
    *w, integer *iblock, integer *isplit, Complex *z__, integer *ldz,
    Real *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(csteqr)(char *compz, integer *n, Real *d__, Real *e,
    Complex *z__, integer *ldz, Real *work, integer *info);

/* Subroutine */ int F77NAME(csycon)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, Real *anorm, Real *rcond, Complex *work, integer *
    info);

/* Subroutine */ int F77NAME(csymv)(char *uplo, integer *n, Complex *alpha, Complex *
    a, integer *lda, Complex *x, integer *incx, Complex *beta, Complex *y,
     integer *incy);

/* Subroutine */ int F77NAME(csyr)(char *uplo, integer *n, Complex *alpha, Complex *x,
     integer *incx, Complex *a, integer *lda);

/* Subroutine */ int F77NAME(csyrfs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, integer *ipiv, Complex *
    b, integer *ldb, Complex *x, integer *ldx, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(csysv)(char *uplo, integer *n, integer *nrhs, Complex *a,
     integer *lda, integer *ipiv, Complex *b, integer *ldb, Complex *work,
     integer *lwork, integer *info);

/* Subroutine */ int F77NAME(csysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, integer *
    ipiv, Complex *b, integer *ldb, Complex *x, integer *ldx, Real *rcond,
     Real *ferr, Real *berr, Complex *work, integer *lwork, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(csytf2)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(csytrf)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(csytri)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, Complex *work, integer *info);

/* Subroutine */ int F77NAME(csytrs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, integer *ipiv, Complex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(ctbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, Complex *ab, integer *ldab, Real *rcond, Complex *work,
    Real *rwork, integer *info);

/* Subroutine */ int F77NAME(ctbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, Complex *ab, integer *ldab, Complex *b,
    integer *ldb, Complex *x, integer *ldx, Real *ferr, Real *berr,
    Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(ctbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, Complex *ab, integer *ldab, Complex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ctgevc)(char *side, char *howmny, logical *select,
    integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, integer *mm,
    integer *m, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(ctgex2)(logical *wantq, logical *wantz, integer *n,
    Complex *a, integer *lda, Complex *b, integer *ldb, Complex *q,
    integer *ldq, Complex *z__, integer *ldz, integer *j1, integer *info);

/* Subroutine */ int F77NAME(ctgexc)(logical *wantq, logical *wantz, integer *n,
    Complex *a, integer *lda, Complex *b, integer *ldb, Complex *q,
    integer *ldq, Complex *z__, integer *ldz, integer *ifst, integer *
    ilst, integer *info);

/* Subroutine */ int F77NAME(ctgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, Complex *a, integer *lda, Complex *b,
    integer *ldb, Complex *alpha, Complex *beta, Complex *q, integer *ldq,
     Complex *z__, integer *ldz, integer *m, Real *pl, Real *pr, Real *
    dif, Complex *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ctgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, Complex *a, integer *
    lda, Complex *b, integer *ldb, Real *tola, Real *tolb, Real *alpha,
    Real *beta, Complex *u, integer *ldu, Complex *v, integer *ldv,
    Complex *q, integer *ldq, Complex *work, integer *ncycle, integer *
    info);

/* Subroutine */ int F77NAME(ctgsna)(char *job, char *howmny, logical *select,
    integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, Real *s, Real
    *dif, integer *mm, integer *m, Complex *work, integer *lwork, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(ctgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, Complex *c__,
    integer *ldc, Complex *d__, integer *ldd, Complex *e, integer *lde,
    Complex *f, integer *ldf, Real *scale, Real *rdsum, Real *rdscal,
    integer *info);

/* Subroutine */ int F77NAME(ctgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, Complex *c__,
    integer *ldc, Complex *d__, integer *ldd, Complex *e, integer *lde,
    Complex *f, integer *ldf, Real *scale, Real *dif, Complex *work,
    integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ctpcon)(char *norm, char *uplo, char *diag, integer *n,
    Complex *ap, Real *rcond, Complex *work, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(ctprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Complex *ap, Complex *b, integer *ldb, Complex *x,
    integer *ldx, Real *ferr, Real *berr, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(ctptri)(char *uplo, char *diag, integer *n, Complex *ap,
    integer *info);

/* Subroutine */ int F77NAME(ctptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Complex *ap, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ctrcon)(char *norm, char *uplo, char *diag, integer *n,
    Complex *a, integer *lda, Real *rcond, Complex *work, Real *rwork,
    integer *info);

/* Subroutine */ int F77NAME(ctrevc)(char *side, char *howmny, logical *select,
    integer *n, Complex *t, integer *ldt, Complex *vl, integer *ldvl,
    Complex *vr, integer *ldvr, integer *mm, integer *m, Complex *work,
    Real *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrexc)(char *compq, integer *n, Complex *t, integer *
    ldt, Complex *q, integer *ldq, integer *ifst, integer *ilst, integer *
    info);

/* Subroutine */ int F77NAME(ctrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *x, integer *ldx, Real *ferr, Real *berr, Complex *work, Real
    *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrsen)(char *job, char *compq, logical *select, integer
    *n, Complex *t, integer *ldt, Complex *q, integer *ldq, Complex *w,
    integer *m, Real *s, Real *sep, Complex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(ctrsna)(char *job, char *howmny, logical *select,
    integer *n, Complex *t, integer *ldt, Complex *vl, integer *ldvl,
    Complex *vr, integer *ldvr, Real *s, Real *sep, integer *mm, integer *
    m, Complex *work, integer *ldwork, Real *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *c__, integer *ldc, Real *scale, integer *info);

/* Subroutine */ int F77NAME(ctrti2)(char *uplo, char *diag, integer *n, Complex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(ctrtri)(char *uplo, char *diag, integer *n, Complex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(ctrtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Complex *a, integer *lda, Complex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(ctzrqf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, integer *info);

/* Subroutine */ int F77NAME(ctzrzf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cung2l)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cung2r)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cungbr)(char *vect, integer *m, integer *n, integer *k,
    Complex *a, integer *lda, Complex *tau, Complex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(cunghr)(integer *n, integer *ilo, integer *ihi, Complex *
    a, integer *lda, Complex *tau, Complex *work, integer *lwork, integer
    *info);

/* Subroutine */ int F77NAME(cungl2)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cunglq)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungql)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungqr)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungr2)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cungrq)(integer *m, integer *n, integer *k, Complex *a,
    integer *lda, Complex *tau, Complex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungtr)(char *uplo, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cunm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cunmbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, Complex *a, integer *lda, Complex *tau,
    Complex *c__, integer *ldc, Complex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cunmhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, Complex *a, integer *lda, Complex *tau,
    Complex *c__, integer *ldc, Complex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cunml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cunmlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmql)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cunmr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, Complex *a, integer *lda, Complex *tau,
    Complex *c__, integer *ldc, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cunmrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, Complex *a, integer *lda, Complex *tau,
    Complex *c__, integer *ldc, Complex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cunmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, Complex *a, integer *lda, Complex *tau, Complex *c__,
    integer *ldc, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cupgtr)(char *uplo, integer *n, Complex *ap, Complex *
    tau, Complex *q, integer *ldq, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cupmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, Complex *ap, Complex *tau, Complex *c__, integer *ldc,
    Complex *work, integer *info);

/* Subroutine */ int F77NAME(dbdsdc)(char *uplo, char *compq, integer *n, doubleReal *
    d__, doubleReal *e, doubleReal *u, integer *ldu, doubleReal *vt,
    integer *ldvt, doubleReal *q, integer *iq, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, doubleReal *d__, doubleReal *e, doubleReal *vt,
    integer *ldvt, doubleReal *u, integer *ldu, doubleReal *c__, integer *
    ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(ddisna)(char *job, integer *m, integer *n, doubleReal *
    d__, doubleReal *sep, integer *info);

/* Subroutine */ int F77NAME(dgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, doubleReal *ab, integer *ldab, doubleReal *
    d__, doubleReal *e, doubleReal *q, integer *ldq, doubleReal *pt,
    integer *ldpt, doubleReal *c__, integer *ldc, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     doubleReal *ab, integer *ldab, integer *ipiv, doubleReal *anorm,
    doubleReal *rcond, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     doubleReal *ab, integer *ldab, doubleReal *r__, doubleReal *c__,
    doubleReal *rowcnd, doubleReal *colcnd, doubleReal *amax, integer *
    info);

/* Subroutine */ int F77NAME(dgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleReal *ab, integer *ldab, doubleReal *afb,
    integer *ldafb, integer *ipiv, doubleReal *b, integer *ldb,
    doubleReal *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, doubleReal *ab, integer *ldab, integer *ipiv, doubleReal *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, doubleReal *ab, integer *ldab,
    doubleReal *afb, integer *ldafb, integer *ipiv, char *equed,
    doubleReal *r__, doubleReal *c__, doubleReal *b, integer *ldb,
    doubleReal *x, integer *ldx, doubleReal *rcond, doubleReal *ferr,
    doubleReal *berr, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     doubleReal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     doubleReal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleReal *ab, integer *ldab, integer *ipiv,
    doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doubleReal *scale, integer *m, doubleReal *v, integer *
    ldv, integer *info);

/* Subroutine */ int F77NAME(dgebal)(char *job, integer *n, doubleReal *a, integer *
    lda, integer *ilo, integer *ihi, doubleReal *scale, integer *info);

/* Subroutine */ int F77NAME(dgebd2)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *d__, doubleReal *e, doubleReal *tauq, doubleReal *
    taup, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dgebrd)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *d__, doubleReal *e, doubleReal *tauq, doubleReal *
    taup, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgecon)(char *norm, integer *n, doubleReal *a, integer *
    lda, doubleReal *anorm, doubleReal *rcond, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgeequ)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *r__, doubleReal *c__, doubleReal *rowcnd, doubleReal
    *colcnd, doubleReal *amax, integer *info);

/* Subroutine */ int F77NAME(dgees)(char *jobvs, char *sort, L_fp select, integer *n,
    doubleReal *a, integer *lda, integer *sdim, doubleReal *wr,
    doubleReal *wi, doubleReal *vs, integer *ldvs, doubleReal *work,
    integer *lwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, doubleReal *a, integer *lda, integer *sdim,
    doubleReal *wr, doubleReal *wi, doubleReal *vs, integer *ldvs,
    doubleReal *rconde, doubleReal *rcondv, doubleReal *work, integer *
    lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dgeev)(char *jobvl, char *jobvr, integer *n, doubleReal *
    a, integer *lda, doubleReal *wr, doubleReal *wi, doubleReal *vl,
    integer *ldvl, doubleReal *vr, integer *ldvr, doubleReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleReal *a, integer *lda, doubleReal *wr,
    doubleReal *wi, doubleReal *vl, integer *ldvl, doubleReal *vr,
    integer *ldvr, integer *ilo, integer *ihi, doubleReal *scale,
    doubleReal *abnrm, doubleReal *rconde, doubleReal *rcondv, doubleReal
    *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgegs)(char *jobvsl, char *jobvsr, integer *n,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, doubleReal *
    alphar, doubleReal *alphai, doubleReal *beta, doubleReal *vsl,
    integer *ldvsl, doubleReal *vsr, integer *ldvsr, doubleReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgegv)(char *jobvl, char *jobvr, integer *n, doubleReal *
    a, integer *lda, doubleReal *b, integer *ldb, doubleReal *alphar,
    doubleReal *alphai, doubleReal *beta, doubleReal *vl, integer *ldvl,
    doubleReal *vr, integer *ldvr, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgehd2)(integer *n, integer *ilo, integer *ihi,
    doubleReal *a, integer *lda, doubleReal *tau, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dgehrd)(integer *n, integer *ilo, integer *ihi,
    doubleReal *a, integer *lda, doubleReal *tau, doubleReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgelq2)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dgelqf)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgels)(char *trans, integer *m, integer *n, integer *
    nrhs, doubleReal *a, integer *lda, doubleReal *b, integer *ldb,
    doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgelsd)(integer *m, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, doubleReal *
    s, doubleReal *rcond, integer *rank, doubleReal *work, integer *lwork,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgelss)(integer *m, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, doubleReal *
    s, doubleReal *rcond, integer *rank, doubleReal *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(dgelsx)(integer *m, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, integer *
    jpvt, doubleReal *rcond, integer *rank, doubleReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dgelsy)(integer *m, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, integer *
    jpvt, doubleReal *rcond, integer *rank, doubleReal *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(dgeql2)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqlf)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgeqp3)(integer *m, integer *n, doubleReal *a, integer *
    lda, integer *jpvt, doubleReal *tau, doubleReal *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(dgeqpf)(integer *m, integer *n, doubleReal *a, integer *
    lda, integer *jpvt, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqr2)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqrf)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgerfs)(char *trans, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *af, integer *ldaf, integer *
    ipiv, doubleReal *b, integer *ldb, doubleReal *x, integer *ldx,
    doubleReal *ferr, doubleReal *berr, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dgerq2)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dgerqf)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgesc2)(integer *n, doubleReal *a, integer *lda,
    doubleReal *rhs, integer *ipiv, integer *jpiv, doubleReal *scale);

/* Subroutine */ int F77NAME(dgesdd)(char *jobz, integer *m, integer *n, doubleReal *
    a, integer *lda, doubleReal *s, doubleReal *u, integer *ldu,
    doubleReal *vt, integer *ldvt, doubleReal *work, integer *lwork,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgesv)(integer *n, integer *nrhs, doubleReal *a, integer
    *lda, integer *ipiv, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgesvd)(char *jobu, char *jobvt, integer *m, integer *n,
    doubleReal *a, integer *lda, doubleReal *s, doubleReal *u, integer *
    ldu, doubleReal *vt, integer *ldvt, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleReal *a, integer *lda, doubleReal *af, integer *ldaf,
    integer *ipiv, char *equed, doubleReal *r__, doubleReal *c__,
    doubleReal *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *
    rcond, doubleReal *ferr, doubleReal *berr, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgetc2)(integer *n, doubleReal *a, integer *lda, integer
    *ipiv, integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(dgetf2)(integer *m, integer *n, doubleReal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgetrf)(integer *m, integer *n, doubleReal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgetri)(integer *n, doubleReal *a, integer *lda, integer
    *ipiv, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgetrs)(char *trans, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, integer *ipiv, doubleReal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doubleReal *lscale, doubleReal *rscale, integer *m,
    doubleReal *v, integer *ldv, integer *info);

/* Subroutine */ int F77NAME(dggbal)(char *job, integer *n, doubleReal *a, integer *
    lda, doubleReal *b, integer *ldb, integer *ilo, integer *ihi,
    doubleReal *lscale, doubleReal *rscale, doubleReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, integer *n, doubleReal *a, integer *lda, doubleReal *b,
    integer *ldb, integer *sdim, doubleReal *alphar, doubleReal *alphai,
    doubleReal *beta, doubleReal *vsl, integer *ldvsl, doubleReal *vsr,
    integer *ldvsr, doubleReal *work, integer *lwork, logical *bwork,
    integer *info);

/* Subroutine */ int F77NAME(dggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, char *sense, integer *n, doubleReal *a, integer *lda,
    doubleReal *b, integer *ldb, integer *sdim, doubleReal *alphar,
    doubleReal *alphai, doubleReal *beta, doubleReal *vsl, integer *ldvsl,
     doubleReal *vsr, integer *ldvsr, doubleReal *rconde, doubleReal *
    rcondv, doubleReal *work, integer *lwork, integer *iwork, integer *
    liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dggev)(char *jobvl, char *jobvr, integer *n, doubleReal *
    a, integer *lda, doubleReal *b, integer *ldb, doubleReal *alphar,
    doubleReal *alphai, doubleReal *beta, doubleReal *vl, integer *ldvl,
    doubleReal *vr, integer *ldvr, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleReal *a, integer *lda, doubleReal *b,
    integer *ldb, doubleReal *alphar, doubleReal *alphai, doubleReal *
    beta, doubleReal *vl, integer *ldvl, doubleReal *vr, integer *ldvr,
    integer *ilo, integer *ihi, doubleReal *lscale, doubleReal *rscale,
    doubleReal *abnrm, doubleReal *bbnrm, doubleReal *rconde, doubleReal *
    rcondv, doubleReal *work, integer *lwork, integer *iwork, logical *
    bwork, integer *info);

/* Subroutine */ int F77NAME(dggglm)(integer *n, integer *m, integer *p, doubleReal *
    a, integer *lda, doubleReal *b, integer *ldb, doubleReal *d__,
    doubleReal *x, doubleReal *y, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, doubleReal *a, integer *lda, doubleReal *b,
    integer *ldb, doubleReal *q, integer *ldq, doubleReal *z__, integer *
    ldz, integer *info);

/* Subroutine */ int F77NAME(dgglse)(integer *m, integer *n, integer *p, doubleReal *
    a, integer *lda, doubleReal *b, integer *ldb, doubleReal *c__,
    doubleReal *d__, doubleReal *x, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dggqrf)(integer *n, integer *m, integer *p, doubleReal *
    a, integer *lda, doubleReal *taua, doubleReal *b, integer *ldb,
    doubleReal *taub, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dggrqf)(integer *m, integer *p, integer *n, doubleReal *
    a, integer *lda, doubleReal *taua, doubleReal *b, integer *ldb,
    doubleReal *taub, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, doubleReal *a,
    integer *lda, doubleReal *b, integer *ldb, doubleReal *alpha,
    doubleReal *beta, doubleReal *u, integer *ldu, doubleReal *v, integer
    *ldv, doubleReal *q, integer *ldq, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, doubleReal *a, integer *lda, doubleReal *b,
    integer *ldb, doubleReal *tola, doubleReal *tolb, integer *k, integer
    *l, doubleReal *u, integer *ldu, doubleReal *v, integer *ldv,
    doubleReal *q, integer *ldq, integer *iwork, doubleReal *tau,
    doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dgtcon)(char *norm, integer *n, doubleReal *dl,
    doubleReal *d__, doubleReal *du, doubleReal *du2, integer *ipiv,
    doubleReal *anorm, doubleReal *rcond, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgtrfs)(char *trans, integer *n, integer *nrhs,
    doubleReal *dl, doubleReal *d__, doubleReal *du, doubleReal *dlf,
    doubleReal *df, doubleReal *duf, doubleReal *du2, integer *ipiv,
    doubleReal *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *
    ferr, doubleReal *berr, doubleReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dgtsv)(integer *n, integer *nrhs, doubleReal *dl,
    doubleReal *d__, doubleReal *du, doubleReal *b, integer *ldb, integer
    *info);

/* Subroutine */ int F77NAME(dgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleReal *dl, doubleReal *d__, doubleReal *du, doubleReal *
    dlf, doubleReal *df, doubleReal *duf, doubleReal *du2, integer *ipiv,
    doubleReal *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *
    rcond, doubleReal *ferr, doubleReal *berr, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgttrf)(integer *n, doubleReal *dl, doubleReal *d__,
    doubleReal *du, doubleReal *du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgttrs)(char *trans, integer *n, integer *nrhs,
    doubleReal *dl, doubleReal *d__, doubleReal *du, doubleReal *du2,
    integer *ipiv, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgtts2)(integer *itrans, integer *n, integer *nrhs,
    doubleReal *dl, doubleReal *d__, doubleReal *du, doubleReal *du2,
    integer *ipiv, doubleReal *b, integer *ldb);

/* Subroutine */ int F77NAME(dhgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, doubleReal *a, integer *lda, doubleReal *
    b, integer *ldb, doubleReal *alphar, doubleReal *alphai, doubleReal *
    beta, doubleReal *q, integer *ldq, doubleReal *z__, integer *ldz,
    doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dhsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, doubleReal *h__, integer *ldh, doubleReal *wr,
    doubleReal *wi, doubleReal *vl, integer *ldvl, doubleReal *vr,
    integer *ldvr, integer *mm, integer *m, doubleReal *work, integer *
    ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(dhseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, doubleReal *h__, integer *ldh, doubleReal *wr,
    doubleReal *wi, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dlabad)(doubleReal *small, doubleReal *large);

/* Subroutine */ int F77NAME(dlabrd)(integer *m, integer *n, integer *nb, doubleReal *
    a, integer *lda, doubleReal *d__, doubleReal *e, doubleReal *tauq,
    doubleReal *taup, doubleReal *x, integer *ldx, doubleReal *y, integer
    *ldy);

/* Subroutine */ int F77NAME(dlacon)(integer *n, doubleReal *v, doubleReal *x,
    integer *isgn, doubleReal *est, integer *kase);

/* Subroutine */ int F77NAME(dlacpy)(char *uplo, integer *m, integer *n, doubleReal *
    a, integer *lda, doubleReal *b, integer *ldb);

/* Subroutine */ int F77NAME(dladiv)(doubleReal *a, doubleReal *b, doubleReal *c__,
    doubleReal *d__, doubleReal *p, doubleReal *q);

/* Subroutine */ int F77NAME(dlae2)(doubleReal *a, doubleReal *b, doubleReal *c__,
    doubleReal *rt1, doubleReal *rt2);

/* Subroutine */ int F77NAME(dlaebz)(integer *ijob, integer *nitmax, integer *n,
    integer *mmax, integer *minp, integer *nbmin, doubleReal *abstol,
    doubleReal *reltol, doubleReal *pivmin, doubleReal *d__, doubleReal *
    e, doubleReal *e2, integer *nval, doubleReal *ab, doubleReal *c__,
    integer *mout, integer *nab, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlaed0)(integer *icompq, integer *qsiz, integer *n,
    doubleReal *d__, doubleReal *e, doubleReal *q, integer *ldq,
    doubleReal *qstore, integer *ldqs, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlaed1)(integer *n, doubleReal *d__, doubleReal *q,
    integer *ldq, integer *indxq, doubleReal *rho, integer *cutpnt,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlaed2)(integer *k, integer *n, integer *n1, doubleReal *
    d__, doubleReal *q, integer *ldq, integer *indxq, doubleReal *rho,
    doubleReal *z__, doubleReal *dlamda, doubleReal *w, doubleReal *q2,
    integer *indx, integer *indxc, integer *indxp, integer *coltyp,
    integer *info);

/* Subroutine */ int F77NAME(dlaed3)(integer *k, integer *n, integer *n1, doubleReal *
    d__, doubleReal *q, integer *ldq, doubleReal *rho, doubleReal *dlamda,
     doubleReal *q2, integer *indx, integer *ctot, doubleReal *w,
    doubleReal *s, integer *info);

/* Subroutine */ int F77NAME(dlaed4)(integer *n, integer *i__, doubleReal *d__,
    doubleReal *z__, doubleReal *delta, doubleReal *rho, doubleReal *dlam,
     integer *info);

/* Subroutine */ int F77NAME(dlaed5)(integer *i__, doubleReal *d__, doubleReal *z__,
    doubleReal *delta, doubleReal *rho, doubleReal *dlam);

/* Subroutine */ int F77NAME(dlaed6)(integer *kniter, logical *orgati, doubleReal *
    rho, doubleReal *d__, doubleReal *z__, doubleReal *finit, doubleReal *
    tau, integer *info);

/* Subroutine */ int F77NAME(dlaed7)(integer *icompq, integer *n, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, doubleReal *d__,
    doubleReal *q, integer *ldq, integer *indxq, doubleReal *rho, integer
    *cutpnt, doubleReal *qstore, integer *qptr, integer *prmptr, integer *
    perm, integer *givptr, integer *givcol, doubleReal *givnum,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlaed8)(integer *icompq, integer *k, integer *n, integer
    *qsiz, doubleReal *d__, doubleReal *q, integer *ldq, integer *indxq,
    doubleReal *rho, integer *cutpnt, doubleReal *z__, doubleReal *dlamda,
     doubleReal *q2, integer *ldq2, doubleReal *w, integer *perm, integer
    *givptr, integer *givcol, doubleReal *givnum, integer *indxp, integer
    *indx, integer *info);

/* Subroutine */ int F77NAME(dlaed9)(integer *k, integer *kstart, integer *kstop,
    integer *n, doubleReal *d__, doubleReal *q, integer *ldq, doubleReal *
    rho, doubleReal *dlamda, doubleReal *w, doubleReal *s, integer *lds,
    integer *info);

/* Subroutine */ int F77NAME(dlaeda)(integer *n, integer *tlvls, integer *curlvl,
    integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
    integer *givcol, doubleReal *givnum, doubleReal *q, integer *qptr,
    doubleReal *z__, doubleReal *ztemp, integer *info);

/* Subroutine */ int F77NAME(dlaein)(logical *rightv, logical *noinit, integer *n,
    doubleReal *h__, integer *ldh, doubleReal *wr, doubleReal *wi,
    doubleReal *vr, doubleReal *vi, doubleReal *b, integer *ldb,
    doubleReal *work, doubleReal *eps3, doubleReal *smlnum, doubleReal *
    bignum, integer *info);

/* Subroutine */ int F77NAME(dlaev2)(doubleReal *a, doubleReal *b, doubleReal *c__,
    doubleReal *rt1, doubleReal *rt2, doubleReal *cs1, doubleReal *sn1);

/* Subroutine */ int F77NAME(dlaexc)(logical *wantq, integer *n, doubleReal *t,
    integer *ldt, doubleReal *q, integer *ldq, integer *j1, integer *n1,
    integer *n2, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dlag2)(doubleReal *a, integer *lda, doubleReal *b,
    integer *ldb, doubleReal *safmin, doubleReal *scale1, doubleReal *
    scale2, doubleReal *wr1, doubleReal *wr2, doubleReal *wi);

/* Subroutine */ int F77NAME(dlags2)(logical *upper, doubleReal *a1, doubleReal *a2,
    doubleReal *a3, doubleReal *b1, doubleReal *b2, doubleReal *b3,
    doubleReal *csu, doubleReal *snu, doubleReal *csv, doubleReal *snv,
    doubleReal *csq, doubleReal *snq);

/* Subroutine */ int F77NAME(dlagtf)(integer *n, doubleReal *a, doubleReal *lambda,
    doubleReal *b, doubleReal *c__, doubleReal *tol, doubleReal *d__,
    integer *in, integer *info);

/* Subroutine */ int F77NAME(dlagtm)(char *trans, integer *n, integer *nrhs,
    doubleReal *alpha, doubleReal *dl, doubleReal *d__, doubleReal *du,
    doubleReal *x, integer *ldx, doubleReal *beta, doubleReal *b, integer
    *ldb);

/* Subroutine */ int F77NAME(dlagts)(integer *job, integer *n, doubleReal *a,
    doubleReal *b, doubleReal *c__, doubleReal *d__, integer *in,
    doubleReal *y, doubleReal *tol, integer *info);

/* Subroutine */ int F77NAME(dlagv2)(doubleReal *a, integer *lda, doubleReal *b,
    integer *ldb, doubleReal *alphar, doubleReal *alphai, doubleReal *
    beta, doubleReal *csl, doubleReal *snl, doubleReal *csr, doubleReal *
    snr);

/* Subroutine */ int F77NAME(dlahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, doubleReal *h__, integer *ldh, doubleReal
    *wr, doubleReal *wi, integer *iloz, integer *ihiz, doubleReal *z__,
    integer *ldz, integer *info);

/* Subroutine */ int F77NAME(dlahrd)(integer *n, integer *k, integer *nb, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *t, integer *ldt,
    doubleReal *y, integer *ldy);

/* Subroutine */ int F77NAME(dlaic1)(integer *job, integer *j, doubleReal *x,
    doubleReal *sest, doubleReal *w, doubleReal *gamma, doubleReal *
    sestpr, doubleReal *s, doubleReal *c__);

/* Subroutine */ int F77NAME(dlaln2)(logical *ltrans, integer *na, integer *nw,
    doubleReal *smin, doubleReal *ca, doubleReal *a, integer *lda,
    doubleReal *d1, doubleReal *d2, doubleReal *b, integer *ldb,
    doubleReal *wr, doubleReal *wi, doubleReal *x, integer *ldx,
    doubleReal *scale, doubleReal *xnorm, integer *info);

/* Subroutine */ int F77NAME(dlals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, doubleReal *b, integer *ldb, doubleReal
    *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, doubleReal *givnum, integer *ldgnum, doubleReal *
    poles, doubleReal *difl, doubleReal *difr, doubleReal *z__, integer *
    k, doubleReal *c__, doubleReal *s, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dlalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, doubleReal *b, integer *ldb, doubleReal *bx, integer *
    ldbx, doubleReal *u, integer *ldu, doubleReal *vt, integer *k,
    doubleReal *difl, doubleReal *difr, doubleReal *z__, doubleReal *
    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
    perm, doubleReal *givnum, doubleReal *c__, doubleReal *s, doubleReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlalsd)(char *uplo, integer *smlsiz, integer *n, integer
    *nrhs, doubleReal *d__, doubleReal *e, doubleReal *b, integer *ldb,
    doubleReal *rcond, integer *rank, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlamc1)(integer *beta, integer *t, logical *rnd, logical
    *ieee1);

/* Subroutine */ int F77NAME(dlamc2)(integer *beta, integer *t, logical *rnd,
    doubleReal *eps, integer *emin, doubleReal *rmin, integer *emax,
    doubleReal *rmax);

/* Subroutine */ int F77NAME(dlamc4)(integer *emin, doubleReal *start, integer *base);

/* Subroutine */ int F77NAME(dlamc5)(integer *beta, integer *p, integer *emin,
    logical *ieee, integer *emax, doubleReal *rmax);

/* Subroutine */ int F77NAME(dlamrg)(integer *n1, integer *n2, doubleReal *a, integer
    *dtrd1, integer *dtrd2, integer *index);

/* Subroutine */ int F77NAME(dlanv2)(doubleReal *a, doubleReal *b, doubleReal *c__,
    doubleReal *d__, doubleReal *rt1r, doubleReal *rt1i, doubleReal *rt2r,
     doubleReal *rt2i, doubleReal *cs, doubleReal *sn);

/* Subroutine */ int F77NAME(dlapll)(integer *n, doubleReal *x, integer *incx,
    doubleReal *y, integer *incy, doubleReal *ssmin);

/* Subroutine */ int F77NAME(dlapmt)(logical *forwrd, integer *m, integer *n,
    doubleReal *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(dlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     doubleReal *ab, integer *ldab, doubleReal *r__, doubleReal *c__,
    doubleReal *rowcnd, doubleReal *colcnd, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqge)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *r__, doubleReal *c__, doubleReal *rowcnd, doubleReal
    *colcnd, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqp2)(integer *m, integer *n, integer *offset,
    doubleReal *a, integer *lda, integer *jpvt, doubleReal *tau,
    doubleReal *vn1, doubleReal *vn2, doubleReal *work);

/* Subroutine */ int F77NAME(dlaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, doubleReal *a, integer *lda, integer *jpvt,
    doubleReal *tau, doubleReal *vn1, doubleReal *vn2, doubleReal *auxv,
    doubleReal *f, integer *ldf);

/* Subroutine */ int F77NAME(dlaqsb)(char *uplo, integer *n, integer *kd, doubleReal *
    ab, integer *ldab, doubleReal *s, doubleReal *scond, doubleReal *amax,
     char *equed);

/* Subroutine */ int F77NAME(dlaqsp)(char *uplo, integer *n, doubleReal *ap,
    doubleReal *s, doubleReal *scond, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqsy)(char *uplo, integer *n, doubleReal *a, integer *
    lda, doubleReal *s, doubleReal *scond, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqtr)(logical *ltran, logical *lReal, integer *n,
    doubleReal *t, integer *ldt, doubleReal *b, doubleReal *w, doubleReal
    *scale, doubleReal *x, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dlar1v)(integer *n, integer *b1, integer *bn, doubleReal
    *sigma, doubleReal *d__, doubleReal *l, doubleReal *ld, doubleReal *
    lld, doubleReal *gersch, doubleReal *z__, doubleReal *ztz, doubleReal
    *mingma, integer *r__, integer *isuppz, doubleReal *work);

/* Subroutine */ int F77NAME(dlar2v)(integer *n, doubleReal *x, doubleReal *y,
    doubleReal *z__, integer *incx, doubleReal *c__, doubleReal *s,
    integer *incc);

/* Subroutine */ int F77NAME(dlarf)(char *side, integer *m, integer *n, doubleReal *v,
     integer *incv, doubleReal *tau, doubleReal *c__, integer *ldc,
    doubleReal *work);

/* Subroutine */ int F77NAME(dlarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, doubleReal *v, integer *
    ldv, doubleReal *t, integer *ldt, doubleReal *c__, integer *ldc,
    doubleReal *work, integer *ldwork);

/* Subroutine */ int F77NAME(dlarfg)(integer *n, doubleReal *alpha, doubleReal *x,
    integer *incx, doubleReal *tau);

/* Subroutine */ int F77NAME(dlarft)(char *direct, char *storev, integer *n, integer *
    k, doubleReal *v, integer *ldv, doubleReal *tau, doubleReal *t,
    integer *ldt);

/* Subroutine */ int F77NAME(dlarfx)(char *side, integer *m, integer *n, doubleReal *
    v, doubleReal *tau, doubleReal *c__, integer *ldc, doubleReal *work);

/* Subroutine */ int F77NAME(dlargv)(integer *n, doubleReal *x, integer *incx,
    doubleReal *y, integer *incy, doubleReal *c__, integer *incc);

/* Subroutine */ int F77NAME(dlarnv)(integer *idist, integer *iseed, integer *n,
    doubleReal *x);

/* Subroutine */ int F77NAME(dlarrb)(integer *n, doubleReal *d__, doubleReal *l,
    doubleReal *ld, doubleReal *lld, integer *ifirst, integer *ilast,
    doubleReal *sigma, doubleReal *reltol, doubleReal *w, doubleReal *
    wgap, doubleReal *werr, doubleReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dlarre)(integer *n, doubleReal *d__, doubleReal *e,
    doubleReal *tol, integer *nsplit, integer *isplit, integer *m,
    doubleReal *w, doubleReal *woff, doubleReal *gersch, doubleReal *work,
     integer *info);

/* Subroutine */ int F77NAME(dlarrf)(integer *n, doubleReal *d__, doubleReal *l,
    doubleReal *ld, doubleReal *lld, integer *ifirst, integer *ilast,
    doubleReal *w, doubleReal *dplus, doubleReal *lplus, doubleReal *work,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlarrv)(integer *n, doubleReal *d__, doubleReal *l,
    integer *isplit, integer *m, doubleReal *w, integer *iblock,
    doubleReal *gersch, doubleReal *tol, doubleReal *z__, integer *ldz,
    integer *isuppz, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlartg)(doubleReal *f, doubleReal *g, doubleReal *cs,
    doubleReal *sn, doubleReal *r__);

/* Subroutine */ int F77NAME(dlartv)(integer *n, doubleReal *x, integer *incx,
    doubleReal *y, integer *incy, doubleReal *c__, doubleReal *s, integer
    *incc);

/* Subroutine */ int F77NAME(dlaruv)(integer *iseed, integer *n, doubleReal *x);

/* Subroutine */ int F77NAME(dlarz)(char *side, integer *m, integer *n, integer *l,
    doubleReal *v, integer *incv, doubleReal *tau, doubleReal *c__,
    integer *ldc, doubleReal *work);

/* Subroutine */ int F77NAME(dlarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, doubleReal *v,
     integer *ldv, doubleReal *t, integer *ldt, doubleReal *c__, integer *
    ldc, doubleReal *work, integer *ldwork);

/* Subroutine */ int F77NAME(dlarzt)(char *direct, char *storev, integer *n, integer *
    k, doubleReal *v, integer *ldv, doubleReal *tau, doubleReal *t,
    integer *ldt);

/* Subroutine */ int F77NAME(dlas2)(doubleReal *f, doubleReal *g, doubleReal *h__,
    doubleReal *ssmin, doubleReal *ssmax);

/* Subroutine */ int F77NAME(dlascl)(char *type__, integer *kl, integer *ku,
    doubleReal *cfrom, doubleReal *cto, integer *m, integer *n,
    doubleReal *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(dlasd0)(integer *n, integer *sqre, doubleReal *d__,
    doubleReal *e, doubleReal *u, integer *ldu, doubleReal *vt, integer *
    ldvt, integer *smlsiz, integer *iwork, doubleReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dlasd1)(integer *nl, integer *nr, integer *sqre,
    doubleReal *d__, doubleReal *alpha, doubleReal *beta, doubleReal *u,
    integer *ldu, doubleReal *vt, integer *ldvt, integer *idxq, integer *
    iwork, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasd2)(integer *nl, integer *nr, integer *sqre, integer
    *k, doubleReal *d__, doubleReal *z__, doubleReal *alpha, doubleReal *
    beta, doubleReal *u, integer *ldu, doubleReal *vt, integer *ldvt,
    doubleReal *dsigma, doubleReal *u2, integer *ldu2, doubleReal *vt2,
    integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
    idxq, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(dlasd3)(integer *nl, integer *nr, integer *sqre, integer
    *k, doubleReal *d__, doubleReal *q, integer *ldq, doubleReal *dsigma,
    doubleReal *u, integer *ldu, doubleReal *u2, integer *ldu2,
    doubleReal *vt, integer *ldvt, doubleReal *vt2, integer *ldvt2,
    integer *idxc, integer *ctot, doubleReal *z__, integer *info);

/* Subroutine */ int F77NAME(dlasd4)(integer *n, integer *i__, doubleReal *d__,
    doubleReal *z__, doubleReal *delta, doubleReal *rho, doubleReal *
    sigma, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasd5)(integer *i__, doubleReal *d__, doubleReal *z__,
    doubleReal *delta, doubleReal *rho, doubleReal *dsigma, doubleReal *
    work);

/* Subroutine */ int F77NAME(dlasd6)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, doubleReal *d__, doubleReal *vf, doubleReal *vl,
    doubleReal *alpha, doubleReal *beta, integer *idxq, integer *perm,
    integer *givptr, integer *givcol, integer *ldgcol, doubleReal *givnum,
     integer *ldgnum, doubleReal *poles, doubleReal *difl, doubleReal *
    difr, doubleReal *z__, integer *k, doubleReal *c__, doubleReal *s,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlasd7)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *k, doubleReal *d__, doubleReal *z__,
    doubleReal *zw, doubleReal *vf, doubleReal *vfw, doubleReal *vl,
    doubleReal *vlw, doubleReal *alpha, doubleReal *beta, doubleReal *
    dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm,
    integer *givptr, integer *givcol, integer *ldgcol, doubleReal *givnum,
     integer *ldgnum, doubleReal *c__, doubleReal *s, integer *info);

/* Subroutine */ int F77NAME(dlasd8)(integer *icompq, integer *k, doubleReal *d__,
    doubleReal *z__, doubleReal *vf, doubleReal *vl, doubleReal *difl,
    doubleReal *difr, integer *lddifr, doubleReal *dsigma, doubleReal *
    work, integer *info);

/* Subroutine */ int F77NAME(dlasd9)(integer *icompq, integer *ldu, integer *k,
    doubleReal *d__, doubleReal *z__, doubleReal *vf, doubleReal *vl,
    doubleReal *difl, doubleReal *difr, doubleReal *dsigma, doubleReal *
    work, integer *info);

/* Subroutine */ int F77NAME(dlasda)(integer *icompq, integer *smlsiz, integer *n,
    integer *sqre, doubleReal *d__, doubleReal *e, doubleReal *u, integer
    *ldu, doubleReal *vt, integer *k, doubleReal *difl, doubleReal *difr,
    doubleReal *z__, doubleReal *poles, integer *givptr, integer *givcol,
    integer *ldgcol, integer *perm, doubleReal *givnum, doubleReal *c__,
    doubleReal *s, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlasdq)(char *uplo, integer *sqre, integer *n, integer *
    ncvt, integer *nru, integer *ncc, doubleReal *d__, doubleReal *e,
    doubleReal *vt, integer *ldvt, doubleReal *u, integer *ldu,
    doubleReal *c__, integer *ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasdt)(integer *n, integer *lvl, integer *nd, integer *
    inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int F77NAME(dlaset)(char *uplo, integer *m, integer *n, doubleReal *
    alpha, doubleReal *beta, doubleReal *a, integer *lda);

/* Subroutine */ int F77NAME(dlasq1)(integer *n, doubleReal *d__, doubleReal *e,
    doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasq2)(integer *n, doubleReal *z__, integer *info);

/* Subroutine */ int F77NAME(dlasq3)(integer *i0, integer *n0, doubleReal *z__,
    integer *pp, doubleReal *dmin__, doubleReal *sigma, doubleReal *desig,
     doubleReal *qmax, integer *nfail, integer *iter, integer *ndiv,
    logical *ieee);

/* Subroutine */ int F77NAME(dlasq4)(integer *i0, integer *n0, doubleReal *z__,
    integer *pp, integer *n0in, doubleReal *dmin__, doubleReal *dmin1,
    doubleReal *dmin2, doubleReal *dn, doubleReal *dn1, doubleReal *dn2,
    doubleReal *tau, integer *ttype);

/* Subroutine */ int F77NAME(dlasq5)(integer *i0, integer *n0, doubleReal *z__,
    integer *pp, doubleReal *tau, doubleReal *dmin__, doubleReal *dmin1,
    doubleReal *dmin2, doubleReal *dn, doubleReal *dnm1, doubleReal *dnm2,
     logical *ieee);

/* Subroutine */ int F77NAME(dlasq6)(integer *i0, integer *n0, doubleReal *z__,
    integer *pp, doubleReal *dmin__, doubleReal *dmin1, doubleReal *dmin2,
     doubleReal *dn, doubleReal *dnm1, doubleReal *dnm2);

/* Subroutine */ int F77NAME(dlasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, doubleReal *c__, doubleReal *s, doubleReal *a, integer *
    lda);

/* Subroutine */ int F77NAME(dlasrt)(char *id, integer *n, doubleReal *d__, integer *
    info);

/* Subroutine */ int F77NAME(dlassq)(integer *n, doubleReal *x, integer *incx,
    doubleReal *scale, doubleReal *sumsq);

/* Subroutine */ int F77NAME(dlasv2)(doubleReal *f, doubleReal *g, doubleReal *h__,
    doubleReal *ssmin, doubleReal *ssmax, doubleReal *snr, doubleReal *
    csr, doubleReal *snl, doubleReal *csl);

/* Subroutine */ int F77NAME(dlaswp)(integer *n, doubleReal *a, integer *lda, integer
    *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(dlasy2)(logical *ltranl, logical *ltranr, integer *isgn,
    integer *n1, integer *n2, doubleReal *tl, integer *ldtl, doubleReal *
    tr, integer *ldtr, doubleReal *b, integer *ldb, doubleReal *scale,
    doubleReal *x, integer *ldx, doubleReal *xnorm, integer *info);

/* Subroutine */ int F77NAME(dlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     doubleReal *a, integer *lda, integer *ipiv, doubleReal *w, integer *
    ldw, integer *info);

/* Subroutine */ int F77NAME(dlatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, doubleReal *ab, integer *ldab,
    doubleReal *x, doubleReal *scale, doubleReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatdf)(integer *ijob, integer *n, doubleReal *z__,
    integer *ldz, doubleReal *rhs, doubleReal *rdsum, doubleReal *rdscal,
    integer *ipiv, integer *jpiv);

/* Subroutine */ int F77NAME(dlatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleReal *ap, doubleReal *x, doubleReal *scale,
    doubleReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatrd)(char *uplo, integer *n, integer *nb, doubleReal *
    a, integer *lda, doubleReal *e, doubleReal *tau, doubleReal *w,
    integer *ldw);

/* Subroutine */ int F77NAME(dlatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleReal *a, integer *lda, doubleReal *x,
    doubleReal *scale, doubleReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatrz)(integer *m, integer *n, integer *l, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work);

/* Subroutine */ int F77NAME(dlatzm)(char *side, integer *m, integer *n, doubleReal *
    v, integer *incv, doubleReal *tau, doubleReal *c1, doubleReal *c2,
    integer *ldc, doubleReal *work);

/* Subroutine */ int F77NAME(dlauu2)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dlauum)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dopgtr)(char *uplo, integer *n, doubleReal *ap,
    doubleReal *tau, doubleReal *q, integer *ldq, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dopmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doubleReal *ap, doubleReal *tau, doubleReal *c__, integer
    *ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dorg2l)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dorg2r)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dorgbr)(char *vect, integer *m, integer *n, integer *k,
    doubleReal *a, integer *lda, doubleReal *tau, doubleReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorghr)(integer *n, integer *ilo, integer *ihi,
    doubleReal *a, integer *lda, doubleReal *tau, doubleReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorgl2)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dorglq)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgql)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgqr)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgr2)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dorgrq)(integer *m, integer *n, integer *k, doubleReal *
    a, integer *lda, doubleReal *tau, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgtr)(char *uplo, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dorm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dormbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, doubleReal *a, integer *lda, doubleReal *tau,
    doubleReal *c__, integer *ldc, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dormhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, doubleReal *a, integer *lda, doubleReal *
    tau, doubleReal *c__, integer *ldc, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dormlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormql)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dormr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doubleReal *a, integer *lda, doubleReal *tau,
    doubleReal *c__, integer *ldc, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dormrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doubleReal *a, integer *lda, doubleReal *tau,
    doubleReal *c__, integer *ldc, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dormtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doubleReal *a, integer *lda, doubleReal *tau, doubleReal *
    c__, integer *ldc, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dpbcon)(char *uplo, integer *n, integer *kd, doubleReal *
    ab, integer *ldab, doubleReal *anorm, doubleReal *rcond, doubleReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpbequ)(char *uplo, integer *n, integer *kd, doubleReal *
    ab, integer *ldab, doubleReal *s, doubleReal *scond, doubleReal *amax,
     integer *info);

/* Subroutine */ int F77NAME(dpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleReal *ab, integer *ldab, doubleReal *afb, integer *ldafb,
    doubleReal *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *
    ferr, doubleReal *berr, doubleReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dpbstf)(char *uplo, integer *n, integer *kd, doubleReal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleReal *ab, integer *ldab, doubleReal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(dpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, doubleReal *ab, integer *ldab, doubleReal *afb,
    integer *ldafb, char *equed, doubleReal *s, doubleReal *b, integer *
    ldb, doubleReal *x, integer *ldx, doubleReal *rcond, doubleReal *ferr,
     doubleReal *berr, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpbtf2)(char *uplo, integer *n, integer *kd, doubleReal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbtrf)(char *uplo, integer *n, integer *kd, doubleReal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleReal *ab, integer *ldab, doubleReal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(dpocon)(char *uplo, integer *n, doubleReal *a, integer *
    lda, doubleReal *anorm, doubleReal *rcond, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dpoequ)(integer *n, doubleReal *a, integer *lda,
    doubleReal *s, doubleReal *scond, doubleReal *amax, integer *info);

/* Subroutine */ int F77NAME(dporfs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *af, integer *ldaf,
    doubleReal *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *
    ferr, doubleReal *berr, doubleReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dposv)(char *uplo, integer *n, integer *nrhs, doubleReal
    *a, integer *lda, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleReal *a, integer *lda, doubleReal *af, integer *ldaf,
    char *equed, doubleReal *s, doubleReal *b, integer *ldb, doubleReal *
    x, integer *ldx, doubleReal *rcond, doubleReal *ferr, doubleReal *
    berr, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpotf2)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotrf)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotri)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotrs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dppcon)(char *uplo, integer *n, doubleReal *ap,
    doubleReal *anorm, doubleReal *rcond, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dppequ)(char *uplo, integer *n, doubleReal *ap,
    doubleReal *s, doubleReal *scond, doubleReal *amax, integer *info);

/* Subroutine */ int F77NAME(dpprfs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *ap, doubleReal *afp, doubleReal *b, integer *ldb,
    doubleReal *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dppsv)(char *uplo, integer *n, integer *nrhs, doubleReal
    *ap, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleReal *ap, doubleReal *afp, char *equed, doubleReal *s,
    doubleReal *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *
    rcond, doubleReal *ferr, doubleReal *berr, doubleReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dpptrf)(char *uplo, integer *n, doubleReal *ap, integer *
    info);

/* Subroutine */ int F77NAME(dpptri)(char *uplo, integer *n, doubleReal *ap, integer *
    info);

/* Subroutine */ int F77NAME(dpptrs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *ap, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptcon)(integer *n, doubleReal *d__, doubleReal *e,
    doubleReal *anorm, doubleReal *rcond, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dpteqr)(char *compz, integer *n, doubleReal *d__,
    doubleReal *e, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dptrfs)(integer *n, integer *nrhs, doubleReal *d__,
    doubleReal *e, doubleReal *df, doubleReal *ef, doubleReal *b, integer
    *ldb, doubleReal *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
     doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dptsv)(integer *n, integer *nrhs, doubleReal *d__,
    doubleReal *e, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptsvx)(char *fact, integer *n, integer *nrhs,
    doubleReal *d__, doubleReal *e, doubleReal *df, doubleReal *ef,
    doubleReal *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *
    rcond, doubleReal *ferr, doubleReal *berr, doubleReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dpttrf)(integer *n, doubleReal *d__, doubleReal *e,
    integer *info);

/* Subroutine */ int F77NAME(dpttrs)(integer *n, integer *nrhs, doubleReal *d__,
    doubleReal *e, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptts2)(integer *n, integer *nrhs, doubleReal *d__,
    doubleReal *e, doubleReal *b, integer *ldb);

/* Subroutine */ int F77NAME(drscl)(integer *n, doubleReal *sa, doubleReal *sx,
    integer *incx);

/* Subroutine */ int F77NAME(dsbev)(char *jobz, char *uplo, integer *n, integer *kd,
    doubleReal *ab, integer *ldab, doubleReal *w, doubleReal *z__,
    integer *ldz, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dsbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    doubleReal *ab, integer *ldab, doubleReal *w, doubleReal *z__,
    integer *ldz, doubleReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, doubleReal *ab, integer *ldab, doubleReal *q, integer *
    ldq, doubleReal *vl, doubleReal *vu, integer *il, integer *iu,
    doubleReal *abstol, integer *m, doubleReal *w, doubleReal *z__,
    integer *ldz, doubleReal *work, integer *iwork, integer *ifail,
    integer *info);

/* Subroutine */ int F77NAME(dsbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, doubleReal *ab, integer *ldab, doubleReal *bb, integer *
    ldbb, doubleReal *x, integer *ldx, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dsbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doubleReal *ab, integer *ldab, doubleReal *bb, integer *
    ldbb, doubleReal *w, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dsbgvd)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doubleReal *ab, integer *ldab, doubleReal *bb, integer *
    ldbb, doubleReal *w, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, doubleReal *ab, integer *ldab, doubleReal *
    bb, integer *ldbb, doubleReal *q, integer *ldq, doubleReal *vl,
    doubleReal *vu, integer *il, integer *iu, doubleReal *abstol, integer
    *m, doubleReal *w, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    doubleReal *ab, integer *ldab, doubleReal *d__, doubleReal *e,
    doubleReal *q, integer *ldq, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dspcon)(char *uplo, integer *n, doubleReal *ap, integer *
    ipiv, doubleReal *anorm, doubleReal *rcond, doubleReal *work, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(dspev)(char *jobz, char *uplo, integer *n, doubleReal *
    ap, doubleReal *w, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dspevd)(char *jobz, char *uplo, integer *n, doubleReal *
    ap, doubleReal *w, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dspevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleReal *ap, doubleReal *vl, doubleReal *vu, integer *il, integer *
    iu, doubleReal *abstol, integer *m, doubleReal *w, doubleReal *z__,
    integer *ldz, doubleReal *work, integer *iwork, integer *ifail,
    integer *info);

/* Subroutine */ int F77NAME(dspgst)(integer *itype, char *uplo, integer *n,
    doubleReal *ap, doubleReal *bp, integer *info);

/* Subroutine */ int F77NAME(dspgv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleReal *ap, doubleReal *bp, doubleReal *w, doubleReal *z__,
    integer *ldz, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dspgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleReal *ap, doubleReal *bp, doubleReal *w, doubleReal *z__,
    integer *ldz, doubleReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dspgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleReal *ap, doubleReal *bp, doubleReal *vl,
    doubleReal *vu, integer *il, integer *iu, doubleReal *abstol, integer
    *m, doubleReal *w, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsprfs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *ap, doubleReal *afp, integer *ipiv, doubleReal *b,
    integer *ldb, doubleReal *x, integer *ldx, doubleReal *ferr,
    doubleReal *berr, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dspsv)(char *uplo, integer *n, integer *nrhs, doubleReal
    *ap, integer *ipiv, doubleReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleReal *ap, doubleReal *afp, integer *ipiv, doubleReal *b,
    integer *ldb, doubleReal *x, integer *ldx, doubleReal *rcond,
    doubleReal *ferr, doubleReal *berr, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dsptrd)(char *uplo, integer *n, doubleReal *ap,
    doubleReal *d__, doubleReal *e, doubleReal *tau, integer *info);

/* Subroutine */ int F77NAME(dsptrf)(char *uplo, integer *n, doubleReal *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(dsptri)(char *uplo, integer *n, doubleReal *ap, integer *
    ipiv, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dsptrs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *ap, integer *ipiv, doubleReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dstebz)(char *range, char *order, integer *n, doubleReal
    *vl, doubleReal *vu, integer *il, integer *iu, doubleReal *abstol,
    doubleReal *d__, doubleReal *e, integer *m, integer *nsplit,
    doubleReal *w, integer *iblock, integer *isplit, doubleReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dstedc)(char *compz, integer *n, doubleReal *d__,
    doubleReal *e, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstegr)(char *jobz, char *range, integer *n, doubleReal *
    d__, doubleReal *e, doubleReal *vl, doubleReal *vu, integer *il,
    integer *iu, doubleReal *abstol, integer *m, doubleReal *w,
    doubleReal *z__, integer *ldz, integer *isuppz, doubleReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstein)(integer *n, doubleReal *d__, doubleReal *e,
    integer *m, doubleReal *w, integer *iblock, integer *isplit,
    doubleReal *z__, integer *ldz, doubleReal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsteqr)(char *compz, integer *n, doubleReal *d__,
    doubleReal *e, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dsterf)(integer *n, doubleReal *d__, doubleReal *e,
    integer *info);

/* Subroutine */ int F77NAME(dstev)(char *jobz, integer *n, doubleReal *d__,
    doubleReal *e, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dstevd)(char *jobz, integer *n, doubleReal *d__,
    doubleReal *e, doubleReal *z__, integer *ldz, doubleReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstevr)(char *jobz, char *range, integer *n, doubleReal *
    d__, doubleReal *e, doubleReal *vl, doubleReal *vu, integer *il,
    integer *iu, doubleReal *abstol, integer *m, doubleReal *w,
    doubleReal *z__, integer *ldz, integer *isuppz, doubleReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstevx)(char *jobz, char *range, integer *n, doubleReal *
    d__, doubleReal *e, doubleReal *vl, doubleReal *vu, integer *il,
    integer *iu, doubleReal *abstol, integer *m, doubleReal *w,
    doubleReal *z__, integer *ldz, doubleReal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsycon)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *ipiv, doubleReal *anorm, doubleReal *rcond, doubleReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dsyev)(char *jobz, char *uplo, integer *n, doubleReal *a,
     integer *lda, doubleReal *w, doubleReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dsyevd)(char *jobz, char *uplo, integer *n, doubleReal *
    a, integer *lda, doubleReal *w, doubleReal *work, integer *lwork,
    integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsyevr)(char *jobz, char *range, char *uplo, integer *n,
    doubleReal *a, integer *lda, doubleReal *vl, doubleReal *vu, integer *
    il, integer *iu, doubleReal *abstol, integer *m, doubleReal *w,
    doubleReal *z__, integer *ldz, integer *isuppz, doubleReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsyevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleReal *a, integer *lda, doubleReal *vl, doubleReal *vu, integer *
    il, integer *iu, doubleReal *abstol, integer *m, doubleReal *w,
    doubleReal *z__, integer *ldz, doubleReal *work, integer *lwork,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsygs2)(integer *itype, char *uplo, integer *n,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dsygst)(integer *itype, char *uplo, integer *n,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dsygv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleReal *a, integer *lda, doubleReal *b, integer *ldb,
    doubleReal *w, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsygvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleReal *a, integer *lda, doubleReal *b, integer *ldb,
    doubleReal *w, doubleReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsygvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleReal *a, integer *lda, doubleReal *b, integer
    *ldb, doubleReal *vl, doubleReal *vu, integer *il, integer *iu,
    doubleReal *abstol, integer *m, doubleReal *w, doubleReal *z__,
    integer *ldz, doubleReal *work, integer *lwork, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsyrfs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, doubleReal *af, integer *ldaf, integer *
    ipiv, doubleReal *b, integer *ldb, doubleReal *x, integer *ldx,
    doubleReal *ferr, doubleReal *berr, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dsysv)(char *uplo, integer *n, integer *nrhs, doubleReal
    *a, integer *lda, integer *ipiv, doubleReal *b, integer *ldb,
    doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleReal *a, integer *lda, doubleReal *af, integer *ldaf,
    integer *ipiv, doubleReal *b, integer *ldb, doubleReal *x, integer *
    ldx, doubleReal *rcond, doubleReal *ferr, doubleReal *berr,
    doubleReal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dsytd2)(char *uplo, integer *n, doubleReal *a, integer *
    lda, doubleReal *d__, doubleReal *e, doubleReal *tau, integer *info);

/* Subroutine */ int F77NAME(dsytf2)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dsytrd)(char *uplo, integer *n, doubleReal *a, integer *
    lda, doubleReal *d__, doubleReal *e, doubleReal *tau, doubleReal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsytrf)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *ipiv, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsytri)(char *uplo, integer *n, doubleReal *a, integer *
    lda, integer *ipiv, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dsytrs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *a, integer *lda, integer *ipiv, doubleReal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dtbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, doubleReal *ab, integer *ldab, doubleReal *rcond,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleReal *ab, integer *ldab, doubleReal
    *b, integer *ldb, doubleReal *x, integer *ldx, doubleReal *ferr,
    doubleReal *berr, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleReal *ab, integer *ldab, doubleReal
    *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dtgevc)(char *side, char *howmny, logical *select,
    integer *n, doubleReal *a, integer *lda, doubleReal *b, integer *ldb,
    doubleReal *vl, integer *ldvl, doubleReal *vr, integer *ldvr, integer
    *mm, integer *m, doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dtgex2)(logical *wantq, logical *wantz, integer *n,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, doubleReal *
    q, integer *ldq, doubleReal *z__, integer *ldz, integer *j1, integer *
    n1, integer *n2, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dtgexc)(logical *wantq, logical *wantz, integer *n,
    doubleReal *a, integer *lda, doubleReal *b, integer *ldb, doubleReal *
    q, integer *ldq, doubleReal *z__, integer *ldz, integer *ifst,
    integer *ilst, doubleReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dtgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, doubleReal *a, integer *lda, doubleReal *
    b, integer *ldb, doubleReal *alphar, doubleReal *alphai, doubleReal *
    beta, doubleReal *q, integer *ldq, doubleReal *z__, integer *ldz,
    integer *m, doubleReal *pl, doubleReal *pr, doubleReal *dif,
    doubleReal *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(dtgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, doubleReal *a,
    integer *lda, doubleReal *b, integer *ldb, doubleReal *tola,
    doubleReal *tolb, doubleReal *alpha, doubleReal *beta, doubleReal *u,
    integer *ldu, doubleReal *v, integer *ldv, doubleReal *q, integer *
    ldq, doubleReal *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(dtgsna)(char *job, char *howmny, logical *select,
    integer *n, doubleReal *a, integer *lda, doubleReal *b, integer *ldb,
    doubleReal *vl, integer *ldvl, doubleReal *vr, integer *ldvr,
    doubleReal *s, doubleReal *dif, integer *mm, integer *m, doubleReal *
    work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, doubleReal *a, integer *lda, doubleReal *b, integer *ldb,
    doubleReal *c__, integer *ldc, doubleReal *d__, integer *ldd,
    doubleReal *e, integer *lde, doubleReal *f, integer *ldf, doubleReal *
    scale, doubleReal *rdsum, doubleReal *rdscal, integer *iwork, integer
    *pq, integer *info);

/* Subroutine */ int F77NAME(dtgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, doubleReal *a, integer *lda, doubleReal *b, integer *ldb,
    doubleReal *c__, integer *ldc, doubleReal *d__, integer *ldd,
    doubleReal *e, integer *lde, doubleReal *f, integer *ldf, doubleReal *
    scale, doubleReal *dif, doubleReal *work, integer *lwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dtpcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleReal *ap, doubleReal *rcond, doubleReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dtprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleReal *ap, doubleReal *b, integer *ldb,
    doubleReal *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtptri)(char *uplo, char *diag, integer *n, doubleReal *
    ap, integer *info);

/* Subroutine */ int F77NAME(dtptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleReal *ap, doubleReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dtrcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleReal *a, integer *lda, doubleReal *rcond, doubleReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtrevc)(char *side, char *howmny, logical *select,
    integer *n, doubleReal *t, integer *ldt, doubleReal *vl, integer *
    ldvl, doubleReal *vr, integer *ldvr, integer *mm, integer *m,
    doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dtrexc)(char *compq, integer *n, doubleReal *t, integer *
    ldt, doubleReal *q, integer *ldq, integer *ifst, integer *ilst,
    doubleReal *work, integer *info);

/* Subroutine */ int F77NAME(dtrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleReal *a, integer *lda, doubleReal *b, integer *
    ldb, doubleReal *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtrsen)(char *job, char *compq, logical *select, integer
    *n, doubleReal *t, integer *ldt, doubleReal *q, integer *ldq,
    doubleReal *wr, doubleReal *wi, integer *m, doubleReal *s, doubleReal
    *sep, doubleReal *work, integer *lwork, integer *iwork, integer *
    liwork, integer *info);

/* Subroutine */ int F77NAME(dtrsna)(char *job, char *howmny, logical *select,
    integer *n, doubleReal *t, integer *ldt, doubleReal *vl, integer *
    ldvl, doubleReal *vr, integer *ldvr, doubleReal *s, doubleReal *sep,
    integer *mm, integer *m, doubleReal *work, integer *ldwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dtrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, doubleReal *a, integer *lda, doubleReal *b, integer *
    ldb, doubleReal *c__, integer *ldc, doubleReal *scale, integer *info);

/* Subroutine */ int F77NAME(dtrti2)(char *uplo, char *diag, integer *n, doubleReal *
    a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(dtrtri)(char *uplo, char *diag, integer *n, doubleReal *
    a, integer *lda, integer *info);

/**
  Lapack routine : DTRTRS solves a triangular system of the form
  A * X = B  or  A**T * X = B,
  where A is a triangular matrix of order N, and B is an N-by-NRHS
  matrix.  A check is made to verify that A is nonsingular.
*/
int F77NAME(dtrtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleReal *a, integer *lda, doubleReal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dtzrqf)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, integer *info);

/* Subroutine */ int F77NAME(dtzrzf)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleReal *tau, doubleReal *work, integer *lwork, integer *info);

integer F77NAME(icmax1)(integer *n, Complex *cx, integer *incx);

integer F77NAME(ieeeck)(integer *ispec, Real *zero, Real *one);

integer F77NAME(ilaenv)(integer *ispec, char *name__, char *opts, integer *n1,
    integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen
    opts_len);

integer F77NAME(izmax1)(integer *n, doubleComplex *cx, integer *incx);

/* Subroutine */ int F77NAME(sbdsdc)(char *uplo, char *compq, integer *n, Real *d__,
    Real *e, Real *u, integer *ldu, Real *vt, integer *ldvt, Real *q,
    integer *iq, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, Real *d__, Real *e, Real *vt, integer *ldvt, Real *
    u, integer *ldu, Real *c__, integer *ldc, Real *work, integer *info);

/* Subroutine */ int F77NAME(sdisna)(char *job, integer *m, integer *n, Real *d__,
    Real *sep, integer *info);

/* Subroutine */ int F77NAME(sgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, Real *ab, integer *ldab, Real *d__, Real *
    e, Real *q, integer *ldq, Real *pt, integer *ldpt, Real *c__, integer
    *ldc, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     Real *ab, integer *ldab, integer *ipiv, Real *anorm, Real *rcond,
    Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     Real *ab, integer *ldab, Real *r__, Real *c__, Real *rowcnd, Real *
    colcnd, Real *amax, integer *info);

/* Subroutine */ int F77NAME(sgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, Real *ab, integer *ldab, Real *afb, integer *ldafb,
     integer *ipiv, Real *b, integer *ldb, Real *x, integer *ldx, Real *
    ferr, Real *berr, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, Real *ab, integer *ldab, integer *ipiv, Real *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(sgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, Real *ab, integer *ldab, Real *afb,
    integer *ldafb, integer *ipiv, char *equed, Real *r__, Real *c__,
    Real *b, integer *ldb, Real *x, integer *ldx, Real *rcond, Real *ferr,
     Real *berr, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     Real *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     Real *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, Real *ab, integer *ldab, integer *ipiv, Real *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, Real *scale, integer *m, Real *v, integer *ldv, integer
    *info);

/* Subroutine */ int F77NAME(sgebal)(char *job, integer *n, Real *a, integer *lda,
    integer *ilo, integer *ihi, Real *scale, integer *info);

/* Subroutine */ int F77NAME(sgebd2)(integer *m, integer *n, Real *a, integer *lda,
    Real *d__, Real *e, Real *tauq, Real *taup, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgebrd)(integer *m, integer *n, Real *a, integer *lda,
    Real *d__, Real *e, Real *tauq, Real *taup, Real *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(sgecon)(char *norm, integer *n, Real *a, integer *lda,
    Real *anorm, Real *rcond, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgeequ)(integer *m, integer *n, Real *a, integer *lda,
    Real *r__, Real *c__, Real *rowcnd, Real *colcnd, Real *amax, integer
    *info);

/* Subroutine */ int F77NAME(sgees)(char *jobvs, char *sort, L_fp select, integer *n,
    Real *a, integer *lda, integer *sdim, Real *wr, Real *wi, Real *vs,
    integer *ldvs, Real *work, integer *lwork, logical *bwork, integer *
    info);

/* Subroutine */ int F77NAME(sgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, Real *a, integer *lda, integer *sdim, Real *wr,
    Real *wi, Real *vs, integer *ldvs, Real *rconde, Real *rcondv, Real *
    work, integer *lwork, integer *iwork, integer *liwork, logical *bwork,
     integer *info);

/* Subroutine */ int F77NAME(sgeev)(char *jobvl, char *jobvr, integer *n, Real *a,
    integer *lda, Real *wr, Real *wi, Real *vl, integer *ldvl, Real *vr,
    integer *ldvr, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, Real *a, integer *lda, Real *wr, Real *wi, Real *
    vl, integer *ldvl, Real *vr, integer *ldvr, integer *ilo, integer *
    ihi, Real *scale, Real *abnrm, Real *rconde, Real *rcondv, Real *work,
     integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgegs)(char *jobvsl, char *jobvsr, integer *n, Real *a,
    integer *lda, Real *b, integer *ldb, Real *alphar, Real *alphai, Real
    *beta, Real *vsl, integer *ldvsl, Real *vsr, integer *ldvsr, Real *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgegv)(char *jobvl, char *jobvr, integer *n, Real *a,
    integer *lda, Real *b, integer *ldb, Real *alphar, Real *alphai, Real
    *beta, Real *vl, integer *ldvl, Real *vr, integer *ldvr, Real *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgehd2)(integer *n, integer *ilo, integer *ihi, Real *a,
    integer *lda, Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgehrd)(integer *n, integer *ilo, integer *ihi, Real *a,
    integer *lda, Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelq2)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgelqf)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgels)(char *trans, integer *m, integer *n, integer *
    nrhs, Real *a, integer *lda, Real *b, integer *ldb, Real *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelsd)(integer *m, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *b, integer *ldb, Real *s, Real *rcond, integer *
    rank, Real *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgelss)(integer *m, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *b, integer *ldb, Real *s, Real *rcond, integer *
    rank, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelsx)(integer *m, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *b, integer *ldb, integer *jpvt, Real *rcond,
    integer *rank, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgelsy)(integer *m, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *b, integer *ldb, integer *jpvt, Real *rcond,
    integer *rank, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeql2)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgeqlf)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeqp3)(integer *m, integer *n, Real *a, integer *lda,
    integer *jpvt, Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeqpf)(integer *m, integer *n, Real *a, integer *lda,
    integer *jpvt, Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgeqr2)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgeqrf)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgerfs)(char *trans, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *af, integer *ldaf, integer *ipiv, Real *b,
    integer *ldb, Real *x, integer *ldx, Real *ferr, Real *berr, Real *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgerq2)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgerqf)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgesc2)(integer *n, Real *a, integer *lda, Real *rhs,
    integer *ipiv, integer *jpiv, Real *scale);

/* Subroutine */ int F77NAME(sgesdd)(char *jobz, integer *m, integer *n, Real *a,
    integer *lda, Real *s, Real *u, integer *ldu, Real *vt, integer *ldvt,
     Real *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgesv)(integer *n, integer *nrhs, Real *a, integer *lda,
    integer *ipiv, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgesvd)(char *jobu, char *jobvt, integer *m, integer *n,
    Real *a, integer *lda, Real *s, Real *u, integer *ldu, Real *vt,
    integer *ldvt, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, Real *a, integer *lda, Real *af, integer *ldaf, integer *ipiv,
    char *equed, Real *r__, Real *c__, Real *b, integer *ldb, Real *x,
    integer *ldx, Real *rcond, Real *ferr, Real *berr, Real *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgetc2)(integer *n, Real *a, integer *lda, integer *ipiv,
     integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(sgetf2)(integer *m, integer *n, Real *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgetrf)(integer *m, integer *n, Real *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgetri)(integer *n, Real *a, integer *lda, integer *ipiv,
     Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgetrs)(char *trans, integer *n, integer *nrhs, Real *a,
    integer *lda, integer *ipiv, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, Real *lscale, Real *rscale, integer *m, Real *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(sggbal)(char *job, integer *n, Real *a, integer *lda,
    Real *b, integer *ldb, integer *ilo, integer *ihi, Real *lscale, Real
    *rscale, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, integer *n, Real *a, integer *lda, Real *b, integer *ldb,
    integer *sdim, Real *alphar, Real *alphai, Real *beta, Real *vsl,
    integer *ldvsl, Real *vsr, integer *ldvsr, Real *work, integer *lwork,
     logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, char *sense, integer *n, Real *a, integer *lda, Real *b,
    integer *ldb, integer *sdim, Real *alphar, Real *alphai, Real *beta,
    Real *vsl, integer *ldvsl, Real *vsr, integer *ldvsr, Real *rconde,
    Real *rcondv, Real *work, integer *lwork, integer *iwork, integer *
    liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggev)(char *jobvl, char *jobvr, integer *n, Real *a,
    integer *lda, Real *b, integer *ldb, Real *alphar, Real *alphai, Real
    *beta, Real *vl, integer *ldvl, Real *vr, integer *ldvr, Real *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, Real *a, integer *lda, Real *b, integer *ldb, Real
    *alphar, Real *alphai, Real *beta, Real *vl, integer *ldvl, Real *vr,
    integer *ldvr, integer *ilo, integer *ihi, Real *lscale, Real *rscale,
     Real *abnrm, Real *bbnrm, Real *rconde, Real *rcondv, Real *work,
    integer *lwork, integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggglm)(integer *n, integer *m, integer *p, Real *a,
    integer *lda, Real *b, integer *ldb, Real *d__, Real *x, Real *y,
    Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, Real *a, integer *lda, Real *b, integer *ldb, Real
    *q, integer *ldq, Real *z__, integer *ldz, integer *info);

/* Subroutine */ int F77NAME(sgglse)(integer *m, integer *n, integer *p, Real *a,
    integer *lda, Real *b, integer *ldb, Real *c__, Real *d__, Real *x,
    Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggqrf)(integer *n, integer *m, integer *p, Real *a,
    integer *lda, Real *taua, Real *b, integer *ldb, Real *taub, Real *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggrqf)(integer *m, integer *p, integer *n, Real *a,
    integer *lda, Real *taua, Real *b, integer *ldb, Real *taub, Real *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, Real *a, integer *lda,
     Real *b, integer *ldb, Real *alpha, Real *beta, Real *u, integer *
    ldu, Real *v, integer *ldv, Real *q, integer *ldq, Real *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, Real *a, integer *lda, Real *b, integer *ldb,
    Real *tola, Real *tolb, integer *k, integer *l, Real *u, integer *ldu,
     Real *v, integer *ldv, Real *q, integer *ldq, integer *iwork, Real *
    tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sgtcon)(char *norm, integer *n, Real *dl, Real *d__,
    Real *du, Real *du2, integer *ipiv, Real *anorm, Real *rcond, Real *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgtrfs)(char *trans, integer *n, integer *nrhs, Real *dl,
     Real *d__, Real *du, Real *dlf, Real *df, Real *duf, Real *du2,
    integer *ipiv, Real *b, integer *ldb, Real *x, integer *ldx, Real *
    ferr, Real *berr, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgtsv)(integer *n, integer *nrhs, Real *dl, Real *d__,
    Real *du, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, Real *dl, Real *d__, Real *du, Real *dlf, Real *df, Real *duf,
    Real *du2, integer *ipiv, Real *b, integer *ldb, Real *x, integer *
    ldx, Real *rcond, Real *ferr, Real *berr, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(sgttrf)(integer *n, Real *dl, Real *d__, Real *du, Real *
    du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgttrs)(char *trans, integer *n, integer *nrhs, Real *dl,
     Real *d__, Real *du, Real *du2, integer *ipiv, Real *b, integer *ldb,
     integer *info);

/* Subroutine */ int F77NAME(sgtts2)(integer *itrans, integer *n, integer *nrhs, Real
    *dl, Real *d__, Real *du, Real *du2, integer *ipiv, Real *b, integer *
    ldb);

/* Subroutine */ int F77NAME(shgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, Real *a, integer *lda, Real *b, integer *
    ldb, Real *alphar, Real *alphai, Real *beta, Real *q, integer *ldq,
    Real *z__, integer *ldz, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(shsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, Real *h__, integer *ldh, Real *wr, Real *wi, Real
    *vl, integer *ldvl, Real *vr, integer *ldvr, integer *mm, integer *m,
    Real *work, integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(shseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, Real *h__, integer *ldh, Real *wr, Real *wi, Real *z__,
     integer *ldz, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(slabad)(Real *small, Real *large);

/* Subroutine */ int F77NAME(slabrd)(integer *m, integer *n, integer *nb, Real *a,
    integer *lda, Real *d__, Real *e, Real *tauq, Real *taup, Real *x,
    integer *ldx, Real *y, integer *ldy);

/* Subroutine */ int F77NAME(slacon)(integer *n, Real *v, Real *x, integer *isgn,
    Real *est, integer *kase);

/* Subroutine */ int F77NAME(slacpy)(char *uplo, integer *m, integer *n, Real *a,
    integer *lda, Real *b, integer *ldb);

/* Subroutine */ int F77NAME(sladiv)(Real *a, Real *b, Real *c__, Real *d__, Real *p,
    Real *q);

/* Subroutine */ int F77NAME(slae2)(Real *a, Real *b, Real *c__, Real *rt1, Real *rt2);

/* Subroutine */ int F77NAME(slaebz)(integer *ijob, integer *nitmax, integer *n,
    integer *mmax, integer *minp, integer *nbmin, Real *abstol, Real *
    reltol, Real *pivmin, Real *d__, Real *e, Real *e2, integer *nval,
    Real *ab, Real *c__, integer *mout, integer *nab, Real *work, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(slaed0)(integer *icompq, integer *qsiz, integer *n, Real
    *d__, Real *e, Real *q, integer *ldq, Real *qstore, integer *ldqs,
    Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slaed1)(integer *n, Real *d__, Real *q, integer *ldq,
    integer *indxq, Real *rho, integer *cutpnt, Real *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(slaed2)(integer *k, integer *n, integer *n1, Real *d__,
    Real *q, integer *ldq, integer *indxq, Real *rho, Real *z__, Real *
    dlamda, Real *w, Real *q2, integer *indx, integer *indxc, integer *
    indxp, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(slaed3)(integer *k, integer *n, integer *n1, Real *d__,
    Real *q, integer *ldq, Real *rho, Real *dlamda, Real *q2, integer *
    indx, integer *ctot, Real *w, Real *s, integer *info);

/* Subroutine */ int F77NAME(slaed4)(integer *n, integer *i__, Real *d__, Real *z__,
    Real *delta, Real *rho, Real *dlam, integer *info);

/* Subroutine */ int F77NAME(slaed5)(integer *i__, Real *d__, Real *z__, Real *delta,
    Real *rho, Real *dlam);

/* Subroutine */ int F77NAME(slaed6)(integer *kniter, logical *orgati, Real *rho,
    Real *d__, Real *z__, Real *finit, Real *tau, integer *info);

/* Subroutine */ int F77NAME(slaed7)(integer *icompq, integer *n, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, Real *d__, Real *q,
    integer *ldq, integer *indxq, Real *rho, integer *cutpnt, Real *
    qstore, integer *qptr, integer *prmptr, integer *perm, integer *
    givptr, integer *givcol, Real *givnum, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(slaed8)(integer *icompq, integer *k, integer *n, integer
    *qsiz, Real *d__, Real *q, integer *ldq, integer *indxq, Real *rho,
    integer *cutpnt, Real *z__, Real *dlamda, Real *q2, integer *ldq2,
    Real *w, integer *perm, integer *givptr, integer *givcol, Real *
    givnum, integer *indxp, integer *indx, integer *info);

/* Subroutine */ int F77NAME(slaed9)(integer *k, integer *kstart, integer *kstop,
    integer *n, Real *d__, Real *q, integer *ldq, Real *rho, Real *dlamda,
     Real *w, Real *s, integer *lds, integer *info);

/* Subroutine */ int F77NAME(slaeda)(integer *n, integer *tlvls, integer *curlvl,
    integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
    integer *givcol, Real *givnum, Real *q, integer *qptr, Real *z__,
    Real *ztemp, integer *info);

/* Subroutine */ int F77NAME(slaein)(logical *rightv, logical *noinit, integer *n,
    Real *h__, integer *ldh, Real *wr, Real *wi, Real *vr, Real *vi, Real
    *b, integer *ldb, Real *work, Real *eps3, Real *smlnum, Real *bignum,
    integer *info);

/* Subroutine */ int F77NAME(slaev2)(Real *a, Real *b, Real *c__, Real *rt1, Real *
    rt2, Real *cs1, Real *sn1);

/* Subroutine */ int F77NAME(slaexc)(logical *wantq, integer *n, Real *t, integer *
    ldt, Real *q, integer *ldq, integer *j1, integer *n1, integer *n2,
    Real *work, integer *info);

/* Subroutine */ int F77NAME(slag2)(Real *a, integer *lda, Real *b, integer *ldb,
    Real *safmin, Real *scale1, Real *scale2, Real *wr1, Real *wr2, Real *
    wi);

/* Subroutine */ int F77NAME(slags2)(logical *upper, Real *a1, Real *a2, Real *a3,
    Real *b1, Real *b2, Real *b3, Real *csu, Real *snu, Real *csv, Real *
    snv, Real *csq, Real *snq);

/* Subroutine */ int F77NAME(slagtf)(integer *n, Real *a, Real *lambda, Real *b, Real
    *c__, Real *tol, Real *d__, integer *in, integer *info);

/* Subroutine */ int F77NAME(slagtm)(char *trans, integer *n, integer *nrhs, Real *
    alpha, Real *dl, Real *d__, Real *du, Real *x, integer *ldx, Real *
    beta, Real *b, integer *ldb);

/* Subroutine */ int F77NAME(slagts)(integer *job, integer *n, Real *a, Real *b, Real
    *c__, Real *d__, integer *in, Real *y, Real *tol, integer *info);

/* Subroutine */ int F77NAME(slagv2)(Real *a, integer *lda, Real *b, integer *ldb,
    Real *alphar, Real *alphai, Real *beta, Real *csl, Real *snl, Real *
    csr, Real *snr);

/* Subroutine */ int F77NAME(slahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, Real *h__, integer *ldh, Real *wr, Real *
    wi, integer *iloz, integer *ihiz, Real *z__, integer *ldz, integer *
    info);

/* Subroutine */ int F77NAME(slahrd)(integer *n, integer *k, integer *nb, Real *a,
    integer *lda, Real *tau, Real *t, integer *ldt, Real *y, integer *ldy);

/* Subroutine */ int F77NAME(slaic1)(integer *job, integer *j, Real *x, Real *sest,
    Real *w, Real *gamma, Real *sestpr, Real *s, Real *c__);

/* Subroutine */ int F77NAME(slaln2)(logical *ltrans, integer *na, integer *nw, Real *
    smin, Real *ca, Real *a, integer *lda, Real *d1, Real *d2, Real *b,
    integer *ldb, Real *wr, Real *wi, Real *x, integer *ldx, Real *scale,
    Real *xnorm, integer *info);

/* Subroutine */ int F77NAME(slals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, Real *b, integer *ldb, Real *bx,
    integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, Real *givnum, integer *ldgnum, Real *poles, Real *
    difl, Real *difr, Real *z__, integer *k, Real *c__, Real *s, Real *
    work, integer *info);

/* Subroutine */ int F77NAME(slalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, Real *b, integer *ldb, Real *bx, integer *ldbx, Real *
    u, integer *ldu, Real *vt, integer *k, Real *difl, Real *difr, Real *
    z__, Real *poles, integer *givptr, integer *givcol, integer *ldgcol,
    integer *perm, Real *givnum, Real *c__, Real *s, Real *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(slalsd)(char *uplo, integer *smlsiz, integer *n, integer
    *nrhs, Real *d__, Real *e, Real *b, integer *ldb, Real *rcond,
    integer *rank, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slamc1)(integer *beta, integer *t, logical *rnd, logical
    *ieee1);

/* Subroutine */ int F77NAME(slamc2)(integer *beta, integer *t, logical *rnd, Real *
    eps, integer *emin, Real *rmin, integer *emax, Real *rmax);

/* Subroutine */ int F77NAME(slamc4)(integer *emin, Real *start, integer *base);

/* Subroutine */ int F77NAME(slamc5)(integer *beta, integer *p, integer *emin,
    logical *ieee, integer *emax, Real *rmax);

/* Subroutine */ int F77NAME(slamrg)(integer *n1, integer *n2, Real *a, integer *
    strd1, integer *strd2, integer *index);

/* Subroutine */ int F77NAME(slanv2)(Real *a, Real *b, Real *c__, Real *d__, Real *
    rt1r, Real *rt1i, Real *rt2r, Real *rt2i, Real *cs, Real *sn);

/* Subroutine */ int F77NAME(slapll)(integer *n, Real *x, integer *incx, Real *y,
    integer *incy, Real *ssmin);

/* Subroutine */ int F77NAME(slapmt)(logical *forwrd, integer *m, integer *n, Real *x,
     integer *ldx, integer *k);

/* Subroutine */ int F77NAME(slaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     Real *ab, integer *ldab, Real *r__, Real *c__, Real *rowcnd, Real *
    colcnd, Real *amax, char *equed);

/* Subroutine */ int F77NAME(slaqge)(integer *m, integer *n, Real *a, integer *lda,
    Real *r__, Real *c__, Real *rowcnd, Real *colcnd, Real *amax, char *
    equed);

/* Subroutine */ int F77NAME(slaqp2)(integer *m, integer *n, integer *offset, Real *a,
     integer *lda, integer *jpvt, Real *tau, Real *vn1, Real *vn2, Real *
    work);

/* Subroutine */ int F77NAME(slaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, Real *a, integer *lda, integer *jpvt, Real *tau,
    Real *vn1, Real *vn2, Real *auxv, Real *f, integer *ldf);

/* Subroutine */ int F77NAME(slaqsb)(char *uplo, integer *n, integer *kd, Real *ab,
    integer *ldab, Real *s, Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(slaqsp)(char *uplo, integer *n, Real *ap, Real *s, Real *
    scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(slaqsy)(char *uplo, integer *n, Real *a, integer *lda,
    Real *s, Real *scond, Real *amax, char *equed);

/* Subroutine */ int F77NAME(slaqtr)(logical *ltran, logical *lReal, integer *n, Real
    *t, integer *ldt, Real *b, Real *w, Real *scale, Real *x, Real *work,
    integer *info);

/* Subroutine */ int F77NAME(slar1v)(integer *n, integer *b1, integer *bn, Real *
    sigma, Real *d__, Real *l, Real *ld, Real *lld, Real *gersch, Real *
    z__, Real *ztz, Real *mingma, integer *r__, integer *isuppz, Real *
    work);

/* Subroutine */ int F77NAME(slar2v)(integer *n, Real *x, Real *y, Real *z__, integer
    *incx, Real *c__, Real *s, integer *incc);

/* Subroutine */ int F77NAME(slarf)(char *side, integer *m, integer *n, Real *v,
    integer *incv, Real *tau, Real *c__, integer *ldc, Real *work);

/* Subroutine */ int F77NAME(slarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, Real *v, integer *ldv,
    Real *t, integer *ldt, Real *c__, integer *ldc, Real *work, integer *
    ldwork);

/* Subroutine */ int F77NAME(slarfg)(integer *n, Real *alpha, Real *x, integer *incx,
    Real *tau);

/* Subroutine */ int F77NAME(slarft)(char *direct, char *storev, integer *n, integer *
    k, Real *v, integer *ldv, Real *tau, Real *t, integer *ldt);

/* Subroutine */ int F77NAME(slarfx)(char *side, integer *m, integer *n, Real *v,
    Real *tau, Real *c__, integer *ldc, Real *work);

/* Subroutine */ int F77NAME(slargv)(integer *n, Real *x, integer *incx, Real *y,
    integer *incy, Real *c__, integer *incc);

/* Subroutine */ int F77NAME(slarnv)(integer *idist, integer *iseed, integer *n, Real
    *x);

/* Subroutine */ int F77NAME(slarrb)(integer *n, Real *d__, Real *l, Real *ld, Real *
    lld, integer *ifirst, integer *ilast, Real *sigma, Real *reltol, Real
    *w, Real *wgap, Real *werr, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slarre)(integer *n, Real *d__, Real *e, Real *tol,
    integer *nsplit, integer *isplit, integer *m, Real *w, Real *woff,
    Real *gersch, Real *work, integer *info);

/* Subroutine */ int F77NAME(slarrf)(integer *n, Real *d__, Real *l, Real *ld, Real *
    lld, integer *ifirst, integer *ilast, Real *w, Real *dplus, Real *
    lplus, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slarrv)(integer *n, Real *d__, Real *l, integer *isplit,
    integer *m, Real *w, integer *iblock, Real *gersch, Real *tol, Real *
    z__, integer *ldz, integer *isuppz, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(slartg)(Real *f, Real *g, Real *cs, Real *sn, Real *r__);

/* Subroutine */ int F77NAME(slartv)(integer *n, Real *x, integer *incx, Real *y,
    integer *incy, Real *c__, Real *s, integer *incc);

/* Subroutine */ int F77NAME(slaruv)(integer *iseed, integer *n, Real *x);

/* Subroutine */ int F77NAME(slarz)(char *side, integer *m, integer *n, integer *l,
    Real *v, integer *incv, Real *tau, Real *c__, integer *ldc, Real *
    work);

/* Subroutine */ int F77NAME(slarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, Real *v,
    integer *ldv, Real *t, integer *ldt, Real *c__, integer *ldc, Real *
    work, integer *ldwork);

/* Subroutine */ int F77NAME(slarzt)(char *direct, char *storev, integer *n, integer *
    k, Real *v, integer *ldv, Real *tau, Real *t, integer *ldt);

/* Subroutine */ int F77NAME(slas2)(Real *f, Real *g, Real *h__, Real *ssmin, Real *
    ssmax);

/* Subroutine */ int F77NAME(slascl)(char *type__, integer *kl, integer *ku, Real *
    cfrom, Real *cto, integer *m, integer *n, Real *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(slasd0)(integer *n, integer *sqre, Real *d__, Real *e,
    Real *u, integer *ldu, Real *vt, integer *ldvt, integer *smlsiz,
    integer *iwork, Real *work, integer *info);

/* Subroutine */ int F77NAME(slasd1)(integer *nl, integer *nr, integer *sqre, Real *
    d__, Real *alpha, Real *beta, Real *u, integer *ldu, Real *vt,
    integer *ldvt, integer *idxq, integer *iwork, Real *work, integer *
    info);

/* Subroutine */ int F77NAME(slasd2)(integer *nl, integer *nr, integer *sqre, integer
    *k, Real *d__, Real *z__, Real *alpha, Real *beta, Real *u, integer *
    ldu, Real *vt, integer *ldvt, Real *dsigma, Real *u2, integer *ldu2,
    Real *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc,
     integer *idxq, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(slasd3)(integer *nl, integer *nr, integer *sqre, integer
    *k, Real *d__, Real *q, integer *ldq, Real *dsigma, Real *u, integer *
    ldu, Real *u2, integer *ldu2, Real *vt, integer *ldvt, Real *vt2,
    integer *ldvt2, integer *idxc, integer *ctot, Real *z__, integer *
    info);

/* Subroutine */ int F77NAME(slasd4)(integer *n, integer *i__, Real *d__, Real *z__,
    Real *delta, Real *rho, Real *sigma, Real *work, integer *info);

/* Subroutine */ int F77NAME(slasd5)(integer *i__, Real *d__, Real *z__, Real *delta,
    Real *rho, Real *dsigma, Real *work);

/* Subroutine */ int F77NAME(slasd6)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, Real *d__, Real *vf, Real *vl, Real *alpha, Real *beta,
     integer *idxq, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, Real *givnum, integer *ldgnum, Real *poles, Real *
    difl, Real *difr, Real *z__, integer *k, Real *c__, Real *s, Real *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slasd7)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *k, Real *d__, Real *z__, Real *zw, Real *vf,
    Real *vfw, Real *vl, Real *vlw, Real *alpha, Real *beta, Real *dsigma,
     integer *idx, integer *idxp, integer *idxq, integer *perm, integer *
    givptr, integer *givcol, integer *ldgcol, Real *givnum, integer *
    ldgnum, Real *c__, Real *s, integer *info);

/* Subroutine */ int F77NAME(slasd8)(integer *icompq, integer *k, Real *d__, Real *
    z__, Real *vf, Real *vl, Real *difl, Real *difr, integer *lddifr,
    Real *dsigma, Real *work, integer *info);

/* Subroutine */ int F77NAME(slasd9)(integer *icompq, integer *ldu, integer *k, Real *
    d__, Real *z__, Real *vf, Real *vl, Real *difl, Real *difr, Real *
    dsigma, Real *work, integer *info);

/* Subroutine */ int F77NAME(slasda)(integer *icompq, integer *smlsiz, integer *n,
    integer *sqre, Real *d__, Real *e, Real *u, integer *ldu, Real *vt,
    integer *k, Real *difl, Real *difr, Real *z__, Real *poles, integer *
    givptr, integer *givcol, integer *ldgcol, integer *perm, Real *givnum,
     Real *c__, Real *s, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slasdq)(char *uplo, integer *sqre, integer *n, integer *
    ncvt, integer *nru, integer *ncc, Real *d__, Real *e, Real *vt,
    integer *ldvt, Real *u, integer *ldu, Real *c__, integer *ldc, Real *
    work, integer *info);

/* Subroutine */ int F77NAME(slasdt)(integer *n, integer *lvl, integer *nd, integer *
    inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int F77NAME(slaset)(char *uplo, integer *m, integer *n, Real *alpha,
    Real *beta, Real *a, integer *lda);

/* Subroutine */ int F77NAME(slasq1)(integer *n, Real *d__, Real *e, Real *work,
    integer *info);

/* Subroutine */ int F77NAME(slasq2)(integer *n, Real *z__, integer *info);

/* Subroutine */ int F77NAME(slasq3)(integer *i0, integer *n0, Real *z__, integer *pp,
     Real *dmin__, Real *sigma, Real *desig, Real *qmax, integer *nfail,
    integer *iter, integer *ndiv, logical *ieee);

/* Subroutine */ int F77NAME(slasq4)(integer *i0, integer *n0, Real *z__, integer *pp,
     integer *n0in, Real *dmin__, Real *dmin1, Real *dmin2, Real *dn,
    Real *dn1, Real *dn2, Real *tau, integer *ttype);

/* Subroutine */ int F77NAME(slasq5)(integer *i0, integer *n0, Real *z__, integer *pp,
     Real *tau, Real *dmin__, Real *dmin1, Real *dmin2, Real *dn, Real *
    dnm1, Real *dnm2, logical *ieee);

/* Subroutine */ int F77NAME(slasq6)(integer *i0, integer *n0, Real *z__, integer *pp,
     Real *dmin__, Real *dmin1, Real *dmin2, Real *dn, Real *dnm1, Real *
    dnm2);

/* Subroutine */ int F77NAME(slasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, Real *c__, Real *s, Real *a, integer *lda);

/* Subroutine */ int F77NAME(slasrt)(char *id, integer *n, Real *d__, integer *info);

/* Subroutine */ int F77NAME(slassq)(integer *n, Real *x, integer *incx, Real *scale,
    Real *sumsq);

/* Subroutine */ int F77NAME(slasv2)(Real *f, Real *g, Real *h__, Real *ssmin, Real *
    ssmax, Real *snr, Real *csr, Real *snl, Real *csl);

/* Subroutine */ int F77NAME(slaswp)(integer *n, Real *a, integer *lda, integer *k1,
    integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(slasy2)(logical *ltranl, logical *ltranr, integer *isgn,
    integer *n1, integer *n2, Real *tl, integer *ldtl, Real *tr, integer *
    ldtr, Real *b, integer *ldb, Real *scale, Real *x, integer *ldx, Real
    *xnorm, integer *info);

/* Subroutine */ int F77NAME(slasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     Real *a, integer *lda, integer *ipiv, Real *w, integer *ldw, integer
    *info);

/* Subroutine */ int F77NAME(slatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, Real *ab, integer *ldab, Real *x,
    Real *scale, Real *cnorm, integer *info);

/* Subroutine */ int F77NAME(slatdf)(integer *ijob, integer *n, Real *z__, integer *
    ldz, Real *rhs, Real *rdsum, Real *rdscal, integer *ipiv, integer *
    jpiv);

/* Subroutine */ int F77NAME(slatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, Real *ap, Real *x, Real *scale, Real *cnorm,
    integer *info);

/* Subroutine */ int F77NAME(slatrd)(char *uplo, integer *n, integer *nb, Real *a,
    integer *lda, Real *e, Real *tau, Real *w, integer *ldw);

/* Subroutine */ int F77NAME(slatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, Real *a, integer *lda, Real *x, Real *scale, Real
    *cnorm, integer *info);

/* Subroutine */ int F77NAME(slatrz)(integer *m, integer *n, integer *l, Real *a,
    integer *lda, Real *tau, Real *work);

/* Subroutine */ int F77NAME(slatzm)(char *side, integer *m, integer *n, Real *v,
    integer *incv, Real *tau, Real *c1, Real *c2, integer *ldc, Real *
    work);

/* Subroutine */ int F77NAME(slauu2)(char *uplo, integer *n, Real *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(slauum)(char *uplo, integer *n, Real *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(sopgtr)(char *uplo, integer *n, Real *ap, Real *tau,
    Real *q, integer *ldq, Real *work, integer *info);

/* Subroutine */ int F77NAME(sopmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, Real *ap, Real *tau, Real *c__, integer *ldc, Real *work,
    integer *info);

/* Subroutine */ int F77NAME(sorg2l)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sorg2r)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sorgbr)(char *vect, integer *m, integer *n, integer *k,
    Real *a, integer *lda, Real *tau, Real *work, integer *lwork, integer
    *info);

/* Subroutine */ int F77NAME(sorghr)(integer *n, integer *ilo, integer *ihi, Real *a,
    integer *lda, Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgl2)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sorglq)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgql)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgqr)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgr2)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *info);

/* Subroutine */ int F77NAME(sorgrq)(integer *m, integer *n, integer *k, Real *a,
    integer *lda, Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgtr)(char *uplo, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *info);

/* Subroutine */ int F77NAME(sorm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *info);

/* Subroutine */ int F77NAME(sormbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, Real *a, integer *lda, Real *tau, Real *c__,
    integer *ldc, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, Real *a, integer *lda, Real *tau, Real *
    c__, integer *ldc, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *info);

/* Subroutine */ int F77NAME(sormlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormql)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *info);

/* Subroutine */ int F77NAME(sormr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, Real *a, integer *lda, Real *tau, Real *c__,
    integer *ldc, Real *work, integer *info);

/* Subroutine */ int F77NAME(sormrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, Real *a, integer *lda, Real *tau, Real *c__,
    integer *ldc, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, Real *a, integer *lda, Real *tau, Real *c__, integer *ldc,
     Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(spbcon)(char *uplo, integer *n, integer *kd, Real *ab,
    integer *ldab, Real *anorm, Real *rcond, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(spbequ)(char *uplo, integer *n, integer *kd, Real *ab,
    integer *ldab, Real *s, Real *scond, Real *amax, integer *info);

/* Subroutine */ int F77NAME(spbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Real *ab, integer *ldab, Real *afb, integer *ldafb, Real *b,
    integer *ldb, Real *x, integer *ldx, Real *ferr, Real *berr, Real *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spbstf)(char *uplo, integer *n, integer *kd, Real *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Real *ab, integer *ldab, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(spbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, Real *ab, integer *ldab, Real *afb, integer *ldafb,
    char *equed, Real *s, Real *b, integer *ldb, Real *x, integer *ldx,
    Real *rcond, Real *ferr, Real *berr, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(spbtf2)(char *uplo, integer *n, integer *kd, Real *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbtrf)(char *uplo, integer *n, integer *kd, Real *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Real *ab, integer *ldab, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(spocon)(char *uplo, integer *n, Real *a, integer *lda,
    Real *anorm, Real *rcond, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spoequ)(integer *n, Real *a, integer *lda, Real *s, Real
    *scond, Real *amax, integer *info);

/* Subroutine */ int F77NAME(sporfs)(char *uplo, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *af, integer *ldaf, Real *b, integer *ldb, Real *x,
     integer *ldx, Real *ferr, Real *berr, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(sposv)(char *uplo, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Real *a, integer *lda, Real *af, integer *ldaf, char *equed,
    Real *s, Real *b, integer *ldb, Real *x, integer *ldx, Real *rcond,
    Real *ferr, Real *berr, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spotf2)(char *uplo, integer *n, Real *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotrf)(char *uplo, integer *n, Real *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotri)(char *uplo, integer *n, Real *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotrs)(char *uplo, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sppcon)(char *uplo, integer *n, Real *ap, Real *anorm,
    Real *rcond, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sppequ)(char *uplo, integer *n, Real *ap, Real *s, Real *
    scond, Real *amax, integer *info);

/* Subroutine */ int F77NAME(spprfs)(char *uplo, integer *n, integer *nrhs, Real *ap,
    Real *afp, Real *b, integer *ldb, Real *x, integer *ldx, Real *ferr,
    Real *berr, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sppsv)(char *uplo, integer *n, integer *nrhs, Real *ap,
    Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Real *ap, Real *afp, char *equed, Real *s, Real *b, integer *
    ldb, Real *x, integer *ldx, Real *rcond, Real *ferr, Real *berr, Real
    *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spptrf)(char *uplo, integer *n, Real *ap, integer *info);

/* Subroutine */ int F77NAME(spptri)(char *uplo, integer *n, Real *ap, integer *info);

/* Subroutine */ int F77NAME(spptrs)(char *uplo, integer *n, integer *nrhs, Real *ap,
    Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptcon)(integer *n, Real *d__, Real *e, Real *anorm,
    Real *rcond, Real *work, integer *info);

/* Subroutine */ int F77NAME(spteqr)(char *compz, integer *n, Real *d__, Real *e,
    Real *z__, integer *ldz, Real *work, integer *info);

/* Subroutine */ int F77NAME(sptrfs)(integer *n, integer *nrhs, Real *d__, Real *e,
    Real *df, Real *ef, Real *b, integer *ldb, Real *x, integer *ldx,
    Real *ferr, Real *berr, Real *work, integer *info);

/* Subroutine */ int F77NAME(sptsv)(integer *n, integer *nrhs, Real *d__, Real *e,
    Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptsvx)(char *fact, integer *n, integer *nrhs, Real *d__,
     Real *e, Real *df, Real *ef, Real *b, integer *ldb, Real *x, integer
    *ldx, Real *rcond, Real *ferr, Real *berr, Real *work, integer *info);

/* Subroutine */ int F77NAME(spttrf)(integer *n, Real *d__, Real *e, integer *info);

/* Subroutine */ int F77NAME(spttrs)(integer *n, integer *nrhs, Real *d__, Real *e,
    Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptts2)(integer *n, integer *nrhs, Real *d__, Real *e,
    Real *b, integer *ldb);

/* Subroutine */ int F77NAME(srscl)(integer *n, Real *sa, Real *sx, integer *incx);

/* Subroutine */ int F77NAME(ssbev)(char *jobz, char *uplo, integer *n, integer *kd,
    Real *ab, integer *ldab, Real *w, Real *z__, integer *ldz, Real *work,
     integer *info);

/* Subroutine */ int F77NAME(ssbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    Real *ab, integer *ldab, Real *w, Real *z__, integer *ldz, Real *work,
     integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, Real *ab, integer *ldab, Real *q, integer *ldq, Real *vl,
     Real *vu, integer *il, integer *iu, Real *abstol, integer *m, Real *
    w, Real *z__, integer *ldz, Real *work, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(ssbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, Real *ab, integer *ldab, Real *bb, integer *ldbb, Real *
    x, integer *ldx, Real *work, integer *info);

/* Subroutine */ int F77NAME(ssbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, Real *ab, integer *ldab, Real *bb, integer *ldbb, Real *
    w, Real *z__, integer *ldz, Real *work, integer *info);

/* Subroutine */ int F77NAME(ssbgvd)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, Real *ab, integer *ldab, Real *bb, integer *ldbb, Real *
    w, Real *z__, integer *ldz, Real *work, integer *lwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, Real *ab, integer *ldab, Real *bb, integer *
    ldbb, Real *q, integer *ldq, Real *vl, Real *vu, integer *il, integer
    *iu, Real *abstol, integer *m, Real *w, Real *z__, integer *ldz, Real
    *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    Real *ab, integer *ldab, Real *d__, Real *e, Real *q, integer *ldq,
    Real *work, integer *info);

/* Subroutine */ int F77NAME(sspcon)(char *uplo, integer *n, Real *ap, integer *ipiv,
    Real *anorm, Real *rcond, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sspev)(char *jobz, char *uplo, integer *n, Real *ap,
    Real *w, Real *z__, integer *ldz, Real *work, integer *info);

/* Subroutine */ int F77NAME(sspevd)(char *jobz, char *uplo, integer *n, Real *ap,
    Real *w, Real *z__, integer *ldz, Real *work, integer *lwork, integer
    *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sspevx)(char *jobz, char *range, char *uplo, integer *n,
    Real *ap, Real *vl, Real *vu, integer *il, integer *iu, Real *abstol,
    integer *m, Real *w, Real *z__, integer *ldz, Real *work, integer *
    iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(sspgst)(integer *itype, char *uplo, integer *n, Real *ap,
     Real *bp, integer *info);

/* Subroutine */ int F77NAME(sspgv)(integer *itype, char *jobz, char *uplo, integer *
    n, Real *ap, Real *bp, Real *w, Real *z__, integer *ldz, Real *work,
    integer *info);

/* Subroutine */ int F77NAME(sspgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, Real *ap, Real *bp, Real *w, Real *z__, integer *ldz, Real *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sspgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, Real *ap, Real *bp, Real *vl, Real *vu, integer *il,
     integer *iu, Real *abstol, integer *m, Real *w, Real *z__, integer *
    ldz, Real *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssprfs)(char *uplo, integer *n, integer *nrhs, Real *ap,
    Real *afp, integer *ipiv, Real *b, integer *ldb, Real *x, integer *
    ldx, Real *ferr, Real *berr, Real *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(sspsv)(char *uplo, integer *n, integer *nrhs, Real *ap,
    integer *ipiv, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Real *ap, Real *afp, integer *ipiv, Real *b, integer *ldb, Real
    *x, integer *ldx, Real *rcond, Real *ferr, Real *berr, Real *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ssptrd)(char *uplo, integer *n, Real *ap, Real *d__,
    Real *e, Real *tau, integer *info);

/* Subroutine */ int F77NAME(ssptrf)(char *uplo, integer *n, Real *ap, integer *ipiv,
    integer *info);

/* Subroutine */ int F77NAME(ssptri)(char *uplo, integer *n, Real *ap, integer *ipiv,
    Real *work, integer *info);

/* Subroutine */ int F77NAME(ssptrs)(char *uplo, integer *n, integer *nrhs, Real *ap,
    integer *ipiv, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sstebz)(char *range, char *order, integer *n, Real *vl,
    Real *vu, integer *il, integer *iu, Real *abstol, Real *d__, Real *e,
    integer *m, integer *nsplit, Real *w, integer *iblock, integer *
    isplit, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sstedc)(char *compz, integer *n, Real *d__, Real *e,
    Real *z__, integer *ldz, Real *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstegr)(char *jobz, char *range, integer *n, Real *d__,
    Real *e, Real *vl, Real *vu, integer *il, integer *iu, Real *abstol,
    integer *m, Real *w, Real *z__, integer *ldz, integer *isuppz, Real *
    work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstein)(integer *n, Real *d__, Real *e, integer *m, Real
    *w, integer *iblock, integer *isplit, Real *z__, integer *ldz, Real *
    work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssteqr)(char *compz, integer *n, Real *d__, Real *e,
    Real *z__, integer *ldz, Real *work, integer *info);

/* Subroutine */ int F77NAME(ssterf)(integer *n, Real *d__, Real *e, integer *info);

/* Subroutine */ int F77NAME(sstev)(char *jobz, integer *n, Real *d__, Real *e, Real *
    z__, integer *ldz, Real *work, integer *info);

/* Subroutine */ int F77NAME(sstevd)(char *jobz, integer *n, Real *d__, Real *e, Real
    *z__, integer *ldz, Real *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstevr)(char *jobz, char *range, integer *n, Real *d__,
    Real *e, Real *vl, Real *vu, integer *il, integer *iu, Real *abstol,
    integer *m, Real *w, Real *z__, integer *ldz, integer *isuppz, Real *
    work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstevx)(char *jobz, char *range, integer *n, Real *d__,
    Real *e, Real *vl, Real *vu, integer *il, integer *iu, Real *abstol,
    integer *m, Real *w, Real *z__, integer *ldz, Real *work, integer *
    iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssycon)(char *uplo, integer *n, Real *a, integer *lda,
    integer *ipiv, Real *anorm, Real *rcond, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(ssyev)(char *jobz, char *uplo, integer *n, Real *a,
    integer *lda, Real *w, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssyevd)(char *jobz, char *uplo, integer *n, Real *a,
    integer *lda, Real *w, Real *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssyevr)(char *jobz, char *range, char *uplo, integer *n,
    Real *a, integer *lda, Real *vl, Real *vu, integer *il, integer *iu,
    Real *abstol, integer *m, Real *w, Real *z__, integer *ldz, integer *
    isuppz, Real *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ssyevx)(char *jobz, char *range, char *uplo, integer *n,
    Real *a, integer *lda, Real *vl, Real *vu, integer *il, integer *iu,
    Real *abstol, integer *m, Real *w, Real *z__, integer *ldz, Real *
    work, integer *lwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssygs2)(integer *itype, char *uplo, integer *n, Real *a,
    integer *lda, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ssygst)(integer *itype, char *uplo, integer *n, Real *a,
    integer *lda, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ssygv)(integer *itype, char *jobz, char *uplo, integer *
    n, Real *a, integer *lda, Real *b, integer *ldb, Real *w, Real *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssygvd)(integer *itype, char *jobz, char *uplo, integer *
    n, Real *a, integer *lda, Real *b, integer *ldb, Real *w, Real *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssygvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, Real *a, integer *lda, Real *b, integer *ldb, Real *
    vl, Real *vu, integer *il, integer *iu, Real *abstol, integer *m,
    Real *w, Real *z__, integer *ldz, Real *work, integer *lwork, integer
    *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssyrfs)(char *uplo, integer *n, integer *nrhs, Real *a,
    integer *lda, Real *af, integer *ldaf, integer *ipiv, Real *b,
    integer *ldb, Real *x, integer *ldx, Real *ferr, Real *berr, Real *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ssysv)(char *uplo, integer *n, integer *nrhs, Real *a,
    integer *lda, integer *ipiv, Real *b, integer *ldb, Real *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Real *a, integer *lda, Real *af, integer *ldaf, integer *ipiv,
    Real *b, integer *ldb, Real *x, integer *ldx, Real *rcond, Real *ferr,
     Real *berr, Real *work, integer *lwork, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(ssytd2)(char *uplo, integer *n, Real *a, integer *lda,
    Real *d__, Real *e, Real *tau, integer *info);

/* Subroutine */ int F77NAME(ssytf2)(char *uplo, integer *n, Real *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(ssytrd)(char *uplo, integer *n, Real *a, integer *lda,
    Real *d__, Real *e, Real *tau, Real *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(ssytrf)(char *uplo, integer *n, Real *a, integer *lda,
    integer *ipiv, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssytri)(char *uplo, integer *n, Real *a, integer *lda,
    integer *ipiv, Real *work, integer *info);

/* Subroutine */ int F77NAME(ssytrs)(char *uplo, integer *n, integer *nrhs, Real *a,
    integer *lda, integer *ipiv, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(stbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, Real *ab, integer *ldab, Real *rcond, Real *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, Real *ab, integer *ldab, Real *b, integer
    *ldb, Real *x, integer *ldx, Real *ferr, Real *berr, Real *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, Real *ab, integer *ldab, Real *b, integer
    *ldb, integer *info);

/* Subroutine */ int F77NAME(stgevc)(char *side, char *howmny, logical *select,
    integer *n, Real *a, integer *lda, Real *b, integer *ldb, Real *vl,
    integer *ldvl, Real *vr, integer *ldvr, integer *mm, integer *m, Real
    *work, integer *info);

/* Subroutine */ int F77NAME(stgex2)(logical *wantq, logical *wantz, integer *n, Real
    *a, integer *lda, Real *b, integer *ldb, Real *q, integer *ldq, Real *
    z__, integer *ldz, integer *j1, integer *n1, integer *n2, Real *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(stgexc)(logical *wantq, logical *wantz, integer *n, Real
    *a, integer *lda, Real *b, integer *ldb, Real *q, integer *ldq, Real *
    z__, integer *ldz, integer *ifst, integer *ilst, Real *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(stgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, Real *a, integer *lda, Real *b, integer *
    ldb, Real *alphar, Real *alphai, Real *beta, Real *q, integer *ldq,
    Real *z__, integer *ldz, integer *m, Real *pl, Real *pr, Real *dif,
    Real *work, integer *lwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(stgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, Real *a, integer *lda,
     Real *b, integer *ldb, Real *tola, Real *tolb, Real *alpha, Real *
    beta, Real *u, integer *ldu, Real *v, integer *ldv, Real *q, integer *
    ldq, Real *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(stgsna)(char *job, char *howmny, logical *select,
    integer *n, Real *a, integer *lda, Real *b, integer *ldb, Real *vl,
    integer *ldvl, Real *vr, integer *ldvr, Real *s, Real *dif, integer *
    mm, integer *m, Real *work, integer *lwork, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(stgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, Real *a, integer *lda, Real *b, integer *ldb, Real *c__, integer *
    ldc, Real *d__, integer *ldd, Real *e, integer *lde, Real *f, integer
    *ldf, Real *scale, Real *rdsum, Real *rdscal, integer *iwork, integer
    *pq, integer *info);

/* Subroutine */ int F77NAME(stgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, Real *a, integer *lda, Real *b, integer *ldb, Real *c__, integer *
    ldc, Real *d__, integer *ldd, Real *e, integer *lde, Real *f, integer
    *ldf, Real *scale, Real *dif, Real *work, integer *lwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(stpcon)(char *norm, char *uplo, char *diag, integer *n,
    Real *ap, Real *rcond, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Real *ap, Real *b, integer *ldb, Real *x, integer *ldx,
     Real *ferr, Real *berr, Real *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stptri)(char *uplo, char *diag, integer *n, Real *ap,
    integer *info);

/* Subroutine */ int F77NAME(stptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Real *ap, Real *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(strcon)(char *norm, char *uplo, char *diag, integer *n,
    Real *a, integer *lda, Real *rcond, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(strevc)(char *side, char *howmny, logical *select,
    integer *n, Real *t, integer *ldt, Real *vl, integer *ldvl, Real *vr,
    integer *ldvr, integer *mm, integer *m, Real *work, integer *info);

/* Subroutine */ int F77NAME(strexc)(char *compq, integer *n, Real *t, integer *ldt,
    Real *q, integer *ldq, integer *ifst, integer *ilst, Real *work,
    integer *info);

/* Subroutine */ int F77NAME(strrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Real *a, integer *lda, Real *b, integer *ldb, Real *x,
    integer *ldx, Real *ferr, Real *berr, Real *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(strsen)(char *job, char *compq, logical *select, integer
    *n, Real *t, integer *ldt, Real *q, integer *ldq, Real *wr, Real *wi,
    integer *m, Real *s, Real *sep, Real *work, integer *lwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(strsna)(char *job, char *howmny, logical *select,
    integer *n, Real *t, integer *ldt, Real *vl, integer *ldvl, Real *vr,
    integer *ldvr, Real *s, Real *sep, integer *mm, integer *m, Real *
    work, integer *ldwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(strsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, Real *a, integer *lda, Real *b, integer *ldb, Real *
    c__, integer *ldc, Real *scale, integer *info);

/* Subroutine */ int F77NAME(strti2)(char *uplo, char *diag, integer *n, Real *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(strtri)(char *uplo, char *diag, integer *n, Real *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(strtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Real *a, integer *lda, Real *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(stzrqf)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, integer *info);

/* Subroutine */ int F77NAME(stzrzf)(integer *m, integer *n, Real *a, integer *lda,
    Real *tau, Real *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(xerbla)(char *srname, integer *info);

/* Subroutine */ int F77NAME(zbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, doubleReal *d__, doubleReal *e, doubleComplex *vt,
    integer *ldvt, doubleComplex *u, integer *ldu, doubleComplex *c__,
    integer *ldc, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zdrot)(integer *n, doubleComplex *cx, integer *incx,
    doubleComplex *cy, integer *incy, doubleReal *c__, doubleReal *s);

/* Subroutine */ int F77NAME(zdrscl)(integer *n, doubleReal *sa, doubleComplex *sx,
    integer *incx);

/* Subroutine */ int F77NAME(zgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, doubleComplex *ab, integer *ldab,
    doubleReal *d__, doubleReal *e, doubleComplex *q, integer *ldq,
    doubleComplex *pt, integer *ldpt, doubleComplex *c__, integer *ldc,
    doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, integer *ipiv, doubleReal *anorm,
    doubleReal *rcond, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, doubleReal *r__, doubleReal *c__,
    doubleReal *rowcnd, doubleReal *colcnd, doubleReal *amax, integer *
    info);

/* Subroutine */ int F77NAME(zgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleComplex *ab, integer *ldab, doubleComplex *
    afb, integer *ldafb, integer *ipiv, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, doubleComplex *ab, integer *ldab, integer *ipiv, doubleComplex *
    b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, doubleComplex *ab, integer *ldab,
    doubleComplex *afb, integer *ldafb, integer *ipiv, char *equed,
    doubleReal *r__, doubleReal *c__, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleReal *rcond, doubleReal *ferr,
    doubleReal *berr, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleComplex *ab, integer *ldab, integer *ipiv,
    doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doubleReal *scale, integer *m, doubleComplex *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(zgebal)(char *job, integer *n, doubleComplex *a, integer
    *lda, integer *ilo, integer *ihi, doubleReal *scale, integer *info);

/* Subroutine */ int F77NAME(zgebd2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleReal *d__, doubleReal *e, doubleComplex *tauq,
    doubleComplex *taup, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgebrd)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleReal *d__, doubleReal *e, doubleComplex *tauq,
    doubleComplex *taup, doubleComplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(zgecon)(char *norm, integer *n, doubleComplex *a,
    integer *lda, doubleReal *anorm, doubleReal *rcond, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeequ)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleReal *r__, doubleReal *c__, doubleReal *rowcnd,
    doubleReal *colcnd, doubleReal *amax, integer *info);

/* Subroutine */ int F77NAME(zgees)(char *jobvs, char *sort, L_fp select, integer *n,
    doubleComplex *a, integer *lda, integer *sdim, doubleComplex *w,
    doubleComplex *vs, integer *ldvs, doubleComplex *work, integer *lwork,
     doubleReal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, doubleComplex *a, integer *lda, integer *sdim,
    doubleComplex *w, doubleComplex *vs, integer *ldvs, doubleReal *
    rconde, doubleReal *rcondv, doubleComplex *work, integer *lwork,
    doubleReal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zgeev)(char *jobvl, char *jobvr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *w, doubleComplex *vl,
    integer *ldvl, doubleComplex *vr, integer *ldvr, doubleComplex *work,
    integer *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleComplex *a, integer *lda, doubleComplex *w,
    doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *ldvr,
    integer *ilo, integer *ihi, doubleReal *scale, doubleReal *abnrm,
    doubleReal *rconde, doubleReal *rcondv, doubleComplex *work, integer *
    lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgegs)(char *jobvsl, char *jobvsr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *vsl,
    integer *ldvsl, doubleComplex *vsr, integer *ldvsr, doubleComplex *
    work, integer *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgegv)(char *jobvl, char *jobvr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *vl, integer
    *ldvl, doubleComplex *vr, integer *ldvr, doubleComplex *work, integer
    *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgehd2)(integer *n, integer *ilo, integer *ihi,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zgehrd)(integer *n, integer *ilo, integer *ihi,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgelq2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgelqf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgels)(char *trans, integer *m, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgelsx)(integer *m, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *jpvt, doubleReal *rcond, integer *rank, doubleComplex *work,
    doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgelsy)(integer *m, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *jpvt, doubleReal *rcond, integer *rank, doubleComplex *work,
    integer *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeql2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgeqlf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgeqp3)(integer *m, integer *n, doubleComplex *a,
    integer *lda, integer *jpvt, doubleComplex *tau, doubleComplex *work,
    integer *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeqpf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, integer *jpvt, doubleComplex *tau, doubleComplex *work,
    doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeqr2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgeqrf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgerfs)(char *trans, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
    integer *ldx, doubleReal *ferr, doubleReal *berr, doubleComplex *work,
     doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgerq2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgerqf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgesc2)(integer *n, doubleComplex *a, integer *lda,
    doubleComplex *rhs, integer *ipiv, integer *jpiv, doubleReal *scale);

/* Subroutine */ int F77NAME(zgesv)(integer *n, integer *nrhs, doubleComplex *a,
    integer *lda, integer *ipiv, doubleComplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(zgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, integer *ipiv, char *equed, doubleReal *r__, doubleReal *c__,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *rcond, doubleReal *ferr, doubleReal *berr, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgetc2)(integer *n, doubleComplex *a, integer *lda,
    integer *ipiv, integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(zgetf2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgetrf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgetri)(integer *n, doubleComplex *a, integer *lda,
    integer *ipiv, doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgetrs)(char *trans, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doubleReal *lscale, doubleReal *rscale, integer *m,
    doubleComplex *v, integer *ldv, integer *info);

/* Subroutine */ int F77NAME(zggbal)(char *job, integer *n, doubleComplex *a, integer
    *lda, doubleComplex *b, integer *ldb, integer *ilo, integer *ihi,
    doubleReal *lscale, doubleReal *rscale, doubleReal *work, integer *
    info);

/* Subroutine */ int F77NAME(zgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, integer *sdim, doubleComplex *alpha, doubleComplex *
    beta, doubleComplex *vsl, integer *ldvsl, doubleComplex *vsr, integer
    *ldvsr, doubleComplex *work, integer *lwork, doubleReal *rwork,
    logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, char *sense, integer *n, doubleComplex *a, integer *lda,
    doubleComplex *b, integer *ldb, integer *sdim, doubleComplex *alpha,
    doubleComplex *beta, doubleComplex *vsl, integer *ldvsl,
    doubleComplex *vsr, integer *ldvsr, doubleReal *rconde, doubleReal *
    rcondv, doubleComplex *work, integer *lwork, doubleReal *rwork,
    integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggev)(char *jobvl, char *jobvr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *vl, integer
    *ldvl, doubleComplex *vr, integer *ldvr, doubleComplex *work, integer
    *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleComplex *alpha, doubleComplex *beta,
    doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *ldvr,
    integer *ilo, integer *ihi, doubleReal *lscale, doubleReal *rscale,
    doubleReal *abnrm, doubleReal *bbnrm, doubleReal *rconde, doubleReal *
    rcondv, doubleComplex *work, integer *lwork, doubleReal *rwork,
    integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggglm)(integer *n, integer *m, integer *p,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *d__, doubleComplex *x, doubleComplex *y, doubleComplex
    *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleComplex *q, integer *ldq, doubleComplex *z__,
    integer *ldz, integer *info);

/* Subroutine */ int F77NAME(zgglse)(integer *m, integer *n, integer *p,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *c__, doubleComplex *d__, doubleComplex *x,
    doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zggqrf)(integer *n, integer *m, integer *p,
    doubleComplex *a, integer *lda, doubleComplex *taua, doubleComplex *b,
     integer *ldb, doubleComplex *taub, doubleComplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zggrqf)(integer *m, integer *p, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *taua, doubleComplex *b,
     integer *ldb, doubleComplex *taub, doubleComplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, doubleComplex *a,
    integer *lda, doubleComplex *b, integer *ldb, doubleReal *alpha,
    doubleReal *beta, doubleComplex *u, integer *ldu, doubleComplex *v,
    integer *ldv, doubleComplex *q, integer *ldq, doubleComplex *work,
    doubleReal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, doubleComplex *a, integer *lda, doubleComplex
    *b, integer *ldb, doubleReal *tola, doubleReal *tolb, integer *k,
    integer *l, doubleComplex *u, integer *ldu, doubleComplex *v, integer
    *ldv, doubleComplex *q, integer *ldq, integer *iwork, doubleReal *
    rwork, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgtcon)(char *norm, integer *n, doubleComplex *dl,
    doubleComplex *d__, doubleComplex *du, doubleComplex *du2, integer *
    ipiv, doubleReal *anorm, doubleReal *rcond, doubleComplex *work,
    integer *info);

/* Subroutine */ int F77NAME(zgtrfs)(char *trans, integer *n, integer *nrhs,
    doubleComplex *dl, doubleComplex *d__, doubleComplex *du,
    doubleComplex *dlf, doubleComplex *df, doubleComplex *duf,
    doubleComplex *du2, integer *ipiv, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgtsv)(integer *n, integer *nrhs, doubleComplex *dl,
    doubleComplex *d__, doubleComplex *du, doubleComplex *b, integer *ldb,
     integer *info);

/* Subroutine */ int F77NAME(zgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleComplex *dl, doubleComplex *d__, doubleComplex *du,
    doubleComplex *dlf, doubleComplex *df, doubleComplex *duf,
    doubleComplex *du2, integer *ipiv, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleReal *rcond, doubleReal *ferr,
    doubleReal *berr, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgttrf)(integer *n, doubleComplex *dl, doubleComplex *
    d__, doubleComplex *du, doubleComplex *du2, integer *ipiv, integer *
    info);

/* Subroutine */ int F77NAME(zgttrs)(char *trans, integer *n, integer *nrhs,
    doubleComplex *dl, doubleComplex *d__, doubleComplex *du,
    doubleComplex *du2, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zgtts2)(integer *itrans, integer *n, integer *nrhs,
    doubleComplex *dl, doubleComplex *d__, doubleComplex *du,
    doubleComplex *du2, integer *ipiv, doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zhbev)(char *jobz, char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleReal *w, doubleComplex *z__,
    integer *ldz, doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleReal *w, doubleComplex *z__,
    integer *ldz, doubleComplex *work, integer *lwork, doubleReal *rwork,
    integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, doubleComplex *ab, integer *ldab, doubleComplex *q,
    integer *ldq, doubleReal *vl, doubleReal *vu, integer *il, integer *
    iu, doubleReal *abstol, integer *m, doubleReal *w, doubleComplex *z__,
     integer *ldz, doubleComplex *work, doubleReal *rwork, integer *iwork,
     integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zhbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, doubleComplex *ab, integer *ldab, doubleComplex *bb,
    integer *ldbb, doubleComplex *x, integer *ldx, doubleComplex *work,
    doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doubleComplex *ab, integer *ldab, doubleComplex *bb,
    integer *ldbb, doubleReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, doubleComplex *ab, integer *ldab,
    doubleComplex *bb, integer *ldbb, doubleComplex *q, integer *ldq,
    doubleReal *vl, doubleReal *vu, integer *il, integer *iu, doubleReal *
    abstol, integer *m, doubleReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, doubleReal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(zhbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleReal *d__, doubleReal *e,
    doubleComplex *q, integer *ldq, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zhecon)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleReal *anorm, doubleReal *rcond,
    doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zheev)(char *jobz, char *uplo, integer *n, doubleComplex
    *a, integer *lda, doubleReal *w, doubleComplex *work, integer *lwork,
    doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zheevd)(char *jobz, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleReal *w, doubleComplex *work,
    integer *lwork, doubleReal *rwork, integer *lrwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zheevr)(char *jobz, char *range, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleReal *vl, doubleReal *vu,
    integer *il, integer *iu, doubleReal *abstol, integer *m, doubleReal *
    w, doubleComplex *z__, integer *ldz, integer *isuppz, doubleComplex *
    work, integer *lwork, doubleReal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zheevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleReal *vl, doubleReal *vu,
    integer *il, integer *iu, doubleReal *abstol, integer *m, doubleReal *
    w, doubleComplex *z__, integer *ldz, doubleComplex *work, integer *
    lwork, doubleReal *rwork, integer *iwork, integer *ifail, integer *
    info);

/* Subroutine */ int F77NAME(zhegs2)(integer *itype, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhegst)(integer *itype, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhegv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleReal *w, doubleComplex *work, integer *lwork, doubleReal *rwork,
     integer *info);

/* Subroutine */ int F77NAME(zhegvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleReal *w, doubleComplex *work, integer *lwork, doubleReal *rwork,
     integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhegvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleReal *vl, doubleReal *vu, integer *il, integer *
    iu, doubleReal *abstol, integer *m, doubleReal *w, doubleComplex *z__,
     integer *ldz, doubleComplex *work, integer *lwork, doubleReal *rwork,
     integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zherfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
    integer *ldx, doubleReal *ferr, doubleReal *berr, doubleComplex *work,
     doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhesv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *b,
    integer *ldb, doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zhesvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
     integer *ldx, doubleReal *rcond, doubleReal *ferr, doubleReal *berr,
    doubleComplex *work, integer *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhetf2)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zhetrd)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleReal *d__, doubleReal *e, doubleComplex *tau,
    doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zhetrf)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleComplex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(zhetri)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zhetrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zhgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, doubleComplex *a, integer *lda,
    doubleComplex *b, integer *ldb, doubleComplex *alpha, doubleComplex *
    beta, doubleComplex *q, integer *ldq, doubleComplex *z__, integer *
    ldz, doubleComplex *work, integer *lwork, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpcon)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, doubleReal *anorm, doubleReal *rcond, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zhpev)(char *jobz, char *uplo, integer *n, doubleComplex
    *ap, doubleReal *w, doubleComplex *z__, integer *ldz, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhpevd)(char *jobz, char *uplo, integer *n,
    doubleComplex *ap, doubleReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, integer *lwork, doubleReal *rwork, integer *
    lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhpevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleComplex *ap, doubleReal *vl, doubleReal *vu, integer *il,
    integer *iu, doubleReal *abstol, integer *m, doubleReal *w,
    doubleComplex *z__, integer *ldz, doubleComplex *work, doubleReal *
    rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zhpgst)(integer *itype, char *uplo, integer *n,
    doubleComplex *ap, doubleComplex *bp, integer *info);

/* Subroutine */ int F77NAME(zhpgv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *ap, doubleComplex *bp, doubleReal *w, doubleComplex
    *z__, integer *ldz, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *ap, doubleComplex *bp, doubleReal *w, doubleComplex
    *z__, integer *ldz, doubleComplex *work, integer *lwork, doubleReal *
    rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleComplex *ap, doubleComplex *bp, doubleReal *
    vl, doubleReal *vu, integer *il, integer *iu, doubleReal *abstol,
    integer *m, doubleReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, doubleReal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(zhprfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *afp, integer *ipiv, doubleComplex *
    b, integer *ldb, doubleComplex *x, integer *ldx, doubleReal *ferr,
    doubleReal *berr, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpsv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhpsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *ap, doubleComplex *afp, integer *ipiv,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *rcond, doubleReal *ferr, doubleReal *berr, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhptrd)(char *uplo, integer *n, doubleComplex *ap,
    doubleReal *d__, doubleReal *e, doubleComplex *tau, integer *info);

/* Subroutine */ int F77NAME(zhptrf)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zhptri)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zhptrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, doubleComplex *h__, integer *ldh, doubleComplex *
    w, doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *ldvr,
     integer *mm, integer *m, doubleComplex *work, doubleReal *rwork,
    integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(zhseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, doubleComplex *h__, integer *ldh, doubleComplex *w,
    doubleComplex *z__, integer *ldz, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zlabrd)(integer *m, integer *n, integer *nb,
    doubleComplex *a, integer *lda, doubleReal *d__, doubleReal *e,
    doubleComplex *tauq, doubleComplex *taup, doubleComplex *x, integer *
    ldx, doubleComplex *y, integer *ldy);

/* Subroutine */ int F77NAME(zlacgv)(integer *n, doubleComplex *x, integer *incx);

/* Subroutine */ int F77NAME(zlacon)(integer *n, doubleComplex *v, doubleComplex *x,
    doubleReal *est, integer *kase);

/* Subroutine */ int F77NAME(zlacp2)(char *uplo, integer *m, integer *n, doubleReal *
    a, integer *lda, doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlacpy)(char *uplo, integer *m, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlacrm)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleReal *b, integer *ldb, doubleComplex *c__,
    integer *ldc, doubleReal *rwork);

/* Subroutine */ int F77NAME(zlacrt)(integer *n, doubleComplex *cx, integer *incx,
    doubleComplex *cy, integer *incy, doubleComplex *c__, doubleComplex *
    s);

/* Subroutine */ int F77NAME(zlaed0)(integer *qsiz, integer *n, doubleReal *d__,
    doubleReal *e, doubleComplex *q, integer *ldq, doubleComplex *qstore,
    integer *ldqs, doubleReal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlaed7)(integer *n, integer *cutpnt, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, doubleReal *d__,
    doubleComplex *q, integer *ldq, doubleReal *rho, integer *indxq,
    doubleReal *qstore, integer *qptr, integer *prmptr, integer *perm,
    integer *givptr, integer *givcol, doubleReal *givnum, doubleComplex *
    work, doubleReal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlaed8)(integer *k, integer *n, integer *qsiz,
    doubleComplex *q, integer *ldq, doubleReal *d__, doubleReal *rho,
    integer *cutpnt, doubleReal *z__, doubleReal *dlamda, doubleComplex *
    q2, integer *ldq2, doubleReal *w, integer *indxp, integer *indx,
    integer *indxq, integer *perm, integer *givptr, integer *givcol,
    doubleReal *givnum, integer *info);

/* Subroutine */ int F77NAME(zlaein)(logical *rightv, logical *noinit, integer *n,
    doubleComplex *h__, integer *ldh, doubleComplex *w, doubleComplex *v,
    doubleComplex *b, integer *ldb, doubleReal *rwork, doubleReal *eps3,
    doubleReal *smlnum, integer *info);

/* Subroutine */ int F77NAME(zlaesy)(doubleComplex *a, doubleComplex *b,
    doubleComplex *c__, doubleComplex *rt1, doubleComplex *rt2,
    doubleComplex *evscal, doubleComplex *cs1, doubleComplex *sn1);

/* Subroutine */ int F77NAME(zlaev2)(doubleComplex *a, doubleComplex *b,
    doubleComplex *c__, doubleReal *rt1, doubleReal *rt2, doubleReal *cs1,
     doubleComplex *sn1);

/* Subroutine */ int F77NAME(zlags2)(logical *upper, doubleReal *a1, doubleComplex *
    a2, doubleReal *a3, doubleReal *b1, doubleComplex *b2, doubleReal *b3,
     doubleReal *csu, doubleComplex *snu, doubleReal *csv, doubleComplex *
    snv, doubleReal *csq, doubleComplex *snq);

/* Subroutine */ int F77NAME(zlagtm)(char *trans, integer *n, integer *nrhs,
    doubleReal *alpha, doubleComplex *dl, doubleComplex *d__,
    doubleComplex *du, doubleComplex *x, integer *ldx, doubleReal *beta,
    doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlahef)(char *uplo, integer *n, integer *nb, integer *kb,
     doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *w,
    integer *ldw, integer *info);

/* Subroutine */ int F77NAME(zlahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, doubleComplex *h__, integer *ldh,
    doubleComplex *w, integer *iloz, integer *ihiz, doubleComplex *z__,
    integer *ldz, integer *info);

/* Subroutine */ int F77NAME(zlahrd)(integer *n, integer *k, integer *nb,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *t,
    integer *ldt, doubleComplex *y, integer *ldy);

/* Subroutine */ int F77NAME(zlaic1)(integer *job, integer *j, doubleComplex *x,
    doubleReal *sest, doubleComplex *w, doubleComplex *gamma, doubleReal *
    sestpr, doubleComplex *s, doubleComplex *c__);

/* Subroutine */ int F77NAME(zlals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, doubleComplex *b, integer *ldb,
    doubleComplex *bx, integer *ldbx, integer *perm, integer *givptr,
    integer *givcol, integer *ldgcol, doubleReal *givnum, integer *ldgnum,
     doubleReal *poles, doubleReal *difl, doubleReal *difr, doubleReal *
    z__, integer *k, doubleReal *c__, doubleReal *s, doubleReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(zlalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, doubleComplex *b, integer *ldb, doubleComplex *bx,
    integer *ldbx, doubleReal *u, integer *ldu, doubleReal *vt, integer *
    k, doubleReal *difl, doubleReal *difr, doubleReal *z__, doubleReal *
    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
    perm, doubleReal *givnum, doubleReal *c__, doubleReal *s, doubleReal *
    rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlapll)(integer *n, doubleComplex *x, integer *incx,
    doubleComplex *y, integer *incy, doubleReal *ssmin);

/* Subroutine */ int F77NAME(zlapmt)(logical *forwrd, integer *m, integer *n,
    doubleComplex *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(zlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, doubleReal *r__, doubleReal *c__,
    doubleReal *rowcnd, doubleReal *colcnd, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqge)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleReal *r__, doubleReal *c__, doubleReal *rowcnd,
    doubleReal *colcnd, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqhb)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleReal *s, doubleReal *scond,
    doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqhe)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleReal *s, doubleReal *scond, doubleReal *amax,
    char *equed);

/* Subroutine */ int F77NAME(zlaqhp)(char *uplo, integer *n, doubleComplex *ap,
    doubleReal *s, doubleReal *scond, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqp2)(integer *m, integer *n, integer *offset,
    doubleComplex *a, integer *lda, integer *jpvt, doubleComplex *tau,
    doubleReal *vn1, doubleReal *vn2, doubleComplex *work);

/* Subroutine */ int F77NAME(zlaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, doubleComplex *a, integer *lda, integer *jpvt,
    doubleComplex *tau, doubleReal *vn1, doubleReal *vn2, doubleComplex *
    auxv, doubleComplex *f, integer *ldf);

/* Subroutine */ int F77NAME(zlaqsb)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleReal *s, doubleReal *scond,
    doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqsp)(char *uplo, integer *n, doubleComplex *ap,
    doubleReal *s, doubleReal *scond, doubleReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqsy)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleReal *s, doubleReal *scond, doubleReal *amax,
    char *equed);

/* Subroutine */ int F77NAME(zlar1v)(integer *n, integer *b1, integer *bn, doubleReal
    *sigma, doubleReal *d__, doubleReal *l, doubleReal *ld, doubleReal *
    lld, doubleReal *gersch, doubleComplex *z__, doubleReal *ztz,
    doubleReal *mingma, integer *r__, integer *isuppz, doubleReal *work);

/* Subroutine */ int F77NAME(zlar2v)(integer *n, doubleComplex *x, doubleComplex *y,
    doubleComplex *z__, integer *incx, doubleReal *c__, doubleComplex *s,
    integer *incc);

/* Subroutine */ int F77NAME(zlarcm)(integer *m, integer *n, doubleReal *a, integer *
    lda, doubleComplex *b, integer *ldb, doubleComplex *c__, integer *ldc,
     doubleReal *rwork);

/* Subroutine */ int F77NAME(zlarf)(char *side, integer *m, integer *n, doubleComplex
    *v, integer *incv, doubleComplex *tau, doubleComplex *c__, integer *
    ldc, doubleComplex *work);

/* Subroutine */ int F77NAME(zlarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, doubleComplex *v, integer
    *ldv, doubleComplex *t, integer *ldt, doubleComplex *c__, integer *
    ldc, doubleComplex *work, integer *ldwork);

/* Subroutine */ int F77NAME(zlarfg)(integer *n, doubleComplex *alpha, doubleComplex *
    x, integer *incx, doubleComplex *tau);

/* Subroutine */ int F77NAME(zlarft)(char *direct, char *storev, integer *n, integer *
    k, doubleComplex *v, integer *ldv, doubleComplex *tau, doubleComplex *
    t, integer *ldt);

/* Subroutine */ int F77NAME(zlarfx)(char *side, integer *m, integer *n,
    doubleComplex *v, doubleComplex *tau, doubleComplex *c__, integer *
    ldc, doubleComplex *work);

/* Subroutine */ int F77NAME(zlargv)(integer *n, doubleComplex *x, integer *incx,
    doubleComplex *y, integer *incy, doubleReal *c__, integer *incc);

/* Subroutine */ int F77NAME(zlarnv)(integer *idist, integer *iseed, integer *n,
    doubleComplex *x);

/* Subroutine */ int F77NAME(zlarrv)(integer *n, doubleReal *d__, doubleReal *l,
    integer *isplit, integer *m, doubleReal *w, integer *iblock,
    doubleReal *gersch, doubleReal *tol, doubleComplex *z__, integer *ldz,
     integer *isuppz, doubleReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlartg)(doubleComplex *f, doubleComplex *g, doubleReal *
    cs, doubleComplex *sn, doubleComplex *r__);

/* Subroutine */ int F77NAME(zlartv)(integer *n, doubleComplex *x, integer *incx,
    doubleComplex *y, integer *incy, doubleReal *c__, doubleComplex *s,
    integer *incc);

/* Subroutine */ int F77NAME(zlarz)(char *side, integer *m, integer *n, integer *l,
    doubleComplex *v, integer *incv, doubleComplex *tau, doubleComplex *
    c__, integer *ldc, doubleComplex *work);

/* Subroutine */ int F77NAME(zlarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, doubleComplex
    *v, integer *ldv, doubleComplex *t, integer *ldt, doubleComplex *c__,
    integer *ldc, doubleComplex *work, integer *ldwork);

/* Subroutine */ int F77NAME(zlarzt)(char *direct, char *storev, integer *n, integer *
    k, doubleComplex *v, integer *ldv, doubleComplex *tau, doubleComplex *
    t, integer *ldt);

/* Subroutine */ int F77NAME(zlascl)(char *type__, integer *kl, integer *ku,
    doubleReal *cfrom, doubleReal *cto, integer *m, integer *n,
    doubleComplex *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(zlaset)(char *uplo, integer *m, integer *n,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *a, integer *
    lda);

/* Subroutine */ int F77NAME(zlasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, doubleReal *c__, doubleReal *s, doubleComplex *a,
    integer *lda);

/* Subroutine */ int F77NAME(zlassq)(integer *n, doubleComplex *x, integer *incx,
    doubleReal *scale, doubleReal *sumsq);

/* Subroutine */ int F77NAME(zlaswp)(integer *n, doubleComplex *a, integer *lda,
    integer *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(zlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *w,
    integer *ldw, integer *info);

/* Subroutine */ int F77NAME(zlatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, doubleComplex *ab, integer *ldab,
    doubleComplex *x, doubleReal *scale, doubleReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatdf)(integer *ijob, integer *n, doubleComplex *z__,
    integer *ldz, doubleComplex *rhs, doubleReal *rdsum, doubleReal *
    rdscal, integer *ipiv, integer *jpiv);

/* Subroutine */ int F77NAME(zlatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleComplex *ap, doubleComplex *x, doubleReal *
    scale, doubleReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatrd)(char *uplo, integer *n, integer *nb,
    doubleComplex *a, integer *lda, doubleReal *e, doubleComplex *tau,
    doubleComplex *w, integer *ldw);

/* Subroutine */ int F77NAME(zlatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleComplex *a, integer *lda, doubleComplex *x,
    doubleReal *scale, doubleReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatrz)(integer *m, integer *n, integer *l,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work);

/* Subroutine */ int F77NAME(zlatzm)(char *side, integer *m, integer *n,
    doubleComplex *v, integer *incv, doubleComplex *tau, doubleComplex *
    c1, doubleComplex *c2, integer *ldc, doubleComplex *work);

/* Subroutine */ int F77NAME(zlauu2)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zlauum)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpbcon)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleReal *anorm, doubleReal *
    rcond, doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpbequ)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleReal *s, doubleReal *scond,
    doubleReal *amax, integer *info);

/* Subroutine */ int F77NAME(zpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleComplex *ab, integer *ldab, doubleComplex *afb, integer *
    ldafb, doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
     doubleReal *ferr, doubleReal *berr, doubleComplex *work, doubleReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zpbstf)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleComplex *ab, integer *ldab, doubleComplex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(zpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, doubleComplex *ab, integer *ldab, doubleComplex *afb,
    integer *ldafb, char *equed, doubleReal *s, doubleComplex *b, integer
    *ldb, doubleComplex *x, integer *ldx, doubleReal *rcond, doubleReal *
    ferr, doubleReal *berr, doubleComplex *work, doubleReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(zpbtf2)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbtrf)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleComplex *ab, integer *ldab, doubleComplex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(zpocon)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleReal *anorm, doubleReal *rcond, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpoequ)(integer *n, doubleComplex *a, integer *lda,
    doubleReal *s, doubleReal *scond, doubleReal *amax, integer *info);

/* Subroutine */ int F77NAME(zporfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *ferr, doubleReal *berr, doubleComplex *work, doubleReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zposv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, char *equed, doubleReal *s, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleReal *rcond, doubleReal *ferr,
    doubleReal *berr, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zpotf2)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpotrf)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpotri)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpotrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zppcon)(char *uplo, integer *n, doubleComplex *ap,
    doubleReal *anorm, doubleReal *rcond, doubleComplex *work, doubleReal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(zppequ)(char *uplo, integer *n, doubleComplex *ap,
    doubleReal *s, doubleReal *scond, doubleReal *amax, integer *info);

/* Subroutine */ int F77NAME(zpprfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *afp, doubleComplex *b, integer *ldb,
     doubleComplex *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zppsv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *ap, doubleComplex *afp, char *equed, doubleReal *
    s, doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *rcond, doubleReal *ferr, doubleReal *berr, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpptrf)(char *uplo, integer *n, doubleComplex *ap,
    integer *info);

/* Subroutine */ int F77NAME(zpptri)(char *uplo, integer *n, doubleComplex *ap,
    integer *info);

/* Subroutine */ int F77NAME(zpptrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zptcon)(integer *n, doubleReal *d__, doubleComplex *e,
    doubleReal *anorm, doubleReal *rcond, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zptrfs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *d__, doubleComplex *e, doubleReal *df, doubleComplex *ef,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *ferr, doubleReal *berr, doubleComplex *work, doubleReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zptsv)(integer *n, integer *nrhs, doubleReal *d__,
    doubleComplex *e, doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zptsvx)(char *fact, integer *n, integer *nrhs,
    doubleReal *d__, doubleComplex *e, doubleReal *df, doubleComplex *ef,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *rcond, doubleReal *ferr, doubleReal *berr, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpttrf)(integer *n, doubleReal *d__, doubleComplex *e,
    integer *info);

/* Subroutine */ int F77NAME(zpttrs)(char *uplo, integer *n, integer *nrhs,
    doubleReal *d__, doubleComplex *e, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zptts2)(integer *iuplo, integer *n, integer *nrhs,
    doubleReal *d__, doubleComplex *e, doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zrot)(integer *n, doubleComplex *cx, integer *incx,
    doubleComplex *cy, integer *incy, doubleReal *c__, doubleComplex *s);

/* Subroutine */ int F77NAME(zspcon)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, doubleReal *anorm, doubleReal *rcond, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zspmv)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *ap, doubleComplex *x, integer *incx, doubleComplex *
    beta, doubleComplex *y, integer *incy);

/* Subroutine */ int F77NAME(zspr)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *x, integer *incx, doubleComplex *ap);

/* Subroutine */ int F77NAME(zsprfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *afp, integer *ipiv, doubleComplex *
    b, integer *ldb, doubleComplex *x, integer *ldx, doubleReal *ferr,
    doubleReal *berr, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zspsv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *ap, doubleComplex *afp, integer *ipiv,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *rcond, doubleReal *ferr, doubleReal *berr, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsptrf)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zsptri)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zsptrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zstedc)(char *compz, integer *n, doubleReal *d__,
    doubleReal *e, doubleComplex *z__, integer *ldz, doubleComplex *work,
    integer *lwork, doubleReal *rwork, integer *lrwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zstein)(integer *n, doubleReal *d__, doubleReal *e,
    integer *m, doubleReal *w, integer *iblock, integer *isplit,
    doubleComplex *z__, integer *ldz, doubleReal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zsteqr)(char *compz, integer *n, doubleReal *d__,
    doubleReal *e, doubleComplex *z__, integer *ldz, doubleReal *work,
    integer *info);

/* Subroutine */ int F77NAME(zsycon)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleReal *anorm, doubleReal *rcond,
    doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zsymv)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *a, integer *lda, doubleComplex *x, integer *incx,
    doubleComplex *beta, doubleComplex *y, integer *incy);

/* Subroutine */ int F77NAME(zsyr)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *x, integer *incx, doubleComplex *a, integer *lda);

/* Subroutine */ int F77NAME(zsyrfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
    integer *ldx, doubleReal *ferr, doubleReal *berr, doubleComplex *work,
     doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsysv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *b,
    integer *ldb, doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zsysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
     integer *ldx, doubleReal *rcond, doubleReal *ferr, doubleReal *berr,
    doubleComplex *work, integer *lwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsytf2)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zsytrf)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleComplex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(zsytri)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zsytrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ztbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, doubleComplex *ab, integer *ldab, doubleReal *rcond,
    doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleComplex *ab, integer *ldab,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleReal *ferr, doubleReal *berr, doubleComplex *work, doubleReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(ztbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleComplex *ab, integer *ldab,
    doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ztgevc)(char *side, char *howmny, logical *select,
    integer *n, doubleComplex *a, integer *lda, doubleComplex *b, integer
    *ldb, doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *
    ldvr, integer *mm, integer *m, doubleComplex *work, doubleReal *rwork,
     integer *info);

/* Subroutine */ int F77NAME(ztgex2)(logical *wantq, logical *wantz, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *q, integer *ldq, doubleComplex *z__, integer *ldz,
    integer *j1, integer *info);

/* Subroutine */ int F77NAME(ztgexc)(logical *wantq, logical *wantz, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *q, integer *ldq, doubleComplex *z__, integer *ldz,
    integer *ifst, integer *ilst, integer *info);

/* Subroutine */ int F77NAME(ztgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, doubleComplex *a, integer *lda,
    doubleComplex *b, integer *ldb, doubleComplex *alpha, doubleComplex *
    beta, doubleComplex *q, integer *ldq, doubleComplex *z__, integer *
    ldz, integer *m, doubleReal *pl, doubleReal *pr, doubleReal *dif,
    doubleComplex *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ztgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, doubleComplex *a,
    integer *lda, doubleComplex *b, integer *ldb, doubleReal *tola,
    doubleReal *tolb, doubleReal *alpha, doubleReal *beta, doubleComplex *
    u, integer *ldu, doubleComplex *v, integer *ldv, doubleComplex *q,
    integer *ldq, doubleComplex *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(ztgsna)(char *job, char *howmny, logical *select,
    integer *n, doubleComplex *a, integer *lda, doubleComplex *b, integer
    *ldb, doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *
    ldvr, doubleReal *s, doubleReal *dif, integer *mm, integer *m,
    doubleComplex *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ztgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *c__, integer *ldc, doubleComplex *d__, integer *ldd,
    doubleComplex *e, integer *lde, doubleComplex *f, integer *ldf,
    doubleReal *scale, doubleReal *rdsum, doubleReal *rdscal, integer *
    info);

/* Subroutine */ int F77NAME(ztgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *c__, integer *ldc, doubleComplex *d__, integer *ldd,
    doubleComplex *e, integer *lde, doubleComplex *f, integer *ldf,
    doubleReal *scale, doubleReal *dif, doubleComplex *work, integer *
    lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ztpcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleComplex *ap, doubleReal *rcond, doubleComplex *work, doubleReal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(ztprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleComplex *ap, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleReal *ferr, doubleReal *berr,
    doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztptri)(char *uplo, char *diag, integer *n,
    doubleComplex *ap, integer *info);

/* Subroutine */ int F77NAME(ztptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleComplex *ap, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(ztrcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleComplex *a, integer *lda, doubleReal *rcond, doubleComplex *
    work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrevc)(char *side, char *howmny, logical *select,
    integer *n, doubleComplex *t, integer *ldt, doubleComplex *vl,
    integer *ldvl, doubleComplex *vr, integer *ldvr, integer *mm, integer
    *m, doubleComplex *work, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrexc)(char *compq, integer *n, doubleComplex *t,
    integer *ldt, doubleComplex *q, integer *ldq, integer *ifst, integer *
    ilst, integer *info);

/* Subroutine */ int F77NAME(ztrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleComplex *x, integer *ldx, doubleReal *ferr,
    doubleReal *berr, doubleComplex *work, doubleReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(ztrsen)(char *job, char *compq, logical *select, integer
    *n, doubleComplex *t, integer *ldt, doubleComplex *q, integer *ldq,
    doubleComplex *w, integer *m, doubleReal *s, doubleReal *sep,
    doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ztrsna)(char *job, char *howmny, logical *select,
    integer *n, doubleComplex *t, integer *ldt, doubleComplex *vl,
    integer *ldvl, doubleComplex *vr, integer *ldvr, doubleReal *s,
    doubleReal *sep, integer *mm, integer *m, doubleComplex *work,
    integer *ldwork, doubleReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleComplex *c__, integer *ldc, doubleReal *scale,
    integer *info);

/* Subroutine */ int F77NAME(ztrti2)(char *uplo, char *diag, integer *n,
    doubleComplex *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(ztrtri)(char *uplo, char *diag, integer *n,
    doubleComplex *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(ztrtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ztzrqf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, integer *info);

/* Subroutine */ int F77NAME(ztzrzf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zung2l)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zung2r)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zungbr)(char *vect, integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zunghr)(integer *n, integer *ilo, integer *ihi,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungl2)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zunglq)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungql)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungqr)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungr2)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zungrq)(integer *m, integer *n, integer *k,
    doubleComplex *a, integer *lda, doubleComplex *tau, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungtr)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zunm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zunmbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, doubleComplex *a, integer *lda, doubleComplex
    *tau, doubleComplex *c__, integer *ldc, doubleComplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zunmhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, doubleComplex *a, integer *lda,
    doubleComplex *tau, doubleComplex *c__, integer *ldc, doubleComplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zunml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zunmlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmql)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zunmr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doubleComplex *a, integer *lda, doubleComplex
    *tau, doubleComplex *c__, integer *ldc, doubleComplex *work, integer *
    info);

/* Subroutine */ int F77NAME(zunmrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doubleComplex *a, integer *lda, doubleComplex
    *tau, doubleComplex *c__, integer *ldc, doubleComplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zunmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doubleComplex *a, integer *lda, doubleComplex *tau,
    doubleComplex *c__, integer *ldc, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zupgtr)(char *uplo, integer *n, doubleComplex *ap,
    doubleComplex *tau, doubleComplex *q, integer *ldq, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zupmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doubleComplex *ap, doubleComplex *tau, doubleComplex *c__,
     integer *ldc, doubleComplex *work, integer *info);

} // extern "C"

#endif /* __CLAPACK_H */
