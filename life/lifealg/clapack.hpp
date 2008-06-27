#ifndef __CLAPACK_H
#define __CLAPACK_H

#include <life/lifealg/cblas.hpp>

typedef int integer;
typedef unsigned long uinteger;
typedef char *address;
typedef short int shortint;
typedef float CReal;
typedef double doubleCReal;
typedef struct { CReal r, i; } Complex;
typedef struct { doubleCReal r, i; } doubleComplex;
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
typedef CReal (*R_fp)(...);
typedef doubleCReal (*D_fp)(...), (*E_fp)(...);
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
typedef CReal (*R_fp)();
typedef doubleCReal (*D_fp)(), (*E_fp)();
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
    nru, integer *ncc, CReal *d__, CReal *e, Complex *vt, integer *ldvt,
    Complex *u, integer *ldu, Complex *c__, integer *ldc, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, Complex *ab, integer *ldab, CReal *d__,
    CReal *e, Complex *q, integer *ldq, Complex *pt, integer *ldpt,
    Complex *c__, integer *ldc, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, integer *ipiv, CReal *anorm, CReal *rcond,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, CReal *r__, CReal *c__, CReal *rowcnd, CReal
    *colcnd, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(cgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, Complex *ab, integer *ldab, Complex *afb, integer *
    ldafb, integer *ipiv, Complex *b, integer *ldb, Complex *x, integer *
    ldx, CReal *ferr, CReal *berr, Complex *work, CReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, Complex *ab, integer *ldab, integer *ipiv, Complex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(cgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, Complex *ab, integer *ldab, Complex *afb,
     integer *ldafb, integer *ipiv, char *equed, CReal *r__, CReal *c__,
    Complex *b, integer *ldb, Complex *x, integer *ldx, CReal *rcond, CReal
    *ferr, CReal *berr, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, Complex *ab, integer *ldab, integer *ipiv, Complex
    *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, CReal *scale, integer *m, Complex *v, integer *ldv,
    integer *info);

/* Subroutine */ int F77NAME(cgebal)(char *job, integer *n, Complex *a, integer *lda,
    integer *ilo, integer *ihi, CReal *scale, integer *info);

/* Subroutine */ int F77NAME(cgebd2)(integer *m, integer *n, Complex *a, integer *lda,
     CReal *d__, CReal *e, Complex *tauq, Complex *taup, Complex *work,
    integer *info);

/* Subroutine */ int F77NAME(cgebrd)(integer *m, integer *n, Complex *a, integer *lda,
     CReal *d__, CReal *e, Complex *tauq, Complex *taup, Complex *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgecon)(char *norm, integer *n, Complex *a, integer *lda,
     CReal *anorm, CReal *rcond, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgeequ)(integer *m, integer *n, Complex *a, integer *lda,
     CReal *r__, CReal *c__, CReal *rowcnd, CReal *colcnd, CReal *amax,
    integer *info);

/* Subroutine */ int F77NAME(cgees)(char *jobvs, char *sort, L_fp select, integer *n,
    Complex *a, integer *lda, integer *sdim, Complex *w, Complex *vs,
    integer *ldvs, Complex *work, integer *lwork, CReal *rwork, logical *
    bwork, integer *info);

/* Subroutine */ int F77NAME(cgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, Complex *a, integer *lda, integer *sdim, Complex *
    w, Complex *vs, integer *ldvs, CReal *rconde, CReal *rcondv, Complex *
    work, integer *lwork, CReal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cgeev)(char *jobvl, char *jobvr, integer *n, Complex *a,
    integer *lda, Complex *w, Complex *vl, integer *ldvl, Complex *vr,
    integer *ldvr, Complex *work, integer *lwork, CReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, Complex *a, integer *lda, Complex *w, Complex *vl,
    integer *ldvl, Complex *vr, integer *ldvr, integer *ilo, integer *ihi,
     CReal *scale, CReal *abnrm, CReal *rconde, CReal *rcondv, Complex *work,
    integer *lwork, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgegs)(char *jobvsl, char *jobvsr, integer *n, Complex *
    a, integer *lda, Complex *b, integer *ldb, Complex *alpha, Complex *
    beta, Complex *vsl, integer *ldvsl, Complex *vsr, integer *ldvsr,
    Complex *work, integer *lwork, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgegv)(char *jobvl, char *jobvr, integer *n, Complex *a,
    integer *lda, Complex *b, integer *ldb, Complex *alpha, Complex *beta,
     Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, Complex *
    work, integer *lwork, CReal *rwork, integer *info);

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
    a, integer *lda, Complex *b, integer *ldb, integer *jpvt, CReal *rcond,
     integer *rank, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgelsy)(integer *m, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *jpvt, CReal *rcond,
     integer *rank, Complex *work, integer *lwork, CReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeql2)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgeqlf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgeqp3)(integer *m, integer *n, Complex *a, integer *lda,
     integer *jpvt, Complex *tau, Complex *work, integer *lwork, CReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(cgeqpf)(integer *m, integer *n, Complex *a, integer *lda,
     integer *jpvt, Complex *tau, Complex *work, CReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeqr2)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgeqrf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgerfs)(char *trans, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, integer *ipiv, Complex *
    b, integer *ldb, Complex *x, integer *ldx, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgerq2)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgerqf)(integer *m, integer *n, Complex *a, integer *lda,
     Complex *tau, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgesc2)(integer *n, Complex *a, integer *lda, Complex *
    rhs, integer *ipiv, integer *jpiv, CReal *scale);

/* Subroutine */ int F77NAME(cgesv)(integer *n, integer *nrhs, Complex *a, integer *
    lda, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, integer *
    ipiv, char *equed, CReal *r__, CReal *c__, Complex *b, integer *ldb,
    Complex *x, integer *ldx, CReal *rcond, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

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
    integer *ihi, CReal *lscale, CReal *rscale, integer *m, Complex *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(cggbal)(char *job, integer *n, Complex *a, integer *lda,
    Complex *b, integer *ldb, integer *ilo, integer *ihi, CReal *lscale,
    CReal *rscale, CReal *work, integer *info);

/* Subroutine */ int F77NAME(cgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, integer *n, Complex *a, integer *lda, Complex *b, integer *
    ldb, integer *sdim, Complex *alpha, Complex *beta, Complex *vsl,
    integer *ldvsl, Complex *vsr, integer *ldvsr, Complex *work, integer *
    lwork, CReal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, char *sense, integer *n, Complex *a, integer *lda, Complex *b,
     integer *ldb, integer *sdim, Complex *alpha, Complex *beta, Complex *
    vsl, integer *ldvsl, Complex *vsr, integer *ldvsr, CReal *rconde, CReal
    *rcondv, Complex *work, integer *lwork, CReal *rwork, integer *iwork,
    integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cggev)(char *jobvl, char *jobvr, integer *n, Complex *a,
    integer *lda, Complex *b, integer *ldb, Complex *alpha, Complex *beta,
     Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, Complex *
    work, integer *lwork, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
     Complex *alpha, Complex *beta, Complex *vl, integer *ldvl, Complex *
    vr, integer *ldvr, integer *ilo, integer *ihi, CReal *lscale, CReal *
    rscale, CReal *abnrm, CReal *bbnrm, CReal *rconde, CReal *rcondv, Complex
    *work, integer *lwork, CReal *rwork, integer *iwork, logical *bwork,
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
    lda, Complex *b, integer *ldb, CReal *alpha, CReal *beta, Complex *u,
    integer *ldu, Complex *v, integer *ldv, Complex *q, integer *ldq,
    Complex *work, CReal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(cggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, Complex *a, integer *lda, Complex *b, integer
    *ldb, CReal *tola, CReal *tolb, integer *k, integer *l, Complex *u,
    integer *ldu, Complex *v, integer *ldv, Complex *q, integer *ldq,
    integer *iwork, CReal *rwork, Complex *tau, Complex *work, integer *
    info);

/* Subroutine */ int F77NAME(cgtcon)(char *norm, integer *n, Complex *dl, Complex *
    d__, Complex *du, Complex *du2, integer *ipiv, CReal *anorm, CReal *
    rcond, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cgtrfs)(char *trans, integer *n, integer *nrhs, Complex *
    dl, Complex *d__, Complex *du, Complex *dlf, Complex *df, Complex *
    duf, Complex *du2, integer *ipiv, Complex *b, integer *ldb, Complex *
    x, integer *ldx, CReal *ferr, CReal *berr, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cgtsv)(integer *n, integer *nrhs, Complex *dl, Complex *
    d__, Complex *du, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, Complex *dl, Complex *d__, Complex *du, Complex *dlf, Complex *
    df, Complex *duf, Complex *du2, integer *ipiv, Complex *b, integer *
    ldb, Complex *x, integer *ldx, CReal *rcond, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgttrf)(integer *n, Complex *dl, Complex *d__, Complex *
    du, Complex *du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgttrs)(char *trans, integer *n, integer *nrhs, Complex *
    dl, Complex *d__, Complex *du, Complex *du2, integer *ipiv, Complex *
    b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgtts2)(integer *itrans, integer *n, integer *nrhs,
    Complex *dl, Complex *d__, Complex *du, Complex *du2, integer *ipiv,
    Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(chbev)(char *jobz, char *uplo, integer *n, integer *kd,
    Complex *ab, integer *ldab, CReal *w, Complex *z__, integer *ldz,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(chbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    Complex *ab, integer *ldab, CReal *w, Complex *z__, integer *ldz,
    Complex *work, integer *lwork, CReal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, Complex *ab, integer *ldab, Complex *q, integer *ldq,
    CReal *vl, CReal *vu, integer *il, integer *iu, CReal *abstol, integer *
    m, CReal *w, Complex *z__, integer *ldz, Complex *work, CReal *rwork,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, Complex *ab, integer *ldab, Complex *bb, integer *ldbb,
    Complex *x, integer *ldx, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(chbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, Complex *ab, integer *ldab, Complex *bb, integer *ldbb,
    CReal *w, Complex *z__, integer *ldz, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, Complex *ab, integer *ldab, Complex *bb,
    integer *ldbb, Complex *q, integer *ldq, CReal *vl, CReal *vu, integer *
    il, integer *iu, CReal *abstol, integer *m, CReal *w, Complex *z__,
    integer *ldz, Complex *work, CReal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(chbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    Complex *ab, integer *ldab, CReal *d__, CReal *e, Complex *q, integer *
    ldq, Complex *work, integer *info);

/* Subroutine */ int F77NAME(checon)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, CReal *anorm, CReal *rcond, Complex *work, integer *
    info);

/* Subroutine */ int F77NAME(cheev)(char *jobz, char *uplo, integer *n, Complex *a,
    integer *lda, CReal *w, Complex *work, integer *lwork, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cheevd)(char *jobz, char *uplo, integer *n, Complex *a,
    integer *lda, CReal *w, Complex *work, integer *lwork, CReal *rwork,
    integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(cheevr)(char *jobz, char *range, char *uplo, integer *n,
    Complex *a, integer *lda, CReal *vl, CReal *vu, integer *il, integer *
    iu, CReal *abstol, integer *m, CReal *w, Complex *z__, integer *ldz,
    integer *isuppz, Complex *work, integer *lwork, CReal *rwork, integer *
    lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(cheevx)(char *jobz, char *range, char *uplo, integer *n,
    Complex *a, integer *lda, CReal *vl, CReal *vu, integer *il, integer *
    iu, CReal *abstol, integer *m, CReal *w, Complex *z__, integer *ldz,
    Complex *work, integer *lwork, CReal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(chegs2)(integer *itype, char *uplo, integer *n, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chegst)(integer *itype, char *uplo, integer *n, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chegv)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, CReal *w,
    Complex *work, integer *lwork, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(chegvd)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, CReal *w,
    Complex *work, integer *lwork, CReal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chegvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    CReal *vl, CReal *vu, integer *il, integer *iu, CReal *abstol, integer *
    m, CReal *w, Complex *z__, integer *ldz, Complex *work, integer *lwork,
     CReal *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(cherfs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, integer *ipiv, Complex *
    b, integer *ldb, Complex *x, integer *ldx, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(chesv)(char *uplo, integer *n, integer *nrhs, Complex *a,
     integer *lda, integer *ipiv, Complex *b, integer *ldb, Complex *work,
     integer *lwork, integer *info);

/* Subroutine */ int F77NAME(chesvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, integer *
    ipiv, Complex *b, integer *ldb, Complex *x, integer *ldx, CReal *rcond,
     CReal *ferr, CReal *berr, Complex *work, integer *lwork, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chetf2)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(chetrd)(char *uplo, integer *n, Complex *a, integer *lda,
     CReal *d__, CReal *e, Complex *tau, Complex *work, integer *lwork,
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
     Complex *z__, integer *ldz, Complex *work, integer *lwork, CReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(chpcon)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, CReal *anorm, CReal *rcond, Complex *work, integer *info);

/* Subroutine */ int F77NAME(chpev)(char *jobz, char *uplo, integer *n, Complex *ap,
    CReal *w, Complex *z__, integer *ldz, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chpevd)(char *jobz, char *uplo, integer *n, Complex *ap,
    CReal *w, Complex *z__, integer *ldz, Complex *work, integer *lwork,
    CReal *rwork, integer *lrwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(chpevx)(char *jobz, char *range, char *uplo, integer *n,
    Complex *ap, CReal *vl, CReal *vu, integer *il, integer *iu, CReal *
    abstol, integer *m, CReal *w, Complex *z__, integer *ldz, Complex *
    work, CReal *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chpgst)(integer *itype, char *uplo, integer *n, Complex *
    ap, Complex *bp, integer *info);

/* Subroutine */ int F77NAME(chpgv)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *ap, Complex *bp, CReal *w, Complex *z__, integer *ldz,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(chpgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, Complex *ap, Complex *bp, CReal *w, Complex *z__, integer *ldz,
    Complex *work, integer *lwork, CReal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chpgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, Complex *ap, Complex *bp, CReal *vl, CReal *vu,
    integer *il, integer *iu, CReal *abstol, integer *m, CReal *w, Complex *
    z__, integer *ldz, Complex *work, CReal *rwork, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chprfs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *afp, integer *ipiv, Complex *b, integer *ldb, Complex *x,
     integer *ldx, CReal *ferr, CReal *berr, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chpsv)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chpsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *ap, Complex *afp, integer *ipiv, Complex *b, integer *
    ldb, Complex *x, integer *ldx, CReal *rcond, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(chptrd)(char *uplo, integer *n, Complex *ap, CReal *d__,
    CReal *e, Complex *tau, integer *info);

/* Subroutine */ int F77NAME(chptrf)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(chptri)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, Complex *work, integer *info);

/* Subroutine */ int F77NAME(chptrs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, Complex *h__, integer *ldh, Complex *w, Complex *
    vl, integer *ldvl, Complex *vr, integer *ldvr, integer *mm, integer *
    m, Complex *work, CReal *rwork, integer *ifaill, integer *ifailr,
    integer *info);

/* Subroutine */ int F77NAME(chseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, Complex *h__, integer *ldh, Complex *w, Complex *z__,
    integer *ldz, Complex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(clabrd)(integer *m, integer *n, integer *nb, Complex *a,
    integer *lda, CReal *d__, CReal *e, Complex *tauq, Complex *taup,
    Complex *x, integer *ldx, Complex *y, integer *ldy);

/* Subroutine */ int F77NAME(clacgv)(integer *n, Complex *x, integer *incx);

/* Subroutine */ int F77NAME(clacon)(integer *n, Complex *v, Complex *x, CReal *est,
    integer *kase);

/* Subroutine */ int F77NAME(clacp2)(char *uplo, integer *m, integer *n, CReal *a,
    integer *lda, Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(clacpy)(char *uplo, integer *m, integer *n, Complex *a,
    integer *lda, Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(clacrm)(integer *m, integer *n, Complex *a, integer *lda,
     CReal *b, integer *ldb, Complex *c__, integer *ldc, CReal *rwork);

/* Subroutine */ int F77NAME(clacrt)(integer *n, Complex *cx, integer *incx, Complex *
    cy, integer *incy, Complex *c__, Complex *s);

/* Subroutine */ int F77NAME(claed0)(integer *qsiz, integer *n, CReal *d__, CReal *e,
    Complex *q, integer *ldq, Complex *qstore, integer *ldqs, CReal *rwork,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(claed7)(integer *n, integer *cutpnt, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, CReal *d__, Complex *
    q, integer *ldq, CReal *rho, integer *indxq, CReal *qstore, integer *
    qptr, integer *prmptr, integer *perm, integer *givptr, integer *
    givcol, CReal *givnum, Complex *work, CReal *rwork, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(claed8)(integer *k, integer *n, integer *qsiz, Complex *
    q, integer *ldq, CReal *d__, CReal *rho, integer *cutpnt, CReal *z__,
    CReal *dlamda, Complex *q2, integer *ldq2, CReal *w, integer *indxp,
    integer *indx, integer *indxq, integer *perm, integer *givptr,
    integer *givcol, CReal *givnum, integer *info);

/* Subroutine */ int F77NAME(claein)(logical *rightv, logical *noinit, integer *n,
    Complex *h__, integer *ldh, Complex *w, Complex *v, Complex *b,
    integer *ldb, CReal *rwork, CReal *eps3, CReal *smlnum, integer *info);

/* Subroutine */ int F77NAME(claesy)(Complex *a, Complex *b, Complex *c__, Complex *
    rt1, Complex *rt2, Complex *evscal, Complex *cs1, Complex *sn1);

/* Subroutine */ int F77NAME(claev2)(Complex *a, Complex *b, Complex *c__, CReal *rt1,
    CReal *rt2, CReal *cs1, Complex *sn1);

/* Subroutine */ int F77NAME(clags2)(logical *upper, CReal *a1, Complex *a2, CReal *a3,
    CReal *b1, Complex *b2, CReal *b3, CReal *csu, Complex *snu, CReal *csv,
    Complex *snv, CReal *csq, Complex *snq);

/* Subroutine */ int F77NAME(clagtm)(char *trans, integer *n, integer *nrhs, CReal *
    alpha, Complex *dl, Complex *d__, Complex *du, Complex *x, integer *
    ldx, CReal *beta, Complex *b, integer *ldb);

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

/* Subroutine */ int F77NAME(claic1)(integer *job, integer *j, Complex *x, CReal *sest,
     Complex *w, Complex *gamma, CReal *sestpr, Complex *s, Complex *c__);

/* Subroutine */ int F77NAME(clals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, Complex *b, integer *ldb, Complex *bx,
    integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, CReal *givnum, integer *ldgnum, CReal *poles, CReal *
    difl, CReal *difr, CReal *z__, integer *k, CReal *c__, CReal *s, CReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(clalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, Complex *b, integer *ldb, Complex *bx, integer *ldbx,
    CReal *u, integer *ldu, CReal *vt, integer *k, CReal *difl, CReal *difr,
    CReal *z__, CReal *poles, integer *givptr, integer *givcol, integer *
    ldgcol, integer *perm, CReal *givnum, CReal *c__, CReal *s, CReal *rwork,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(clapll)(integer *n, Complex *x, integer *incx, Complex *
    y, integer *incy, CReal *ssmin);

/* Subroutine */ int F77NAME(clapmt)(logical *forwrd, integer *m, integer *n, Complex
    *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(claqgb)(integer *m, integer *n, integer *kl, integer *ku,
     Complex *ab, integer *ldab, CReal *r__, CReal *c__, CReal *rowcnd, CReal
    *colcnd, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(claqge)(integer *m, integer *n, Complex *a, integer *lda,
     CReal *r__, CReal *c__, CReal *rowcnd, CReal *colcnd, CReal *amax, char *
    equed);

/* Subroutine */ int F77NAME(claqhb)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, CReal *s, CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(claqhe)(char *uplo, integer *n, Complex *a, integer *lda,
     CReal *s, CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(claqhp)(char *uplo, integer *n, Complex *ap, CReal *s,
    CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(claqp2)(integer *m, integer *n, integer *offset, Complex
    *a, integer *lda, integer *jpvt, Complex *tau, CReal *vn1, CReal *vn2,
    Complex *work);

/* Subroutine */ int F77NAME(claqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, Complex *a, integer *lda, integer *jpvt, Complex *
    tau, CReal *vn1, CReal *vn2, Complex *auxv, Complex *f, integer *ldf);

/* Subroutine */ int F77NAME(claqsb)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, CReal *s, CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(claqsp)(char *uplo, integer *n, Complex *ap, CReal *s,
    CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(claqsy)(char *uplo, integer *n, Complex *a, integer *lda,
     CReal *s, CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(clar1v)(integer *n, integer *b1, integer *bn, CReal *
    sigma, CReal *d__, CReal *l, CReal *ld, CReal *lld, CReal *gersch, Complex
    *z__, CReal *ztz, CReal *mingma, integer *r__, integer *isuppz, CReal *
    work);

/* Subroutine */ int F77NAME(clar2v)(integer *n, Complex *x, Complex *y, Complex *z__,
     integer *incx, CReal *c__, Complex *s, integer *incc);

/* Subroutine */ int F77NAME(clarcm)(integer *m, integer *n, CReal *a, integer *lda,
    Complex *b, integer *ldb, Complex *c__, integer *ldc, CReal *rwork);

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
    y, integer *incy, CReal *c__, integer *incc);

/* Subroutine */ int F77NAME(clarnv)(integer *idist, integer *iseed, integer *n,
    Complex *x);

/* Subroutine */ int F77NAME(clarrv)(integer *n, CReal *d__, CReal *l, integer *isplit,
    integer *m, CReal *w, integer *iblock, CReal *gersch, CReal *tol,
    Complex *z__, integer *ldz, integer *isuppz, CReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(clartg)(Complex *f, Complex *g, CReal *cs, Complex *sn,
    Complex *r__);

/* Subroutine */ int F77NAME(clartv)(integer *n, Complex *x, integer *incx, Complex *
    y, integer *incy, CReal *c__, Complex *s, integer *incc);

/* Subroutine */ int F77NAME(clarz)(char *side, integer *m, integer *n, integer *l,
    Complex *v, integer *incv, Complex *tau, Complex *c__, integer *ldc,
    Complex *work);

/* Subroutine */ int F77NAME(clarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, Complex *v,
    integer *ldv, Complex *t, integer *ldt, Complex *c__, integer *ldc,
    Complex *work, integer *ldwork);

/* Subroutine */ int F77NAME(clarzt)(char *direct, char *storev, integer *n, integer *
    k, Complex *v, integer *ldv, Complex *tau, Complex *t, integer *ldt);

/* Subroutine */ int F77NAME(clascl)(char *type__, integer *kl, integer *ku, CReal *
    cfrom, CReal *cto, integer *m, integer *n, Complex *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(claset)(char *uplo, integer *m, integer *n, Complex *
    alpha, Complex *beta, Complex *a, integer *lda);

/* Subroutine */ int F77NAME(clasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, CReal *c__, CReal *s, Complex *a, integer *lda);

/* Subroutine */ int F77NAME(classq)(integer *n, Complex *x, integer *incx, CReal *
    scale, CReal *sumsq);

/* Subroutine */ int F77NAME(claswp)(integer *n, Complex *a, integer *lda, integer *
    k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(clasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     Complex *a, integer *lda, integer *ipiv, Complex *w, integer *ldw,
    integer *info);

/* Subroutine */ int F77NAME(clatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, Complex *ab, integer *ldab, Complex *
    x, CReal *scale, CReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(clatdf)(integer *ijob, integer *n, Complex *z__, integer
    *ldz, Complex *rhs, CReal *rdsum, CReal *rdscal, integer *ipiv, integer
    *jpiv);

/* Subroutine */ int F77NAME(clatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, Complex *ap, Complex *x, CReal *scale, CReal *cnorm,
     integer *info);

/* Subroutine */ int F77NAME(clatrd)(char *uplo, integer *n, integer *nb, Complex *a,
    integer *lda, CReal *e, Complex *tau, Complex *w, integer *ldw);

/* Subroutine */ int F77NAME(clatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, Complex *a, integer *lda, Complex *x, CReal *scale,
     CReal *cnorm, integer *info);

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
     integer *ldab, CReal *anorm, CReal *rcond, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cpbequ)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, CReal *s, CReal *scond, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(cpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Complex *ab, integer *ldab, Complex *afb, integer *ldafb,
    Complex *b, integer *ldb, Complex *x, integer *ldx, CReal *ferr, CReal *
    berr, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpbstf)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Complex *ab, integer *ldab, Complex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, Complex *ab, integer *ldab, Complex *afb, integer *
    ldafb, char *equed, CReal *s, Complex *b, integer *ldb, Complex *x,
    integer *ldx, CReal *rcond, CReal *ferr, CReal *berr, Complex *work,
    CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpbtf2)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbtrf)(char *uplo, integer *n, integer *kd, Complex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, Complex *ab, integer *ldab, Complex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cpocon)(char *uplo, integer *n, Complex *a, integer *lda,
     CReal *anorm, CReal *rcond, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpoequ)(integer *n, Complex *a, integer *lda, CReal *s,
    CReal *scond, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(cporfs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, Complex *b, integer *ldb,
     Complex *x, integer *ldx, CReal *ferr, CReal *berr, Complex *work,
    CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cposv)(char *uplo, integer *n, integer *nrhs, Complex *a,
     integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, char *
    equed, CReal *s, Complex *b, integer *ldb, Complex *x, integer *ldx,
    CReal *rcond, CReal *ferr, CReal *berr, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cpotf2)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotrf)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotri)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotrs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cppcon)(char *uplo, integer *n, Complex *ap, CReal *anorm,
     CReal *rcond, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cppequ)(char *uplo, integer *n, Complex *ap, CReal *s,
    CReal *scond, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(cpprfs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *afp, Complex *b, integer *ldb, Complex *x, integer *ldx,
    CReal *ferr, CReal *berr, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cppsv)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *ap, Complex *afp, char *equed, CReal *s, Complex *b,
    integer *ldb, Complex *x, integer *ldx, CReal *rcond, CReal *ferr, CReal
    *berr, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpptrf)(char *uplo, integer *n, Complex *ap, integer *
    info);

/* Subroutine */ int F77NAME(cpptri)(char *uplo, integer *n, Complex *ap, integer *
    info);

/* Subroutine */ int F77NAME(cpptrs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptcon)(integer *n, CReal *d__, Complex *e, CReal *anorm,
    CReal *rcond, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cptrfs)(char *uplo, integer *n, integer *nrhs, CReal *d__,
     Complex *e, CReal *df, Complex *ef, Complex *b, integer *ldb, Complex
    *x, integer *ldx, CReal *ferr, CReal *berr, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cptsv)(integer *n, integer *nrhs, CReal *d__, Complex *e,
    Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptsvx)(char *fact, integer *n, integer *nrhs, CReal *d__,
     Complex *e, CReal *df, Complex *ef, Complex *b, integer *ldb, Complex
    *x, integer *ldx, CReal *rcond, CReal *ferr, CReal *berr, Complex *work,
    CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpttrf)(integer *n, CReal *d__, Complex *e, integer *info);

/* Subroutine */ int F77NAME(cpttrs)(char *uplo, integer *n, integer *nrhs, CReal *d__,
     Complex *e, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptts2)(integer *iuplo, integer *n, integer *nrhs, CReal *
    d__, Complex *e, Complex *b, integer *ldb);

/* Subroutine */ int F77NAME(crot)(integer *n, Complex *cx, integer *incx, Complex *
    cy, integer *incy, CReal *c__, Complex *s);

/* Subroutine */ int F77NAME(cspcon)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, CReal *anorm, CReal *rcond, Complex *work, integer *info);

/* Subroutine */ int F77NAME(cspmv)(char *uplo, integer *n, Complex *alpha, Complex *
    ap, Complex *x, integer *incx, Complex *beta, Complex *y, integer *
    incy);

/* Subroutine */ int F77NAME(cspr)(char *uplo, integer *n, Complex *alpha, Complex *x,
     integer *incx, Complex *ap);

/* Subroutine */ int F77NAME(csprfs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, Complex *afp, integer *ipiv, Complex *b, integer *ldb, Complex *x,
     integer *ldx, CReal *ferr, CReal *berr, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cspsv)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *ap, Complex *afp, integer *ipiv, Complex *b, integer *
    ldb, Complex *x, integer *ldx, CReal *rcond, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(csptrf)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(csptri)(char *uplo, integer *n, Complex *ap, integer *
    ipiv, Complex *work, integer *info);

/* Subroutine */ int F77NAME(csptrs)(char *uplo, integer *n, integer *nrhs, Complex *
    ap, integer *ipiv, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(csrot)(integer *n, Complex *cx, integer *incx, Complex *
    cy, integer *incy, CReal *c__, CReal *s);

/* Subroutine */ int F77NAME(csrscl)(integer *n, CReal *sa, Complex *sx, integer *incx);

/* Subroutine */ int F77NAME(cstedc)(char *compz, integer *n, CReal *d__, CReal *e,
    Complex *z__, integer *ldz, Complex *work, integer *lwork, CReal *
    rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(cstein)(integer *n, CReal *d__, CReal *e, integer *m, CReal
    *w, integer *iblock, integer *isplit, Complex *z__, integer *ldz,
    CReal *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(csteqr)(char *compz, integer *n, CReal *d__, CReal *e,
    Complex *z__, integer *ldz, CReal *work, integer *info);

/* Subroutine */ int F77NAME(csycon)(char *uplo, integer *n, Complex *a, integer *lda,
     integer *ipiv, CReal *anorm, CReal *rcond, Complex *work, integer *
    info);

/* Subroutine */ int F77NAME(csymv)(char *uplo, integer *n, Complex *alpha, Complex *
    a, integer *lda, Complex *x, integer *incx, Complex *beta, Complex *y,
     integer *incy);

/* Subroutine */ int F77NAME(csyr)(char *uplo, integer *n, Complex *alpha, Complex *x,
     integer *incx, Complex *a, integer *lda);

/* Subroutine */ int F77NAME(csyrfs)(char *uplo, integer *n, integer *nrhs, Complex *
    a, integer *lda, Complex *af, integer *ldaf, integer *ipiv, Complex *
    b, integer *ldb, Complex *x, integer *ldx, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(csysv)(char *uplo, integer *n, integer *nrhs, Complex *a,
     integer *lda, integer *ipiv, Complex *b, integer *ldb, Complex *work,
     integer *lwork, integer *info);

/* Subroutine */ int F77NAME(csysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, Complex *a, integer *lda, Complex *af, integer *ldaf, integer *
    ipiv, Complex *b, integer *ldb, Complex *x, integer *ldx, CReal *rcond,
     CReal *ferr, CReal *berr, Complex *work, integer *lwork, CReal *rwork,
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
    integer *kd, Complex *ab, integer *ldab, CReal *rcond, Complex *work,
    CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, Complex *ab, integer *ldab, Complex *b,
    integer *ldb, Complex *x, integer *ldx, CReal *ferr, CReal *berr,
    Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, Complex *ab, integer *ldab, Complex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ctgevc)(char *side, char *howmny, logical *select,
    integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, integer *mm,
    integer *m, Complex *work, CReal *rwork, integer *info);

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
     Complex *z__, integer *ldz, integer *m, CReal *pl, CReal *pr, CReal *
    dif, Complex *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ctgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, Complex *a, integer *
    lda, Complex *b, integer *ldb, CReal *tola, CReal *tolb, CReal *alpha,
    CReal *beta, Complex *u, integer *ldu, Complex *v, integer *ldv,
    Complex *q, integer *ldq, Complex *work, integer *ncycle, integer *
    info);

/* Subroutine */ int F77NAME(ctgsna)(char *job, char *howmny, logical *select,
    integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *vl, integer *ldvl, Complex *vr, integer *ldvr, CReal *s, CReal
    *dif, integer *mm, integer *m, Complex *work, integer *lwork, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(ctgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, Complex *c__,
    integer *ldc, Complex *d__, integer *ldd, Complex *e, integer *lde,
    Complex *f, integer *ldf, CReal *scale, CReal *rdsum, CReal *rdscal,
    integer *info);

/* Subroutine */ int F77NAME(ctgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, Complex *a, integer *lda, Complex *b, integer *ldb, Complex *c__,
    integer *ldc, Complex *d__, integer *ldd, Complex *e, integer *lde,
    Complex *f, integer *ldf, CReal *scale, CReal *dif, Complex *work,
    integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ctpcon)(char *norm, char *uplo, char *diag, integer *n,
    Complex *ap, CReal *rcond, Complex *work, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Complex *ap, Complex *b, integer *ldb, Complex *x,
    integer *ldx, CReal *ferr, CReal *berr, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(ctptri)(char *uplo, char *diag, integer *n, Complex *ap,
    integer *info);

/* Subroutine */ int F77NAME(ctptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Complex *ap, Complex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ctrcon)(char *norm, char *uplo, char *diag, integer *n,
    Complex *a, integer *lda, CReal *rcond, Complex *work, CReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(ctrevc)(char *side, char *howmny, logical *select,
    integer *n, Complex *t, integer *ldt, Complex *vl, integer *ldvl,
    Complex *vr, integer *ldvr, integer *mm, integer *m, Complex *work,
    CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrexc)(char *compq, integer *n, Complex *t, integer *
    ldt, Complex *q, integer *ldq, integer *ifst, integer *ilst, integer *
    info);

/* Subroutine */ int F77NAME(ctrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *x, integer *ldx, CReal *ferr, CReal *berr, Complex *work, CReal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrsen)(char *job, char *compq, logical *select, integer
    *n, Complex *t, integer *ldt, Complex *q, integer *ldq, Complex *w,
    integer *m, CReal *s, CReal *sep, Complex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(ctrsna)(char *job, char *howmny, logical *select,
    integer *n, Complex *t, integer *ldt, Complex *vl, integer *ldvl,
    Complex *vr, integer *ldvr, CReal *s, CReal *sep, integer *mm, integer *
    m, Complex *work, integer *ldwork, CReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, Complex *a, integer *lda, Complex *b, integer *ldb,
    Complex *c__, integer *ldc, CReal *scale, integer *info);

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

/* Subroutine */ int F77NAME(dbdsdc)(char *uplo, char *compq, integer *n, doubleCReal *
    d__, doubleCReal *e, doubleCReal *u, integer *ldu, doubleCReal *vt,
    integer *ldvt, doubleCReal *q, integer *iq, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, doubleCReal *d__, doubleCReal *e, doubleCReal *vt,
    integer *ldvt, doubleCReal *u, integer *ldu, doubleCReal *c__, integer *
    ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(ddisna)(char *job, integer *m, integer *n, doubleCReal *
    d__, doubleCReal *sep, integer *info);

/* Subroutine */ int F77NAME(dgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, doubleCReal *ab, integer *ldab, doubleCReal *
    d__, doubleCReal *e, doubleCReal *q, integer *ldq, doubleCReal *pt,
    integer *ldpt, doubleCReal *c__, integer *ldc, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     doubleCReal *ab, integer *ldab, integer *ipiv, doubleCReal *anorm,
    doubleCReal *rcond, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     doubleCReal *ab, integer *ldab, doubleCReal *r__, doubleCReal *c__,
    doubleCReal *rowcnd, doubleCReal *colcnd, doubleCReal *amax, integer *
    info);

/* Subroutine */ int F77NAME(dgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleCReal *ab, integer *ldab, doubleCReal *afb,
    integer *ldafb, integer *ipiv, doubleCReal *b, integer *ldb,
    doubleCReal *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, doubleCReal *ab, integer *ldab, integer *ipiv, doubleCReal *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, doubleCReal *ab, integer *ldab,
    doubleCReal *afb, integer *ldafb, integer *ipiv, char *equed,
    doubleCReal *r__, doubleCReal *c__, doubleCReal *b, integer *ldb,
    doubleCReal *x, integer *ldx, doubleCReal *rcond, doubleCReal *ferr,
    doubleCReal *berr, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     doubleCReal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     doubleCReal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleCReal *ab, integer *ldab, integer *ipiv,
    doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doubleCReal *scale, integer *m, doubleCReal *v, integer *
    ldv, integer *info);

/* Subroutine */ int F77NAME(dgebal)(char *job, integer *n, doubleCReal *a, integer *
    lda, integer *ilo, integer *ihi, doubleCReal *scale, integer *info);

/* Subroutine */ int F77NAME(dgebd2)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *d__, doubleCReal *e, doubleCReal *tauq, doubleCReal *
    taup, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dgebrd)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *d__, doubleCReal *e, doubleCReal *tauq, doubleCReal *
    taup, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgecon)(char *norm, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *anorm, doubleCReal *rcond, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgeequ)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *r__, doubleCReal *c__, doubleCReal *rowcnd, doubleCReal
    *colcnd, doubleCReal *amax, integer *info);

/* Subroutine */ int F77NAME(dgees)(char *jobvs, char *sort, L_fp select, integer *n,
    doubleCReal *a, integer *lda, integer *sdim, doubleCReal *wr,
    doubleCReal *wi, doubleCReal *vs, integer *ldvs, doubleCReal *work,
    integer *lwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, doubleCReal *a, integer *lda, integer *sdim,
    doubleCReal *wr, doubleCReal *wi, doubleCReal *vs, integer *ldvs,
    doubleCReal *rconde, doubleCReal *rcondv, doubleCReal *work, integer *
    lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dgeev)(char *jobvl, char *jobvr, integer *n, doubleCReal *
    a, integer *lda, doubleCReal *wr, doubleCReal *wi, doubleCReal *vl,
    integer *ldvl, doubleCReal *vr, integer *ldvr, doubleCReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleCReal *a, integer *lda, doubleCReal *wr,
    doubleCReal *wi, doubleCReal *vl, integer *ldvl, doubleCReal *vr,
    integer *ldvr, integer *ilo, integer *ihi, doubleCReal *scale,
    doubleCReal *abnrm, doubleCReal *rconde, doubleCReal *rcondv, doubleCReal
    *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgegs)(char *jobvsl, char *jobvsr, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *
    alphar, doubleCReal *alphai, doubleCReal *beta, doubleCReal *vsl,
    integer *ldvsl, doubleCReal *vsr, integer *ldvsr, doubleCReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgegv)(char *jobvl, char *jobvr, integer *n, doubleCReal *
    a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *alphar,
    doubleCReal *alphai, doubleCReal *beta, doubleCReal *vl, integer *ldvl,
    doubleCReal *vr, integer *ldvr, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgehd2)(integer *n, integer *ilo, integer *ihi,
    doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dgehrd)(integer *n, integer *ilo, integer *ihi,
    doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgelq2)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dgelqf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgels)(char *trans, integer *m, integer *n, integer *
    nrhs, doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb,
    doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgelsd)(integer *m, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *
    s, doubleCReal *rcond, integer *rank, doubleCReal *work, integer *lwork,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgelss)(integer *m, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *
    s, doubleCReal *rcond, integer *rank, doubleCReal *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(dgelsx)(integer *m, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, integer *
    jpvt, doubleCReal *rcond, integer *rank, doubleCReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dgelsy)(integer *m, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, integer *
    jpvt, doubleCReal *rcond, integer *rank, doubleCReal *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(dgeql2)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqlf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgeqp3)(integer *m, integer *n, doubleCReal *a, integer *
    lda, integer *jpvt, doubleCReal *tau, doubleCReal *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(dgeqpf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, integer *jpvt, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqr2)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqrf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgerfs)(char *trans, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *af, integer *ldaf, integer *
    ipiv, doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx,
    doubleCReal *ferr, doubleCReal *berr, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dgerq2)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dgerqf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgesc2)(integer *n, doubleCReal *a, integer *lda,
    doubleCReal *rhs, integer *ipiv, integer *jpiv, doubleCReal *scale);

/* Subroutine */ int F77NAME(dgesdd)(char *jobz, integer *m, integer *n, doubleCReal *
    a, integer *lda, doubleCReal *s, doubleCReal *u, integer *ldu,
    doubleCReal *vt, integer *ldvt, doubleCReal *work, integer *lwork,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgesv)(integer *n, integer *nrhs, doubleCReal *a, integer
    *lda, integer *ipiv, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgesvd)(char *jobu, char *jobvt, integer *m, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *s, doubleCReal *u, integer *
    ldu, doubleCReal *vt, integer *ldvt, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleCReal *a, integer *lda, doubleCReal *af, integer *ldaf,
    integer *ipiv, char *equed, doubleCReal *r__, doubleCReal *c__,
    doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *
    rcond, doubleCReal *ferr, doubleCReal *berr, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgetc2)(integer *n, doubleCReal *a, integer *lda, integer
    *ipiv, integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(dgetf2)(integer *m, integer *n, doubleCReal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgetrf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgetri)(integer *n, doubleCReal *a, integer *lda, integer
    *ipiv, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgetrs)(char *trans, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, integer *ipiv, doubleCReal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doubleCReal *lscale, doubleCReal *rscale, integer *m,
    doubleCReal *v, integer *ldv, integer *info);

/* Subroutine */ int F77NAME(dggbal)(char *job, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *b, integer *ldb, integer *ilo, integer *ihi,
    doubleCReal *lscale, doubleCReal *rscale, doubleCReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, integer *n, doubleCReal *a, integer *lda, doubleCReal *b,
    integer *ldb, integer *sdim, doubleCReal *alphar, doubleCReal *alphai,
    doubleCReal *beta, doubleCReal *vsl, integer *ldvsl, doubleCReal *vsr,
    integer *ldvsr, doubleCReal *work, integer *lwork, logical *bwork,
    integer *info);

/* Subroutine */ int F77NAME(dggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, char *sense, integer *n, doubleCReal *a, integer *lda,
    doubleCReal *b, integer *ldb, integer *sdim, doubleCReal *alphar,
    doubleCReal *alphai, doubleCReal *beta, doubleCReal *vsl, integer *ldvsl,
     doubleCReal *vsr, integer *ldvsr, doubleCReal *rconde, doubleCReal *
    rcondv, doubleCReal *work, integer *lwork, integer *iwork, integer *
    liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dggev)(char *jobvl, char *jobvr, integer *n, doubleCReal *
    a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *alphar,
    doubleCReal *alphai, doubleCReal *beta, doubleCReal *vl, integer *ldvl,
    doubleCReal *vr, integer *ldvr, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleCReal *a, integer *lda, doubleCReal *b,
    integer *ldb, doubleCReal *alphar, doubleCReal *alphai, doubleCReal *
    beta, doubleCReal *vl, integer *ldvl, doubleCReal *vr, integer *ldvr,
    integer *ilo, integer *ihi, doubleCReal *lscale, doubleCReal *rscale,
    doubleCReal *abnrm, doubleCReal *bbnrm, doubleCReal *rconde, doubleCReal *
    rcondv, doubleCReal *work, integer *lwork, integer *iwork, logical *
    bwork, integer *info);

/* Subroutine */ int F77NAME(dggglm)(integer *n, integer *m, integer *p, doubleCReal *
    a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *d__,
    doubleCReal *x, doubleCReal *y, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, doubleCReal *a, integer *lda, doubleCReal *b,
    integer *ldb, doubleCReal *q, integer *ldq, doubleCReal *z__, integer *
    ldz, integer *info);

/* Subroutine */ int F77NAME(dgglse)(integer *m, integer *n, integer *p, doubleCReal *
    a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *c__,
    doubleCReal *d__, doubleCReal *x, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dggqrf)(integer *n, integer *m, integer *p, doubleCReal *
    a, integer *lda, doubleCReal *taua, doubleCReal *b, integer *ldb,
    doubleCReal *taub, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dggrqf)(integer *m, integer *p, integer *n, doubleCReal *
    a, integer *lda, doubleCReal *taua, doubleCReal *b, integer *ldb,
    doubleCReal *taub, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, doubleCReal *a,
    integer *lda, doubleCReal *b, integer *ldb, doubleCReal *alpha,
    doubleCReal *beta, doubleCReal *u, integer *ldu, doubleCReal *v, integer
    *ldv, doubleCReal *q, integer *ldq, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, doubleCReal *a, integer *lda, doubleCReal *b,
    integer *ldb, doubleCReal *tola, doubleCReal *tolb, integer *k, integer
    *l, doubleCReal *u, integer *ldu, doubleCReal *v, integer *ldv,
    doubleCReal *q, integer *ldq, integer *iwork, doubleCReal *tau,
    doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dgtcon)(char *norm, integer *n, doubleCReal *dl,
    doubleCReal *d__, doubleCReal *du, doubleCReal *du2, integer *ipiv,
    doubleCReal *anorm, doubleCReal *rcond, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgtrfs)(char *trans, integer *n, integer *nrhs,
    doubleCReal *dl, doubleCReal *d__, doubleCReal *du, doubleCReal *dlf,
    doubleCReal *df, doubleCReal *duf, doubleCReal *du2, integer *ipiv,
    doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *
    ferr, doubleCReal *berr, doubleCReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dgtsv)(integer *n, integer *nrhs, doubleCReal *dl,
    doubleCReal *d__, doubleCReal *du, doubleCReal *b, integer *ldb, integer
    *info);

/* Subroutine */ int F77NAME(dgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleCReal *dl, doubleCReal *d__, doubleCReal *du, doubleCReal *
    dlf, doubleCReal *df, doubleCReal *duf, doubleCReal *du2, integer *ipiv,
    doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *
    rcond, doubleCReal *ferr, doubleCReal *berr, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgttrf)(integer *n, doubleCReal *dl, doubleCReal *d__,
    doubleCReal *du, doubleCReal *du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgttrs)(char *trans, integer *n, integer *nrhs,
    doubleCReal *dl, doubleCReal *d__, doubleCReal *du, doubleCReal *du2,
    integer *ipiv, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgtts2)(integer *itrans, integer *n, integer *nrhs,
    doubleCReal *dl, doubleCReal *d__, doubleCReal *du, doubleCReal *du2,
    integer *ipiv, doubleCReal *b, integer *ldb);

/* Subroutine */ int F77NAME(dhgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, doubleCReal *a, integer *lda, doubleCReal *
    b, integer *ldb, doubleCReal *alphar, doubleCReal *alphai, doubleCReal *
    beta, doubleCReal *q, integer *ldq, doubleCReal *z__, integer *ldz,
    doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dhsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, doubleCReal *h__, integer *ldh, doubleCReal *wr,
    doubleCReal *wi, doubleCReal *vl, integer *ldvl, doubleCReal *vr,
    integer *ldvr, integer *mm, integer *m, doubleCReal *work, integer *
    ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(dhseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, doubleCReal *h__, integer *ldh, doubleCReal *wr,
    doubleCReal *wi, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dlabad)(doubleCReal *small, doubleCReal *large);

/* Subroutine */ int F77NAME(dlabrd)(integer *m, integer *n, integer *nb, doubleCReal *
    a, integer *lda, doubleCReal *d__, doubleCReal *e, doubleCReal *tauq,
    doubleCReal *taup, doubleCReal *x, integer *ldx, doubleCReal *y, integer
    *ldy);

/* Subroutine */ int F77NAME(dlacon)(integer *n, doubleCReal *v, doubleCReal *x,
    integer *isgn, doubleCReal *est, integer *kase);

/* Subroutine */ int F77NAME(dlacpy)(char *uplo, integer *m, integer *n, doubleCReal *
    a, integer *lda, doubleCReal *b, integer *ldb);

/* Subroutine */ int F77NAME(dladiv)(doubleCReal *a, doubleCReal *b, doubleCReal *c__,
    doubleCReal *d__, doubleCReal *p, doubleCReal *q);

/* Subroutine */ int F77NAME(dlae2)(doubleCReal *a, doubleCReal *b, doubleCReal *c__,
    doubleCReal *rt1, doubleCReal *rt2);

/* Subroutine */ int F77NAME(dlaebz)(integer *ijob, integer *nitmax, integer *n,
    integer *mmax, integer *minp, integer *nbmin, doubleCReal *abstol,
    doubleCReal *reltol, doubleCReal *pivmin, doubleCReal *d__, doubleCReal *
    e, doubleCReal *e2, integer *nval, doubleCReal *ab, doubleCReal *c__,
    integer *mout, integer *nab, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlaed0)(integer *icompq, integer *qsiz, integer *n,
    doubleCReal *d__, doubleCReal *e, doubleCReal *q, integer *ldq,
    doubleCReal *qstore, integer *ldqs, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlaed1)(integer *n, doubleCReal *d__, doubleCReal *q,
    integer *ldq, integer *indxq, doubleCReal *rho, integer *cutpnt,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlaed2)(integer *k, integer *n, integer *n1, doubleCReal *
    d__, doubleCReal *q, integer *ldq, integer *indxq, doubleCReal *rho,
    doubleCReal *z__, doubleCReal *dlamda, doubleCReal *w, doubleCReal *q2,
    integer *indx, integer *indxc, integer *indxp, integer *coltyp,
    integer *info);

/* Subroutine */ int F77NAME(dlaed3)(integer *k, integer *n, integer *n1, doubleCReal *
    d__, doubleCReal *q, integer *ldq, doubleCReal *rho, doubleCReal *dlamda,
     doubleCReal *q2, integer *indx, integer *ctot, doubleCReal *w,
    doubleCReal *s, integer *info);

/* Subroutine */ int F77NAME(dlaed4)(integer *n, integer *i__, doubleCReal *d__,
    doubleCReal *z__, doubleCReal *delta, doubleCReal *rho, doubleCReal *dlam,
     integer *info);

/* Subroutine */ int F77NAME(dlaed5)(integer *i__, doubleCReal *d__, doubleCReal *z__,
    doubleCReal *delta, doubleCReal *rho, doubleCReal *dlam);

/* Subroutine */ int F77NAME(dlaed6)(integer *kniter, logical *orgati, doubleCReal *
    rho, doubleCReal *d__, doubleCReal *z__, doubleCReal *finit, doubleCReal *
    tau, integer *info);

/* Subroutine */ int F77NAME(dlaed7)(integer *icompq, integer *n, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, doubleCReal *d__,
    doubleCReal *q, integer *ldq, integer *indxq, doubleCReal *rho, integer
    *cutpnt, doubleCReal *qstore, integer *qptr, integer *prmptr, integer *
    perm, integer *givptr, integer *givcol, doubleCReal *givnum,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlaed8)(integer *icompq, integer *k, integer *n, integer
    *qsiz, doubleCReal *d__, doubleCReal *q, integer *ldq, integer *indxq,
    doubleCReal *rho, integer *cutpnt, doubleCReal *z__, doubleCReal *dlamda,
     doubleCReal *q2, integer *ldq2, doubleCReal *w, integer *perm, integer
    *givptr, integer *givcol, doubleCReal *givnum, integer *indxp, integer
    *indx, integer *info);

/* Subroutine */ int F77NAME(dlaed9)(integer *k, integer *kstart, integer *kstop,
    integer *n, doubleCReal *d__, doubleCReal *q, integer *ldq, doubleCReal *
    rho, doubleCReal *dlamda, doubleCReal *w, doubleCReal *s, integer *lds,
    integer *info);

/* Subroutine */ int F77NAME(dlaeda)(integer *n, integer *tlvls, integer *curlvl,
    integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
    integer *givcol, doubleCReal *givnum, doubleCReal *q, integer *qptr,
    doubleCReal *z__, doubleCReal *ztemp, integer *info);

/* Subroutine */ int F77NAME(dlaein)(logical *rightv, logical *noinit, integer *n,
    doubleCReal *h__, integer *ldh, doubleCReal *wr, doubleCReal *wi,
    doubleCReal *vr, doubleCReal *vi, doubleCReal *b, integer *ldb,
    doubleCReal *work, doubleCReal *eps3, doubleCReal *smlnum, doubleCReal *
    bignum, integer *info);

/* Subroutine */ int F77NAME(dlaev2)(doubleCReal *a, doubleCReal *b, doubleCReal *c__,
    doubleCReal *rt1, doubleCReal *rt2, doubleCReal *cs1, doubleCReal *sn1);

/* Subroutine */ int F77NAME(dlaexc)(logical *wantq, integer *n, doubleCReal *t,
    integer *ldt, doubleCReal *q, integer *ldq, integer *j1, integer *n1,
    integer *n2, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dlag2)(doubleCReal *a, integer *lda, doubleCReal *b,
    integer *ldb, doubleCReal *safmin, doubleCReal *scale1, doubleCReal *
    scale2, doubleCReal *wr1, doubleCReal *wr2, doubleCReal *wi);

/* Subroutine */ int F77NAME(dlags2)(logical *upper, doubleCReal *a1, doubleCReal *a2,
    doubleCReal *a3, doubleCReal *b1, doubleCReal *b2, doubleCReal *b3,
    doubleCReal *csu, doubleCReal *snu, doubleCReal *csv, doubleCReal *snv,
    doubleCReal *csq, doubleCReal *snq);

/* Subroutine */ int F77NAME(dlagtf)(integer *n, doubleCReal *a, doubleCReal *lambda,
    doubleCReal *b, doubleCReal *c__, doubleCReal *tol, doubleCReal *d__,
    integer *in, integer *info);

/* Subroutine */ int F77NAME(dlagtm)(char *trans, integer *n, integer *nrhs,
    doubleCReal *alpha, doubleCReal *dl, doubleCReal *d__, doubleCReal *du,
    doubleCReal *x, integer *ldx, doubleCReal *beta, doubleCReal *b, integer
    *ldb);

/* Subroutine */ int F77NAME(dlagts)(integer *job, integer *n, doubleCReal *a,
    doubleCReal *b, doubleCReal *c__, doubleCReal *d__, integer *in,
    doubleCReal *y, doubleCReal *tol, integer *info);

/* Subroutine */ int F77NAME(dlagv2)(doubleCReal *a, integer *lda, doubleCReal *b,
    integer *ldb, doubleCReal *alphar, doubleCReal *alphai, doubleCReal *
    beta, doubleCReal *csl, doubleCReal *snl, doubleCReal *csr, doubleCReal *
    snr);

/* Subroutine */ int F77NAME(dlahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, doubleCReal *h__, integer *ldh, doubleCReal
    *wr, doubleCReal *wi, integer *iloz, integer *ihiz, doubleCReal *z__,
    integer *ldz, integer *info);

/* Subroutine */ int F77NAME(dlahrd)(integer *n, integer *k, integer *nb, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *t, integer *ldt,
    doubleCReal *y, integer *ldy);

/* Subroutine */ int F77NAME(dlaic1)(integer *job, integer *j, doubleCReal *x,
    doubleCReal *sest, doubleCReal *w, doubleCReal *gamma, doubleCReal *
    sestpr, doubleCReal *s, doubleCReal *c__);

/* Subroutine */ int F77NAME(dlaln2)(logical *ltrans, integer *na, integer *nw,
    doubleCReal *smin, doubleCReal *ca, doubleCReal *a, integer *lda,
    doubleCReal *d1, doubleCReal *d2, doubleCReal *b, integer *ldb,
    doubleCReal *wr, doubleCReal *wi, doubleCReal *x, integer *ldx,
    doubleCReal *scale, doubleCReal *xnorm, integer *info);

/* Subroutine */ int F77NAME(dlals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, doubleCReal *b, integer *ldb, doubleCReal
    *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, doubleCReal *givnum, integer *ldgnum, doubleCReal *
    poles, doubleCReal *difl, doubleCReal *difr, doubleCReal *z__, integer *
    k, doubleCReal *c__, doubleCReal *s, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dlalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, doubleCReal *b, integer *ldb, doubleCReal *bx, integer *
    ldbx, doubleCReal *u, integer *ldu, doubleCReal *vt, integer *k,
    doubleCReal *difl, doubleCReal *difr, doubleCReal *z__, doubleCReal *
    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
    perm, doubleCReal *givnum, doubleCReal *c__, doubleCReal *s, doubleCReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlalsd)(char *uplo, integer *smlsiz, integer *n, integer
    *nrhs, doubleCReal *d__, doubleCReal *e, doubleCReal *b, integer *ldb,
    doubleCReal *rcond, integer *rank, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlamc1)(integer *beta, integer *t, logical *rnd, logical
    *ieee1);

/* Subroutine */ int F77NAME(dlamc2)(integer *beta, integer *t, logical *rnd,
    doubleCReal *eps, integer *emin, doubleCReal *rmin, integer *emax,
    doubleCReal *rmax);

/* Subroutine */ int F77NAME(dlamc4)(integer *emin, doubleCReal *start, integer *base);

/* Subroutine */ int F77NAME(dlamc5)(integer *beta, integer *p, integer *emin,
    logical *ieee, integer *emax, doubleCReal *rmax);

/* Subroutine */ int F77NAME(dlamrg)(integer *n1, integer *n2, doubleCReal *a, integer
    *dtrd1, integer *dtrd2, integer *index);

/* Subroutine */ int F77NAME(dlanv2)(doubleCReal *a, doubleCReal *b, doubleCReal *c__,
    doubleCReal *d__, doubleCReal *rt1r, doubleCReal *rt1i, doubleCReal *rt2r,
     doubleCReal *rt2i, doubleCReal *cs, doubleCReal *sn);

/* Subroutine */ int F77NAME(dlapll)(integer *n, doubleCReal *x, integer *incx,
    doubleCReal *y, integer *incy, doubleCReal *ssmin);

/* Subroutine */ int F77NAME(dlapmt)(logical *forwrd, integer *m, integer *n,
    doubleCReal *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(dlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     doubleCReal *ab, integer *ldab, doubleCReal *r__, doubleCReal *c__,
    doubleCReal *rowcnd, doubleCReal *colcnd, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqge)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *r__, doubleCReal *c__, doubleCReal *rowcnd, doubleCReal
    *colcnd, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqp2)(integer *m, integer *n, integer *offset,
    doubleCReal *a, integer *lda, integer *jpvt, doubleCReal *tau,
    doubleCReal *vn1, doubleCReal *vn2, doubleCReal *work);

/* Subroutine */ int F77NAME(dlaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, doubleCReal *a, integer *lda, integer *jpvt,
    doubleCReal *tau, doubleCReal *vn1, doubleCReal *vn2, doubleCReal *auxv,
    doubleCReal *f, integer *ldf);

/* Subroutine */ int F77NAME(dlaqsb)(char *uplo, integer *n, integer *kd, doubleCReal *
    ab, integer *ldab, doubleCReal *s, doubleCReal *scond, doubleCReal *amax,
     char *equed);

/* Subroutine */ int F77NAME(dlaqsp)(char *uplo, integer *n, doubleCReal *ap,
    doubleCReal *s, doubleCReal *scond, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqsy)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *s, doubleCReal *scond, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqtr)(logical *ltran, logical *lCReal, integer *n,
    doubleCReal *t, integer *ldt, doubleCReal *b, doubleCReal *w, doubleCReal
    *scale, doubleCReal *x, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dlar1v)(integer *n, integer *b1, integer *bn, doubleCReal
    *sigma, doubleCReal *d__, doubleCReal *l, doubleCReal *ld, doubleCReal *
    lld, doubleCReal *gersch, doubleCReal *z__, doubleCReal *ztz, doubleCReal
    *mingma, integer *r__, integer *isuppz, doubleCReal *work);

/* Subroutine */ int F77NAME(dlar2v)(integer *n, doubleCReal *x, doubleCReal *y,
    doubleCReal *z__, integer *incx, doubleCReal *c__, doubleCReal *s,
    integer *incc);

/* Subroutine */ int F77NAME(dlarf)(char *side, integer *m, integer *n, doubleCReal *v,
     integer *incv, doubleCReal *tau, doubleCReal *c__, integer *ldc,
    doubleCReal *work);

/* Subroutine */ int F77NAME(dlarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, doubleCReal *v, integer *
    ldv, doubleCReal *t, integer *ldt, doubleCReal *c__, integer *ldc,
    doubleCReal *work, integer *ldwork);

/* Subroutine */ int F77NAME(dlarfg)(integer *n, doubleCReal *alpha, doubleCReal *x,
    integer *incx, doubleCReal *tau);

/* Subroutine */ int F77NAME(dlarft)(char *direct, char *storev, integer *n, integer *
    k, doubleCReal *v, integer *ldv, doubleCReal *tau, doubleCReal *t,
    integer *ldt);

/* Subroutine */ int F77NAME(dlarfx)(char *side, integer *m, integer *n, doubleCReal *
    v, doubleCReal *tau, doubleCReal *c__, integer *ldc, doubleCReal *work);

/* Subroutine */ int F77NAME(dlargv)(integer *n, doubleCReal *x, integer *incx,
    doubleCReal *y, integer *incy, doubleCReal *c__, integer *incc);

/* Subroutine */ int F77NAME(dlarnv)(integer *idist, integer *iseed, integer *n,
    doubleCReal *x);

/* Subroutine */ int F77NAME(dlarrb)(integer *n, doubleCReal *d__, doubleCReal *l,
    doubleCReal *ld, doubleCReal *lld, integer *ifirst, integer *ilast,
    doubleCReal *sigma, doubleCReal *reltol, doubleCReal *w, doubleCReal *
    wgap, doubleCReal *werr, doubleCReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dlarre)(integer *n, doubleCReal *d__, doubleCReal *e,
    doubleCReal *tol, integer *nsplit, integer *isplit, integer *m,
    doubleCReal *w, doubleCReal *woff, doubleCReal *gersch, doubleCReal *work,
     integer *info);

/* Subroutine */ int F77NAME(dlarrf)(integer *n, doubleCReal *d__, doubleCReal *l,
    doubleCReal *ld, doubleCReal *lld, integer *ifirst, integer *ilast,
    doubleCReal *w, doubleCReal *dplus, doubleCReal *lplus, doubleCReal *work,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlarrv)(integer *n, doubleCReal *d__, doubleCReal *l,
    integer *isplit, integer *m, doubleCReal *w, integer *iblock,
    doubleCReal *gersch, doubleCReal *tol, doubleCReal *z__, integer *ldz,
    integer *isuppz, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlartg)(doubleCReal *f, doubleCReal *g, doubleCReal *cs,
    doubleCReal *sn, doubleCReal *r__);

/* Subroutine */ int F77NAME(dlartv)(integer *n, doubleCReal *x, integer *incx,
    doubleCReal *y, integer *incy, doubleCReal *c__, doubleCReal *s, integer
    *incc);

/* Subroutine */ int F77NAME(dlaruv)(integer *iseed, integer *n, doubleCReal *x);

/* Subroutine */ int F77NAME(dlarz)(char *side, integer *m, integer *n, integer *l,
    doubleCReal *v, integer *incv, doubleCReal *tau, doubleCReal *c__,
    integer *ldc, doubleCReal *work);

/* Subroutine */ int F77NAME(dlarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, doubleCReal *v,
     integer *ldv, doubleCReal *t, integer *ldt, doubleCReal *c__, integer *
    ldc, doubleCReal *work, integer *ldwork);

/* Subroutine */ int F77NAME(dlarzt)(char *direct, char *storev, integer *n, integer *
    k, doubleCReal *v, integer *ldv, doubleCReal *tau, doubleCReal *t,
    integer *ldt);

/* Subroutine */ int F77NAME(dlas2)(doubleCReal *f, doubleCReal *g, doubleCReal *h__,
    doubleCReal *ssmin, doubleCReal *ssmax);

/* Subroutine */ int F77NAME(dlascl)(char *type__, integer *kl, integer *ku,
    doubleCReal *cfrom, doubleCReal *cto, integer *m, integer *n,
    doubleCReal *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(dlasd0)(integer *n, integer *sqre, doubleCReal *d__,
    doubleCReal *e, doubleCReal *u, integer *ldu, doubleCReal *vt, integer *
    ldvt, integer *smlsiz, integer *iwork, doubleCReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dlasd1)(integer *nl, integer *nr, integer *sqre,
    doubleCReal *d__, doubleCReal *alpha, doubleCReal *beta, doubleCReal *u,
    integer *ldu, doubleCReal *vt, integer *ldvt, integer *idxq, integer *
    iwork, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasd2)(integer *nl, integer *nr, integer *sqre, integer
    *k, doubleCReal *d__, doubleCReal *z__, doubleCReal *alpha, doubleCReal *
    beta, doubleCReal *u, integer *ldu, doubleCReal *vt, integer *ldvt,
    doubleCReal *dsigma, doubleCReal *u2, integer *ldu2, doubleCReal *vt2,
    integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
    idxq, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(dlasd3)(integer *nl, integer *nr, integer *sqre, integer
    *k, doubleCReal *d__, doubleCReal *q, integer *ldq, doubleCReal *dsigma,
    doubleCReal *u, integer *ldu, doubleCReal *u2, integer *ldu2,
    doubleCReal *vt, integer *ldvt, doubleCReal *vt2, integer *ldvt2,
    integer *idxc, integer *ctot, doubleCReal *z__, integer *info);

/* Subroutine */ int F77NAME(dlasd4)(integer *n, integer *i__, doubleCReal *d__,
    doubleCReal *z__, doubleCReal *delta, doubleCReal *rho, doubleCReal *
    sigma, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasd5)(integer *i__, doubleCReal *d__, doubleCReal *z__,
    doubleCReal *delta, doubleCReal *rho, doubleCReal *dsigma, doubleCReal *
    work);

/* Subroutine */ int F77NAME(dlasd6)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, doubleCReal *d__, doubleCReal *vf, doubleCReal *vl,
    doubleCReal *alpha, doubleCReal *beta, integer *idxq, integer *perm,
    integer *givptr, integer *givcol, integer *ldgcol, doubleCReal *givnum,
     integer *ldgnum, doubleCReal *poles, doubleCReal *difl, doubleCReal *
    difr, doubleCReal *z__, integer *k, doubleCReal *c__, doubleCReal *s,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlasd7)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *k, doubleCReal *d__, doubleCReal *z__,
    doubleCReal *zw, doubleCReal *vf, doubleCReal *vfw, doubleCReal *vl,
    doubleCReal *vlw, doubleCReal *alpha, doubleCReal *beta, doubleCReal *
    dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm,
    integer *givptr, integer *givcol, integer *ldgcol, doubleCReal *givnum,
     integer *ldgnum, doubleCReal *c__, doubleCReal *s, integer *info);

/* Subroutine */ int F77NAME(dlasd8)(integer *icompq, integer *k, doubleCReal *d__,
    doubleCReal *z__, doubleCReal *vf, doubleCReal *vl, doubleCReal *difl,
    doubleCReal *difr, integer *lddifr, doubleCReal *dsigma, doubleCReal *
    work, integer *info);

/* Subroutine */ int F77NAME(dlasd9)(integer *icompq, integer *ldu, integer *k,
    doubleCReal *d__, doubleCReal *z__, doubleCReal *vf, doubleCReal *vl,
    doubleCReal *difl, doubleCReal *difr, doubleCReal *dsigma, doubleCReal *
    work, integer *info);

/* Subroutine */ int F77NAME(dlasda)(integer *icompq, integer *smlsiz, integer *n,
    integer *sqre, doubleCReal *d__, doubleCReal *e, doubleCReal *u, integer
    *ldu, doubleCReal *vt, integer *k, doubleCReal *difl, doubleCReal *difr,
    doubleCReal *z__, doubleCReal *poles, integer *givptr, integer *givcol,
    integer *ldgcol, integer *perm, doubleCReal *givnum, doubleCReal *c__,
    doubleCReal *s, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlasdq)(char *uplo, integer *sqre, integer *n, integer *
    ncvt, integer *nru, integer *ncc, doubleCReal *d__, doubleCReal *e,
    doubleCReal *vt, integer *ldvt, doubleCReal *u, integer *ldu,
    doubleCReal *c__, integer *ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasdt)(integer *n, integer *lvl, integer *nd, integer *
    inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int F77NAME(dlaset)(char *uplo, integer *m, integer *n, doubleCReal *
    alpha, doubleCReal *beta, doubleCReal *a, integer *lda);

/* Subroutine */ int F77NAME(dlasq1)(integer *n, doubleCReal *d__, doubleCReal *e,
    doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dlasq2)(integer *n, doubleCReal *z__, integer *info);

/* Subroutine */ int F77NAME(dlasq3)(integer *i0, integer *n0, doubleCReal *z__,
    integer *pp, doubleCReal *dmin__, doubleCReal *sigma, doubleCReal *desig,
     doubleCReal *qmax, integer *nfail, integer *iter, integer *ndiv,
    logical *ieee);

/* Subroutine */ int F77NAME(dlasq4)(integer *i0, integer *n0, doubleCReal *z__,
    integer *pp, integer *n0in, doubleCReal *dmin__, doubleCReal *dmin1,
    doubleCReal *dmin2, doubleCReal *dn, doubleCReal *dn1, doubleCReal *dn2,
    doubleCReal *tau, integer *ttype);

/* Subroutine */ int F77NAME(dlasq5)(integer *i0, integer *n0, doubleCReal *z__,
    integer *pp, doubleCReal *tau, doubleCReal *dmin__, doubleCReal *dmin1,
    doubleCReal *dmin2, doubleCReal *dn, doubleCReal *dnm1, doubleCReal *dnm2,
     logical *ieee);

/* Subroutine */ int F77NAME(dlasq6)(integer *i0, integer *n0, doubleCReal *z__,
    integer *pp, doubleCReal *dmin__, doubleCReal *dmin1, doubleCReal *dmin2,
     doubleCReal *dn, doubleCReal *dnm1, doubleCReal *dnm2);

/* Subroutine */ int F77NAME(dlasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, doubleCReal *c__, doubleCReal *s, doubleCReal *a, integer *
    lda);

/* Subroutine */ int F77NAME(dlasrt)(char *id, integer *n, doubleCReal *d__, integer *
    info);

/* Subroutine */ int F77NAME(dlassq)(integer *n, doubleCReal *x, integer *incx,
    doubleCReal *scale, doubleCReal *sumsq);

/* Subroutine */ int F77NAME(dlasv2)(doubleCReal *f, doubleCReal *g, doubleCReal *h__,
    doubleCReal *ssmin, doubleCReal *ssmax, doubleCReal *snr, doubleCReal *
    csr, doubleCReal *snl, doubleCReal *csl);

/* Subroutine */ int F77NAME(dlaswp)(integer *n, doubleCReal *a, integer *lda, integer
    *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(dlasy2)(logical *ltranl, logical *ltranr, integer *isgn,
    integer *n1, integer *n2, doubleCReal *tl, integer *ldtl, doubleCReal *
    tr, integer *ldtr, doubleCReal *b, integer *ldb, doubleCReal *scale,
    doubleCReal *x, integer *ldx, doubleCReal *xnorm, integer *info);

/* Subroutine */ int F77NAME(dlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     doubleCReal *a, integer *lda, integer *ipiv, doubleCReal *w, integer *
    ldw, integer *info);

/* Subroutine */ int F77NAME(dlatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, doubleCReal *ab, integer *ldab,
    doubleCReal *x, doubleCReal *scale, doubleCReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatdf)(integer *ijob, integer *n, doubleCReal *z__,
    integer *ldz, doubleCReal *rhs, doubleCReal *rdsum, doubleCReal *rdscal,
    integer *ipiv, integer *jpiv);

/* Subroutine */ int F77NAME(dlatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleCReal *ap, doubleCReal *x, doubleCReal *scale,
    doubleCReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatrd)(char *uplo, integer *n, integer *nb, doubleCReal *
    a, integer *lda, doubleCReal *e, doubleCReal *tau, doubleCReal *w,
    integer *ldw);

/* Subroutine */ int F77NAME(dlatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleCReal *a, integer *lda, doubleCReal *x,
    doubleCReal *scale, doubleCReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatrz)(integer *m, integer *n, integer *l, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work);

/* Subroutine */ int F77NAME(dlatzm)(char *side, integer *m, integer *n, doubleCReal *
    v, integer *incv, doubleCReal *tau, doubleCReal *c1, doubleCReal *c2,
    integer *ldc, doubleCReal *work);

/* Subroutine */ int F77NAME(dlauu2)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dlauum)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dopgtr)(char *uplo, integer *n, doubleCReal *ap,
    doubleCReal *tau, doubleCReal *q, integer *ldq, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dopmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doubleCReal *ap, doubleCReal *tau, doubleCReal *c__, integer
    *ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dorg2l)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dorg2r)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dorgbr)(char *vect, integer *m, integer *n, integer *k,
    doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorghr)(integer *n, integer *ilo, integer *ihi,
    doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorgl2)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dorglq)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgql)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgqr)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgr2)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dorgrq)(integer *m, integer *n, integer *k, doubleCReal *
    a, integer *lda, doubleCReal *tau, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgtr)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dorm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dormbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, doubleCReal *a, integer *lda, doubleCReal *tau,
    doubleCReal *c__, integer *ldc, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dormhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, doubleCReal *a, integer *lda, doubleCReal *
    tau, doubleCReal *c__, integer *ldc, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dormlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormql)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dormr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doubleCReal *a, integer *lda, doubleCReal *tau,
    doubleCReal *c__, integer *ldc, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dormrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doubleCReal *a, integer *lda, doubleCReal *tau,
    doubleCReal *c__, integer *ldc, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dormtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doubleCReal *a, integer *lda, doubleCReal *tau, doubleCReal *
    c__, integer *ldc, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dpbcon)(char *uplo, integer *n, integer *kd, doubleCReal *
    ab, integer *ldab, doubleCReal *anorm, doubleCReal *rcond, doubleCReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpbequ)(char *uplo, integer *n, integer *kd, doubleCReal *
    ab, integer *ldab, doubleCReal *s, doubleCReal *scond, doubleCReal *amax,
     integer *info);

/* Subroutine */ int F77NAME(dpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleCReal *ab, integer *ldab, doubleCReal *afb, integer *ldafb,
    doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *
    ferr, doubleCReal *berr, doubleCReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dpbstf)(char *uplo, integer *n, integer *kd, doubleCReal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleCReal *ab, integer *ldab, doubleCReal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(dpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, doubleCReal *ab, integer *ldab, doubleCReal *afb,
    integer *ldafb, char *equed, doubleCReal *s, doubleCReal *b, integer *
    ldb, doubleCReal *x, integer *ldx, doubleCReal *rcond, doubleCReal *ferr,
     doubleCReal *berr, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpbtf2)(char *uplo, integer *n, integer *kd, doubleCReal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbtrf)(char *uplo, integer *n, integer *kd, doubleCReal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleCReal *ab, integer *ldab, doubleCReal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(dpocon)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *anorm, doubleCReal *rcond, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dpoequ)(integer *n, doubleCReal *a, integer *lda,
    doubleCReal *s, doubleCReal *scond, doubleCReal *amax, integer *info);

/* Subroutine */ int F77NAME(dporfs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *af, integer *ldaf,
    doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *
    ferr, doubleCReal *berr, doubleCReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dposv)(char *uplo, integer *n, integer *nrhs, doubleCReal
    *a, integer *lda, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleCReal *a, integer *lda, doubleCReal *af, integer *ldaf,
    char *equed, doubleCReal *s, doubleCReal *b, integer *ldb, doubleCReal *
    x, integer *ldx, doubleCReal *rcond, doubleCReal *ferr, doubleCReal *
    berr, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpotf2)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotrf)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotri)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotrs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dppcon)(char *uplo, integer *n, doubleCReal *ap,
    doubleCReal *anorm, doubleCReal *rcond, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dppequ)(char *uplo, integer *n, doubleCReal *ap,
    doubleCReal *s, doubleCReal *scond, doubleCReal *amax, integer *info);

/* Subroutine */ int F77NAME(dpprfs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *ap, doubleCReal *afp, doubleCReal *b, integer *ldb,
    doubleCReal *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dppsv)(char *uplo, integer *n, integer *nrhs, doubleCReal
    *ap, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleCReal *ap, doubleCReal *afp, char *equed, doubleCReal *s,
    doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *
    rcond, doubleCReal *ferr, doubleCReal *berr, doubleCReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dpptrf)(char *uplo, integer *n, doubleCReal *ap, integer *
    info);

/* Subroutine */ int F77NAME(dpptri)(char *uplo, integer *n, doubleCReal *ap, integer *
    info);

/* Subroutine */ int F77NAME(dpptrs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *ap, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptcon)(integer *n, doubleCReal *d__, doubleCReal *e,
    doubleCReal *anorm, doubleCReal *rcond, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dpteqr)(char *compz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dptrfs)(integer *n, integer *nrhs, doubleCReal *d__,
    doubleCReal *e, doubleCReal *df, doubleCReal *ef, doubleCReal *b, integer
    *ldb, doubleCReal *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
     doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dptsv)(integer *n, integer *nrhs, doubleCReal *d__,
    doubleCReal *e, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptsvx)(char *fact, integer *n, integer *nrhs,
    doubleCReal *d__, doubleCReal *e, doubleCReal *df, doubleCReal *ef,
    doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *
    rcond, doubleCReal *ferr, doubleCReal *berr, doubleCReal *work, integer *
    info);

/* Subroutine */ int F77NAME(dpttrf)(integer *n, doubleCReal *d__, doubleCReal *e,
    integer *info);

/* Subroutine */ int F77NAME(dpttrs)(integer *n, integer *nrhs, doubleCReal *d__,
    doubleCReal *e, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptts2)(integer *n, integer *nrhs, doubleCReal *d__,
    doubleCReal *e, doubleCReal *b, integer *ldb);

/* Subroutine */ int F77NAME(drscl)(integer *n, doubleCReal *sa, doubleCReal *sx,
    integer *incx);

/* Subroutine */ int F77NAME(dsbev)(char *jobz, char *uplo, integer *n, integer *kd,
    doubleCReal *ab, integer *ldab, doubleCReal *w, doubleCReal *z__,
    integer *ldz, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dsbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    doubleCReal *ab, integer *ldab, doubleCReal *w, doubleCReal *z__,
    integer *ldz, doubleCReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, doubleCReal *ab, integer *ldab, doubleCReal *q, integer *
    ldq, doubleCReal *vl, doubleCReal *vu, integer *il, integer *iu,
    doubleCReal *abstol, integer *m, doubleCReal *w, doubleCReal *z__,
    integer *ldz, doubleCReal *work, integer *iwork, integer *ifail,
    integer *info);

/* Subroutine */ int F77NAME(dsbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, doubleCReal *ab, integer *ldab, doubleCReal *bb, integer *
    ldbb, doubleCReal *x, integer *ldx, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dsbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doubleCReal *ab, integer *ldab, doubleCReal *bb, integer *
    ldbb, doubleCReal *w, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dsbgvd)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doubleCReal *ab, integer *ldab, doubleCReal *bb, integer *
    ldbb, doubleCReal *w, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, doubleCReal *ab, integer *ldab, doubleCReal *
    bb, integer *ldbb, doubleCReal *q, integer *ldq, doubleCReal *vl,
    doubleCReal *vu, integer *il, integer *iu, doubleCReal *abstol, integer
    *m, doubleCReal *w, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    doubleCReal *ab, integer *ldab, doubleCReal *d__, doubleCReal *e,
    doubleCReal *q, integer *ldq, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dspcon)(char *uplo, integer *n, doubleCReal *ap, integer *
    ipiv, doubleCReal *anorm, doubleCReal *rcond, doubleCReal *work, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(dspev)(char *jobz, char *uplo, integer *n, doubleCReal *
    ap, doubleCReal *w, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dspevd)(char *jobz, char *uplo, integer *n, doubleCReal *
    ap, doubleCReal *w, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dspevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleCReal *ap, doubleCReal *vl, doubleCReal *vu, integer *il, integer *
    iu, doubleCReal *abstol, integer *m, doubleCReal *w, doubleCReal *z__,
    integer *ldz, doubleCReal *work, integer *iwork, integer *ifail,
    integer *info);

/* Subroutine */ int F77NAME(dspgst)(integer *itype, char *uplo, integer *n,
    doubleCReal *ap, doubleCReal *bp, integer *info);

/* Subroutine */ int F77NAME(dspgv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleCReal *ap, doubleCReal *bp, doubleCReal *w, doubleCReal *z__,
    integer *ldz, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dspgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleCReal *ap, doubleCReal *bp, doubleCReal *w, doubleCReal *z__,
    integer *ldz, doubleCReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dspgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleCReal *ap, doubleCReal *bp, doubleCReal *vl,
    doubleCReal *vu, integer *il, integer *iu, doubleCReal *abstol, integer
    *m, doubleCReal *w, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsprfs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *ap, doubleCReal *afp, integer *ipiv, doubleCReal *b,
    integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *ferr,
    doubleCReal *berr, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dspsv)(char *uplo, integer *n, integer *nrhs, doubleCReal
    *ap, integer *ipiv, doubleCReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleCReal *ap, doubleCReal *afp, integer *ipiv, doubleCReal *b,
    integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *rcond,
    doubleCReal *ferr, doubleCReal *berr, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dsptrd)(char *uplo, integer *n, doubleCReal *ap,
    doubleCReal *d__, doubleCReal *e, doubleCReal *tau, integer *info);

/* Subroutine */ int F77NAME(dsptrf)(char *uplo, integer *n, doubleCReal *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(dsptri)(char *uplo, integer *n, doubleCReal *ap, integer *
    ipiv, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dsptrs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *ap, integer *ipiv, doubleCReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dstebz)(char *range, char *order, integer *n, doubleCReal
    *vl, doubleCReal *vu, integer *il, integer *iu, doubleCReal *abstol,
    doubleCReal *d__, doubleCReal *e, integer *m, integer *nsplit,
    doubleCReal *w, integer *iblock, integer *isplit, doubleCReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dstedc)(char *compz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstegr)(char *jobz, char *range, integer *n, doubleCReal *
    d__, doubleCReal *e, doubleCReal *vl, doubleCReal *vu, integer *il,
    integer *iu, doubleCReal *abstol, integer *m, doubleCReal *w,
    doubleCReal *z__, integer *ldz, integer *isuppz, doubleCReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstein)(integer *n, doubleCReal *d__, doubleCReal *e,
    integer *m, doubleCReal *w, integer *iblock, integer *isplit,
    doubleCReal *z__, integer *ldz, doubleCReal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsteqr)(char *compz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dsterf)(integer *n, doubleCReal *d__, doubleCReal *e,
    integer *info);

/* Subroutine */ int F77NAME(dstev)(char *jobz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(dstevd)(char *jobz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleCReal *z__, integer *ldz, doubleCReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstevr)(char *jobz, char *range, integer *n, doubleCReal *
    d__, doubleCReal *e, doubleCReal *vl, doubleCReal *vu, integer *il,
    integer *iu, doubleCReal *abstol, integer *m, doubleCReal *w,
    doubleCReal *z__, integer *ldz, integer *isuppz, doubleCReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstevx)(char *jobz, char *range, integer *n, doubleCReal *
    d__, doubleCReal *e, doubleCReal *vl, doubleCReal *vu, integer *il,
    integer *iu, doubleCReal *abstol, integer *m, doubleCReal *w,
    doubleCReal *z__, integer *ldz, doubleCReal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsycon)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *ipiv, doubleCReal *anorm, doubleCReal *rcond, doubleCReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dsyev)(char *jobz, char *uplo, integer *n, doubleCReal *a,
     integer *lda, doubleCReal *w, doubleCReal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dsyevd)(char *jobz, char *uplo, integer *n, doubleCReal *
    a, integer *lda, doubleCReal *w, doubleCReal *work, integer *lwork,
    integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsyevr)(char *jobz, char *range, char *uplo, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *vl, doubleCReal *vu, integer *
    il, integer *iu, doubleCReal *abstol, integer *m, doubleCReal *w,
    doubleCReal *z__, integer *ldz, integer *isuppz, doubleCReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsyevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *vl, doubleCReal *vu, integer *
    il, integer *iu, doubleCReal *abstol, integer *m, doubleCReal *w,
    doubleCReal *z__, integer *ldz, doubleCReal *work, integer *lwork,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsygs2)(integer *itype, char *uplo, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dsygst)(integer *itype, char *uplo, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dsygv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb,
    doubleCReal *w, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsygvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb,
    doubleCReal *w, doubleCReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsygvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleCReal *a, integer *lda, doubleCReal *b, integer
    *ldb, doubleCReal *vl, doubleCReal *vu, integer *il, integer *iu,
    doubleCReal *abstol, integer *m, doubleCReal *w, doubleCReal *z__,
    integer *ldz, doubleCReal *work, integer *lwork, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsyrfs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, doubleCReal *af, integer *ldaf, integer *
    ipiv, doubleCReal *b, integer *ldb, doubleCReal *x, integer *ldx,
    doubleCReal *ferr, doubleCReal *berr, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dsysv)(char *uplo, integer *n, integer *nrhs, doubleCReal
    *a, integer *lda, integer *ipiv, doubleCReal *b, integer *ldb,
    doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleCReal *a, integer *lda, doubleCReal *af, integer *ldaf,
    integer *ipiv, doubleCReal *b, integer *ldb, doubleCReal *x, integer *
    ldx, doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr,
    doubleCReal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dsytd2)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *d__, doubleCReal *e, doubleCReal *tau, integer *info);

/* Subroutine */ int F77NAME(dsytf2)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dsytrd)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *d__, doubleCReal *e, doubleCReal *tau, doubleCReal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsytrf)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *ipiv, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsytri)(char *uplo, integer *n, doubleCReal *a, integer *
    lda, integer *ipiv, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dsytrs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *a, integer *lda, integer *ipiv, doubleCReal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dtbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, doubleCReal *ab, integer *ldab, doubleCReal *rcond,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleCReal *ab, integer *ldab, doubleCReal
    *b, integer *ldb, doubleCReal *x, integer *ldx, doubleCReal *ferr,
    doubleCReal *berr, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleCReal *ab, integer *ldab, doubleCReal
    *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dtgevc)(char *side, char *howmny, logical *select,
    integer *n, doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb,
    doubleCReal *vl, integer *ldvl, doubleCReal *vr, integer *ldvr, integer
    *mm, integer *m, doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dtgex2)(logical *wantq, logical *wantz, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *
    q, integer *ldq, doubleCReal *z__, integer *ldz, integer *j1, integer *
    n1, integer *n2, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dtgexc)(logical *wantq, logical *wantz, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb, doubleCReal *
    q, integer *ldq, doubleCReal *z__, integer *ldz, integer *ifst,
    integer *ilst, doubleCReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dtgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, doubleCReal *a, integer *lda, doubleCReal *
    b, integer *ldb, doubleCReal *alphar, doubleCReal *alphai, doubleCReal *
    beta, doubleCReal *q, integer *ldq, doubleCReal *z__, integer *ldz,
    integer *m, doubleCReal *pl, doubleCReal *pr, doubleCReal *dif,
    doubleCReal *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(dtgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, doubleCReal *a,
    integer *lda, doubleCReal *b, integer *ldb, doubleCReal *tola,
    doubleCReal *tolb, doubleCReal *alpha, doubleCReal *beta, doubleCReal *u,
    integer *ldu, doubleCReal *v, integer *ldv, doubleCReal *q, integer *
    ldq, doubleCReal *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(dtgsna)(char *job, char *howmny, logical *select,
    integer *n, doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb,
    doubleCReal *vl, integer *ldvl, doubleCReal *vr, integer *ldvr,
    doubleCReal *s, doubleCReal *dif, integer *mm, integer *m, doubleCReal *
    work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb,
    doubleCReal *c__, integer *ldc, doubleCReal *d__, integer *ldd,
    doubleCReal *e, integer *lde, doubleCReal *f, integer *ldf, doubleCReal *
    scale, doubleCReal *rdsum, doubleCReal *rdscal, integer *iwork, integer
    *pq, integer *info);

/* Subroutine */ int F77NAME(dtgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, doubleCReal *a, integer *lda, doubleCReal *b, integer *ldb,
    doubleCReal *c__, integer *ldc, doubleCReal *d__, integer *ldd,
    doubleCReal *e, integer *lde, doubleCReal *f, integer *ldf, doubleCReal *
    scale, doubleCReal *dif, doubleCReal *work, integer *lwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dtpcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleCReal *ap, doubleCReal *rcond, doubleCReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dtprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleCReal *ap, doubleCReal *b, integer *ldb,
    doubleCReal *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtptri)(char *uplo, char *diag, integer *n, doubleCReal *
    ap, integer *info);

/* Subroutine */ int F77NAME(dtptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleCReal *ap, doubleCReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dtrcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleCReal *a, integer *lda, doubleCReal *rcond, doubleCReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtrevc)(char *side, char *howmny, logical *select,
    integer *n, doubleCReal *t, integer *ldt, doubleCReal *vl, integer *
    ldvl, doubleCReal *vr, integer *ldvr, integer *mm, integer *m,
    doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dtrexc)(char *compq, integer *n, doubleCReal *t, integer *
    ldt, doubleCReal *q, integer *ldq, integer *ifst, integer *ilst,
    doubleCReal *work, integer *info);

/* Subroutine */ int F77NAME(dtrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleCReal *a, integer *lda, doubleCReal *b, integer *
    ldb, doubleCReal *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtrsen)(char *job, char *compq, logical *select, integer
    *n, doubleCReal *t, integer *ldt, doubleCReal *q, integer *ldq,
    doubleCReal *wr, doubleCReal *wi, integer *m, doubleCReal *s, doubleCReal
    *sep, doubleCReal *work, integer *lwork, integer *iwork, integer *
    liwork, integer *info);

/* Subroutine */ int F77NAME(dtrsna)(char *job, char *howmny, logical *select,
    integer *n, doubleCReal *t, integer *ldt, doubleCReal *vl, integer *
    ldvl, doubleCReal *vr, integer *ldvr, doubleCReal *s, doubleCReal *sep,
    integer *mm, integer *m, doubleCReal *work, integer *ldwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dtrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, doubleCReal *a, integer *lda, doubleCReal *b, integer *
    ldb, doubleCReal *c__, integer *ldc, doubleCReal *scale, integer *info);

/* Subroutine */ int F77NAME(dtrti2)(char *uplo, char *diag, integer *n, doubleCReal *
    a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(dtrtri)(char *uplo, char *diag, integer *n, doubleCReal *
    a, integer *lda, integer *info);

/**
  Lapack routine : DTRTRS solves a triangular system of the form
  A * X = B  or  A**T * X = B,
  where A is a triangular matrix of order N, and B is an N-by-NRHS
  matrix.  A check is made to verify that A is nonsingular.
*/
int F77NAME(dtrtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleCReal *a, integer *lda, doubleCReal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dtzrqf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, integer *info);

/* Subroutine */ int F77NAME(dtzrzf)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleCReal *tau, doubleCReal *work, integer *lwork, integer *info);

integer F77NAME(icmax1)(integer *n, Complex *cx, integer *incx);

integer F77NAME(ieeeck)(integer *ispec, CReal *zero, CReal *one);

integer F77NAME(ilaenv)(integer *ispec, char *name__, char *opts, integer *n1,
    integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen
    opts_len);

integer F77NAME(izmax1)(integer *n, doubleComplex *cx, integer *incx);

/* Subroutine */ int F77NAME(sbdsdc)(char *uplo, char *compq, integer *n, CReal *d__,
    CReal *e, CReal *u, integer *ldu, CReal *vt, integer *ldvt, CReal *q,
    integer *iq, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, CReal *d__, CReal *e, CReal *vt, integer *ldvt, CReal *
    u, integer *ldu, CReal *c__, integer *ldc, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sdisna)(char *job, integer *m, integer *n, CReal *d__,
    CReal *sep, integer *info);

/* Subroutine */ int F77NAME(sgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, CReal *ab, integer *ldab, CReal *d__, CReal *
    e, CReal *q, integer *ldq, CReal *pt, integer *ldpt, CReal *c__, integer
    *ldc, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     CReal *ab, integer *ldab, integer *ipiv, CReal *anorm, CReal *rcond,
    CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     CReal *ab, integer *ldab, CReal *r__, CReal *c__, CReal *rowcnd, CReal *
    colcnd, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(sgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, CReal *ab, integer *ldab, CReal *afb, integer *ldafb,
     integer *ipiv, CReal *b, integer *ldb, CReal *x, integer *ldx, CReal *
    ferr, CReal *berr, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, CReal *ab, integer *ldab, integer *ipiv, CReal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(sgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, CReal *ab, integer *ldab, CReal *afb,
    integer *ldafb, integer *ipiv, char *equed, CReal *r__, CReal *c__,
    CReal *b, integer *ldb, CReal *x, integer *ldx, CReal *rcond, CReal *ferr,
     CReal *berr, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     CReal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     CReal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, CReal *ab, integer *ldab, integer *ipiv, CReal *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, CReal *scale, integer *m, CReal *v, integer *ldv, integer
    *info);

/* Subroutine */ int F77NAME(sgebal)(char *job, integer *n, CReal *a, integer *lda,
    integer *ilo, integer *ihi, CReal *scale, integer *info);

/* Subroutine */ int F77NAME(sgebd2)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *d__, CReal *e, CReal *tauq, CReal *taup, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgebrd)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *d__, CReal *e, CReal *tauq, CReal *taup, CReal *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(sgecon)(char *norm, integer *n, CReal *a, integer *lda,
    CReal *anorm, CReal *rcond, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgeequ)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *r__, CReal *c__, CReal *rowcnd, CReal *colcnd, CReal *amax, integer
    *info);

/* Subroutine */ int F77NAME(sgees)(char *jobvs, char *sort, L_fp select, integer *n,
    CReal *a, integer *lda, integer *sdim, CReal *wr, CReal *wi, CReal *vs,
    integer *ldvs, CReal *work, integer *lwork, logical *bwork, integer *
    info);

/* Subroutine */ int F77NAME(sgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, CReal *a, integer *lda, integer *sdim, CReal *wr,
    CReal *wi, CReal *vs, integer *ldvs, CReal *rconde, CReal *rcondv, CReal *
    work, integer *lwork, integer *iwork, integer *liwork, logical *bwork,
     integer *info);

/* Subroutine */ int F77NAME(sgeev)(char *jobvl, char *jobvr, integer *n, CReal *a,
    integer *lda, CReal *wr, CReal *wi, CReal *vl, integer *ldvl, CReal *vr,
    integer *ldvr, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, CReal *a, integer *lda, CReal *wr, CReal *wi, CReal *
    vl, integer *ldvl, CReal *vr, integer *ldvr, integer *ilo, integer *
    ihi, CReal *scale, CReal *abnrm, CReal *rconde, CReal *rcondv, CReal *work,
     integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgegs)(char *jobvsl, char *jobvsr, integer *n, CReal *a,
    integer *lda, CReal *b, integer *ldb, CReal *alphar, CReal *alphai, CReal
    *beta, CReal *vsl, integer *ldvsl, CReal *vsr, integer *ldvsr, CReal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgegv)(char *jobvl, char *jobvr, integer *n, CReal *a,
    integer *lda, CReal *b, integer *ldb, CReal *alphar, CReal *alphai, CReal
    *beta, CReal *vl, integer *ldvl, CReal *vr, integer *ldvr, CReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgehd2)(integer *n, integer *ilo, integer *ihi, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgehrd)(integer *n, integer *ilo, integer *ihi, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelq2)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgelqf)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgels)(char *trans, integer *m, integer *n, integer *
    nrhs, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelsd)(integer *m, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *b, integer *ldb, CReal *s, CReal *rcond, integer *
    rank, CReal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgelss)(integer *m, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *b, integer *ldb, CReal *s, CReal *rcond, integer *
    rank, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelsx)(integer *m, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *b, integer *ldb, integer *jpvt, CReal *rcond,
    integer *rank, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgelsy)(integer *m, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *b, integer *ldb, integer *jpvt, CReal *rcond,
    integer *rank, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeql2)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgeqlf)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeqp3)(integer *m, integer *n, CReal *a, integer *lda,
    integer *jpvt, CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeqpf)(integer *m, integer *n, CReal *a, integer *lda,
    integer *jpvt, CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgeqr2)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgeqrf)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgerfs)(char *trans, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *af, integer *ldaf, integer *ipiv, CReal *b,
    integer *ldb, CReal *x, integer *ldx, CReal *ferr, CReal *berr, CReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgerq2)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgerqf)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgesc2)(integer *n, CReal *a, integer *lda, CReal *rhs,
    integer *ipiv, integer *jpiv, CReal *scale);

/* Subroutine */ int F77NAME(sgesdd)(char *jobz, integer *m, integer *n, CReal *a,
    integer *lda, CReal *s, CReal *u, integer *ldu, CReal *vt, integer *ldvt,
     CReal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgesv)(integer *n, integer *nrhs, CReal *a, integer *lda,
    integer *ipiv, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgesvd)(char *jobu, char *jobvt, integer *m, integer *n,
    CReal *a, integer *lda, CReal *s, CReal *u, integer *ldu, CReal *vt,
    integer *ldvt, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, CReal *a, integer *lda, CReal *af, integer *ldaf, integer *ipiv,
    char *equed, CReal *r__, CReal *c__, CReal *b, integer *ldb, CReal *x,
    integer *ldx, CReal *rcond, CReal *ferr, CReal *berr, CReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgetc2)(integer *n, CReal *a, integer *lda, integer *ipiv,
     integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(sgetf2)(integer *m, integer *n, CReal *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgetrf)(integer *m, integer *n, CReal *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgetri)(integer *n, CReal *a, integer *lda, integer *ipiv,
     CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgetrs)(char *trans, integer *n, integer *nrhs, CReal *a,
    integer *lda, integer *ipiv, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, CReal *lscale, CReal *rscale, integer *m, CReal *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(sggbal)(char *job, integer *n, CReal *a, integer *lda,
    CReal *b, integer *ldb, integer *ilo, integer *ihi, CReal *lscale, CReal
    *rscale, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, integer *n, CReal *a, integer *lda, CReal *b, integer *ldb,
    integer *sdim, CReal *alphar, CReal *alphai, CReal *beta, CReal *vsl,
    integer *ldvsl, CReal *vsr, integer *ldvsr, CReal *work, integer *lwork,
     logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, char *sense, integer *n, CReal *a, integer *lda, CReal *b,
    integer *ldb, integer *sdim, CReal *alphar, CReal *alphai, CReal *beta,
    CReal *vsl, integer *ldvsl, CReal *vsr, integer *ldvsr, CReal *rconde,
    CReal *rcondv, CReal *work, integer *lwork, integer *iwork, integer *
    liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggev)(char *jobvl, char *jobvr, integer *n, CReal *a,
    integer *lda, CReal *b, integer *ldb, CReal *alphar, CReal *alphai, CReal
    *beta, CReal *vl, integer *ldvl, CReal *vr, integer *ldvr, CReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal
    *alphar, CReal *alphai, CReal *beta, CReal *vl, integer *ldvl, CReal *vr,
    integer *ldvr, integer *ilo, integer *ihi, CReal *lscale, CReal *rscale,
     CReal *abnrm, CReal *bbnrm, CReal *rconde, CReal *rcondv, CReal *work,
    integer *lwork, integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggglm)(integer *n, integer *m, integer *p, CReal *a,
    integer *lda, CReal *b, integer *ldb, CReal *d__, CReal *x, CReal *y,
    CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, CReal *a, integer *lda, CReal *b, integer *ldb, CReal
    *q, integer *ldq, CReal *z__, integer *ldz, integer *info);

/* Subroutine */ int F77NAME(sgglse)(integer *m, integer *n, integer *p, CReal *a,
    integer *lda, CReal *b, integer *ldb, CReal *c__, CReal *d__, CReal *x,
    CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggqrf)(integer *n, integer *m, integer *p, CReal *a,
    integer *lda, CReal *taua, CReal *b, integer *ldb, CReal *taub, CReal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggrqf)(integer *m, integer *p, integer *n, CReal *a,
    integer *lda, CReal *taua, CReal *b, integer *ldb, CReal *taub, CReal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, CReal *a, integer *lda,
     CReal *b, integer *ldb, CReal *alpha, CReal *beta, CReal *u, integer *
    ldu, CReal *v, integer *ldv, CReal *q, integer *ldq, CReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, CReal *a, integer *lda, CReal *b, integer *ldb,
    CReal *tola, CReal *tolb, integer *k, integer *l, CReal *u, integer *ldu,
     CReal *v, integer *ldv, CReal *q, integer *ldq, integer *iwork, CReal *
    tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sgtcon)(char *norm, integer *n, CReal *dl, CReal *d__,
    CReal *du, CReal *du2, integer *ipiv, CReal *anorm, CReal *rcond, CReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgtrfs)(char *trans, integer *n, integer *nrhs, CReal *dl,
     CReal *d__, CReal *du, CReal *dlf, CReal *df, CReal *duf, CReal *du2,
    integer *ipiv, CReal *b, integer *ldb, CReal *x, integer *ldx, CReal *
    ferr, CReal *berr, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgtsv)(integer *n, integer *nrhs, CReal *dl, CReal *d__,
    CReal *du, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, CReal *dl, CReal *d__, CReal *du, CReal *dlf, CReal *df, CReal *duf,
    CReal *du2, integer *ipiv, CReal *b, integer *ldb, CReal *x, integer *
    ldx, CReal *rcond, CReal *ferr, CReal *berr, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(sgttrf)(integer *n, CReal *dl, CReal *d__, CReal *du, CReal *
    du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgttrs)(char *trans, integer *n, integer *nrhs, CReal *dl,
     CReal *d__, CReal *du, CReal *du2, integer *ipiv, CReal *b, integer *ldb,
     integer *info);

/* Subroutine */ int F77NAME(sgtts2)(integer *itrans, integer *n, integer *nrhs, CReal
    *dl, CReal *d__, CReal *du, CReal *du2, integer *ipiv, CReal *b, integer *
    ldb);

/* Subroutine */ int F77NAME(shgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, CReal *a, integer *lda, CReal *b, integer *
    ldb, CReal *alphar, CReal *alphai, CReal *beta, CReal *q, integer *ldq,
    CReal *z__, integer *ldz, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(shsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, CReal *h__, integer *ldh, CReal *wr, CReal *wi, CReal
    *vl, integer *ldvl, CReal *vr, integer *ldvr, integer *mm, integer *m,
    CReal *work, integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(shseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, CReal *h__, integer *ldh, CReal *wr, CReal *wi, CReal *z__,
     integer *ldz, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(slabad)(CReal *small, CReal *large);

/* Subroutine */ int F77NAME(slabrd)(integer *m, integer *n, integer *nb, CReal *a,
    integer *lda, CReal *d__, CReal *e, CReal *tauq, CReal *taup, CReal *x,
    integer *ldx, CReal *y, integer *ldy);

/* Subroutine */ int F77NAME(slacon)(integer *n, CReal *v, CReal *x, integer *isgn,
    CReal *est, integer *kase);

/* Subroutine */ int F77NAME(slacpy)(char *uplo, integer *m, integer *n, CReal *a,
    integer *lda, CReal *b, integer *ldb);

/* Subroutine */ int F77NAME(sladiv)(CReal *a, CReal *b, CReal *c__, CReal *d__, CReal *p,
    CReal *q);

/* Subroutine */ int F77NAME(slae2)(CReal *a, CReal *b, CReal *c__, CReal *rt1, CReal *rt2);

/* Subroutine */ int F77NAME(slaebz)(integer *ijob, integer *nitmax, integer *n,
    integer *mmax, integer *minp, integer *nbmin, CReal *abstol, CReal *
    reltol, CReal *pivmin, CReal *d__, CReal *e, CReal *e2, integer *nval,
    CReal *ab, CReal *c__, integer *mout, integer *nab, CReal *work, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(slaed0)(integer *icompq, integer *qsiz, integer *n, CReal
    *d__, CReal *e, CReal *q, integer *ldq, CReal *qstore, integer *ldqs,
    CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slaed1)(integer *n, CReal *d__, CReal *q, integer *ldq,
    integer *indxq, CReal *rho, integer *cutpnt, CReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(slaed2)(integer *k, integer *n, integer *n1, CReal *d__,
    CReal *q, integer *ldq, integer *indxq, CReal *rho, CReal *z__, CReal *
    dlamda, CReal *w, CReal *q2, integer *indx, integer *indxc, integer *
    indxp, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(slaed3)(integer *k, integer *n, integer *n1, CReal *d__,
    CReal *q, integer *ldq, CReal *rho, CReal *dlamda, CReal *q2, integer *
    indx, integer *ctot, CReal *w, CReal *s, integer *info);

/* Subroutine */ int F77NAME(slaed4)(integer *n, integer *i__, CReal *d__, CReal *z__,
    CReal *delta, CReal *rho, CReal *dlam, integer *info);

/* Subroutine */ int F77NAME(slaed5)(integer *i__, CReal *d__, CReal *z__, CReal *delta,
    CReal *rho, CReal *dlam);

/* Subroutine */ int F77NAME(slaed6)(integer *kniter, logical *orgati, CReal *rho,
    CReal *d__, CReal *z__, CReal *finit, CReal *tau, integer *info);

/* Subroutine */ int F77NAME(slaed7)(integer *icompq, integer *n, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, CReal *d__, CReal *q,
    integer *ldq, integer *indxq, CReal *rho, integer *cutpnt, CReal *
    qstore, integer *qptr, integer *prmptr, integer *perm, integer *
    givptr, integer *givcol, CReal *givnum, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(slaed8)(integer *icompq, integer *k, integer *n, integer
    *qsiz, CReal *d__, CReal *q, integer *ldq, integer *indxq, CReal *rho,
    integer *cutpnt, CReal *z__, CReal *dlamda, CReal *q2, integer *ldq2,
    CReal *w, integer *perm, integer *givptr, integer *givcol, CReal *
    givnum, integer *indxp, integer *indx, integer *info);

/* Subroutine */ int F77NAME(slaed9)(integer *k, integer *kstart, integer *kstop,
    integer *n, CReal *d__, CReal *q, integer *ldq, CReal *rho, CReal *dlamda,
     CReal *w, CReal *s, integer *lds, integer *info);

/* Subroutine */ int F77NAME(slaeda)(integer *n, integer *tlvls, integer *curlvl,
    integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
    integer *givcol, CReal *givnum, CReal *q, integer *qptr, CReal *z__,
    CReal *ztemp, integer *info);

/* Subroutine */ int F77NAME(slaein)(logical *rightv, logical *noinit, integer *n,
    CReal *h__, integer *ldh, CReal *wr, CReal *wi, CReal *vr, CReal *vi, CReal
    *b, integer *ldb, CReal *work, CReal *eps3, CReal *smlnum, CReal *bignum,
    integer *info);

/* Subroutine */ int F77NAME(slaev2)(CReal *a, CReal *b, CReal *c__, CReal *rt1, CReal *
    rt2, CReal *cs1, CReal *sn1);

/* Subroutine */ int F77NAME(slaexc)(logical *wantq, integer *n, CReal *t, integer *
    ldt, CReal *q, integer *ldq, integer *j1, integer *n1, integer *n2,
    CReal *work, integer *info);

/* Subroutine */ int F77NAME(slag2)(CReal *a, integer *lda, CReal *b, integer *ldb,
    CReal *safmin, CReal *scale1, CReal *scale2, CReal *wr1, CReal *wr2, CReal *
    wi);

/* Subroutine */ int F77NAME(slags2)(logical *upper, CReal *a1, CReal *a2, CReal *a3,
    CReal *b1, CReal *b2, CReal *b3, CReal *csu, CReal *snu, CReal *csv, CReal *
    snv, CReal *csq, CReal *snq);

/* Subroutine */ int F77NAME(slagtf)(integer *n, CReal *a, CReal *lambda, CReal *b, CReal
    *c__, CReal *tol, CReal *d__, integer *in, integer *info);

/* Subroutine */ int F77NAME(slagtm)(char *trans, integer *n, integer *nrhs, CReal *
    alpha, CReal *dl, CReal *d__, CReal *du, CReal *x, integer *ldx, CReal *
    beta, CReal *b, integer *ldb);

/* Subroutine */ int F77NAME(slagts)(integer *job, integer *n, CReal *a, CReal *b, CReal
    *c__, CReal *d__, integer *in, CReal *y, CReal *tol, integer *info);

/* Subroutine */ int F77NAME(slagv2)(CReal *a, integer *lda, CReal *b, integer *ldb,
    CReal *alphar, CReal *alphai, CReal *beta, CReal *csl, CReal *snl, CReal *
    csr, CReal *snr);

/* Subroutine */ int F77NAME(slahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, CReal *h__, integer *ldh, CReal *wr, CReal *
    wi, integer *iloz, integer *ihiz, CReal *z__, integer *ldz, integer *
    info);

/* Subroutine */ int F77NAME(slahrd)(integer *n, integer *k, integer *nb, CReal *a,
    integer *lda, CReal *tau, CReal *t, integer *ldt, CReal *y, integer *ldy);

/* Subroutine */ int F77NAME(slaic1)(integer *job, integer *j, CReal *x, CReal *sest,
    CReal *w, CReal *gamma, CReal *sestpr, CReal *s, CReal *c__);

/* Subroutine */ int F77NAME(slaln2)(logical *ltrans, integer *na, integer *nw, CReal *
    smin, CReal *ca, CReal *a, integer *lda, CReal *d1, CReal *d2, CReal *b,
    integer *ldb, CReal *wr, CReal *wi, CReal *x, integer *ldx, CReal *scale,
    CReal *xnorm, integer *info);

/* Subroutine */ int F77NAME(slals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, CReal *b, integer *ldb, CReal *bx,
    integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, CReal *givnum, integer *ldgnum, CReal *poles, CReal *
    difl, CReal *difr, CReal *z__, integer *k, CReal *c__, CReal *s, CReal *
    work, integer *info);

/* Subroutine */ int F77NAME(slalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, CReal *b, integer *ldb, CReal *bx, integer *ldbx, CReal *
    u, integer *ldu, CReal *vt, integer *k, CReal *difl, CReal *difr, CReal *
    z__, CReal *poles, integer *givptr, integer *givcol, integer *ldgcol,
    integer *perm, CReal *givnum, CReal *c__, CReal *s, CReal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(slalsd)(char *uplo, integer *smlsiz, integer *n, integer
    *nrhs, CReal *d__, CReal *e, CReal *b, integer *ldb, CReal *rcond,
    integer *rank, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slamc1)(integer *beta, integer *t, logical *rnd, logical
    *ieee1);

/* Subroutine */ int F77NAME(slamc2)(integer *beta, integer *t, logical *rnd, CReal *
    eps, integer *emin, CReal *rmin, integer *emax, CReal *rmax);

/* Subroutine */ int F77NAME(slamc4)(integer *emin, CReal *start, integer *base);

/* Subroutine */ int F77NAME(slamc5)(integer *beta, integer *p, integer *emin,
    logical *ieee, integer *emax, CReal *rmax);

/* Subroutine */ int F77NAME(slamrg)(integer *n1, integer *n2, CReal *a, integer *
    strd1, integer *strd2, integer *index);

/* Subroutine */ int F77NAME(slanv2)(CReal *a, CReal *b, CReal *c__, CReal *d__, CReal *
    rt1r, CReal *rt1i, CReal *rt2r, CReal *rt2i, CReal *cs, CReal *sn);

/* Subroutine */ int F77NAME(slapll)(integer *n, CReal *x, integer *incx, CReal *y,
    integer *incy, CReal *ssmin);

/* Subroutine */ int F77NAME(slapmt)(logical *forwrd, integer *m, integer *n, CReal *x,
     integer *ldx, integer *k);

/* Subroutine */ int F77NAME(slaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     CReal *ab, integer *ldab, CReal *r__, CReal *c__, CReal *rowcnd, CReal *
    colcnd, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqge)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *r__, CReal *c__, CReal *rowcnd, CReal *colcnd, CReal *amax, char *
    equed);

/* Subroutine */ int F77NAME(slaqp2)(integer *m, integer *n, integer *offset, CReal *a,
     integer *lda, integer *jpvt, CReal *tau, CReal *vn1, CReal *vn2, CReal *
    work);

/* Subroutine */ int F77NAME(slaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, CReal *a, integer *lda, integer *jpvt, CReal *tau,
    CReal *vn1, CReal *vn2, CReal *auxv, CReal *f, integer *ldf);

/* Subroutine */ int F77NAME(slaqsb)(char *uplo, integer *n, integer *kd, CReal *ab,
    integer *ldab, CReal *s, CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqsp)(char *uplo, integer *n, CReal *ap, CReal *s, CReal *
    scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqsy)(char *uplo, integer *n, CReal *a, integer *lda,
    CReal *s, CReal *scond, CReal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqtr)(logical *ltran, logical *lCReal, integer *n, CReal
    *t, integer *ldt, CReal *b, CReal *w, CReal *scale, CReal *x, CReal *work,
    integer *info);

/* Subroutine */ int F77NAME(slar1v)(integer *n, integer *b1, integer *bn, CReal *
    sigma, CReal *d__, CReal *l, CReal *ld, CReal *lld, CReal *gersch, CReal *
    z__, CReal *ztz, CReal *mingma, integer *r__, integer *isuppz, CReal *
    work);

/* Subroutine */ int F77NAME(slar2v)(integer *n, CReal *x, CReal *y, CReal *z__, integer
    *incx, CReal *c__, CReal *s, integer *incc);

/* Subroutine */ int F77NAME(slarf)(char *side, integer *m, integer *n, CReal *v,
    integer *incv, CReal *tau, CReal *c__, integer *ldc, CReal *work);

/* Subroutine */ int F77NAME(slarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, CReal *v, integer *ldv,
    CReal *t, integer *ldt, CReal *c__, integer *ldc, CReal *work, integer *
    ldwork);

/* Subroutine */ int F77NAME(slarfg)(integer *n, CReal *alpha, CReal *x, integer *incx,
    CReal *tau);

/* Subroutine */ int F77NAME(slarft)(char *direct, char *storev, integer *n, integer *
    k, CReal *v, integer *ldv, CReal *tau, CReal *t, integer *ldt);

/* Subroutine */ int F77NAME(slarfx)(char *side, integer *m, integer *n, CReal *v,
    CReal *tau, CReal *c__, integer *ldc, CReal *work);

/* Subroutine */ int F77NAME(slargv)(integer *n, CReal *x, integer *incx, CReal *y,
    integer *incy, CReal *c__, integer *incc);

/* Subroutine */ int F77NAME(slarnv)(integer *idist, integer *iseed, integer *n, CReal
    *x);

/* Subroutine */ int F77NAME(slarrb)(integer *n, CReal *d__, CReal *l, CReal *ld, CReal *
    lld, integer *ifirst, integer *ilast, CReal *sigma, CReal *reltol, CReal
    *w, CReal *wgap, CReal *werr, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slarre)(integer *n, CReal *d__, CReal *e, CReal *tol,
    integer *nsplit, integer *isplit, integer *m, CReal *w, CReal *woff,
    CReal *gersch, CReal *work, integer *info);

/* Subroutine */ int F77NAME(slarrf)(integer *n, CReal *d__, CReal *l, CReal *ld, CReal *
    lld, integer *ifirst, integer *ilast, CReal *w, CReal *dplus, CReal *
    lplus, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slarrv)(integer *n, CReal *d__, CReal *l, integer *isplit,
    integer *m, CReal *w, integer *iblock, CReal *gersch, CReal *tol, CReal *
    z__, integer *ldz, integer *isuppz, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(slartg)(CReal *f, CReal *g, CReal *cs, CReal *sn, CReal *r__);

/* Subroutine */ int F77NAME(slartv)(integer *n, CReal *x, integer *incx, CReal *y,
    integer *incy, CReal *c__, CReal *s, integer *incc);

/* Subroutine */ int F77NAME(slaruv)(integer *iseed, integer *n, CReal *x);

/* Subroutine */ int F77NAME(slarz)(char *side, integer *m, integer *n, integer *l,
    CReal *v, integer *incv, CReal *tau, CReal *c__, integer *ldc, CReal *
    work);

/* Subroutine */ int F77NAME(slarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, CReal *v,
    integer *ldv, CReal *t, integer *ldt, CReal *c__, integer *ldc, CReal *
    work, integer *ldwork);

/* Subroutine */ int F77NAME(slarzt)(char *direct, char *storev, integer *n, integer *
    k, CReal *v, integer *ldv, CReal *tau, CReal *t, integer *ldt);

/* Subroutine */ int F77NAME(slas2)(CReal *f, CReal *g, CReal *h__, CReal *ssmin, CReal *
    ssmax);

/* Subroutine */ int F77NAME(slascl)(char *type__, integer *kl, integer *ku, CReal *
    cfrom, CReal *cto, integer *m, integer *n, CReal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(slasd0)(integer *n, integer *sqre, CReal *d__, CReal *e,
    CReal *u, integer *ldu, CReal *vt, integer *ldvt, integer *smlsiz,
    integer *iwork, CReal *work, integer *info);

/* Subroutine */ int F77NAME(slasd1)(integer *nl, integer *nr, integer *sqre, CReal *
    d__, CReal *alpha, CReal *beta, CReal *u, integer *ldu, CReal *vt,
    integer *ldvt, integer *idxq, integer *iwork, CReal *work, integer *
    info);

/* Subroutine */ int F77NAME(slasd2)(integer *nl, integer *nr, integer *sqre, integer
    *k, CReal *d__, CReal *z__, CReal *alpha, CReal *beta, CReal *u, integer *
    ldu, CReal *vt, integer *ldvt, CReal *dsigma, CReal *u2, integer *ldu2,
    CReal *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc,
     integer *idxq, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(slasd3)(integer *nl, integer *nr, integer *sqre, integer
    *k, CReal *d__, CReal *q, integer *ldq, CReal *dsigma, CReal *u, integer *
    ldu, CReal *u2, integer *ldu2, CReal *vt, integer *ldvt, CReal *vt2,
    integer *ldvt2, integer *idxc, integer *ctot, CReal *z__, integer *
    info);

/* Subroutine */ int F77NAME(slasd4)(integer *n, integer *i__, CReal *d__, CReal *z__,
    CReal *delta, CReal *rho, CReal *sigma, CReal *work, integer *info);

/* Subroutine */ int F77NAME(slasd5)(integer *i__, CReal *d__, CReal *z__, CReal *delta,
    CReal *rho, CReal *dsigma, CReal *work);

/* Subroutine */ int F77NAME(slasd6)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, CReal *d__, CReal *vf, CReal *vl, CReal *alpha, CReal *beta,
     integer *idxq, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, CReal *givnum, integer *ldgnum, CReal *poles, CReal *
    difl, CReal *difr, CReal *z__, integer *k, CReal *c__, CReal *s, CReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slasd7)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *k, CReal *d__, CReal *z__, CReal *zw, CReal *vf,
    CReal *vfw, CReal *vl, CReal *vlw, CReal *alpha, CReal *beta, CReal *dsigma,
     integer *idx, integer *idxp, integer *idxq, integer *perm, integer *
    givptr, integer *givcol, integer *ldgcol, CReal *givnum, integer *
    ldgnum, CReal *c__, CReal *s, integer *info);

/* Subroutine */ int F77NAME(slasd8)(integer *icompq, integer *k, CReal *d__, CReal *
    z__, CReal *vf, CReal *vl, CReal *difl, CReal *difr, integer *lddifr,
    CReal *dsigma, CReal *work, integer *info);

/* Subroutine */ int F77NAME(slasd9)(integer *icompq, integer *ldu, integer *k, CReal *
    d__, CReal *z__, CReal *vf, CReal *vl, CReal *difl, CReal *difr, CReal *
    dsigma, CReal *work, integer *info);

/* Subroutine */ int F77NAME(slasda)(integer *icompq, integer *smlsiz, integer *n,
    integer *sqre, CReal *d__, CReal *e, CReal *u, integer *ldu, CReal *vt,
    integer *k, CReal *difl, CReal *difr, CReal *z__, CReal *poles, integer *
    givptr, integer *givcol, integer *ldgcol, integer *perm, CReal *givnum,
     CReal *c__, CReal *s, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slasdq)(char *uplo, integer *sqre, integer *n, integer *
    ncvt, integer *nru, integer *ncc, CReal *d__, CReal *e, CReal *vt,
    integer *ldvt, CReal *u, integer *ldu, CReal *c__, integer *ldc, CReal *
    work, integer *info);

/* Subroutine */ int F77NAME(slasdt)(integer *n, integer *lvl, integer *nd, integer *
    inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int F77NAME(slaset)(char *uplo, integer *m, integer *n, CReal *alpha,
    CReal *beta, CReal *a, integer *lda);

/* Subroutine */ int F77NAME(slasq1)(integer *n, CReal *d__, CReal *e, CReal *work,
    integer *info);

/* Subroutine */ int F77NAME(slasq2)(integer *n, CReal *z__, integer *info);

/* Subroutine */ int F77NAME(slasq3)(integer *i0, integer *n0, CReal *z__, integer *pp,
     CReal *dmin__, CReal *sigma, CReal *desig, CReal *qmax, integer *nfail,
    integer *iter, integer *ndiv, logical *ieee);

/* Subroutine */ int F77NAME(slasq4)(integer *i0, integer *n0, CReal *z__, integer *pp,
     integer *n0in, CReal *dmin__, CReal *dmin1, CReal *dmin2, CReal *dn,
    CReal *dn1, CReal *dn2, CReal *tau, integer *ttype);

/* Subroutine */ int F77NAME(slasq5)(integer *i0, integer *n0, CReal *z__, integer *pp,
     CReal *tau, CReal *dmin__, CReal *dmin1, CReal *dmin2, CReal *dn, CReal *
    dnm1, CReal *dnm2, logical *ieee);

/* Subroutine */ int F77NAME(slasq6)(integer *i0, integer *n0, CReal *z__, integer *pp,
     CReal *dmin__, CReal *dmin1, CReal *dmin2, CReal *dn, CReal *dnm1, CReal *
    dnm2);

/* Subroutine */ int F77NAME(slasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, CReal *c__, CReal *s, CReal *a, integer *lda);

/* Subroutine */ int F77NAME(slasrt)(char *id, integer *n, CReal *d__, integer *info);

/* Subroutine */ int F77NAME(slassq)(integer *n, CReal *x, integer *incx, CReal *scale,
    CReal *sumsq);

/* Subroutine */ int F77NAME(slasv2)(CReal *f, CReal *g, CReal *h__, CReal *ssmin, CReal *
    ssmax, CReal *snr, CReal *csr, CReal *snl, CReal *csl);

/* Subroutine */ int F77NAME(slaswp)(integer *n, CReal *a, integer *lda, integer *k1,
    integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(slasy2)(logical *ltranl, logical *ltranr, integer *isgn,
    integer *n1, integer *n2, CReal *tl, integer *ldtl, CReal *tr, integer *
    ldtr, CReal *b, integer *ldb, CReal *scale, CReal *x, integer *ldx, CReal
    *xnorm, integer *info);

/* Subroutine */ int F77NAME(slasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     CReal *a, integer *lda, integer *ipiv, CReal *w, integer *ldw, integer
    *info);

/* Subroutine */ int F77NAME(slatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, CReal *ab, integer *ldab, CReal *x,
    CReal *scale, CReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(slatdf)(integer *ijob, integer *n, CReal *z__, integer *
    ldz, CReal *rhs, CReal *rdsum, CReal *rdscal, integer *ipiv, integer *
    jpiv);

/* Subroutine */ int F77NAME(slatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, CReal *ap, CReal *x, CReal *scale, CReal *cnorm,
    integer *info);

/* Subroutine */ int F77NAME(slatrd)(char *uplo, integer *n, integer *nb, CReal *a,
    integer *lda, CReal *e, CReal *tau, CReal *w, integer *ldw);

/* Subroutine */ int F77NAME(slatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, CReal *a, integer *lda, CReal *x, CReal *scale, CReal
    *cnorm, integer *info);

/* Subroutine */ int F77NAME(slatrz)(integer *m, integer *n, integer *l, CReal *a,
    integer *lda, CReal *tau, CReal *work);

/* Subroutine */ int F77NAME(slatzm)(char *side, integer *m, integer *n, CReal *v,
    integer *incv, CReal *tau, CReal *c1, CReal *c2, integer *ldc, CReal *
    work);

/* Subroutine */ int F77NAME(slauu2)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(slauum)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(sopgtr)(char *uplo, integer *n, CReal *ap, CReal *tau,
    CReal *q, integer *ldq, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sopmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, CReal *ap, CReal *tau, CReal *c__, integer *ldc, CReal *work,
    integer *info);

/* Subroutine */ int F77NAME(sorg2l)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sorg2r)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sorgbr)(char *vect, integer *m, integer *n, integer *k,
    CReal *a, integer *lda, CReal *tau, CReal *work, integer *lwork, integer
    *info);

/* Subroutine */ int F77NAME(sorghr)(integer *n, integer *ilo, integer *ihi, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgl2)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sorglq)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgql)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgqr)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgr2)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sorgrq)(integer *m, integer *n, integer *k, CReal *a,
    integer *lda, CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgtr)(char *uplo, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *info);

/* Subroutine */ int F77NAME(sorm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *info);

/* Subroutine */ int F77NAME(sormbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__,
    integer *ldc, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, CReal *a, integer *lda, CReal *tau, CReal *
    c__, integer *ldc, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *info);

/* Subroutine */ int F77NAME(sormlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormql)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *info);

/* Subroutine */ int F77NAME(sormr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, CReal *a, integer *lda, CReal *tau, CReal *c__,
    integer *ldc, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sormrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, CReal *a, integer *lda, CReal *tau, CReal *c__,
    integer *ldc, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, CReal *a, integer *lda, CReal *tau, CReal *c__, integer *ldc,
     CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(spbcon)(char *uplo, integer *n, integer *kd, CReal *ab,
    integer *ldab, CReal *anorm, CReal *rcond, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(spbequ)(char *uplo, integer *n, integer *kd, CReal *ab,
    integer *ldab, CReal *s, CReal *scond, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(spbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, CReal *ab, integer *ldab, CReal *afb, integer *ldafb, CReal *b,
    integer *ldb, CReal *x, integer *ldx, CReal *ferr, CReal *berr, CReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spbstf)(char *uplo, integer *n, integer *kd, CReal *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, CReal *ab, integer *ldab, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(spbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, CReal *ab, integer *ldab, CReal *afb, integer *ldafb,
    char *equed, CReal *s, CReal *b, integer *ldb, CReal *x, integer *ldx,
    CReal *rcond, CReal *ferr, CReal *berr, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(spbtf2)(char *uplo, integer *n, integer *kd, CReal *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbtrf)(char *uplo, integer *n, integer *kd, CReal *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, CReal *ab, integer *ldab, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(spocon)(char *uplo, integer *n, CReal *a, integer *lda,
    CReal *anorm, CReal *rcond, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spoequ)(integer *n, CReal *a, integer *lda, CReal *s, CReal
    *scond, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(sporfs)(char *uplo, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *af, integer *ldaf, CReal *b, integer *ldb, CReal *x,
     integer *ldx, CReal *ferr, CReal *berr, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(sposv)(char *uplo, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, CReal *a, integer *lda, CReal *af, integer *ldaf, char *equed,
    CReal *s, CReal *b, integer *ldb, CReal *x, integer *ldx, CReal *rcond,
    CReal *ferr, CReal *berr, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spotf2)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotrf)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotri)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotrs)(char *uplo, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sppcon)(char *uplo, integer *n, CReal *ap, CReal *anorm,
    CReal *rcond, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sppequ)(char *uplo, integer *n, CReal *ap, CReal *s, CReal *
    scond, CReal *amax, integer *info);

/* Subroutine */ int F77NAME(spprfs)(char *uplo, integer *n, integer *nrhs, CReal *ap,
    CReal *afp, CReal *b, integer *ldb, CReal *x, integer *ldx, CReal *ferr,
    CReal *berr, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sppsv)(char *uplo, integer *n, integer *nrhs, CReal *ap,
    CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, CReal *ap, CReal *afp, char *equed, CReal *s, CReal *b, integer *
    ldb, CReal *x, integer *ldx, CReal *rcond, CReal *ferr, CReal *berr, CReal
    *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spptrf)(char *uplo, integer *n, CReal *ap, integer *info);

/* Subroutine */ int F77NAME(spptri)(char *uplo, integer *n, CReal *ap, integer *info);

/* Subroutine */ int F77NAME(spptrs)(char *uplo, integer *n, integer *nrhs, CReal *ap,
    CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptcon)(integer *n, CReal *d__, CReal *e, CReal *anorm,
    CReal *rcond, CReal *work, integer *info);

/* Subroutine */ int F77NAME(spteqr)(char *compz, integer *n, CReal *d__, CReal *e,
    CReal *z__, integer *ldz, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sptrfs)(integer *n, integer *nrhs, CReal *d__, CReal *e,
    CReal *df, CReal *ef, CReal *b, integer *ldb, CReal *x, integer *ldx,
    CReal *ferr, CReal *berr, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sptsv)(integer *n, integer *nrhs, CReal *d__, CReal *e,
    CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptsvx)(char *fact, integer *n, integer *nrhs, CReal *d__,
     CReal *e, CReal *df, CReal *ef, CReal *b, integer *ldb, CReal *x, integer
    *ldx, CReal *rcond, CReal *ferr, CReal *berr, CReal *work, integer *info);

/* Subroutine */ int F77NAME(spttrf)(integer *n, CReal *d__, CReal *e, integer *info);

/* Subroutine */ int F77NAME(spttrs)(integer *n, integer *nrhs, CReal *d__, CReal *e,
    CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptts2)(integer *n, integer *nrhs, CReal *d__, CReal *e,
    CReal *b, integer *ldb);

/* Subroutine */ int F77NAME(srscl)(integer *n, CReal *sa, CReal *sx, integer *incx);

/* Subroutine */ int F77NAME(ssbev)(char *jobz, char *uplo, integer *n, integer *kd,
    CReal *ab, integer *ldab, CReal *w, CReal *z__, integer *ldz, CReal *work,
     integer *info);

/* Subroutine */ int F77NAME(ssbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    CReal *ab, integer *ldab, CReal *w, CReal *z__, integer *ldz, CReal *work,
     integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, CReal *ab, integer *ldab, CReal *q, integer *ldq, CReal *vl,
     CReal *vu, integer *il, integer *iu, CReal *abstol, integer *m, CReal *
    w, CReal *z__, integer *ldz, CReal *work, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(ssbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, CReal *ab, integer *ldab, CReal *bb, integer *ldbb, CReal *
    x, integer *ldx, CReal *work, integer *info);

/* Subroutine */ int F77NAME(ssbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, CReal *ab, integer *ldab, CReal *bb, integer *ldbb, CReal *
    w, CReal *z__, integer *ldz, CReal *work, integer *info);

/* Subroutine */ int F77NAME(ssbgvd)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, CReal *ab, integer *ldab, CReal *bb, integer *ldbb, CReal *
    w, CReal *z__, integer *ldz, CReal *work, integer *lwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, CReal *ab, integer *ldab, CReal *bb, integer *
    ldbb, CReal *q, integer *ldq, CReal *vl, CReal *vu, integer *il, integer
    *iu, CReal *abstol, integer *m, CReal *w, CReal *z__, integer *ldz, CReal
    *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    CReal *ab, integer *ldab, CReal *d__, CReal *e, CReal *q, integer *ldq,
    CReal *work, integer *info);

/* Subroutine */ int F77NAME(sspcon)(char *uplo, integer *n, CReal *ap, integer *ipiv,
    CReal *anorm, CReal *rcond, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sspev)(char *jobz, char *uplo, integer *n, CReal *ap,
    CReal *w, CReal *z__, integer *ldz, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sspevd)(char *jobz, char *uplo, integer *n, CReal *ap,
    CReal *w, CReal *z__, integer *ldz, CReal *work, integer *lwork, integer
    *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sspevx)(char *jobz, char *range, char *uplo, integer *n,
    CReal *ap, CReal *vl, CReal *vu, integer *il, integer *iu, CReal *abstol,
    integer *m, CReal *w, CReal *z__, integer *ldz, CReal *work, integer *
    iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(sspgst)(integer *itype, char *uplo, integer *n, CReal *ap,
     CReal *bp, integer *info);

/* Subroutine */ int F77NAME(sspgv)(integer *itype, char *jobz, char *uplo, integer *
    n, CReal *ap, CReal *bp, CReal *w, CReal *z__, integer *ldz, CReal *work,
    integer *info);

/* Subroutine */ int F77NAME(sspgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, CReal *ap, CReal *bp, CReal *w, CReal *z__, integer *ldz, CReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sspgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, CReal *ap, CReal *bp, CReal *vl, CReal *vu, integer *il,
     integer *iu, CReal *abstol, integer *m, CReal *w, CReal *z__, integer *
    ldz, CReal *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssprfs)(char *uplo, integer *n, integer *nrhs, CReal *ap,
    CReal *afp, integer *ipiv, CReal *b, integer *ldb, CReal *x, integer *
    ldx, CReal *ferr, CReal *berr, CReal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(sspsv)(char *uplo, integer *n, integer *nrhs, CReal *ap,
    integer *ipiv, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, CReal *ap, CReal *afp, integer *ipiv, CReal *b, integer *ldb, CReal
    *x, integer *ldx, CReal *rcond, CReal *ferr, CReal *berr, CReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ssptrd)(char *uplo, integer *n, CReal *ap, CReal *d__,
    CReal *e, CReal *tau, integer *info);

/* Subroutine */ int F77NAME(ssptrf)(char *uplo, integer *n, CReal *ap, integer *ipiv,
    integer *info);

/* Subroutine */ int F77NAME(ssptri)(char *uplo, integer *n, CReal *ap, integer *ipiv,
    CReal *work, integer *info);

/* Subroutine */ int F77NAME(ssptrs)(char *uplo, integer *n, integer *nrhs, CReal *ap,
    integer *ipiv, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sstebz)(char *range, char *order, integer *n, CReal *vl,
    CReal *vu, integer *il, integer *iu, CReal *abstol, CReal *d__, CReal *e,
    integer *m, integer *nsplit, CReal *w, integer *iblock, integer *
    isplit, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sstedc)(char *compz, integer *n, CReal *d__, CReal *e,
    CReal *z__, integer *ldz, CReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstegr)(char *jobz, char *range, integer *n, CReal *d__,
    CReal *e, CReal *vl, CReal *vu, integer *il, integer *iu, CReal *abstol,
    integer *m, CReal *w, CReal *z__, integer *ldz, integer *isuppz, CReal *
    work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstein)(integer *n, CReal *d__, CReal *e, integer *m, CReal
    *w, integer *iblock, integer *isplit, CReal *z__, integer *ldz, CReal *
    work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssteqr)(char *compz, integer *n, CReal *d__, CReal *e,
    CReal *z__, integer *ldz, CReal *work, integer *info);

/* Subroutine */ int F77NAME(ssterf)(integer *n, CReal *d__, CReal *e, integer *info);

/* Subroutine */ int F77NAME(sstev)(char *jobz, integer *n, CReal *d__, CReal *e, CReal *
    z__, integer *ldz, CReal *work, integer *info);

/* Subroutine */ int F77NAME(sstevd)(char *jobz, integer *n, CReal *d__, CReal *e, CReal
    *z__, integer *ldz, CReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstevr)(char *jobz, char *range, integer *n, CReal *d__,
    CReal *e, CReal *vl, CReal *vu, integer *il, integer *iu, CReal *abstol,
    integer *m, CReal *w, CReal *z__, integer *ldz, integer *isuppz, CReal *
    work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstevx)(char *jobz, char *range, integer *n, CReal *d__,
    CReal *e, CReal *vl, CReal *vu, integer *il, integer *iu, CReal *abstol,
    integer *m, CReal *w, CReal *z__, integer *ldz, CReal *work, integer *
    iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssycon)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *ipiv, CReal *anorm, CReal *rcond, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(ssyev)(char *jobz, char *uplo, integer *n, CReal *a,
    integer *lda, CReal *w, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssyevd)(char *jobz, char *uplo, integer *n, CReal *a,
    integer *lda, CReal *w, CReal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssyevr)(char *jobz, char *range, char *uplo, integer *n,
    CReal *a, integer *lda, CReal *vl, CReal *vu, integer *il, integer *iu,
    CReal *abstol, integer *m, CReal *w, CReal *z__, integer *ldz, integer *
    isuppz, CReal *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ssyevx)(char *jobz, char *range, char *uplo, integer *n,
    CReal *a, integer *lda, CReal *vl, CReal *vu, integer *il, integer *iu,
    CReal *abstol, integer *m, CReal *w, CReal *z__, integer *ldz, CReal *
    work, integer *lwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssygs2)(integer *itype, char *uplo, integer *n, CReal *a,
    integer *lda, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ssygst)(integer *itype, char *uplo, integer *n, CReal *a,
    integer *lda, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ssygv)(integer *itype, char *jobz, char *uplo, integer *
    n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *w, CReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssygvd)(integer *itype, char *jobz, char *uplo, integer *
    n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *w, CReal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssygvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *
    vl, CReal *vu, integer *il, integer *iu, CReal *abstol, integer *m,
    CReal *w, CReal *z__, integer *ldz, CReal *work, integer *lwork, integer
    *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssyrfs)(char *uplo, integer *n, integer *nrhs, CReal *a,
    integer *lda, CReal *af, integer *ldaf, integer *ipiv, CReal *b,
    integer *ldb, CReal *x, integer *ldx, CReal *ferr, CReal *berr, CReal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ssysv)(char *uplo, integer *n, integer *nrhs, CReal *a,
    integer *lda, integer *ipiv, CReal *b, integer *ldb, CReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, CReal *a, integer *lda, CReal *af, integer *ldaf, integer *ipiv,
    CReal *b, integer *ldb, CReal *x, integer *ldx, CReal *rcond, CReal *ferr,
     CReal *berr, CReal *work, integer *lwork, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(ssytd2)(char *uplo, integer *n, CReal *a, integer *lda,
    CReal *d__, CReal *e, CReal *tau, integer *info);

/* Subroutine */ int F77NAME(ssytf2)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(ssytrd)(char *uplo, integer *n, CReal *a, integer *lda,
    CReal *d__, CReal *e, CReal *tau, CReal *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(ssytrf)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *ipiv, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssytri)(char *uplo, integer *n, CReal *a, integer *lda,
    integer *ipiv, CReal *work, integer *info);

/* Subroutine */ int F77NAME(ssytrs)(char *uplo, integer *n, integer *nrhs, CReal *a,
    integer *lda, integer *ipiv, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(stbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, CReal *ab, integer *ldab, CReal *rcond, CReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, CReal *ab, integer *ldab, CReal *b, integer
    *ldb, CReal *x, integer *ldx, CReal *ferr, CReal *berr, CReal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, CReal *ab, integer *ldab, CReal *b, integer
    *ldb, integer *info);

/* Subroutine */ int F77NAME(stgevc)(char *side, char *howmny, logical *select,
    integer *n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *vl,
    integer *ldvl, CReal *vr, integer *ldvr, integer *mm, integer *m, CReal
    *work, integer *info);

/* Subroutine */ int F77NAME(stgex2)(logical *wantq, logical *wantz, integer *n, CReal
    *a, integer *lda, CReal *b, integer *ldb, CReal *q, integer *ldq, CReal *
    z__, integer *ldz, integer *j1, integer *n1, integer *n2, CReal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(stgexc)(logical *wantq, logical *wantz, integer *n, CReal
    *a, integer *lda, CReal *b, integer *ldb, CReal *q, integer *ldq, CReal *
    z__, integer *ldz, integer *ifst, integer *ilst, CReal *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(stgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, CReal *a, integer *lda, CReal *b, integer *
    ldb, CReal *alphar, CReal *alphai, CReal *beta, CReal *q, integer *ldq,
    CReal *z__, integer *ldz, integer *m, CReal *pl, CReal *pr, CReal *dif,
    CReal *work, integer *lwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(stgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, CReal *a, integer *lda,
     CReal *b, integer *ldb, CReal *tola, CReal *tolb, CReal *alpha, CReal *
    beta, CReal *u, integer *ldu, CReal *v, integer *ldv, CReal *q, integer *
    ldq, CReal *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(stgsna)(char *job, char *howmny, logical *select,
    integer *n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *vl,
    integer *ldvl, CReal *vr, integer *ldvr, CReal *s, CReal *dif, integer *
    mm, integer *m, CReal *work, integer *lwork, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(stgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *c__, integer *
    ldc, CReal *d__, integer *ldd, CReal *e, integer *lde, CReal *f, integer
    *ldf, CReal *scale, CReal *rdsum, CReal *rdscal, integer *iwork, integer
    *pq, integer *info);

/* Subroutine */ int F77NAME(stgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *c__, integer *
    ldc, CReal *d__, integer *ldd, CReal *e, integer *lde, CReal *f, integer
    *ldf, CReal *scale, CReal *dif, CReal *work, integer *lwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(stpcon)(char *norm, char *uplo, char *diag, integer *n,
    CReal *ap, CReal *rcond, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, CReal *ap, CReal *b, integer *ldb, CReal *x, integer *ldx,
     CReal *ferr, CReal *berr, CReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stptri)(char *uplo, char *diag, integer *n, CReal *ap,
    integer *info);

/* Subroutine */ int F77NAME(stptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, CReal *ap, CReal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(strcon)(char *norm, char *uplo, char *diag, integer *n,
    CReal *a, integer *lda, CReal *rcond, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(strevc)(char *side, char *howmny, logical *select,
    integer *n, CReal *t, integer *ldt, CReal *vl, integer *ldvl, CReal *vr,
    integer *ldvr, integer *mm, integer *m, CReal *work, integer *info);

/* Subroutine */ int F77NAME(strexc)(char *compq, integer *n, CReal *t, integer *ldt,
    CReal *q, integer *ldq, integer *ifst, integer *ilst, CReal *work,
    integer *info);

/* Subroutine */ int F77NAME(strrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *x,
    integer *ldx, CReal *ferr, CReal *berr, CReal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(strsen)(char *job, char *compq, logical *select, integer
    *n, CReal *t, integer *ldt, CReal *q, integer *ldq, CReal *wr, CReal *wi,
    integer *m, CReal *s, CReal *sep, CReal *work, integer *lwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(strsna)(char *job, char *howmny, logical *select,
    integer *n, CReal *t, integer *ldt, CReal *vl, integer *ldvl, CReal *vr,
    integer *ldvr, CReal *s, CReal *sep, integer *mm, integer *m, CReal *
    work, integer *ldwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(strsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, CReal *a, integer *lda, CReal *b, integer *ldb, CReal *
    c__, integer *ldc, CReal *scale, integer *info);

/* Subroutine */ int F77NAME(strti2)(char *uplo, char *diag, integer *n, CReal *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(strtri)(char *uplo, char *diag, integer *n, CReal *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(strtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, CReal *a, integer *lda, CReal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(stzrqf)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, integer *info);

/* Subroutine */ int F77NAME(stzrzf)(integer *m, integer *n, CReal *a, integer *lda,
    CReal *tau, CReal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(xerbla)(char *srname, integer *info);

/* Subroutine */ int F77NAME(zbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, doubleCReal *d__, doubleCReal *e, doubleComplex *vt,
    integer *ldvt, doubleComplex *u, integer *ldu, doubleComplex *c__,
    integer *ldc, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zdrot)(integer *n, doubleComplex *cx, integer *incx,
    doubleComplex *cy, integer *incy, doubleCReal *c__, doubleCReal *s);

/* Subroutine */ int F77NAME(zdrscl)(integer *n, doubleCReal *sa, doubleComplex *sx,
    integer *incx);

/* Subroutine */ int F77NAME(zgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, doubleComplex *ab, integer *ldab,
    doubleCReal *d__, doubleCReal *e, doubleComplex *q, integer *ldq,
    doubleComplex *pt, integer *ldpt, doubleComplex *c__, integer *ldc,
    doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, integer *ipiv, doubleCReal *anorm,
    doubleCReal *rcond, doubleComplex *work, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, doubleCReal *r__, doubleCReal *c__,
    doubleCReal *rowcnd, doubleCReal *colcnd, doubleCReal *amax, integer *
    info);

/* Subroutine */ int F77NAME(zgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleComplex *ab, integer *ldab, doubleComplex *
    afb, integer *ldafb, integer *ipiv, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, doubleComplex *ab, integer *ldab, integer *ipiv, doubleComplex *
    b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, doubleComplex *ab, integer *ldab,
    doubleComplex *afb, integer *ldafb, integer *ipiv, char *equed,
    doubleCReal *r__, doubleCReal *c__, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleCReal *rcond, doubleCReal *ferr,
    doubleCReal *berr, doubleComplex *work, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doubleComplex *ab, integer *ldab, integer *ipiv,
    doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doubleCReal *scale, integer *m, doubleComplex *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(zgebal)(char *job, integer *n, doubleComplex *a, integer
    *lda, integer *ilo, integer *ihi, doubleCReal *scale, integer *info);

/* Subroutine */ int F77NAME(zgebd2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *d__, doubleCReal *e, doubleComplex *tauq,
    doubleComplex *taup, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgebrd)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *d__, doubleCReal *e, doubleComplex *tauq,
    doubleComplex *taup, doubleComplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(zgecon)(char *norm, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *anorm, doubleCReal *rcond, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeequ)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *r__, doubleCReal *c__, doubleCReal *rowcnd,
    doubleCReal *colcnd, doubleCReal *amax, integer *info);

/* Subroutine */ int F77NAME(zgees)(char *jobvs, char *sort, L_fp select, integer *n,
    doubleComplex *a, integer *lda, integer *sdim, doubleComplex *w,
    doubleComplex *vs, integer *ldvs, doubleComplex *work, integer *lwork,
     doubleCReal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, doubleComplex *a, integer *lda, integer *sdim,
    doubleComplex *w, doubleComplex *vs, integer *ldvs, doubleCReal *
    rconde, doubleCReal *rcondv, doubleComplex *work, integer *lwork,
    doubleCReal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zgeev)(char *jobvl, char *jobvr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *w, doubleComplex *vl,
    integer *ldvl, doubleComplex *vr, integer *ldvr, doubleComplex *work,
    integer *lwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleComplex *a, integer *lda, doubleComplex *w,
    doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *ldvr,
    integer *ilo, integer *ihi, doubleCReal *scale, doubleCReal *abnrm,
    doubleCReal *rconde, doubleCReal *rcondv, doubleComplex *work, integer *
    lwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgegs)(char *jobvsl, char *jobvsr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *vsl,
    integer *ldvsl, doubleComplex *vsr, integer *ldvsr, doubleComplex *
    work, integer *lwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgegv)(char *jobvl, char *jobvr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *vl, integer
    *ldvl, doubleComplex *vr, integer *ldvr, doubleComplex *work, integer
    *lwork, doubleCReal *rwork, integer *info);

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
    integer *jpvt, doubleCReal *rcond, integer *rank, doubleComplex *work,
    doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgelsy)(integer *m, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *jpvt, doubleCReal *rcond, integer *rank, doubleComplex *work,
    integer *lwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeql2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgeqlf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgeqp3)(integer *m, integer *n, doubleComplex *a,
    integer *lda, integer *jpvt, doubleComplex *tau, doubleComplex *work,
    integer *lwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeqpf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, integer *jpvt, doubleComplex *tau, doubleComplex *work,
    doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeqr2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgeqrf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgerfs)(char *trans, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
    integer *ldx, doubleCReal *ferr, doubleCReal *berr, doubleComplex *work,
     doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgerq2)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgerqf)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleComplex *tau, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgesc2)(integer *n, doubleComplex *a, integer *lda,
    doubleComplex *rhs, integer *ipiv, integer *jpiv, doubleCReal *scale);

/* Subroutine */ int F77NAME(zgesv)(integer *n, integer *nrhs, doubleComplex *a,
    integer *lda, integer *ipiv, doubleComplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(zgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, integer *ipiv, char *equed, doubleCReal *r__, doubleCReal *c__,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr, doubleComplex *
    work, doubleCReal *rwork, integer *info);

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
    integer *ihi, doubleCReal *lscale, doubleCReal *rscale, integer *m,
    doubleComplex *v, integer *ldv, integer *info);

/* Subroutine */ int F77NAME(zggbal)(char *job, integer *n, doubleComplex *a, integer
    *lda, doubleComplex *b, integer *ldb, integer *ilo, integer *ihi,
    doubleCReal *lscale, doubleCReal *rscale, doubleCReal *work, integer *
    info);

/* Subroutine */ int F77NAME(zgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, integer *sdim, doubleComplex *alpha, doubleComplex *
    beta, doubleComplex *vsl, integer *ldvsl, doubleComplex *vsr, integer
    *ldvsr, doubleComplex *work, integer *lwork, doubleCReal *rwork,
    logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, char *sense, integer *n, doubleComplex *a, integer *lda,
    doubleComplex *b, integer *ldb, integer *sdim, doubleComplex *alpha,
    doubleComplex *beta, doubleComplex *vsl, integer *ldvsl,
    doubleComplex *vsr, integer *ldvsr, doubleCReal *rconde, doubleCReal *
    rcondv, doubleComplex *work, integer *lwork, doubleCReal *rwork,
    integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggev)(char *jobvl, char *jobvr, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *vl, integer
    *ldvl, doubleComplex *vr, integer *ldvr, doubleComplex *work, integer
    *lwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleComplex *alpha, doubleComplex *beta,
    doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *ldvr,
    integer *ilo, integer *ihi, doubleCReal *lscale, doubleCReal *rscale,
    doubleCReal *abnrm, doubleCReal *bbnrm, doubleCReal *rconde, doubleCReal *
    rcondv, doubleComplex *work, integer *lwork, doubleCReal *rwork,
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
    integer *lda, doubleComplex *b, integer *ldb, doubleCReal *alpha,
    doubleCReal *beta, doubleComplex *u, integer *ldu, doubleComplex *v,
    integer *ldv, doubleComplex *q, integer *ldq, doubleComplex *work,
    doubleCReal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, doubleComplex *a, integer *lda, doubleComplex
    *b, integer *ldb, doubleCReal *tola, doubleCReal *tolb, integer *k,
    integer *l, doubleComplex *u, integer *ldu, doubleComplex *v, integer
    *ldv, doubleComplex *q, integer *ldq, integer *iwork, doubleCReal *
    rwork, doubleComplex *tau, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zgtcon)(char *norm, integer *n, doubleComplex *dl,
    doubleComplex *d__, doubleComplex *du, doubleComplex *du2, integer *
    ipiv, doubleCReal *anorm, doubleCReal *rcond, doubleComplex *work,
    integer *info);

/* Subroutine */ int F77NAME(zgtrfs)(char *trans, integer *n, integer *nrhs,
    doubleComplex *dl, doubleComplex *d__, doubleComplex *du,
    doubleComplex *dlf, doubleComplex *df, doubleComplex *duf,
    doubleComplex *du2, integer *ipiv, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgtsv)(integer *n, integer *nrhs, doubleComplex *dl,
    doubleComplex *d__, doubleComplex *du, doubleComplex *b, integer *ldb,
     integer *info);

/* Subroutine */ int F77NAME(zgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doubleComplex *dl, doubleComplex *d__, doubleComplex *du,
    doubleComplex *dlf, doubleComplex *df, doubleComplex *duf,
    doubleComplex *du2, integer *ipiv, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleCReal *rcond, doubleCReal *ferr,
    doubleCReal *berr, doubleComplex *work, doubleCReal *rwork, integer *
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
    doubleComplex *ab, integer *ldab, doubleCReal *w, doubleComplex *z__,
    integer *ldz, doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleCReal *w, doubleComplex *z__,
    integer *ldz, doubleComplex *work, integer *lwork, doubleCReal *rwork,
    integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, doubleComplex *ab, integer *ldab, doubleComplex *q,
    integer *ldq, doubleCReal *vl, doubleCReal *vu, integer *il, integer *
    iu, doubleCReal *abstol, integer *m, doubleCReal *w, doubleComplex *z__,
     integer *ldz, doubleComplex *work, doubleCReal *rwork, integer *iwork,
     integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zhbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, doubleComplex *ab, integer *ldab, doubleComplex *bb,
    integer *ldbb, doubleComplex *x, integer *ldx, doubleComplex *work,
    doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doubleComplex *ab, integer *ldab, doubleComplex *bb,
    integer *ldbb, doubleCReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, doubleComplex *ab, integer *ldab,
    doubleComplex *bb, integer *ldbb, doubleComplex *q, integer *ldq,
    doubleCReal *vl, doubleCReal *vu, integer *il, integer *iu, doubleCReal *
    abstol, integer *m, doubleCReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, doubleCReal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(zhbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleCReal *d__, doubleCReal *e,
    doubleComplex *q, integer *ldq, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zhecon)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleCReal *anorm, doubleCReal *rcond,
    doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zheev)(char *jobz, char *uplo, integer *n, doubleComplex
    *a, integer *lda, doubleCReal *w, doubleComplex *work, integer *lwork,
    doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zheevd)(char *jobz, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleCReal *w, doubleComplex *work,
    integer *lwork, doubleCReal *rwork, integer *lrwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zheevr)(char *jobz, char *range, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleCReal *vl, doubleCReal *vu,
    integer *il, integer *iu, doubleCReal *abstol, integer *m, doubleCReal *
    w, doubleComplex *z__, integer *ldz, integer *isuppz, doubleComplex *
    work, integer *lwork, doubleCReal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zheevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleCReal *vl, doubleCReal *vu,
    integer *il, integer *iu, doubleCReal *abstol, integer *m, doubleCReal *
    w, doubleComplex *z__, integer *ldz, doubleComplex *work, integer *
    lwork, doubleCReal *rwork, integer *iwork, integer *ifail, integer *
    info);

/* Subroutine */ int F77NAME(zhegs2)(integer *itype, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhegst)(integer *itype, char *uplo, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhegv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleCReal *w, doubleComplex *work, integer *lwork, doubleCReal *rwork,
     integer *info);

/* Subroutine */ int F77NAME(zhegvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleCReal *w, doubleComplex *work, integer *lwork, doubleCReal *rwork,
     integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhegvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleCReal *vl, doubleCReal *vu, integer *il, integer *
    iu, doubleCReal *abstol, integer *m, doubleCReal *w, doubleComplex *z__,
     integer *ldz, doubleComplex *work, integer *lwork, doubleCReal *rwork,
     integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zherfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
    integer *ldx, doubleCReal *ferr, doubleCReal *berr, doubleComplex *work,
     doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhesv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *b,
    integer *ldb, doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zhesvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
     integer *ldx, doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr,
    doubleComplex *work, integer *lwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhetf2)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zhetrd)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *d__, doubleCReal *e, doubleComplex *tau,
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
    ldz, doubleComplex *work, integer *lwork, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpcon)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, doubleCReal *anorm, doubleCReal *rcond, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zhpev)(char *jobz, char *uplo, integer *n, doubleComplex
    *ap, doubleCReal *w, doubleComplex *z__, integer *ldz, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhpevd)(char *jobz, char *uplo, integer *n,
    doubleComplex *ap, doubleCReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, integer *lwork, doubleCReal *rwork, integer *
    lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhpevx)(char *jobz, char *range, char *uplo, integer *n,
    doubleComplex *ap, doubleCReal *vl, doubleCReal *vu, integer *il,
    integer *iu, doubleCReal *abstol, integer *m, doubleCReal *w,
    doubleComplex *z__, integer *ldz, doubleComplex *work, doubleCReal *
    rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zhpgst)(integer *itype, char *uplo, integer *n,
    doubleComplex *ap, doubleComplex *bp, integer *info);

/* Subroutine */ int F77NAME(zhpgv)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *ap, doubleComplex *bp, doubleCReal *w, doubleComplex
    *z__, integer *ldz, doubleComplex *work, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doubleComplex *ap, doubleComplex *bp, doubleCReal *w, doubleComplex
    *z__, integer *ldz, doubleComplex *work, integer *lwork, doubleCReal *
    rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doubleComplex *ap, doubleComplex *bp, doubleCReal *
    vl, doubleCReal *vu, integer *il, integer *iu, doubleCReal *abstol,
    integer *m, doubleCReal *w, doubleComplex *z__, integer *ldz,
    doubleComplex *work, doubleCReal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(zhprfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *afp, integer *ipiv, doubleComplex *
    b, integer *ldb, doubleComplex *x, integer *ldx, doubleCReal *ferr,
    doubleCReal *berr, doubleComplex *work, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpsv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhpsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *ap, doubleComplex *afp, integer *ipiv,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhptrd)(char *uplo, integer *n, doubleComplex *ap,
    doubleCReal *d__, doubleCReal *e, doubleComplex *tau, integer *info);

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
     integer *mm, integer *m, doubleComplex *work, doubleCReal *rwork,
    integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(zhseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, doubleComplex *h__, integer *ldh, doubleComplex *w,
    doubleComplex *z__, integer *ldz, doubleComplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zlabrd)(integer *m, integer *n, integer *nb,
    doubleComplex *a, integer *lda, doubleCReal *d__, doubleCReal *e,
    doubleComplex *tauq, doubleComplex *taup, doubleComplex *x, integer *
    ldx, doubleComplex *y, integer *ldy);

/* Subroutine */ int F77NAME(zlacgv)(integer *n, doubleComplex *x, integer *incx);

/* Subroutine */ int F77NAME(zlacon)(integer *n, doubleComplex *v, doubleComplex *x,
    doubleCReal *est, integer *kase);

/* Subroutine */ int F77NAME(zlacp2)(char *uplo, integer *m, integer *n, doubleCReal *
    a, integer *lda, doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlacpy)(char *uplo, integer *m, integer *n,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlacrm)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *b, integer *ldb, doubleComplex *c__,
    integer *ldc, doubleCReal *rwork);

/* Subroutine */ int F77NAME(zlacrt)(integer *n, doubleComplex *cx, integer *incx,
    doubleComplex *cy, integer *incy, doubleComplex *c__, doubleComplex *
    s);

/* Subroutine */ int F77NAME(zlaed0)(integer *qsiz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleComplex *q, integer *ldq, doubleComplex *qstore,
    integer *ldqs, doubleCReal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlaed7)(integer *n, integer *cutpnt, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, doubleCReal *d__,
    doubleComplex *q, integer *ldq, doubleCReal *rho, integer *indxq,
    doubleCReal *qstore, integer *qptr, integer *prmptr, integer *perm,
    integer *givptr, integer *givcol, doubleCReal *givnum, doubleComplex *
    work, doubleCReal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlaed8)(integer *k, integer *n, integer *qsiz,
    doubleComplex *q, integer *ldq, doubleCReal *d__, doubleCReal *rho,
    integer *cutpnt, doubleCReal *z__, doubleCReal *dlamda, doubleComplex *
    q2, integer *ldq2, doubleCReal *w, integer *indxp, integer *indx,
    integer *indxq, integer *perm, integer *givptr, integer *givcol,
    doubleCReal *givnum, integer *info);

/* Subroutine */ int F77NAME(zlaein)(logical *rightv, logical *noinit, integer *n,
    doubleComplex *h__, integer *ldh, doubleComplex *w, doubleComplex *v,
    doubleComplex *b, integer *ldb, doubleCReal *rwork, doubleCReal *eps3,
    doubleCReal *smlnum, integer *info);

/* Subroutine */ int F77NAME(zlaesy)(doubleComplex *a, doubleComplex *b,
    doubleComplex *c__, doubleComplex *rt1, doubleComplex *rt2,
    doubleComplex *evscal, doubleComplex *cs1, doubleComplex *sn1);

/* Subroutine */ int F77NAME(zlaev2)(doubleComplex *a, doubleComplex *b,
    doubleComplex *c__, doubleCReal *rt1, doubleCReal *rt2, doubleCReal *cs1,
     doubleComplex *sn1);

/* Subroutine */ int F77NAME(zlags2)(logical *upper, doubleCReal *a1, doubleComplex *
    a2, doubleCReal *a3, doubleCReal *b1, doubleComplex *b2, doubleCReal *b3,
     doubleCReal *csu, doubleComplex *snu, doubleCReal *csv, doubleComplex *
    snv, doubleCReal *csq, doubleComplex *snq);

/* Subroutine */ int F77NAME(zlagtm)(char *trans, integer *n, integer *nrhs,
    doubleCReal *alpha, doubleComplex *dl, doubleComplex *d__,
    doubleComplex *du, doubleComplex *x, integer *ldx, doubleCReal *beta,
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
    doubleCReal *sest, doubleComplex *w, doubleComplex *gamma, doubleCReal *
    sestpr, doubleComplex *s, doubleComplex *c__);

/* Subroutine */ int F77NAME(zlals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, doubleComplex *b, integer *ldb,
    doubleComplex *bx, integer *ldbx, integer *perm, integer *givptr,
    integer *givcol, integer *ldgcol, doubleCReal *givnum, integer *ldgnum,
     doubleCReal *poles, doubleCReal *difl, doubleCReal *difr, doubleCReal *
    z__, integer *k, doubleCReal *c__, doubleCReal *s, doubleCReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(zlalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, doubleComplex *b, integer *ldb, doubleComplex *bx,
    integer *ldbx, doubleCReal *u, integer *ldu, doubleCReal *vt, integer *
    k, doubleCReal *difl, doubleCReal *difr, doubleCReal *z__, doubleCReal *
    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
    perm, doubleCReal *givnum, doubleCReal *c__, doubleCReal *s, doubleCReal *
    rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlapll)(integer *n, doubleComplex *x, integer *incx,
    doubleComplex *y, integer *incy, doubleCReal *ssmin);

/* Subroutine */ int F77NAME(zlapmt)(logical *forwrd, integer *m, integer *n,
    doubleComplex *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(zlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     doubleComplex *ab, integer *ldab, doubleCReal *r__, doubleCReal *c__,
    doubleCReal *rowcnd, doubleCReal *colcnd, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqge)(integer *m, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *r__, doubleCReal *c__, doubleCReal *rowcnd,
    doubleCReal *colcnd, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqhb)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleCReal *s, doubleCReal *scond,
    doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqhe)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *s, doubleCReal *scond, doubleCReal *amax,
    char *equed);

/* Subroutine */ int F77NAME(zlaqhp)(char *uplo, integer *n, doubleComplex *ap,
    doubleCReal *s, doubleCReal *scond, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqp2)(integer *m, integer *n, integer *offset,
    doubleComplex *a, integer *lda, integer *jpvt, doubleComplex *tau,
    doubleCReal *vn1, doubleCReal *vn2, doubleComplex *work);

/* Subroutine */ int F77NAME(zlaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, doubleComplex *a, integer *lda, integer *jpvt,
    doubleComplex *tau, doubleCReal *vn1, doubleCReal *vn2, doubleComplex *
    auxv, doubleComplex *f, integer *ldf);

/* Subroutine */ int F77NAME(zlaqsb)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleCReal *s, doubleCReal *scond,
    doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqsp)(char *uplo, integer *n, doubleComplex *ap,
    doubleCReal *s, doubleCReal *scond, doubleCReal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqsy)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *s, doubleCReal *scond, doubleCReal *amax,
    char *equed);

/* Subroutine */ int F77NAME(zlar1v)(integer *n, integer *b1, integer *bn, doubleCReal
    *sigma, doubleCReal *d__, doubleCReal *l, doubleCReal *ld, doubleCReal *
    lld, doubleCReal *gersch, doubleComplex *z__, doubleCReal *ztz,
    doubleCReal *mingma, integer *r__, integer *isuppz, doubleCReal *work);

/* Subroutine */ int F77NAME(zlar2v)(integer *n, doubleComplex *x, doubleComplex *y,
    doubleComplex *z__, integer *incx, doubleCReal *c__, doubleComplex *s,
    integer *incc);

/* Subroutine */ int F77NAME(zlarcm)(integer *m, integer *n, doubleCReal *a, integer *
    lda, doubleComplex *b, integer *ldb, doubleComplex *c__, integer *ldc,
     doubleCReal *rwork);

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
    doubleComplex *y, integer *incy, doubleCReal *c__, integer *incc);

/* Subroutine */ int F77NAME(zlarnv)(integer *idist, integer *iseed, integer *n,
    doubleComplex *x);

/* Subroutine */ int F77NAME(zlarrv)(integer *n, doubleCReal *d__, doubleCReal *l,
    integer *isplit, integer *m, doubleCReal *w, integer *iblock,
    doubleCReal *gersch, doubleCReal *tol, doubleComplex *z__, integer *ldz,
     integer *isuppz, doubleCReal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlartg)(doubleComplex *f, doubleComplex *g, doubleCReal *
    cs, doubleComplex *sn, doubleComplex *r__);

/* Subroutine */ int F77NAME(zlartv)(integer *n, doubleComplex *x, integer *incx,
    doubleComplex *y, integer *incy, doubleCReal *c__, doubleComplex *s,
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
    doubleCReal *cfrom, doubleCReal *cto, integer *m, integer *n,
    doubleComplex *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(zlaset)(char *uplo, integer *m, integer *n,
    doubleComplex *alpha, doubleComplex *beta, doubleComplex *a, integer *
    lda);

/* Subroutine */ int F77NAME(zlasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, doubleCReal *c__, doubleCReal *s, doubleComplex *a,
    integer *lda);

/* Subroutine */ int F77NAME(zlassq)(integer *n, doubleComplex *x, integer *incx,
    doubleCReal *scale, doubleCReal *sumsq);

/* Subroutine */ int F77NAME(zlaswp)(integer *n, doubleComplex *a, integer *lda,
    integer *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(zlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *w,
    integer *ldw, integer *info);

/* Subroutine */ int F77NAME(zlatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, doubleComplex *ab, integer *ldab,
    doubleComplex *x, doubleCReal *scale, doubleCReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatdf)(integer *ijob, integer *n, doubleComplex *z__,
    integer *ldz, doubleComplex *rhs, doubleCReal *rdsum, doubleCReal *
    rdscal, integer *ipiv, integer *jpiv);

/* Subroutine */ int F77NAME(zlatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleComplex *ap, doubleComplex *x, doubleCReal *
    scale, doubleCReal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatrd)(char *uplo, integer *n, integer *nb,
    doubleComplex *a, integer *lda, doubleCReal *e, doubleComplex *tau,
    doubleComplex *w, integer *ldw);

/* Subroutine */ int F77NAME(zlatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doubleComplex *a, integer *lda, doubleComplex *x,
    doubleCReal *scale, doubleCReal *cnorm, integer *info);

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
    doubleComplex *ab, integer *ldab, doubleCReal *anorm, doubleCReal *
    rcond, doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpbequ)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, doubleCReal *s, doubleCReal *scond,
    doubleCReal *amax, integer *info);

/* Subroutine */ int F77NAME(zpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleComplex *ab, integer *ldab, doubleComplex *afb, integer *
    ldafb, doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
     doubleCReal *ferr, doubleCReal *berr, doubleComplex *work, doubleCReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zpbstf)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleComplex *ab, integer *ldab, doubleComplex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(zpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, doubleComplex *ab, integer *ldab, doubleComplex *afb,
    integer *ldafb, char *equed, doubleCReal *s, doubleComplex *b, integer
    *ldb, doubleComplex *x, integer *ldx, doubleCReal *rcond, doubleCReal *
    ferr, doubleCReal *berr, doubleComplex *work, doubleCReal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(zpbtf2)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbtrf)(char *uplo, integer *n, integer *kd,
    doubleComplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doubleComplex *ab, integer *ldab, doubleComplex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(zpocon)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, doubleCReal *anorm, doubleCReal *rcond, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpoequ)(integer *n, doubleComplex *a, integer *lda,
    doubleCReal *s, doubleCReal *scond, doubleCReal *amax, integer *info);

/* Subroutine */ int F77NAME(zporfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *ferr, doubleCReal *berr, doubleComplex *work, doubleCReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zposv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, char *equed, doubleCReal *s, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleCReal *rcond, doubleCReal *ferr,
    doubleCReal *berr, doubleComplex *work, doubleCReal *rwork, integer *
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
    doubleCReal *anorm, doubleCReal *rcond, doubleComplex *work, doubleCReal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(zppequ)(char *uplo, integer *n, doubleComplex *ap,
    doubleCReal *s, doubleCReal *scond, doubleCReal *amax, integer *info);

/* Subroutine */ int F77NAME(zpprfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *afp, doubleComplex *b, integer *ldb,
     doubleComplex *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zppsv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *ap, doubleComplex *afp, char *equed, doubleCReal *
    s, doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpptrf)(char *uplo, integer *n, doubleComplex *ap,
    integer *info);

/* Subroutine */ int F77NAME(zpptri)(char *uplo, integer *n, doubleComplex *ap,
    integer *info);

/* Subroutine */ int F77NAME(zpptrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zptcon)(integer *n, doubleCReal *d__, doubleComplex *e,
    doubleCReal *anorm, doubleCReal *rcond, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zptrfs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *d__, doubleComplex *e, doubleCReal *df, doubleComplex *ef,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *ferr, doubleCReal *berr, doubleComplex *work, doubleCReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zptsv)(integer *n, integer *nrhs, doubleCReal *d__,
    doubleComplex *e, doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zptsvx)(char *fact, integer *n, integer *nrhs,
    doubleCReal *d__, doubleComplex *e, doubleCReal *df, doubleComplex *ef,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpttrf)(integer *n, doubleCReal *d__, doubleComplex *e,
    integer *info);

/* Subroutine */ int F77NAME(zpttrs)(char *uplo, integer *n, integer *nrhs,
    doubleCReal *d__, doubleComplex *e, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zptts2)(integer *iuplo, integer *n, integer *nrhs,
    doubleCReal *d__, doubleComplex *e, doubleComplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zrot)(integer *n, doubleComplex *cx, integer *incx,
    doubleComplex *cy, integer *incy, doubleCReal *c__, doubleComplex *s);

/* Subroutine */ int F77NAME(zspcon)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, doubleCReal *anorm, doubleCReal *rcond, doubleComplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zspmv)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *ap, doubleComplex *x, integer *incx, doubleComplex *
    beta, doubleComplex *y, integer *incy);

/* Subroutine */ int F77NAME(zspr)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *x, integer *incx, doubleComplex *ap);

/* Subroutine */ int F77NAME(zsprfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, doubleComplex *afp, integer *ipiv, doubleComplex *
    b, integer *ldb, doubleComplex *x, integer *ldx, doubleCReal *ferr,
    doubleCReal *berr, doubleComplex *work, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zspsv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *ap, doubleComplex *afp, integer *ipiv,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsptrf)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zsptri)(char *uplo, integer *n, doubleComplex *ap,
    integer *ipiv, doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zsptrs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *ap, integer *ipiv, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zstedc)(char *compz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleComplex *z__, integer *ldz, doubleComplex *work,
    integer *lwork, doubleCReal *rwork, integer *lrwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zstein)(integer *n, doubleCReal *d__, doubleCReal *e,
    integer *m, doubleCReal *w, integer *iblock, integer *isplit,
    doubleComplex *z__, integer *ldz, doubleCReal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zsteqr)(char *compz, integer *n, doubleCReal *d__,
    doubleCReal *e, doubleComplex *z__, integer *ldz, doubleCReal *work,
    integer *info);

/* Subroutine */ int F77NAME(zsycon)(char *uplo, integer *n, doubleComplex *a,
    integer *lda, integer *ipiv, doubleCReal *anorm, doubleCReal *rcond,
    doubleComplex *work, integer *info);

/* Subroutine */ int F77NAME(zsymv)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *a, integer *lda, doubleComplex *x, integer *incx,
    doubleComplex *beta, doubleComplex *y, integer *incy);

/* Subroutine */ int F77NAME(zsyr)(char *uplo, integer *n, doubleComplex *alpha,
    doubleComplex *x, integer *incx, doubleComplex *a, integer *lda);

/* Subroutine */ int F77NAME(zsyrfs)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, doubleComplex *af, integer *ldaf,
    integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
    integer *ldx, doubleCReal *ferr, doubleCReal *berr, doubleComplex *work,
     doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsysv)(char *uplo, integer *n, integer *nrhs,
    doubleComplex *a, integer *lda, integer *ipiv, doubleComplex *b,
    integer *ldb, doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zsysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doubleComplex *a, integer *lda, doubleComplex *af, integer *
    ldaf, integer *ipiv, doubleComplex *b, integer *ldb, doubleComplex *x,
     integer *ldx, doubleCReal *rcond, doubleCReal *ferr, doubleCReal *berr,
    doubleComplex *work, integer *lwork, doubleCReal *rwork, integer *info);

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
    integer *kd, doubleComplex *ab, integer *ldab, doubleCReal *rcond,
    doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleComplex *ab, integer *ldab,
    doubleComplex *b, integer *ldb, doubleComplex *x, integer *ldx,
    doubleCReal *ferr, doubleCReal *berr, doubleComplex *work, doubleCReal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(ztbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doubleComplex *ab, integer *ldab,
    doubleComplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ztgevc)(char *side, char *howmny, logical *select,
    integer *n, doubleComplex *a, integer *lda, doubleComplex *b, integer
    *ldb, doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *
    ldvr, integer *mm, integer *m, doubleComplex *work, doubleCReal *rwork,
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
    ldz, integer *m, doubleCReal *pl, doubleCReal *pr, doubleCReal *dif,
    doubleComplex *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ztgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, doubleComplex *a,
    integer *lda, doubleComplex *b, integer *ldb, doubleCReal *tola,
    doubleCReal *tolb, doubleCReal *alpha, doubleCReal *beta, doubleComplex *
    u, integer *ldu, doubleComplex *v, integer *ldv, doubleComplex *q,
    integer *ldq, doubleComplex *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(ztgsna)(char *job, char *howmny, logical *select,
    integer *n, doubleComplex *a, integer *lda, doubleComplex *b, integer
    *ldb, doubleComplex *vl, integer *ldvl, doubleComplex *vr, integer *
    ldvr, doubleCReal *s, doubleCReal *dif, integer *mm, integer *m,
    doubleComplex *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ztgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *c__, integer *ldc, doubleComplex *d__, integer *ldd,
    doubleComplex *e, integer *lde, doubleComplex *f, integer *ldf,
    doubleCReal *scale, doubleCReal *rdsum, doubleCReal *rdscal, integer *
    info);

/* Subroutine */ int F77NAME(ztgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, doubleComplex *a, integer *lda, doubleComplex *b, integer *ldb,
    doubleComplex *c__, integer *ldc, doubleComplex *d__, integer *ldd,
    doubleComplex *e, integer *lde, doubleComplex *f, integer *ldf,
    doubleCReal *scale, doubleCReal *dif, doubleComplex *work, integer *
    lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ztpcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleComplex *ap, doubleCReal *rcond, doubleComplex *work, doubleCReal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(ztprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleComplex *ap, doubleComplex *b, integer *ldb,
    doubleComplex *x, integer *ldx, doubleCReal *ferr, doubleCReal *berr,
    doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztptri)(char *uplo, char *diag, integer *n,
    doubleComplex *ap, integer *info);

/* Subroutine */ int F77NAME(ztptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleComplex *ap, doubleComplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(ztrcon)(char *norm, char *uplo, char *diag, integer *n,
    doubleComplex *a, integer *lda, doubleCReal *rcond, doubleComplex *
    work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrevc)(char *side, char *howmny, logical *select,
    integer *n, doubleComplex *t, integer *ldt, doubleComplex *vl,
    integer *ldvl, doubleComplex *vr, integer *ldvr, integer *mm, integer
    *m, doubleComplex *work, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrexc)(char *compq, integer *n, doubleComplex *t,
    integer *ldt, doubleComplex *q, integer *ldq, integer *ifst, integer *
    ilst, integer *info);

/* Subroutine */ int F77NAME(ztrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleComplex *x, integer *ldx, doubleCReal *ferr,
    doubleCReal *berr, doubleComplex *work, doubleCReal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(ztrsen)(char *job, char *compq, logical *select, integer
    *n, doubleComplex *t, integer *ldt, doubleComplex *q, integer *ldq,
    doubleComplex *w, integer *m, doubleCReal *s, doubleCReal *sep,
    doubleComplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ztrsna)(char *job, char *howmny, logical *select,
    integer *n, doubleComplex *t, integer *ldt, doubleComplex *vl,
    integer *ldvl, doubleComplex *vr, integer *ldvr, doubleCReal *s,
    doubleCReal *sep, integer *mm, integer *m, doubleComplex *work,
    integer *ldwork, doubleCReal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, doubleComplex *a, integer *lda, doubleComplex *b,
    integer *ldb, doubleComplex *c__, integer *ldc, doubleCReal *scale,
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
