#ifndef __CLAPACK_H
#define __CLAPACK_H

#include <life/lifealg/cblas.hpp>

typedef int integer;
typedef unsigned long uinteger;
typedef char *address;
typedef short int shortint;
typedef float floatreal;
typedef double doublereal;
typedef struct { floatreal r, i; } floatcomplex;
typedef struct { doublereal r, i; } doublecomplex;
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
typedef floatreal (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* floatcomplex */ VOID (*C_fp)(...);
typedef /* Double floatcomplex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef floatreal (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* floatcomplex */ VOID (*C_fp)();
typedef /* Double floatcomplex */ VOID (*Z_fp)();
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
    nru, integer *ncc, floatreal *d__, floatreal *e, floatcomplex *vt, integer *ldvt,
    floatcomplex *u, integer *ldu, floatcomplex *c__, integer *ldc, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, floatcomplex *ab, integer *ldab, floatreal *d__,
    floatreal *e, floatcomplex *q, integer *ldq, floatcomplex *pt, integer *ldpt,
    floatcomplex *c__, integer *ldc, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     floatcomplex *ab, integer *ldab, integer *ipiv, floatreal *anorm, floatreal *rcond,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     floatcomplex *ab, integer *ldab, floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal
    *colcnd, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(cgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, floatcomplex *ab, integer *ldab, floatcomplex *afb, integer *
    ldafb, integer *ipiv, floatcomplex *b, integer *ldb, floatcomplex *x, integer *
    ldx, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, floatcomplex *ab, integer *ldab, integer *ipiv, floatcomplex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(cgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, floatcomplex *ab, integer *ldab, floatcomplex *afb,
     integer *ldafb, integer *ipiv, char *equed, floatreal *r__, floatreal *c__,
    floatcomplex *b, integer *ldb, floatcomplex *x, integer *ldx, floatreal *rcond, floatreal
    *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     floatcomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     floatcomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, floatcomplex *ab, integer *ldab, integer *ipiv, floatcomplex
    *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, floatreal *scale, integer *m, floatcomplex *v, integer *ldv,
    integer *info);

/* Subroutine */ int F77NAME(cgebal)(char *job, integer *n, floatcomplex *a, integer *lda,
    integer *ilo, integer *ihi, floatreal *scale, integer *info);

/* Subroutine */ int F77NAME(cgebd2)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatreal *d__, floatreal *e, floatcomplex *tauq, floatcomplex *taup, floatcomplex *work,
    integer *info);

/* Subroutine */ int F77NAME(cgebrd)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatreal *d__, floatreal *e, floatcomplex *tauq, floatcomplex *taup, floatcomplex *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgecon)(char *norm, integer *n, floatcomplex *a, integer *lda,
     floatreal *anorm, floatreal *rcond, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgeequ)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal *colcnd, floatreal *amax,
    integer *info);

/* Subroutine */ int F77NAME(cgees)(char *jobvs, char *sort, L_fp select, integer *n,
    floatcomplex *a, integer *lda, integer *sdim, floatcomplex *w, floatcomplex *vs,
    integer *ldvs, floatcomplex *work, integer *lwork, floatreal *rwork, logical *
    bwork, integer *info);

/* Subroutine */ int F77NAME(cgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, floatcomplex *a, integer *lda, integer *sdim, floatcomplex *
    w, floatcomplex *vs, integer *ldvs, floatreal *rconde, floatreal *rcondv, floatcomplex *
    work, integer *lwork, floatreal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cgeev)(char *jobvl, char *jobvr, integer *n, floatcomplex *a,
    integer *lda, floatcomplex *w, floatcomplex *vl, integer *ldvl, floatcomplex *vr,
    integer *ldvr, floatcomplex *work, integer *lwork, floatreal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, floatcomplex *a, integer *lda, floatcomplex *w, floatcomplex *vl,
    integer *ldvl, floatcomplex *vr, integer *ldvr, integer *ilo, integer *ihi,
     floatreal *scale, floatreal *abnrm, floatreal *rconde, floatreal *rcondv, floatcomplex *work,
    integer *lwork, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgegs)(char *jobvsl, char *jobvsr, integer *n, floatcomplex *
    a, integer *lda, floatcomplex *b, integer *ldb, floatcomplex *alpha, floatcomplex *
    beta, floatcomplex *vsl, integer *ldvsl, floatcomplex *vsr, integer *ldvsr,
    floatcomplex *work, integer *lwork, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgegv)(char *jobvl, char *jobvr, integer *n, floatcomplex *a,
    integer *lda, floatcomplex *b, integer *ldb, floatcomplex *alpha, floatcomplex *beta,
     floatcomplex *vl, integer *ldvl, floatcomplex *vr, integer *ldvr, floatcomplex *
    work, integer *lwork, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgehd2)(integer *n, integer *ilo, integer *ihi, floatcomplex *
    a, integer *lda, floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cgehrd)(integer *n, integer *ilo, integer *ihi, floatcomplex *
    a, integer *lda, floatcomplex *tau, floatcomplex *work, integer *lwork, integer
    *info);

/* Subroutine */ int F77NAME(cgelq2)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cgelqf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgels)(char *trans, integer *m, integer *n, integer *
    nrhs, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb, floatcomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgelsx)(integer *m, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, floatcomplex *b, integer *ldb, integer *jpvt, floatreal *rcond,
     integer *rank, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgelsy)(integer *m, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, floatcomplex *b, integer *ldb, integer *jpvt, floatreal *rcond,
     integer *rank, floatcomplex *work, integer *lwork, floatreal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeql2)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cgeqlf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgeqp3)(integer *m, integer *n, floatcomplex *a, integer *lda,
     integer *jpvt, floatcomplex *tau, floatcomplex *work, integer *lwork, floatreal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(cgeqpf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     integer *jpvt, floatcomplex *tau, floatcomplex *work, floatreal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(cgeqr2)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cgeqrf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgerfs)(char *trans, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, floatcomplex *af, integer *ldaf, integer *ipiv, floatcomplex *
    b, integer *ldb, floatcomplex *x, integer *ldx, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgerq2)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cgerqf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgesc2)(integer *n, floatcomplex *a, integer *lda, floatcomplex *
    rhs, integer *ipiv, integer *jpiv, floatreal *scale);

/* Subroutine */ int F77NAME(cgesv)(integer *n, integer *nrhs, floatcomplex *a, integer *
    lda, integer *ipiv, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, floatcomplex *a, integer *lda, floatcomplex *af, integer *ldaf, integer *
    ipiv, char *equed, floatreal *r__, floatreal *c__, floatcomplex *b, integer *ldb,
    floatcomplex *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgetc2)(integer *n, floatcomplex *a, integer *lda, integer *
    ipiv, integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(cgetf2)(integer *m, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgetrf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgetri)(integer *n, floatcomplex *a, integer *lda, integer *
    ipiv, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgetrs)(char *trans, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, integer *ipiv, floatcomplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, floatreal *lscale, floatreal *rscale, integer *m, floatcomplex *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(cggbal)(char *job, integer *n, floatcomplex *a, integer *lda,
    floatcomplex *b, integer *ldb, integer *ilo, integer *ihi, floatreal *lscale,
    floatreal *rscale, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(cgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, integer *n, floatcomplex *a, integer *lda, floatcomplex *b, integer *
    ldb, integer *sdim, floatcomplex *alpha, floatcomplex *beta, floatcomplex *vsl,
    integer *ldvsl, floatcomplex *vsr, integer *ldvsr, floatcomplex *work, integer *
    lwork, floatreal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, char *sense, integer *n, floatcomplex *a, integer *lda, floatcomplex *b,
     integer *ldb, integer *sdim, floatcomplex *alpha, floatcomplex *beta, floatcomplex *
    vsl, integer *ldvsl, floatcomplex *vsr, integer *ldvsr, floatreal *rconde, floatreal
    *rcondv, floatcomplex *work, integer *lwork, floatreal *rwork, integer *iwork,
    integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(cggev)(char *jobvl, char *jobvr, integer *n, floatcomplex *a,
    integer *lda, floatcomplex *b, integer *ldb, floatcomplex *alpha, floatcomplex *beta,
     floatcomplex *vl, integer *ldvl, floatcomplex *vr, integer *ldvr, floatcomplex *
    work, integer *lwork, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
     floatcomplex *alpha, floatcomplex *beta, floatcomplex *vl, integer *ldvl, floatcomplex *
    vr, integer *ldvr, integer *ilo, integer *ihi, floatreal *lscale, floatreal *
    rscale, floatreal *abnrm, floatreal *bbnrm, floatreal *rconde, floatreal *rcondv, floatcomplex
    *work, integer *lwork, floatreal *rwork, integer *iwork, logical *bwork,
    integer *info);

/* Subroutine */ int F77NAME(cggglm)(integer *n, integer *m, integer *p, floatcomplex *a,
    integer *lda, floatcomplex *b, integer *ldb, floatcomplex *d__, floatcomplex *x,
    floatcomplex *y, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
     floatcomplex *q, integer *ldq, floatcomplex *z__, integer *ldz, integer *info);

/* Subroutine */ int F77NAME(cgglse)(integer *m, integer *n, integer *p, floatcomplex *a,
    integer *lda, floatcomplex *b, integer *ldb, floatcomplex *c__, floatcomplex *d__,
    floatcomplex *x, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cggqrf)(integer *n, integer *m, integer *p, floatcomplex *a,
    integer *lda, floatcomplex *taua, floatcomplex *b, integer *ldb, floatcomplex *taub,
    floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cggrqf)(integer *m, integer *p, integer *n, floatcomplex *a,
    integer *lda, floatcomplex *taua, floatcomplex *b, integer *ldb, floatcomplex *taub,
    floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, floatcomplex *a, integer *
    lda, floatcomplex *b, integer *ldb, floatreal *alpha, floatreal *beta, floatcomplex *u,
    integer *ldu, floatcomplex *v, integer *ldv, floatcomplex *q, integer *ldq,
    floatcomplex *work, floatreal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(cggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, floatcomplex *a, integer *lda, floatcomplex *b, integer
    *ldb, floatreal *tola, floatreal *tolb, integer *k, integer *l, floatcomplex *u,
    integer *ldu, floatcomplex *v, integer *ldv, floatcomplex *q, integer *ldq,
    integer *iwork, floatreal *rwork, floatcomplex *tau, floatcomplex *work, integer *
    info);

/* Subroutine */ int F77NAME(cgtcon)(char *norm, integer *n, floatcomplex *dl, floatcomplex *
    d__, floatcomplex *du, floatcomplex *du2, integer *ipiv, floatreal *anorm, floatreal *
    rcond, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cgtrfs)(char *trans, integer *n, integer *nrhs, floatcomplex *
    dl, floatcomplex *d__, floatcomplex *du, floatcomplex *dlf, floatcomplex *df, floatcomplex *
    duf, floatcomplex *du2, integer *ipiv, floatcomplex *b, integer *ldb, floatcomplex *
    x, integer *ldx, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cgtsv)(integer *n, integer *nrhs, floatcomplex *dl, floatcomplex *
    d__, floatcomplex *du, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, floatcomplex *dl, floatcomplex *d__, floatcomplex *du, floatcomplex *dlf, floatcomplex *
    df, floatcomplex *duf, floatcomplex *du2, integer *ipiv, floatcomplex *b, integer *
    ldb, floatcomplex *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cgttrf)(integer *n, floatcomplex *dl, floatcomplex *d__, floatcomplex *
    du, floatcomplex *du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(cgttrs)(char *trans, integer *n, integer *nrhs, floatcomplex *
    dl, floatcomplex *d__, floatcomplex *du, floatcomplex *du2, integer *ipiv, floatcomplex *
    b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cgtts2)(integer *itrans, integer *n, integer *nrhs,
    floatcomplex *dl, floatcomplex *d__, floatcomplex *du, floatcomplex *du2, integer *ipiv,
    floatcomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(chbev)(char *jobz, char *uplo, integer *n, integer *kd,
    floatcomplex *ab, integer *ldab, floatreal *w, floatcomplex *z__, integer *ldz,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(chbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    floatcomplex *ab, integer *ldab, floatreal *w, floatcomplex *z__, integer *ldz,
    floatcomplex *work, integer *lwork, floatreal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, floatcomplex *ab, integer *ldab, floatcomplex *q, integer *ldq,
    floatreal *vl, floatreal *vu, integer *il, integer *iu, floatreal *abstol, integer *
    m, floatreal *w, floatcomplex *z__, integer *ldz, floatcomplex *work, floatreal *rwork,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, floatcomplex *ab, integer *ldab, floatcomplex *bb, integer *ldbb,
    floatcomplex *x, integer *ldx, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(chbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, floatcomplex *ab, integer *ldab, floatcomplex *bb, integer *ldbb,
    floatreal *w, floatcomplex *z__, integer *ldz, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, floatcomplex *ab, integer *ldab, floatcomplex *bb,
    integer *ldbb, floatcomplex *q, integer *ldq, floatreal *vl, floatreal *vu, integer *
    il, integer *iu, floatreal *abstol, integer *m, floatreal *w, floatcomplex *z__,
    integer *ldz, floatcomplex *work, floatreal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(chbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    floatcomplex *ab, integer *ldab, floatreal *d__, floatreal *e, floatcomplex *q, integer *
    ldq, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(checon)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, floatreal *anorm, floatreal *rcond, floatcomplex *work, integer *
    info);

/* Subroutine */ int F77NAME(cheev)(char *jobz, char *uplo, integer *n, floatcomplex *a,
    integer *lda, floatreal *w, floatcomplex *work, integer *lwork, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cheevd)(char *jobz, char *uplo, integer *n, floatcomplex *a,
    integer *lda, floatreal *w, floatcomplex *work, integer *lwork, floatreal *rwork,
    integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(cheevr)(char *jobz, char *range, char *uplo, integer *n,
    floatcomplex *a, integer *lda, floatreal *vl, floatreal *vu, integer *il, integer *
    iu, floatreal *abstol, integer *m, floatreal *w, floatcomplex *z__, integer *ldz,
    integer *isuppz, floatcomplex *work, integer *lwork, floatreal *rwork, integer *
    lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(cheevx)(char *jobz, char *range, char *uplo, integer *n,
    floatcomplex *a, integer *lda, floatreal *vl, floatreal *vu, integer *il, integer *
    iu, floatreal *abstol, integer *m, floatreal *w, floatcomplex *z__, integer *ldz,
    floatcomplex *work, integer *lwork, floatreal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(chegs2)(integer *itype, char *uplo, integer *n, floatcomplex *
    a, integer *lda, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chegst)(integer *itype, char *uplo, integer *n, floatcomplex *
    a, integer *lda, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chegv)(integer *itype, char *jobz, char *uplo, integer *
    n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb, floatreal *w,
    floatcomplex *work, integer *lwork, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(chegvd)(integer *itype, char *jobz, char *uplo, integer *
    n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb, floatreal *w,
    floatcomplex *work, integer *lwork, floatreal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chegvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
    floatreal *vl, floatreal *vu, integer *il, integer *iu, floatreal *abstol, integer *
    m, floatreal *w, floatcomplex *z__, integer *ldz, floatcomplex *work, integer *lwork,
     floatreal *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(cherfs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, floatcomplex *af, integer *ldaf, integer *ipiv, floatcomplex *
    b, integer *ldb, floatcomplex *x, integer *ldx, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(chesv)(char *uplo, integer *n, integer *nrhs, floatcomplex *a,
     integer *lda, integer *ipiv, floatcomplex *b, integer *ldb, floatcomplex *work,
     integer *lwork, integer *info);

/* Subroutine */ int F77NAME(chesvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatcomplex *a, integer *lda, floatcomplex *af, integer *ldaf, integer *
    ipiv, floatcomplex *b, integer *ldb, floatcomplex *x, integer *ldx, floatreal *rcond,
     floatreal *ferr, floatreal *berr, floatcomplex *work, integer *lwork, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chetf2)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(chetrd)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     floatreal *d__, floatreal *e, floatcomplex *tau, floatcomplex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(chetrf)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(chetri)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(chetrs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, integer *ipiv, floatcomplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(chgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, floatcomplex *a, integer *lda, floatcomplex *b,
    integer *ldb, floatcomplex *alpha, floatcomplex *beta, floatcomplex *q, integer *ldq,
     floatcomplex *z__, integer *ldz, floatcomplex *work, integer *lwork, floatreal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(chpcon)(char *uplo, integer *n, floatcomplex *ap, integer *
    ipiv, floatreal *anorm, floatreal *rcond, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(chpev)(char *jobz, char *uplo, integer *n, floatcomplex *ap,
    floatreal *w, floatcomplex *z__, integer *ldz, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chpevd)(char *jobz, char *uplo, integer *n, floatcomplex *ap,
    floatreal *w, floatcomplex *z__, integer *ldz, floatcomplex *work, integer *lwork,
    floatreal *rwork, integer *lrwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(chpevx)(char *jobz, char *range, char *uplo, integer *n,
    floatcomplex *ap, floatreal *vl, floatreal *vu, integer *il, integer *iu, floatreal *
    abstol, integer *m, floatreal *w, floatcomplex *z__, integer *ldz, floatcomplex *
    work, floatreal *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chpgst)(integer *itype, char *uplo, integer *n, floatcomplex *
    ap, floatcomplex *bp, integer *info);

/* Subroutine */ int F77NAME(chpgv)(integer *itype, char *jobz, char *uplo, integer *
    n, floatcomplex *ap, floatcomplex *bp, floatreal *w, floatcomplex *z__, integer *ldz,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(chpgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, floatcomplex *ap, floatcomplex *bp, floatreal *w, floatcomplex *z__, integer *ldz,
    floatcomplex *work, integer *lwork, floatreal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(chpgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, floatcomplex *ap, floatcomplex *bp, floatreal *vl, floatreal *vu,
    integer *il, integer *iu, floatreal *abstol, integer *m, floatreal *w, floatcomplex *
    z__, integer *ldz, floatcomplex *work, floatreal *rwork, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(chprfs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, floatcomplex *afp, integer *ipiv, floatcomplex *b, integer *ldb, floatcomplex *x,
     integer *ldx, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(chpsv)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, integer *ipiv, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chpsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatcomplex *ap, floatcomplex *afp, integer *ipiv, floatcomplex *b, integer *
    ldb, floatcomplex *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(chptrd)(char *uplo, integer *n, floatcomplex *ap, floatreal *d__,
    floatreal *e, floatcomplex *tau, integer *info);

/* Subroutine */ int F77NAME(chptrf)(char *uplo, integer *n, floatcomplex *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(chptri)(char *uplo, integer *n, floatcomplex *ap, integer *
    ipiv, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(chptrs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, integer *ipiv, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(chsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, floatcomplex *h__, integer *ldh, floatcomplex *w, floatcomplex *
    vl, integer *ldvl, floatcomplex *vr, integer *ldvr, integer *mm, integer *
    m, floatcomplex *work, floatreal *rwork, integer *ifaill, integer *ifailr,
    integer *info);

/* Subroutine */ int F77NAME(chseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, floatcomplex *h__, integer *ldh, floatcomplex *w, floatcomplex *z__,
    integer *ldz, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(clabrd)(integer *m, integer *n, integer *nb, floatcomplex *a,
    integer *lda, floatreal *d__, floatreal *e, floatcomplex *tauq, floatcomplex *taup,
    floatcomplex *x, integer *ldx, floatcomplex *y, integer *ldy);

/* Subroutine */ int F77NAME(clacgv)(integer *n, floatcomplex *x, integer *incx);

/* Subroutine */ int F77NAME(clacon)(integer *n, floatcomplex *v, floatcomplex *x, floatreal *est,
    integer *kase);

/* Subroutine */ int F77NAME(clacp2)(char *uplo, integer *m, integer *n, floatreal *a,
    integer *lda, floatcomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(clacpy)(char *uplo, integer *m, integer *n, floatcomplex *a,
    integer *lda, floatcomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(clacrm)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatreal *b, integer *ldb, floatcomplex *c__, integer *ldc, floatreal *rwork);

/* Subroutine */ int F77NAME(clacrt)(integer *n, floatcomplex *cx, integer *incx, floatcomplex *
    cy, integer *incy, floatcomplex *c__, floatcomplex *s);

/* Subroutine */ int F77NAME(claed0)(integer *qsiz, integer *n, floatreal *d__, floatreal *e,
    floatcomplex *q, integer *ldq, floatcomplex *qstore, integer *ldqs, floatreal *rwork,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(claed7)(integer *n, integer *cutpnt, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, floatreal *d__, floatcomplex *
    q, integer *ldq, floatreal *rho, integer *indxq, floatreal *qstore, integer *
    qptr, integer *prmptr, integer *perm, integer *givptr, integer *
    givcol, floatreal *givnum, floatcomplex *work, floatreal *rwork, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(claed8)(integer *k, integer *n, integer *qsiz, floatcomplex *
    q, integer *ldq, floatreal *d__, floatreal *rho, integer *cutpnt, floatreal *z__,
    floatreal *dlamda, floatcomplex *q2, integer *ldq2, floatreal *w, integer *indxp,
    integer *indx, integer *indxq, integer *perm, integer *givptr,
    integer *givcol, floatreal *givnum, integer *info);

/* Subroutine */ int F77NAME(claein)(logical *rightv, logical *noinit, integer *n,
    floatcomplex *h__, integer *ldh, floatcomplex *w, floatcomplex *v, floatcomplex *b,
    integer *ldb, floatreal *rwork, floatreal *eps3, floatreal *smlnum, integer *info);

/* Subroutine */ int F77NAME(claesy)(floatcomplex *a, floatcomplex *b, floatcomplex *c__, floatcomplex *
    rt1, floatcomplex *rt2, floatcomplex *evscal, floatcomplex *cs1, floatcomplex *sn1);

/* Subroutine */ int F77NAME(claev2)(floatcomplex *a, floatcomplex *b, floatcomplex *c__, floatreal *rt1,
    floatreal *rt2, floatreal *cs1, floatcomplex *sn1);

/* Subroutine */ int F77NAME(clags2)(logical *upper, floatreal *a1, floatcomplex *a2, floatreal *a3,
    floatreal *b1, floatcomplex *b2, floatreal *b3, floatreal *csu, floatcomplex *snu, floatreal *csv,
    floatcomplex *snv, floatreal *csq, floatcomplex *snq);

/* Subroutine */ int F77NAME(clagtm)(char *trans, integer *n, integer *nrhs, floatreal *
    alpha, floatcomplex *dl, floatcomplex *d__, floatcomplex *du, floatcomplex *x, integer *
    ldx, floatreal *beta, floatcomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(clahef)(char *uplo, integer *n, integer *nb, integer *kb,
     floatcomplex *a, integer *lda, integer *ipiv, floatcomplex *w, integer *ldw,
    integer *info);

/* Subroutine */ int F77NAME(clahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, floatcomplex *h__, integer *ldh, floatcomplex *w,
    integer *iloz, integer *ihiz, floatcomplex *z__, integer *ldz, integer *
    info);

/* Subroutine */ int F77NAME(clahrd)(integer *n, integer *k, integer *nb, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *t, integer *ldt, floatcomplex *y,
    integer *ldy);

/* Subroutine */ int F77NAME(claic1)(integer *job, integer *j, floatcomplex *x, floatreal *sest,
     floatcomplex *w, floatcomplex *gamma, floatreal *sestpr, floatcomplex *s, floatcomplex *c__);

/* Subroutine */ int F77NAME(clals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, floatcomplex *b, integer *ldb, floatcomplex *bx,
    integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, floatreal *givnum, integer *ldgnum, floatreal *poles, floatreal *
    difl, floatreal *difr, floatreal *z__, integer *k, floatreal *c__, floatreal *s, floatreal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(clalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, floatcomplex *b, integer *ldb, floatcomplex *bx, integer *ldbx,
    floatreal *u, integer *ldu, floatreal *vt, integer *k, floatreal *difl, floatreal *difr,
    floatreal *z__, floatreal *poles, integer *givptr, integer *givcol, integer *
    ldgcol, integer *perm, floatreal *givnum, floatreal *c__, floatreal *s, floatreal *rwork,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(clapll)(integer *n, floatcomplex *x, integer *incx, floatcomplex *
    y, integer *incy, floatreal *ssmin);

/* Subroutine */ int F77NAME(clapmt)(logical *forwrd, integer *m, integer *n, floatcomplex
    *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(claqgb)(integer *m, integer *n, integer *kl, integer *ku,
     floatcomplex *ab, integer *ldab, floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal
    *colcnd, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(claqge)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal *colcnd, floatreal *amax, char *
    equed);

/* Subroutine */ int F77NAME(claqhb)(char *uplo, integer *n, integer *kd, floatcomplex *ab,
     integer *ldab, floatreal *s, floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(claqhe)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     floatreal *s, floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(claqhp)(char *uplo, integer *n, floatcomplex *ap, floatreal *s,
    floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(claqp2)(integer *m, integer *n, integer *offset, floatcomplex
    *a, integer *lda, integer *jpvt, floatcomplex *tau, floatreal *vn1, floatreal *vn2,
    floatcomplex *work);

/* Subroutine */ int F77NAME(claqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, floatcomplex *a, integer *lda, integer *jpvt, floatcomplex *
    tau, floatreal *vn1, floatreal *vn2, floatcomplex *auxv, floatcomplex *f, integer *ldf);

/* Subroutine */ int F77NAME(claqsb)(char *uplo, integer *n, integer *kd, floatcomplex *ab,
     integer *ldab, floatreal *s, floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(claqsp)(char *uplo, integer *n, floatcomplex *ap, floatreal *s,
    floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(claqsy)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     floatreal *s, floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(clar1v)(integer *n, integer *b1, integer *bn, floatreal *
    sigma, floatreal *d__, floatreal *l, floatreal *ld, floatreal *lld, floatreal *gersch, floatcomplex
    *z__, floatreal *ztz, floatreal *mingma, integer *r__, integer *isuppz, floatreal *
    work);

/* Subroutine */ int F77NAME(clar2v)(integer *n, floatcomplex *x, floatcomplex *y, floatcomplex *z__,
     integer *incx, floatreal *c__, floatcomplex *s, integer *incc);

/* Subroutine */ int F77NAME(clarcm)(integer *m, integer *n, floatreal *a, integer *lda,
    floatcomplex *b, integer *ldb, floatcomplex *c__, integer *ldc, floatreal *rwork);

/* Subroutine */ int F77NAME(clarf)(char *side, integer *m, integer *n, floatcomplex *v,
    integer *incv, floatcomplex *tau, floatcomplex *c__, integer *ldc, floatcomplex *
    work);

/* Subroutine */ int F77NAME(clarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, floatcomplex *v, integer *ldv,
    floatcomplex *t, integer *ldt, floatcomplex *c__, integer *ldc, floatcomplex *work,
    integer *ldwork);

/* Subroutine */ int F77NAME(clarfg)(integer *n, floatcomplex *alpha, floatcomplex *x, integer *
    incx, floatcomplex *tau);

/* Subroutine */ int F77NAME(clarft)(char *direct, char *storev, integer *n, integer *
    k, floatcomplex *v, integer *ldv, floatcomplex *tau, floatcomplex *t, integer *ldt);

/* Subroutine */ int F77NAME(clarfx)(char *side, integer *m, integer *n, floatcomplex *v,
    floatcomplex *tau, floatcomplex *c__, integer *ldc, floatcomplex *work);

/* Subroutine */ int F77NAME(clargv)(integer *n, floatcomplex *x, integer *incx, floatcomplex *
    y, integer *incy, floatreal *c__, integer *incc);

/* Subroutine */ int F77NAME(clarnv)(integer *idist, integer *iseed, integer *n,
    floatcomplex *x);

/* Subroutine */ int F77NAME(clarrv)(integer *n, floatreal *d__, floatreal *l, integer *isplit,
    integer *m, floatreal *w, integer *iblock, floatreal *gersch, floatreal *tol,
    floatcomplex *z__, integer *ldz, integer *isuppz, floatreal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(clartg)(floatcomplex *f, floatcomplex *g, floatreal *cs, floatcomplex *sn,
    floatcomplex *r__);

/* Subroutine */ int F77NAME(clartv)(integer *n, floatcomplex *x, integer *incx, floatcomplex *
    y, integer *incy, floatreal *c__, floatcomplex *s, integer *incc);

/* Subroutine */ int F77NAME(clarz)(char *side, integer *m, integer *n, integer *l,
    floatcomplex *v, integer *incv, floatcomplex *tau, floatcomplex *c__, integer *ldc,
    floatcomplex *work);

/* Subroutine */ int F77NAME(clarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, floatcomplex *v,
    integer *ldv, floatcomplex *t, integer *ldt, floatcomplex *c__, integer *ldc,
    floatcomplex *work, integer *ldwork);

/* Subroutine */ int F77NAME(clarzt)(char *direct, char *storev, integer *n, integer *
    k, floatcomplex *v, integer *ldv, floatcomplex *tau, floatcomplex *t, integer *ldt);

/* Subroutine */ int F77NAME(clascl)(char *type__, integer *kl, integer *ku, floatreal *
    cfrom, floatreal *cto, integer *m, integer *n, floatcomplex *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(claset)(char *uplo, integer *m, integer *n, floatcomplex *
    alpha, floatcomplex *beta, floatcomplex *a, integer *lda);

/* Subroutine */ int F77NAME(clasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, floatreal *c__, floatreal *s, floatcomplex *a, integer *lda);

/* Subroutine */ int F77NAME(classq)(integer *n, floatcomplex *x, integer *incx, floatreal *
    scale, floatreal *sumsq);

/* Subroutine */ int F77NAME(claswp)(integer *n, floatcomplex *a, integer *lda, integer *
    k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(clasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     floatcomplex *a, integer *lda, integer *ipiv, floatcomplex *w, integer *ldw,
    integer *info);

/* Subroutine */ int F77NAME(clatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, floatcomplex *ab, integer *ldab, floatcomplex *
    x, floatreal *scale, floatreal *cnorm, integer *info);

/* Subroutine */ int F77NAME(clatdf)(integer *ijob, integer *n, floatcomplex *z__, integer
    *ldz, floatcomplex *rhs, floatreal *rdsum, floatreal *rdscal, integer *ipiv, integer
    *jpiv);

/* Subroutine */ int F77NAME(clatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, floatcomplex *ap, floatcomplex *x, floatreal *scale, floatreal *cnorm,
     integer *info);

/* Subroutine */ int F77NAME(clatrd)(char *uplo, integer *n, integer *nb, floatcomplex *a,
    integer *lda, floatreal *e, floatcomplex *tau, floatcomplex *w, integer *ldw);

/* Subroutine */ int F77NAME(clatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, floatcomplex *a, integer *lda, floatcomplex *x, floatreal *scale,
     floatreal *cnorm, integer *info);

/* Subroutine */ int F77NAME(clatrz)(integer *m, integer *n, integer *l, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work);

/* Subroutine */ int F77NAME(clatzm)(char *side, integer *m, integer *n, floatcomplex *v,
    integer *incv, floatcomplex *tau, floatcomplex *c1, floatcomplex *c2, integer *ldc,
    floatcomplex *work);

/* Subroutine */ int F77NAME(clauu2)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(clauum)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpbcon)(char *uplo, integer *n, integer *kd, floatcomplex *ab,
     integer *ldab, floatreal *anorm, floatreal *rcond, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cpbequ)(char *uplo, integer *n, integer *kd, floatcomplex *ab,
     integer *ldab, floatreal *s, floatreal *scond, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(cpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, floatcomplex *ab, integer *ldab, floatcomplex *afb, integer *ldafb,
    floatcomplex *b, integer *ldb, floatcomplex *x, integer *ldx, floatreal *ferr, floatreal *
    berr, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpbstf)(char *uplo, integer *n, integer *kd, floatcomplex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, floatcomplex *ab, integer *ldab, floatcomplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, floatcomplex *ab, integer *ldab, floatcomplex *afb, integer *
    ldafb, char *equed, floatreal *s, floatcomplex *b, integer *ldb, floatcomplex *x,
    integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr, floatcomplex *work,
    floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpbtf2)(char *uplo, integer *n, integer *kd, floatcomplex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbtrf)(char *uplo, integer *n, integer *kd, floatcomplex *ab,
     integer *ldab, integer *info);

/* Subroutine */ int F77NAME(cpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, floatcomplex *ab, integer *ldab, floatcomplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(cpocon)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     floatreal *anorm, floatreal *rcond, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpoequ)(integer *n, floatcomplex *a, integer *lda, floatreal *s,
    floatreal *scond, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(cporfs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, floatcomplex *af, integer *ldaf, floatcomplex *b, integer *ldb,
     floatcomplex *x, integer *ldx, floatreal *ferr, floatreal *berr, floatcomplex *work,
    floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cposv)(char *uplo, integer *n, integer *nrhs, floatcomplex *a,
     integer *lda, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatcomplex *a, integer *lda, floatcomplex *af, integer *ldaf, char *
    equed, floatreal *s, floatcomplex *b, integer *ldb, floatcomplex *x, integer *ldx,
    floatreal *rcond, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cpotf2)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotrf)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotri)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *info);

/* Subroutine */ int F77NAME(cpotrs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cppcon)(char *uplo, integer *n, floatcomplex *ap, floatreal *anorm,
     floatreal *rcond, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cppequ)(char *uplo, integer *n, floatcomplex *ap, floatreal *s,
    floatreal *scond, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(cpprfs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, floatcomplex *afp, floatcomplex *b, integer *ldb, floatcomplex *x, integer *ldx,
    floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cppsv)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatcomplex *ap, floatcomplex *afp, char *equed, floatreal *s, floatcomplex *b,
    integer *ldb, floatcomplex *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal
    *berr, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpptrf)(char *uplo, integer *n, floatcomplex *ap, integer *
    info);

/* Subroutine */ int F77NAME(cpptri)(char *uplo, integer *n, floatcomplex *ap, integer *
    info);

/* Subroutine */ int F77NAME(cpptrs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptcon)(integer *n, floatreal *d__, floatcomplex *e, floatreal *anorm,
    floatreal *rcond, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cptrfs)(char *uplo, integer *n, integer *nrhs, floatreal *d__,
     floatcomplex *e, floatreal *df, floatcomplex *ef, floatcomplex *b, integer *ldb, floatcomplex
    *x, integer *ldx, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cptsv)(integer *n, integer *nrhs, floatreal *d__, floatcomplex *e,
    floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptsvx)(char *fact, integer *n, integer *nrhs, floatreal *d__,
     floatcomplex *e, floatreal *df, floatcomplex *ef, floatcomplex *b, integer *ldb, floatcomplex
    *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr, floatcomplex *work,
    floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(cpttrf)(integer *n, floatreal *d__, floatcomplex *e, integer *info);

/* Subroutine */ int F77NAME(cpttrs)(char *uplo, integer *n, integer *nrhs, floatreal *d__,
     floatcomplex *e, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cptts2)(integer *iuplo, integer *n, integer *nrhs, floatreal *
    d__, floatcomplex *e, floatcomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(crot)(integer *n, floatcomplex *cx, integer *incx, floatcomplex *
    cy, integer *incy, floatreal *c__, floatcomplex *s);

/* Subroutine */ int F77NAME(cspcon)(char *uplo, integer *n, floatcomplex *ap, integer *
    ipiv, floatreal *anorm, floatreal *rcond, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cspmv)(char *uplo, integer *n, floatcomplex *alpha, floatcomplex *
    ap, floatcomplex *x, integer *incx, floatcomplex *beta, floatcomplex *y, integer *
    incy);

/* Subroutine */ int F77NAME(cspr)(char *uplo, integer *n, floatcomplex *alpha, floatcomplex *x,
     integer *incx, floatcomplex *ap);

/* Subroutine */ int F77NAME(csprfs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, floatcomplex *afp, integer *ipiv, floatcomplex *b, integer *ldb, floatcomplex *x,
     integer *ldx, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(cspsv)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, integer *ipiv, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(cspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatcomplex *ap, floatcomplex *afp, integer *ipiv, floatcomplex *b, integer *
    ldb, floatcomplex *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(csptrf)(char *uplo, integer *n, floatcomplex *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(csptri)(char *uplo, integer *n, floatcomplex *ap, integer *
    ipiv, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(csptrs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    ap, integer *ipiv, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(csrot)(integer *n, floatcomplex *cx, integer *incx, floatcomplex *
    cy, integer *incy, floatreal *c__, floatreal *s);

/* Subroutine */ int F77NAME(csrscl)(integer *n, floatreal *sa, floatcomplex *sx, integer *incx);

/* Subroutine */ int F77NAME(cstedc)(char *compz, integer *n, floatreal *d__, floatreal *e,
    floatcomplex *z__, integer *ldz, floatcomplex *work, integer *lwork, floatreal *
    rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(cstein)(integer *n, floatreal *d__, floatreal *e, integer *m, floatreal
    *w, integer *iblock, integer *isplit, floatcomplex *z__, integer *ldz,
    floatreal *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(csteqr)(char *compz, integer *n, floatreal *d__, floatreal *e,
    floatcomplex *z__, integer *ldz, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(csycon)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, floatreal *anorm, floatreal *rcond, floatcomplex *work, integer *
    info);

/* Subroutine */ int F77NAME(csymv)(char *uplo, integer *n, floatcomplex *alpha, floatcomplex *
    a, integer *lda, floatcomplex *x, integer *incx, floatcomplex *beta, floatcomplex *y,
     integer *incy);

/* Subroutine */ int F77NAME(csyr)(char *uplo, integer *n, floatcomplex *alpha, floatcomplex *x,
     integer *incx, floatcomplex *a, integer *lda);

/* Subroutine */ int F77NAME(csyrfs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, floatcomplex *af, integer *ldaf, integer *ipiv, floatcomplex *
    b, integer *ldb, floatcomplex *x, integer *ldx, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(csysv)(char *uplo, integer *n, integer *nrhs, floatcomplex *a,
     integer *lda, integer *ipiv, floatcomplex *b, integer *ldb, floatcomplex *work,
     integer *lwork, integer *info);

/* Subroutine */ int F77NAME(csysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatcomplex *a, integer *lda, floatcomplex *af, integer *ldaf, integer *
    ipiv, floatcomplex *b, integer *ldb, floatcomplex *x, integer *ldx, floatreal *rcond,
     floatreal *ferr, floatreal *berr, floatcomplex *work, integer *lwork, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(csytf2)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(csytrf)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(csytri)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     integer *ipiv, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(csytrs)(char *uplo, integer *n, integer *nrhs, floatcomplex *
    a, integer *lda, integer *ipiv, floatcomplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(ctbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, floatcomplex *ab, integer *ldab, floatreal *rcond, floatcomplex *work,
    floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, floatcomplex *ab, integer *ldab, floatcomplex *b,
    integer *ldb, floatcomplex *x, integer *ldx, floatreal *ferr, floatreal *berr,
    floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, floatcomplex *ab, integer *ldab, floatcomplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ctgevc)(char *side, char *howmny, logical *select,
    integer *n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
    floatcomplex *vl, integer *ldvl, floatcomplex *vr, integer *ldvr, integer *mm,
    integer *m, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctgex2)(logical *wantq, logical *wantz, integer *n,
    floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb, floatcomplex *q,
    integer *ldq, floatcomplex *z__, integer *ldz, integer *j1, integer *info);

/* Subroutine */ int F77NAME(ctgexc)(logical *wantq, logical *wantz, integer *n,
    floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb, floatcomplex *q,
    integer *ldq, floatcomplex *z__, integer *ldz, integer *ifst, integer *
    ilst, integer *info);

/* Subroutine */ int F77NAME(ctgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, floatcomplex *a, integer *lda, floatcomplex *b,
    integer *ldb, floatcomplex *alpha, floatcomplex *beta, floatcomplex *q, integer *ldq,
     floatcomplex *z__, integer *ldz, integer *m, floatreal *pl, floatreal *pr, floatreal *
    dif, floatcomplex *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ctgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, floatcomplex *a, integer *
    lda, floatcomplex *b, integer *ldb, floatreal *tola, floatreal *tolb, floatreal *alpha,
    floatreal *beta, floatcomplex *u, integer *ldu, floatcomplex *v, integer *ldv,
    floatcomplex *q, integer *ldq, floatcomplex *work, integer *ncycle, integer *
    info);

/* Subroutine */ int F77NAME(ctgsna)(char *job, char *howmny, logical *select,
    integer *n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
    floatcomplex *vl, integer *ldvl, floatcomplex *vr, integer *ldvr, floatreal *s, floatreal
    *dif, integer *mm, integer *m, floatcomplex *work, integer *lwork, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(ctgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb, floatcomplex *c__,
    integer *ldc, floatcomplex *d__, integer *ldd, floatcomplex *e, integer *lde,
    floatcomplex *f, integer *ldf, floatreal *scale, floatreal *rdsum, floatreal *rdscal,
    integer *info);

/* Subroutine */ int F77NAME(ctgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb, floatcomplex *c__,
    integer *ldc, floatcomplex *d__, integer *ldd, floatcomplex *e, integer *lde,
    floatcomplex *f, integer *ldf, floatreal *scale, floatreal *dif, floatcomplex *work,
    integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ctpcon)(char *norm, char *uplo, char *diag, integer *n,
    floatcomplex *ap, floatreal *rcond, floatcomplex *work, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatcomplex *ap, floatcomplex *b, integer *ldb, floatcomplex *x,
    integer *ldx, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(ctptri)(char *uplo, char *diag, integer *n, floatcomplex *ap,
    integer *info);

/* Subroutine */ int F77NAME(ctptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatcomplex *ap, floatcomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ctrcon)(char *norm, char *uplo, char *diag, integer *n,
    floatcomplex *a, integer *lda, floatreal *rcond, floatcomplex *work, floatreal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(ctrevc)(char *side, char *howmny, logical *select,
    integer *n, floatcomplex *t, integer *ldt, floatcomplex *vl, integer *ldvl,
    floatcomplex *vr, integer *ldvr, integer *mm, integer *m, floatcomplex *work,
    floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrexc)(char *compq, integer *n, floatcomplex *t, integer *
    ldt, floatcomplex *q, integer *ldq, integer *ifst, integer *ilst, integer *
    info);

/* Subroutine */ int F77NAME(ctrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
    floatcomplex *x, integer *ldx, floatreal *ferr, floatreal *berr, floatcomplex *work, floatreal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrsen)(char *job, char *compq, logical *select, integer
    *n, floatcomplex *t, integer *ldt, floatcomplex *q, integer *ldq, floatcomplex *w,
    integer *m, floatreal *s, floatreal *sep, floatcomplex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(ctrsna)(char *job, char *howmny, logical *select,
    integer *n, floatcomplex *t, integer *ldt, floatcomplex *vl, integer *ldvl,
    floatcomplex *vr, integer *ldvr, floatreal *s, floatreal *sep, integer *mm, integer *
    m, floatcomplex *work, integer *ldwork, floatreal *rwork, integer *info);

/* Subroutine */ int F77NAME(ctrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
    floatcomplex *c__, integer *ldc, floatreal *scale, integer *info);

/* Subroutine */ int F77NAME(ctrti2)(char *uplo, char *diag, integer *n, floatcomplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(ctrtri)(char *uplo, char *diag, integer *n, floatcomplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(ctrtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatcomplex *a, integer *lda, floatcomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(ctzrqf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, integer *info);

/* Subroutine */ int F77NAME(ctzrzf)(integer *m, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cung2l)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cung2r)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cungbr)(char *vect, integer *m, integer *n, integer *k,
    floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(cunghr)(integer *n, integer *ilo, integer *ihi, floatcomplex *
    a, integer *lda, floatcomplex *tau, floatcomplex *work, integer *lwork, integer
    *info);

/* Subroutine */ int F77NAME(cungl2)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cunglq)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungql)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungqr)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungr2)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cungrq)(integer *m, integer *n, integer *k, floatcomplex *a,
    integer *lda, floatcomplex *tau, floatcomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cungtr)(char *uplo, integer *n, floatcomplex *a, integer *lda,
     floatcomplex *tau, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cunm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cunmbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, floatcomplex *a, integer *lda, floatcomplex *tau,
    floatcomplex *c__, integer *ldc, floatcomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cunmhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, floatcomplex *a, integer *lda, floatcomplex *tau,
    floatcomplex *c__, integer *ldc, floatcomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cunml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cunmlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmql)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cunmr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, floatcomplex *a, integer *lda, floatcomplex *tau,
    floatcomplex *c__, integer *ldc, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cunmrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cunmrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, floatcomplex *a, integer *lda, floatcomplex *tau,
    floatcomplex *c__, integer *ldc, floatcomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(cunmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, floatcomplex *a, integer *lda, floatcomplex *tau, floatcomplex *c__,
    integer *ldc, floatcomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(cupgtr)(char *uplo, integer *n, floatcomplex *ap, floatcomplex *
    tau, floatcomplex *q, integer *ldq, floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(cupmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, floatcomplex *ap, floatcomplex *tau, floatcomplex *c__, integer *ldc,
    floatcomplex *work, integer *info);

/* Subroutine */ int F77NAME(dbdsdc)(char *uplo, char *compq, integer *n, doublereal *
    d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt,
    integer *ldvt, doublereal *q, integer *iq, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt,
    integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *
    ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(ddisna)(char *job, integer *m, integer *n, doublereal *
    d__, doublereal *sep, integer *info);

/* Subroutine */ int F77NAME(dgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *
    d__, doublereal *e, doublereal *q, integer *ldq, doublereal *pt,
    integer *ldpt, doublereal *c__, integer *ldc, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     doublereal *ab, integer *ldab, integer *ipiv, doublereal *anorm,
    doublereal *rcond, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__,
    doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
    info);

/* Subroutine */ int F77NAME(dgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb,
    integer *ldafb, integer *ipiv, doublereal *b, integer *ldb,
    doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, doublereal *ab, integer *ldab,
    doublereal *afb, integer *ldafb, integer *ipiv, char *equed,
    doublereal *r__, doublereal *c__, doublereal *b, integer *ldb,
    doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
    doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     doublereal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     doublereal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv,
    doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
    ldv, integer *info);

/* Subroutine */ int F77NAME(dgebal)(char *job, integer *n, doublereal *a, integer *
    lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);

/* Subroutine */ int F77NAME(dgebd2)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
    taup, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dgebrd)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
    taup, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgecon)(char *norm, integer *n, doublereal *a, integer *
    lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgeequ)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal
    *colcnd, doublereal *amax, integer *info);

/* Subroutine */ int F77NAME(dgees)(char *jobvs, char *sort, L_fp select, integer *n,
    doublereal *a, integer *lda, integer *sdim, doublereal *wr,
    doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work,
    integer *lwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, doublereal *a, integer *lda, integer *sdim,
    doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs,
    doublereal *rconde, doublereal *rcondv, doublereal *work, integer *
    lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dgeev)(char *jobvl, char *jobvr, integer *n, doublereal *
    a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl,
    integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doublereal *a, integer *lda, doublereal *wr,
    doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr,
    integer *ldvr, integer *ilo, integer *ihi, doublereal *scale,
    doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal
    *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgegs)(char *jobvsl, char *jobvsr, integer *n,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
    alphar, doublereal *alphai, doublereal *beta, doublereal *vsl,
    integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgegv)(char *jobvl, char *jobvr, integer *n, doublereal *
    a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar,
    doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl,
    doublereal *vr, integer *ldvr, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgehd2)(integer *n, integer *ilo, integer *ihi,
    doublereal *a, integer *lda, doublereal *tau, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dgehrd)(integer *n, integer *ilo, integer *ihi,
    doublereal *a, integer *lda, doublereal *tau, doublereal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgelq2)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dgelqf)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgels)(char *trans, integer *m, integer *n, integer *
    nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb,
    doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgelsd)(integer *m, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
    s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgelss)(integer *m, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
    s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(dgelsx)(integer *m, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
    jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
    info);

/* Subroutine */ int F77NAME(dgelsy)(integer *m, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
    jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(dgeql2)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqlf)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgeqp3)(integer *m, integer *n, doublereal *a, integer *
    lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(dgeqpf)(integer *m, integer *n, doublereal *a, integer *
    lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqr2)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dgeqrf)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgerfs)(char *trans, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
    ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx,
    doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dgerq2)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dgerqf)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgesc2)(integer *n, doublereal *a, integer *lda,
    doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

/* Subroutine */ int F77NAME(dgesdd)(char *jobz, integer *m, integer *n, doublereal *
    a, integer *lda, doublereal *s, doublereal *u, integer *ldu,
    doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dgesv)(integer *n, integer *nrhs, doublereal *a, integer
    *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgesvd)(char *jobu, char *jobvt, integer *m, integer *n,
    doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
    ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf,
    integer *ipiv, char *equed, doublereal *r__, doublereal *c__,
    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgetc2)(integer *n, doublereal *a, integer *lda, integer
    *ipiv, integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(dgetf2)(integer *m, integer *n, doublereal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgetrf)(integer *m, integer *n, doublereal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgetri)(integer *n, doublereal *a, integer *lda, integer
    *ipiv, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dgetrs)(char *trans, integer *n, integer *nrhs,
    doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doublereal *lscale, doublereal *rscale, integer *m,
    doublereal *v, integer *ldv, integer *info);

/* Subroutine */ int F77NAME(dggbal)(char *job, integer *n, doublereal *a, integer *
    lda, doublereal *b, integer *ldb, integer *ilo, integer *ihi,
    doublereal *lscale, doublereal *rscale, doublereal *work, integer *
    info);

/* Subroutine */ int F77NAME(dgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, integer *n, doublereal *a, integer *lda, doublereal *b,
    integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai,
    doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr,
    integer *ldvsr, doublereal *work, integer *lwork, logical *bwork,
    integer *info);

/* Subroutine */ int F77NAME(dggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, char *sense, integer *n, doublereal *a, integer *lda,
    doublereal *b, integer *ldb, integer *sdim, doublereal *alphar,
    doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl,
     doublereal *vsr, integer *ldvsr, doublereal *rconde, doublereal *
    rcondv, doublereal *work, integer *lwork, integer *iwork, integer *
    liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(dggev)(char *jobvl, char *jobvr, integer *n, doublereal *
    a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar,
    doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl,
    doublereal *vr, integer *ldvr, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doublereal *a, integer *lda, doublereal *b,
    integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
    beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr,
    integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale,
    doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
    rcondv, doublereal *work, integer *lwork, integer *iwork, logical *
    bwork, integer *info);

/* Subroutine */ int F77NAME(dggglm)(integer *n, integer *m, integer *p, doublereal *
    a, integer *lda, doublereal *b, integer *ldb, doublereal *d__,
    doublereal *x, doublereal *y, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b,
    integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
    ldz, integer *info);

/* Subroutine */ int F77NAME(dgglse)(integer *m, integer *n, integer *p, doublereal *
    a, integer *lda, doublereal *b, integer *ldb, doublereal *c__,
    doublereal *d__, doublereal *x, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dggqrf)(integer *n, integer *m, integer *p, doublereal *
    a, integer *lda, doublereal *taua, doublereal *b, integer *ldb,
    doublereal *taub, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dggrqf)(integer *m, integer *p, integer *n, doublereal *
    a, integer *lda, doublereal *taua, doublereal *b, integer *ldb,
    doublereal *taub, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, doublereal *a,
    integer *lda, doublereal *b, integer *ldb, doublereal *alpha,
    doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer
    *ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, doublereal *a, integer *lda, doublereal *b,
    integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer
    *l, doublereal *u, integer *ldu, doublereal *v, integer *ldv,
    doublereal *q, integer *ldq, integer *iwork, doublereal *tau,
    doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dgtcon)(char *norm, integer *n, doublereal *dl,
    doublereal *d__, doublereal *du, doublereal *du2, integer *ipiv,
    doublereal *anorm, doublereal *rcond, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgtrfs)(char *trans, integer *n, integer *nrhs,
    doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf,
    doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv,
    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
    ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dgtsv)(integer *n, integer *nrhs, doublereal *dl,
    doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer
    *info);

/* Subroutine */ int F77NAME(dgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *
    dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv,
    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dgttrf)(integer *n, doublereal *dl, doublereal *d__,
    doublereal *du, doublereal *du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dgttrs)(char *trans, integer *n, integer *nrhs,
    doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2,
    integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dgtts2)(integer *itrans, integer *n, integer *nrhs,
    doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2,
    integer *ipiv, doublereal *b, integer *ldb);

/* Subroutine */ int F77NAME(dhgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
    b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
    beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz,
    doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dhsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, doublereal *h__, integer *ldh, doublereal *wr,
    doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr,
    integer *ldvr, integer *mm, integer *m, doublereal *work, integer *
    ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(dhseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, doublereal *h__, integer *ldh, doublereal *wr,
    doublereal *wi, doublereal *z__, integer *ldz, doublereal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dlabad)(doublereal *small, doublereal *large);

/* Subroutine */ int F77NAME(dlabrd)(integer *m, integer *n, integer *nb, doublereal *
    a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq,
    doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer
    *ldy);

/* Subroutine */ int F77NAME(dlacon)(integer *n, doublereal *v, doublereal *x,
    integer *isgn, doublereal *est, integer *kase);

/* Subroutine */ int F77NAME(dlacpy)(char *uplo, integer *m, integer *n, doublereal *
    a, integer *lda, doublereal *b, integer *ldb);

/* Subroutine */ int F77NAME(dladiv)(doublereal *a, doublereal *b, doublereal *c__,
    doublereal *d__, doublereal *p, doublereal *q);

/* Subroutine */ int F77NAME(dlae2)(doublereal *a, doublereal *b, doublereal *c__,
    doublereal *rt1, doublereal *rt2);

/* Subroutine */ int F77NAME(dlaebz)(integer *ijob, integer *nitmax, integer *n,
    integer *mmax, integer *minp, integer *nbmin, doublereal *abstol,
    doublereal *reltol, doublereal *pivmin, doublereal *d__, doublereal *
    e, doublereal *e2, integer *nval, doublereal *ab, doublereal *c__,
    integer *mout, integer *nab, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlaed0)(integer *icompq, integer *qsiz, integer *n,
    doublereal *d__, doublereal *e, doublereal *q, integer *ldq,
    doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlaed1)(integer *n, doublereal *d__, doublereal *q,
    integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlaed2)(integer *k, integer *n, integer *n1, doublereal *
    d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho,
    doublereal *z__, doublereal *dlamda, doublereal *w, doublereal *q2,
    integer *indx, integer *indxc, integer *indxp, integer *coltyp,
    integer *info);

/* Subroutine */ int F77NAME(dlaed3)(integer *k, integer *n, integer *n1, doublereal *
    d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda,
     doublereal *q2, integer *indx, integer *ctot, doublereal *w,
    doublereal *s, integer *info);

/* Subroutine */ int F77NAME(dlaed4)(integer *n, integer *i__, doublereal *d__,
    doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam,
     integer *info);

/* Subroutine */ int F77NAME(dlaed5)(integer *i__, doublereal *d__, doublereal *z__,
    doublereal *delta, doublereal *rho, doublereal *dlam);

/* Subroutine */ int F77NAME(dlaed6)(integer *kniter, logical *orgati, doublereal *
    rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal *
    tau, integer *info);

/* Subroutine */ int F77NAME(dlaed7)(integer *icompq, integer *n, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__,
    doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer
    *cutpnt, doublereal *qstore, integer *qptr, integer *prmptr, integer *
    perm, integer *givptr, integer *givcol, doublereal *givnum,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlaed8)(integer *icompq, integer *k, integer *n, integer
    *qsiz, doublereal *d__, doublereal *q, integer *ldq, integer *indxq,
    doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda,
     doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer
    *givptr, integer *givcol, doublereal *givnum, integer *indxp, integer
    *indx, integer *info);

/* Subroutine */ int F77NAME(dlaed9)(integer *k, integer *kstart, integer *kstop,
    integer *n, doublereal *d__, doublereal *q, integer *ldq, doublereal *
    rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds,
    integer *info);

/* Subroutine */ int F77NAME(dlaeda)(integer *n, integer *tlvls, integer *curlvl,
    integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
    integer *givcol, doublereal *givnum, doublereal *q, integer *qptr,
    doublereal *z__, doublereal *ztemp, integer *info);

/* Subroutine */ int F77NAME(dlaein)(logical *rightv, logical *noinit, integer *n,
    doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi,
    doublereal *vr, doublereal *vi, doublereal *b, integer *ldb,
    doublereal *work, doublereal *eps3, doublereal *smlnum, doublereal *
    bignum, integer *info);

/* Subroutine */ int F77NAME(dlaev2)(doublereal *a, doublereal *b, doublereal *c__,
    doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);

/* Subroutine */ int F77NAME(dlaexc)(logical *wantq, integer *n, doublereal *t,
    integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1,
    integer *n2, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dlag2)(doublereal *a, integer *lda, doublereal *b,
    integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *
    scale2, doublereal *wr1, doublereal *wr2, doublereal *wi);

/* Subroutine */ int F77NAME(dlags2)(logical *upper, doublereal *a1, doublereal *a2,
    doublereal *a3, doublereal *b1, doublereal *b2, doublereal *b3,
    doublereal *csu, doublereal *snu, doublereal *csv, doublereal *snv,
    doublereal *csq, doublereal *snq);

/* Subroutine */ int F77NAME(dlagtf)(integer *n, doublereal *a, doublereal *lambda,
    doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__,
    integer *in, integer *info);

/* Subroutine */ int F77NAME(dlagtm)(char *trans, integer *n, integer *nrhs,
    doublereal *alpha, doublereal *dl, doublereal *d__, doublereal *du,
    doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer
    *ldb);

/* Subroutine */ int F77NAME(dlagts)(integer *job, integer *n, doublereal *a,
    doublereal *b, doublereal *c__, doublereal *d__, integer *in,
    doublereal *y, doublereal *tol, integer *info);

/* Subroutine */ int F77NAME(dlagv2)(doublereal *a, integer *lda, doublereal *b,
    integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
    beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
    snr);

/* Subroutine */ int F77NAME(dlahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal
    *wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__,
    integer *ldz, integer *info);

/* Subroutine */ int F77NAME(dlahrd)(integer *n, integer *k, integer *nb, doublereal *
    a, integer *lda, doublereal *tau, doublereal *t, integer *ldt,
    doublereal *y, integer *ldy);

/* Subroutine */ int F77NAME(dlaic1)(integer *job, integer *j, doublereal *x,
    doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
    sestpr, doublereal *s, doublereal *c__);

/* Subroutine */ int F77NAME(dlaln2)(logical *ltrans, integer *na, integer *nw,
    doublereal *smin, doublereal *ca, doublereal *a, integer *lda,
    doublereal *d1, doublereal *d2, doublereal *b, integer *ldb,
    doublereal *wr, doublereal *wi, doublereal *x, integer *ldx,
    doublereal *scale, doublereal *xnorm, integer *info);

/* Subroutine */ int F77NAME(dlals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, doublereal *b, integer *ldb, doublereal
    *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *
    poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *
    k, doublereal *c__, doublereal *s, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dlalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *
    ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k,
    doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
    perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlalsd)(char *uplo, integer *smlsiz, integer *n, integer
    *nrhs, doublereal *d__, doublereal *e, doublereal *b, integer *ldb,
    doublereal *rcond, integer *rank, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dlamc1)(integer *beta, integer *t, logical *rnd, logical
    *ieee1);

/* Subroutine */ int F77NAME(dlamc2)(integer *beta, integer *t, logical *rnd,
    doublereal *eps, integer *emin, doublereal *rmin, integer *emax,
    doublereal *rmax);

/* Subroutine */ int F77NAME(dlamc4)(integer *emin, doublereal *start, integer *base);

/* Subroutine */ int F77NAME(dlamc5)(integer *beta, integer *p, integer *emin,
    logical *ieee, integer *emax, doublereal *rmax);

/* Subroutine */ int F77NAME(dlamrg)(integer *n1, integer *n2, doublereal *a, integer
    *dtrd1, integer *dtrd2, integer *index);

/* Subroutine */ int F77NAME(dlanv2)(doublereal *a, doublereal *b, doublereal *c__,
    doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r,
     doublereal *rt2i, doublereal *cs, doublereal *sn);

/* Subroutine */ int F77NAME(dlapll)(integer *n, doublereal *x, integer *incx,
    doublereal *y, integer *incy, doublereal *ssmin);

/* Subroutine */ int F77NAME(dlapmt)(logical *forwrd, integer *m, integer *n,
    doublereal *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(dlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__,
    doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqge)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal
    *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqp2)(integer *m, integer *n, integer *offset,
    doublereal *a, integer *lda, integer *jpvt, doublereal *tau,
    doublereal *vn1, doublereal *vn2, doublereal *work);

/* Subroutine */ int F77NAME(dlaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, doublereal *a, integer *lda, integer *jpvt,
    doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *auxv,
    doublereal *f, integer *ldf);

/* Subroutine */ int F77NAME(dlaqsb)(char *uplo, integer *n, integer *kd, doublereal *
    ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
     char *equed);

/* Subroutine */ int F77NAME(dlaqsp)(char *uplo, integer *n, doublereal *ap,
    doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqsy)(char *uplo, integer *n, doublereal *a, integer *
    lda, doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(dlaqtr)(logical *ltran, logical *lfloatreal, integer *n,
    doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal
    *scale, doublereal *x, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dlar1v)(integer *n, integer *b1, integer *bn, doublereal
    *sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
    lld, doublereal *gersch, doublereal *z__, doublereal *ztz, doublereal
    *mingma, integer *r__, integer *isuppz, doublereal *work);

/* Subroutine */ int F77NAME(dlar2v)(integer *n, doublereal *x, doublereal *y,
    doublereal *z__, integer *incx, doublereal *c__, doublereal *s,
    integer *incc);

/* Subroutine */ int F77NAME(dlarf)(char *side, integer *m, integer *n, doublereal *v,
     integer *incv, doublereal *tau, doublereal *c__, integer *ldc,
    doublereal *work);

/* Subroutine */ int F77NAME(dlarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, doublereal *v, integer *
    ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc,
    doublereal *work, integer *ldwork);

/* Subroutine */ int F77NAME(dlarfg)(integer *n, doublereal *alpha, doublereal *x,
    integer *incx, doublereal *tau);

/* Subroutine */ int F77NAME(dlarft)(char *direct, char *storev, integer *n, integer *
    k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t,
    integer *ldt);

/* Subroutine */ int F77NAME(dlarfx)(char *side, integer *m, integer *n, doublereal *
    v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);

/* Subroutine */ int F77NAME(dlargv)(integer *n, doublereal *x, integer *incx,
    doublereal *y, integer *incy, doublereal *c__, integer *incc);

/* Subroutine */ int F77NAME(dlarnv)(integer *idist, integer *iseed, integer *n,
    doublereal *x);

/* Subroutine */ int F77NAME(dlarrb)(integer *n, doublereal *d__, doublereal *l,
    doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast,
    doublereal *sigma, doublereal *reltol, doublereal *w, doublereal *
    wgap, doublereal *werr, doublereal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dlarre)(integer *n, doublereal *d__, doublereal *e,
    doublereal *tol, integer *nsplit, integer *isplit, integer *m,
    doublereal *w, doublereal *woff, doublereal *gersch, doublereal *work,
     integer *info);

/* Subroutine */ int F77NAME(dlarrf)(integer *n, doublereal *d__, doublereal *l,
    doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast,
    doublereal *w, doublereal *dplus, doublereal *lplus, doublereal *work,
     integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlarrv)(integer *n, doublereal *d__, doublereal *l,
    integer *isplit, integer *m, doublereal *w, integer *iblock,
    doublereal *gersch, doublereal *tol, doublereal *z__, integer *ldz,
    integer *isuppz, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlartg)(doublereal *f, doublereal *g, doublereal *cs,
    doublereal *sn, doublereal *r__);

/* Subroutine */ int F77NAME(dlartv)(integer *n, doublereal *x, integer *incx,
    doublereal *y, integer *incy, doublereal *c__, doublereal *s, integer
    *incc);

/* Subroutine */ int F77NAME(dlaruv)(integer *iseed, integer *n, doublereal *x);

/* Subroutine */ int F77NAME(dlarz)(char *side, integer *m, integer *n, integer *l,
    doublereal *v, integer *incv, doublereal *tau, doublereal *c__,
    integer *ldc, doublereal *work);

/* Subroutine */ int F77NAME(dlarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, doublereal *v,
     integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *
    ldc, doublereal *work, integer *ldwork);

/* Subroutine */ int F77NAME(dlarzt)(char *direct, char *storev, integer *n, integer *
    k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t,
    integer *ldt);

/* Subroutine */ int F77NAME(dlas2)(doublereal *f, doublereal *g, doublereal *h__,
    doublereal *ssmin, doublereal *ssmax);

/* Subroutine */ int F77NAME(dlascl)(char *type__, integer *kl, integer *ku,
    doublereal *cfrom, doublereal *cto, integer *m, integer *n,
    doublereal *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(dlasd0)(integer *n, integer *sqre, doublereal *d__,
    doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *
    ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *
    info);

/* Subroutine */ int F77NAME(dlasd1)(integer *nl, integer *nr, integer *sqre,
    doublereal *d__, doublereal *alpha, doublereal *beta, doublereal *u,
    integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *
    iwork, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dlasd2)(integer *nl, integer *nr, integer *sqre, integer
    *k, doublereal *d__, doublereal *z__, doublereal *alpha, doublereal *
    beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
    doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2,
    integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
    idxq, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(dlasd3)(integer *nl, integer *nr, integer *sqre, integer
    *k, doublereal *d__, doublereal *q, integer *ldq, doublereal *dsigma,
    doublereal *u, integer *ldu, doublereal *u2, integer *ldu2,
    doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2,
    integer *idxc, integer *ctot, doublereal *z__, integer *info);

/* Subroutine */ int F77NAME(dlasd4)(integer *n, integer *i__, doublereal *d__,
    doublereal *z__, doublereal *delta, doublereal *rho, doublereal *
    sigma, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dlasd5)(integer *i__, doublereal *d__, doublereal *z__,
    doublereal *delta, doublereal *rho, doublereal *dsigma, doublereal *
    work);

/* Subroutine */ int F77NAME(dlasd6)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, doublereal *d__, doublereal *vf, doublereal *vl,
    doublereal *alpha, doublereal *beta, integer *idxq, integer *perm,
    integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
     integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *
    difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlasd7)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *k, doublereal *d__, doublereal *z__,
    doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl,
    doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *
    dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm,
    integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
     integer *ldgnum, doublereal *c__, doublereal *s, integer *info);

/* Subroutine */ int F77NAME(dlasd8)(integer *icompq, integer *k, doublereal *d__,
    doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl,
    doublereal *difr, integer *lddifr, doublereal *dsigma, doublereal *
    work, integer *info);

/* Subroutine */ int F77NAME(dlasd9)(integer *icompq, integer *ldu, integer *k,
    doublereal *d__, doublereal *z__, doublereal *vf, doublereal *vl,
    doublereal *difl, doublereal *difr, doublereal *dsigma, doublereal *
    work, integer *info);

/* Subroutine */ int F77NAME(dlasda)(integer *icompq, integer *smlsiz, integer *n,
    integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer
    *ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr,
    doublereal *z__, doublereal *poles, integer *givptr, integer *givcol,
    integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__,
    doublereal *s, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dlasdq)(char *uplo, integer *sqre, integer *n, integer *
    ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e,
    doublereal *vt, integer *ldvt, doublereal *u, integer *ldu,
    doublereal *c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dlasdt)(integer *n, integer *lvl, integer *nd, integer *
    inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int F77NAME(dlaset)(char *uplo, integer *m, integer *n, doublereal *
    alpha, doublereal *beta, doublereal *a, integer *lda);

/* Subroutine */ int F77NAME(dlasq1)(integer *n, doublereal *d__, doublereal *e,
    doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dlasq2)(integer *n, doublereal *z__, integer *info);

/* Subroutine */ int F77NAME(dlasq3)(integer *i0, integer *n0, doublereal *z__,
    integer *pp, doublereal *dmin__, doublereal *sigma, doublereal *desig,
     doublereal *qmax, integer *nfail, integer *iter, integer *ndiv,
    logical *ieee);

/* Subroutine */ int F77NAME(dlasq4)(integer *i0, integer *n0, doublereal *z__,
    integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1,
    doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2,
    doublereal *tau, integer *ttype);

/* Subroutine */ int F77NAME(dlasq5)(integer *i0, integer *n0, doublereal *z__,
    integer *pp, doublereal *tau, doublereal *dmin__, doublereal *dmin1,
    doublereal *dmin2, doublereal *dn, doublereal *dnm1, doublereal *dnm2,
     logical *ieee);

/* Subroutine */ int F77NAME(dlasq6)(integer *i0, integer *n0, doublereal *z__,
    integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
     doublereal *dn, doublereal *dnm1, doublereal *dnm2);

/* Subroutine */ int F77NAME(dlasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, doublereal *c__, doublereal *s, doublereal *a, integer *
    lda);

/* Subroutine */ int F77NAME(dlasrt)(char *id, integer *n, doublereal *d__, integer *
    info);

/* Subroutine */ int F77NAME(dlassq)(integer *n, doublereal *x, integer *incx,
    doublereal *scale, doublereal *sumsq);

/* Subroutine */ int F77NAME(dlasv2)(doublereal *f, doublereal *g, doublereal *h__,
    doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *
    csr, doublereal *snl, doublereal *csl);

/* Subroutine */ int F77NAME(dlaswp)(integer *n, doublereal *a, integer *lda, integer
    *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(dlasy2)(logical *ltranl, logical *ltranr, integer *isgn,
    integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
    tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale,
    doublereal *x, integer *ldx, doublereal *xnorm, integer *info);

/* Subroutine */ int F77NAME(dlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     doublereal *a, integer *lda, integer *ipiv, doublereal *w, integer *
    ldw, integer *info);

/* Subroutine */ int F77NAME(dlatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, doublereal *ab, integer *ldab,
    doublereal *x, doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatdf)(integer *ijob, integer *n, doublereal *z__,
    integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal,
    integer *ipiv, integer *jpiv);

/* Subroutine */ int F77NAME(dlatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doublereal *ap, doublereal *x, doublereal *scale,
    doublereal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatrd)(char *uplo, integer *n, integer *nb, doublereal *
    a, integer *lda, doublereal *e, doublereal *tau, doublereal *w,
    integer *ldw);

/* Subroutine */ int F77NAME(dlatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doublereal *a, integer *lda, doublereal *x,
    doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int F77NAME(dlatrz)(integer *m, integer *n, integer *l, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work);

/* Subroutine */ int F77NAME(dlatzm)(char *side, integer *m, integer *n, doublereal *
    v, integer *incv, doublereal *tau, doublereal *c1, doublereal *c2,
    integer *ldc, doublereal *work);

/* Subroutine */ int F77NAME(dlauu2)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dlauum)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dopgtr)(char *uplo, integer *n, doublereal *ap,
    doublereal *tau, doublereal *q, integer *ldq, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dopmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doublereal *ap, doublereal *tau, doublereal *c__, integer
    *ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dorg2l)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dorg2r)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dorgbr)(char *vect, integer *m, integer *n, integer *k,
    doublereal *a, integer *lda, doublereal *tau, doublereal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorghr)(integer *n, integer *ilo, integer *ihi,
    doublereal *a, integer *lda, doublereal *tau, doublereal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorgl2)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dorglq)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgql)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgqr)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgr2)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dorgrq)(integer *m, integer *n, integer *k, doublereal *
    a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorgtr)(char *uplo, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dorm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dorm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dormbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau,
    doublereal *c__, integer *ldc, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dormhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
    tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dorml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dormlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormql)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dormr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau,
    doublereal *c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dormrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dormrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau,
    doublereal *c__, integer *ldc, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dormtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *
    c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dpbcon)(char *uplo, integer *n, integer *kd, doublereal *
    ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpbequ)(char *uplo, integer *n, integer *kd, doublereal *
    ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
     integer *info);

/* Subroutine */ int F77NAME(dpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb,
    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
    ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dpbstf)(char *uplo, integer *n, integer *kd, doublereal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(dpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb,
    integer *ldafb, char *equed, doublereal *s, doublereal *b, integer *
    ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
     doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpbtf2)(char *uplo, integer *n, integer *kd, doublereal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbtrf)(char *uplo, integer *n, integer *kd, doublereal *
    ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(dpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(dpocon)(char *uplo, integer *n, doublereal *a, integer *
    lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dpoequ)(integer *n, doublereal *a, integer *lda,
    doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int F77NAME(dporfs)(char *uplo, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *af, integer *ldaf,
    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
    ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(dposv)(char *uplo, integer *n, integer *nrhs, doublereal
    *a, integer *lda, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf,
    char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *
    x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *
    berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dpotf2)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotrf)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotri)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *info);

/* Subroutine */ int F77NAME(dpotrs)(char *uplo, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dppcon)(char *uplo, integer *n, doublereal *ap,
    doublereal *anorm, doublereal *rcond, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dppequ)(char *uplo, integer *n, doublereal *ap,
    doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int F77NAME(dpprfs)(char *uplo, integer *n, integer *nrhs,
    doublereal *ap, doublereal *afp, doublereal *b, integer *ldb,
    doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dppsv)(char *uplo, integer *n, integer *nrhs, doublereal
    *ap, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublereal *ap, doublereal *afp, char *equed, doublereal *s,
    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dpptrf)(char *uplo, integer *n, doublereal *ap, integer *
    info);

/* Subroutine */ int F77NAME(dpptri)(char *uplo, integer *n, doublereal *ap, integer *
    info);

/* Subroutine */ int F77NAME(dpptrs)(char *uplo, integer *n, integer *nrhs,
    doublereal *ap, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptcon)(integer *n, doublereal *d__, doublereal *e,
    doublereal *anorm, doublereal *rcond, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dpteqr)(char *compz, integer *n, doublereal *d__,
    doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dptrfs)(integer *n, integer *nrhs, doublereal *d__,
    doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer
    *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
     doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dptsv)(integer *n, integer *nrhs, doublereal *d__,
    doublereal *e, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptsvx)(char *fact, integer *n, integer *nrhs,
    doublereal *d__, doublereal *e, doublereal *df, doublereal *ef,
    doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
    rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
    info);

/* Subroutine */ int F77NAME(dpttrf)(integer *n, doublereal *d__, doublereal *e,
    integer *info);

/* Subroutine */ int F77NAME(dpttrs)(integer *n, integer *nrhs, doublereal *d__,
    doublereal *e, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dptts2)(integer *n, integer *nrhs, doublereal *d__,
    doublereal *e, doublereal *b, integer *ldb);

/* Subroutine */ int F77NAME(drscl)(integer *n, doublereal *sa, doublereal *sx,
    integer *incx);

/* Subroutine */ int F77NAME(dsbev)(char *jobz, char *uplo, integer *n, integer *kd,
    doublereal *ab, integer *ldab, doublereal *w, doublereal *z__,
    integer *ldz, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dsbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    doublereal *ab, integer *ldab, doublereal *w, doublereal *z__,
    integer *ldz, doublereal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *
    ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu,
    doublereal *abstol, integer *m, doublereal *w, doublereal *z__,
    integer *ldz, doublereal *work, integer *iwork, integer *ifail,
    integer *info);

/* Subroutine */ int F77NAME(dsbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
    ldbb, doublereal *x, integer *ldx, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dsbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
    ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dsbgvd)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
    ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *
    bb, integer *ldbb, doublereal *q, integer *ldq, doublereal *vl,
    doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer
    *m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    doublereal *ab, integer *ldab, doublereal *d__, doublereal *e,
    doublereal *q, integer *ldq, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dspcon)(char *uplo, integer *n, doublereal *ap, integer *
    ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(dspev)(char *jobz, char *uplo, integer *n, doublereal *
    ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dspevd)(char *jobz, char *uplo, integer *n, doublereal *
    ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dspevx)(char *jobz, char *range, char *uplo, integer *n,
    doublereal *ap, doublereal *vl, doublereal *vu, integer *il, integer *
    iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__,
    integer *ldz, doublereal *work, integer *iwork, integer *ifail,
    integer *info);

/* Subroutine */ int F77NAME(dspgst)(integer *itype, char *uplo, integer *n,
    doublereal *ap, doublereal *bp, integer *info);

/* Subroutine */ int F77NAME(dspgv)(integer *itype, char *jobz, char *uplo, integer *
    n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__,
    integer *ldz, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dspgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__,
    integer *ldz, doublereal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dspgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doublereal *ap, doublereal *bp, doublereal *vl,
    doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer
    *m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsprfs)(char *uplo, integer *n, integer *nrhs,
    doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b,
    integer *ldb, doublereal *x, integer *ldx, doublereal *ferr,
    doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dspsv)(char *uplo, integer *n, integer *nrhs, doublereal
    *ap, integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b,
    integer *ldb, doublereal *x, integer *ldx, doublereal *rcond,
    doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dsptrd)(char *uplo, integer *n, doublereal *ap,
    doublereal *d__, doublereal *e, doublereal *tau, integer *info);

/* Subroutine */ int F77NAME(dsptrf)(char *uplo, integer *n, doublereal *ap, integer *
    ipiv, integer *info);

/* Subroutine */ int F77NAME(dsptri)(char *uplo, integer *n, doublereal *ap, integer *
    ipiv, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dsptrs)(char *uplo, integer *n, integer *nrhs,
    doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dstebz)(char *range, char *order, integer *n, doublereal
    *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol,
    doublereal *d__, doublereal *e, integer *m, integer *nsplit,
    doublereal *w, integer *iblock, integer *isplit, doublereal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dstedc)(char *compz, integer *n, doublereal *d__,
    doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstegr)(char *jobz, char *range, integer *n, doublereal *
    d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il,
    integer *iu, doublereal *abstol, integer *m, doublereal *w,
    doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstein)(integer *n, doublereal *d__, doublereal *e,
    integer *m, doublereal *w, integer *iblock, integer *isplit,
    doublereal *z__, integer *ldz, doublereal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsteqr)(char *compz, integer *n, doublereal *d__,
    doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dsterf)(integer *n, doublereal *d__, doublereal *e,
    integer *info);

/* Subroutine */ int F77NAME(dstev)(char *jobz, integer *n, doublereal *d__,
    doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(dstevd)(char *jobz, integer *n, doublereal *d__,
    doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstevr)(char *jobz, char *range, integer *n, doublereal *
    d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il,
    integer *iu, doublereal *abstol, integer *m, doublereal *w,
    doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dstevx)(char *jobz, char *range, integer *n, doublereal *
    d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il,
    integer *iu, doublereal *abstol, integer *m, doublereal *w,
    doublereal *z__, integer *ldz, doublereal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsycon)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dsyev)(char *jobz, char *uplo, integer *n, doublereal *a,
     integer *lda, doublereal *w, doublereal *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(dsyevd)(char *jobz, char *uplo, integer *n, doublereal *
    a, integer *lda, doublereal *w, doublereal *work, integer *lwork,
    integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsyevr)(char *jobz, char *range, char *uplo, integer *n,
    doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
    il, integer *iu, doublereal *abstol, integer *m, doublereal *w,
    doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsyevx)(char *jobz, char *range, char *uplo, integer *n,
    doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
    il, integer *iu, doublereal *abstol, integer *m, doublereal *w,
    doublereal *z__, integer *ldz, doublereal *work, integer *lwork,
    integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsygs2)(integer *itype, char *uplo, integer *n,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dsygst)(integer *itype, char *uplo, integer *n,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dsygv)(integer *itype, char *jobz, char *uplo, integer *
    n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
    doublereal *w, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsygvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
    doublereal *w, doublereal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(dsygvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doublereal *a, integer *lda, doublereal *b, integer
    *ldb, doublereal *vl, doublereal *vu, integer *il, integer *iu,
    doublereal *abstol, integer *m, doublereal *w, doublereal *z__,
    integer *ldz, doublereal *work, integer *lwork, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(dsyrfs)(char *uplo, integer *n, integer *nrhs,
    doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
    ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx,
    doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dsysv)(char *uplo, integer *n, integer *nrhs, doublereal
    *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb,
    doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf,
    integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *
    ldx, doublereal *rcond, doublereal *ferr, doublereal *berr,
    doublereal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dsytd2)(char *uplo, integer *n, doublereal *a, integer *
    lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info);

/* Subroutine */ int F77NAME(dsytf2)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(dsytrd)(char *uplo, integer *n, doublereal *a, integer *
    lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsytrf)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dsytri)(char *uplo, integer *n, doublereal *a, integer *
    lda, integer *ipiv, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dsytrs)(char *uplo, integer *n, integer *nrhs,
    doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dtbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, doublereal *ab, integer *ldab, doublereal *rcond,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal
    *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr,
    doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal
    *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(dtgevc)(char *side, char *howmny, logical *select,
    integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
    doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer
    *mm, integer *m, doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dtgex2)(logical *wantq, logical *wantz, integer *n,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
    q, integer *ldq, doublereal *z__, integer *ldz, integer *j1, integer *
    n1, integer *n2, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dtgexc)(logical *wantq, logical *wantz, integer *n,
    doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
    q, integer *ldq, doublereal *z__, integer *ldz, integer *ifst,
    integer *ilst, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(dtgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, doublereal *a, integer *lda, doublereal *
    b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
    beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz,
    integer *m, doublereal *pl, doublereal *pr, doublereal *dif,
    doublereal *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(dtgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, doublereal *a,
    integer *lda, doublereal *b, integer *ldb, doublereal *tola,
    doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u,
    integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *
    ldq, doublereal *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(dtgsna)(char *job, char *howmny, logical *select,
    integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
    doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr,
    doublereal *s, doublereal *dif, integer *mm, integer *m, doublereal *
    work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
    doublereal *c__, integer *ldc, doublereal *d__, integer *ldd,
    doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
    scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer
    *pq, integer *info);

/* Subroutine */ int F77NAME(dtgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
    doublereal *c__, integer *ldc, doublereal *d__, integer *ldd,
    doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
    scale, doublereal *dif, doublereal *work, integer *lwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dtpcon)(char *norm, char *uplo, char *diag, integer *n,
    doublereal *ap, doublereal *rcond, doublereal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(dtprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublereal *ap, doublereal *b, integer *ldb,
    doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtptri)(char *uplo, char *diag, integer *n, doublereal *
    ap, integer *info);

/* Subroutine */ int F77NAME(dtptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(dtrcon)(char *norm, char *uplo, char *diag, integer *n,
    doublereal *a, integer *lda, doublereal *rcond, doublereal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtrevc)(char *side, char *howmny, logical *select,
    integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
    ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m,
    doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dtrexc)(char *compq, integer *n, doublereal *t, integer *
    ldt, doublereal *q, integer *ldq, integer *ifst, integer *ilst,
    doublereal *work, integer *info);

/* Subroutine */ int F77NAME(dtrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
    ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(dtrsen)(char *job, char *compq, logical *select, integer
    *n, doublereal *t, integer *ldt, doublereal *q, integer *ldq,
    doublereal *wr, doublereal *wi, integer *m, doublereal *s, doublereal
    *sep, doublereal *work, integer *lwork, integer *iwork, integer *
    liwork, integer *info);

/* Subroutine */ int F77NAME(dtrsna)(char *job, char *howmny, logical *select,
    integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
    ldvl, doublereal *vr, integer *ldvr, doublereal *s, doublereal *sep,
    integer *mm, integer *m, doublereal *work, integer *ldwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(dtrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
    ldb, doublereal *c__, integer *ldc, doublereal *scale, integer *info);

/* Subroutine */ int F77NAME(dtrti2)(char *uplo, char *diag, integer *n, doublereal *
    a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(dtrtri)(char *uplo, char *diag, integer *n, doublereal *
    a, integer *lda, integer *info);

/**
  Lapack routine : DTRTRS solves a triangular system of the form
  A * X = B  or  A**T * X = B,
  where A is a triangular matrix of order N, and B is an N-by-NRHS
  matrix.  A check is made to verify that A is nonsingular.
*/
int F77NAME(dtrtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(dtzrqf)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, integer *info);

/* Subroutine */ int F77NAME(dtzrzf)(integer *m, integer *n, doublereal *a, integer *
    lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

integer F77NAME(icmax1)(integer *n, floatcomplex *cx, integer *incx);

integer F77NAME(ieeeck)(integer *ispec, floatreal *zero, floatreal *one);

integer F77NAME(ilaenv)(integer *ispec, char *name__, char *opts, integer *n1,
    integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen
    opts_len);

integer F77NAME(izmax1)(integer *n, doublecomplex *cx, integer *incx);

/* Subroutine */ int F77NAME(sbdsdc)(char *uplo, char *compq, integer *n, floatreal *d__,
    floatreal *e, floatreal *u, integer *ldu, floatreal *vt, integer *ldvt, floatreal *q,
    integer *iq, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, floatreal *d__, floatreal *e, floatreal *vt, integer *ldvt, floatreal *
    u, integer *ldu, floatreal *c__, integer *ldc, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sdisna)(char *job, integer *m, integer *n, floatreal *d__,
    floatreal *sep, integer *info);

/* Subroutine */ int F77NAME(sgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, floatreal *ab, integer *ldab, floatreal *d__, floatreal *
    e, floatreal *q, integer *ldq, floatreal *pt, integer *ldpt, floatreal *c__, integer
    *ldc, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     floatreal *ab, integer *ldab, integer *ipiv, floatreal *anorm, floatreal *rcond,
    floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     floatreal *ab, integer *ldab, floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal *
    colcnd, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(sgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, floatreal *ab, integer *ldab, floatreal *afb, integer *ldafb,
     integer *ipiv, floatreal *b, integer *ldb, floatreal *x, integer *ldx, floatreal *
    ferr, floatreal *berr, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, floatreal *ab, integer *ldab, integer *ipiv, floatreal *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(sgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, floatreal *ab, integer *ldab, floatreal *afb,
    integer *ldafb, integer *ipiv, char *equed, floatreal *r__, floatreal *c__,
    floatreal *b, integer *ldb, floatreal *x, integer *ldx, floatreal *rcond, floatreal *ferr,
     floatreal *berr, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     floatreal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     floatreal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, floatreal *ab, integer *ldab, integer *ipiv, floatreal *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, floatreal *scale, integer *m, floatreal *v, integer *ldv, integer
    *info);

/* Subroutine */ int F77NAME(sgebal)(char *job, integer *n, floatreal *a, integer *lda,
    integer *ilo, integer *ihi, floatreal *scale, integer *info);

/* Subroutine */ int F77NAME(sgebd2)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *d__, floatreal *e, floatreal *tauq, floatreal *taup, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgebrd)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *d__, floatreal *e, floatreal *tauq, floatreal *taup, floatreal *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(sgecon)(char *norm, integer *n, floatreal *a, integer *lda,
    floatreal *anorm, floatreal *rcond, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgeequ)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal *colcnd, floatreal *amax, integer
    *info);

/* Subroutine */ int F77NAME(sgees)(char *jobvs, char *sort, L_fp select, integer *n,
    floatreal *a, integer *lda, integer *sdim, floatreal *wr, floatreal *wi, floatreal *vs,
    integer *ldvs, floatreal *work, integer *lwork, logical *bwork, integer *
    info);

/* Subroutine */ int F77NAME(sgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, floatreal *a, integer *lda, integer *sdim, floatreal *wr,
    floatreal *wi, floatreal *vs, integer *ldvs, floatreal *rconde, floatreal *rcondv, floatreal *
    work, integer *lwork, integer *iwork, integer *liwork, logical *bwork,
     integer *info);

/* Subroutine */ int F77NAME(sgeev)(char *jobvl, char *jobvr, integer *n, floatreal *a,
    integer *lda, floatreal *wr, floatreal *wi, floatreal *vl, integer *ldvl, floatreal *vr,
    integer *ldvr, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, floatreal *a, integer *lda, floatreal *wr, floatreal *wi, floatreal *
    vl, integer *ldvl, floatreal *vr, integer *ldvr, integer *ilo, integer *
    ihi, floatreal *scale, floatreal *abnrm, floatreal *rconde, floatreal *rcondv, floatreal *work,
     integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgegs)(char *jobvsl, char *jobvsr, integer *n, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, floatreal *alphar, floatreal *alphai, floatreal
    *beta, floatreal *vsl, integer *ldvsl, floatreal *vsr, integer *ldvsr, floatreal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgegv)(char *jobvl, char *jobvr, integer *n, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, floatreal *alphar, floatreal *alphai, floatreal
    *beta, floatreal *vl, integer *ldvl, floatreal *vr, integer *ldvr, floatreal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgehd2)(integer *n, integer *ilo, integer *ihi, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgehrd)(integer *n, integer *ilo, integer *ihi, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelq2)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgelqf)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgels)(char *trans, integer *m, integer *n, integer *
    nrhs, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelsd)(integer *m, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, floatreal *s, floatreal *rcond, integer *
    rank, floatreal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgelss)(integer *m, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, floatreal *s, floatreal *rcond, integer *
    rank, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgelsx)(integer *m, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, integer *jpvt, floatreal *rcond,
    integer *rank, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgelsy)(integer *m, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, integer *jpvt, floatreal *rcond,
    integer *rank, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeql2)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgeqlf)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeqp3)(integer *m, integer *n, floatreal *a, integer *lda,
    integer *jpvt, floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgeqpf)(integer *m, integer *n, floatreal *a, integer *lda,
    integer *jpvt, floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgeqr2)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgeqrf)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgerfs)(char *trans, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *af, integer *ldaf, integer *ipiv, floatreal *b,
    integer *ldb, floatreal *x, integer *ldx, floatreal *ferr, floatreal *berr, floatreal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgerq2)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgerqf)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgesc2)(integer *n, floatreal *a, integer *lda, floatreal *rhs,
    integer *ipiv, integer *jpiv, floatreal *scale);

/* Subroutine */ int F77NAME(sgesdd)(char *jobz, integer *m, integer *n, floatreal *a,
    integer *lda, floatreal *s, floatreal *u, integer *ldu, floatreal *vt, integer *ldvt,
     floatreal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgesv)(integer *n, integer *nrhs, floatreal *a, integer *lda,
    integer *ipiv, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgesvd)(char *jobu, char *jobvt, integer *m, integer *n,
    floatreal *a, integer *lda, floatreal *s, floatreal *u, integer *ldu, floatreal *vt,
    integer *ldvt, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, floatreal *a, integer *lda, floatreal *af, integer *ldaf, integer *ipiv,
    char *equed, floatreal *r__, floatreal *c__, floatreal *b, integer *ldb, floatreal *x,
    integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr, floatreal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgetc2)(integer *n, floatreal *a, integer *lda, integer *ipiv,
     integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(sgetf2)(integer *m, integer *n, floatreal *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgetrf)(integer *m, integer *n, floatreal *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgetri)(integer *n, floatreal *a, integer *lda, integer *ipiv,
     floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgetrs)(char *trans, integer *n, integer *nrhs, floatreal *a,
    integer *lda, integer *ipiv, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, floatreal *lscale, floatreal *rscale, integer *m, floatreal *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(sggbal)(char *job, integer *n, floatreal *a, integer *lda,
    floatreal *b, integer *ldb, integer *ilo, integer *ihi, floatreal *lscale, floatreal
    *rscale, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, integer *n, floatreal *a, integer *lda, floatreal *b, integer *ldb,
    integer *sdim, floatreal *alphar, floatreal *alphai, floatreal *beta, floatreal *vsl,
    integer *ldvsl, floatreal *vsr, integer *ldvsr, floatreal *work, integer *lwork,
     logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    selctg, char *sense, integer *n, floatreal *a, integer *lda, floatreal *b,
    integer *ldb, integer *sdim, floatreal *alphar, floatreal *alphai, floatreal *beta,
    floatreal *vsl, integer *ldvsl, floatreal *vsr, integer *ldvsr, floatreal *rconde,
    floatreal *rcondv, floatreal *work, integer *lwork, integer *iwork, integer *
    liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggev)(char *jobvl, char *jobvr, integer *n, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, floatreal *alphar, floatreal *alphai, floatreal
    *beta, floatreal *vl, integer *ldvl, floatreal *vr, integer *ldvr, floatreal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal
    *alphar, floatreal *alphai, floatreal *beta, floatreal *vl, integer *ldvl, floatreal *vr,
    integer *ldvr, integer *ilo, integer *ihi, floatreal *lscale, floatreal *rscale,
     floatreal *abnrm, floatreal *bbnrm, floatreal *rconde, floatreal *rcondv, floatreal *work,
    integer *lwork, integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(sggglm)(integer *n, integer *m, integer *p, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, floatreal *d__, floatreal *x, floatreal *y,
    floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal
    *q, integer *ldq, floatreal *z__, integer *ldz, integer *info);

/* Subroutine */ int F77NAME(sgglse)(integer *m, integer *n, integer *p, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, floatreal *c__, floatreal *d__, floatreal *x,
    floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggqrf)(integer *n, integer *m, integer *p, floatreal *a,
    integer *lda, floatreal *taua, floatreal *b, integer *ldb, floatreal *taub, floatreal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggrqf)(integer *m, integer *p, integer *n, floatreal *a,
    integer *lda, floatreal *taua, floatreal *b, integer *ldb, floatreal *taub, floatreal *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, floatreal *a, integer *lda,
     floatreal *b, integer *ldb, floatreal *alpha, floatreal *beta, floatreal *u, integer *
    ldu, floatreal *v, integer *ldv, floatreal *q, integer *ldq, floatreal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, floatreal *a, integer *lda, floatreal *b, integer *ldb,
    floatreal *tola, floatreal *tolb, integer *k, integer *l, floatreal *u, integer *ldu,
     floatreal *v, integer *ldv, floatreal *q, integer *ldq, integer *iwork, floatreal *
    tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sgtcon)(char *norm, integer *n, floatreal *dl, floatreal *d__,
    floatreal *du, floatreal *du2, integer *ipiv, floatreal *anorm, floatreal *rcond, floatreal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgtrfs)(char *trans, integer *n, integer *nrhs, floatreal *dl,
     floatreal *d__, floatreal *du, floatreal *dlf, floatreal *df, floatreal *duf, floatreal *du2,
    integer *ipiv, floatreal *b, integer *ldb, floatreal *x, integer *ldx, floatreal *
    ferr, floatreal *berr, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sgtsv)(integer *n, integer *nrhs, floatreal *dl, floatreal *d__,
    floatreal *du, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, floatreal *dl, floatreal *d__, floatreal *du, floatreal *dlf, floatreal *df, floatreal *duf,
    floatreal *du2, integer *ipiv, floatreal *b, integer *ldb, floatreal *x, integer *
    ldx, floatreal *rcond, floatreal *ferr, floatreal *berr, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(sgttrf)(integer *n, floatreal *dl, floatreal *d__, floatreal *du, floatreal *
    du2, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(sgttrs)(char *trans, integer *n, integer *nrhs, floatreal *dl,
     floatreal *d__, floatreal *du, floatreal *du2, integer *ipiv, floatreal *b, integer *ldb,
     integer *info);

/* Subroutine */ int F77NAME(sgtts2)(integer *itrans, integer *n, integer *nrhs, floatreal
    *dl, floatreal *d__, floatreal *du, floatreal *du2, integer *ipiv, floatreal *b, integer *
    ldb);

/* Subroutine */ int F77NAME(shgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, floatreal *a, integer *lda, floatreal *b, integer *
    ldb, floatreal *alphar, floatreal *alphai, floatreal *beta, floatreal *q, integer *ldq,
    floatreal *z__, integer *ldz, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(shsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, floatreal *h__, integer *ldh, floatreal *wr, floatreal *wi, floatreal
    *vl, integer *ldvl, floatreal *vr, integer *ldvr, integer *mm, integer *m,
    floatreal *work, integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(shseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, floatreal *h__, integer *ldh, floatreal *wr, floatreal *wi, floatreal *z__,
     integer *ldz, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(slabad)(floatreal *small, floatreal *large);

/* Subroutine */ int F77NAME(slabrd)(integer *m, integer *n, integer *nb, floatreal *a,
    integer *lda, floatreal *d__, floatreal *e, floatreal *tauq, floatreal *taup, floatreal *x,
    integer *ldx, floatreal *y, integer *ldy);

/* Subroutine */ int F77NAME(slacon)(integer *n, floatreal *v, floatreal *x, integer *isgn,
    floatreal *est, integer *kase);

/* Subroutine */ int F77NAME(slacpy)(char *uplo, integer *m, integer *n, floatreal *a,
    integer *lda, floatreal *b, integer *ldb);

/* Subroutine */ int F77NAME(sladiv)(floatreal *a, floatreal *b, floatreal *c__, floatreal *d__, floatreal *p,
    floatreal *q);

/* Subroutine */ int F77NAME(slae2)(floatreal *a, floatreal *b, floatreal *c__, floatreal *rt1, floatreal *rt2);

/* Subroutine */ int F77NAME(slaebz)(integer *ijob, integer *nitmax, integer *n,
    integer *mmax, integer *minp, integer *nbmin, floatreal *abstol, floatreal *
    reltol, floatreal *pivmin, floatreal *d__, floatreal *e, floatreal *e2, integer *nval,
    floatreal *ab, floatreal *c__, integer *mout, integer *nab, floatreal *work, integer
    *iwork, integer *info);

/* Subroutine */ int F77NAME(slaed0)(integer *icompq, integer *qsiz, integer *n, floatreal
    *d__, floatreal *e, floatreal *q, integer *ldq, floatreal *qstore, integer *ldqs,
    floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slaed1)(integer *n, floatreal *d__, floatreal *q, integer *ldq,
    integer *indxq, floatreal *rho, integer *cutpnt, floatreal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(slaed2)(integer *k, integer *n, integer *n1, floatreal *d__,
    floatreal *q, integer *ldq, integer *indxq, floatreal *rho, floatreal *z__, floatreal *
    dlamda, floatreal *w, floatreal *q2, integer *indx, integer *indxc, integer *
    indxp, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(slaed3)(integer *k, integer *n, integer *n1, floatreal *d__,
    floatreal *q, integer *ldq, floatreal *rho, floatreal *dlamda, floatreal *q2, integer *
    indx, integer *ctot, floatreal *w, floatreal *s, integer *info);

/* Subroutine */ int F77NAME(slaed4)(integer *n, integer *i__, floatreal *d__, floatreal *z__,
    floatreal *delta, floatreal *rho, floatreal *dlam, integer *info);

/* Subroutine */ int F77NAME(slaed5)(integer *i__, floatreal *d__, floatreal *z__, floatreal *delta,
    floatreal *rho, floatreal *dlam);

/* Subroutine */ int F77NAME(slaed6)(integer *kniter, logical *orgati, floatreal *rho,
    floatreal *d__, floatreal *z__, floatreal *finit, floatreal *tau, integer *info);

/* Subroutine */ int F77NAME(slaed7)(integer *icompq, integer *n, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, floatreal *d__, floatreal *q,
    integer *ldq, integer *indxq, floatreal *rho, integer *cutpnt, floatreal *
    qstore, integer *qptr, integer *prmptr, integer *perm, integer *
    givptr, integer *givcol, floatreal *givnum, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(slaed8)(integer *icompq, integer *k, integer *n, integer
    *qsiz, floatreal *d__, floatreal *q, integer *ldq, integer *indxq, floatreal *rho,
    integer *cutpnt, floatreal *z__, floatreal *dlamda, floatreal *q2, integer *ldq2,
    floatreal *w, integer *perm, integer *givptr, integer *givcol, floatreal *
    givnum, integer *indxp, integer *indx, integer *info);

/* Subroutine */ int F77NAME(slaed9)(integer *k, integer *kstart, integer *kstop,
    integer *n, floatreal *d__, floatreal *q, integer *ldq, floatreal *rho, floatreal *dlamda,
     floatreal *w, floatreal *s, integer *lds, integer *info);

/* Subroutine */ int F77NAME(slaeda)(integer *n, integer *tlvls, integer *curlvl,
    integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
    integer *givcol, floatreal *givnum, floatreal *q, integer *qptr, floatreal *z__,
    floatreal *ztemp, integer *info);

/* Subroutine */ int F77NAME(slaein)(logical *rightv, logical *noinit, integer *n,
    floatreal *h__, integer *ldh, floatreal *wr, floatreal *wi, floatreal *vr, floatreal *vi, floatreal
    *b, integer *ldb, floatreal *work, floatreal *eps3, floatreal *smlnum, floatreal *bignum,
    integer *info);

/* Subroutine */ int F77NAME(slaev2)(floatreal *a, floatreal *b, floatreal *c__, floatreal *rt1, floatreal *
    rt2, floatreal *cs1, floatreal *sn1);

/* Subroutine */ int F77NAME(slaexc)(logical *wantq, integer *n, floatreal *t, integer *
    ldt, floatreal *q, integer *ldq, integer *j1, integer *n1, integer *n2,
    floatreal *work, integer *info);

/* Subroutine */ int F77NAME(slag2)(floatreal *a, integer *lda, floatreal *b, integer *ldb,
    floatreal *safmin, floatreal *scale1, floatreal *scale2, floatreal *wr1, floatreal *wr2, floatreal *
    wi);

/* Subroutine */ int F77NAME(slags2)(logical *upper, floatreal *a1, floatreal *a2, floatreal *a3,
    floatreal *b1, floatreal *b2, floatreal *b3, floatreal *csu, floatreal *snu, floatreal *csv, floatreal *
    snv, floatreal *csq, floatreal *snq);

/* Subroutine */ int F77NAME(slagtf)(integer *n, floatreal *a, floatreal *lambda, floatreal *b, floatreal
    *c__, floatreal *tol, floatreal *d__, integer *in, integer *info);

/* Subroutine */ int F77NAME(slagtm)(char *trans, integer *n, integer *nrhs, floatreal *
    alpha, floatreal *dl, floatreal *d__, floatreal *du, floatreal *x, integer *ldx, floatreal *
    beta, floatreal *b, integer *ldb);

/* Subroutine */ int F77NAME(slagts)(integer *job, integer *n, floatreal *a, floatreal *b, floatreal
    *c__, floatreal *d__, integer *in, floatreal *y, floatreal *tol, integer *info);

/* Subroutine */ int F77NAME(slagv2)(floatreal *a, integer *lda, floatreal *b, integer *ldb,
    floatreal *alphar, floatreal *alphai, floatreal *beta, floatreal *csl, floatreal *snl, floatreal *
    csr, floatreal *snr);

/* Subroutine */ int F77NAME(slahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, floatreal *h__, integer *ldh, floatreal *wr, floatreal *
    wi, integer *iloz, integer *ihiz, floatreal *z__, integer *ldz, integer *
    info);

/* Subroutine */ int F77NAME(slahrd)(integer *n, integer *k, integer *nb, floatreal *a,
    integer *lda, floatreal *tau, floatreal *t, integer *ldt, floatreal *y, integer *ldy);

/* Subroutine */ int F77NAME(slaic1)(integer *job, integer *j, floatreal *x, floatreal *sest,
    floatreal *w, floatreal *gamma, floatreal *sestpr, floatreal *s, floatreal *c__);

/* Subroutine */ int F77NAME(slaln2)(logical *ltrans, integer *na, integer *nw, floatreal *
    smin, floatreal *ca, floatreal *a, integer *lda, floatreal *d1, floatreal *d2, floatreal *b,
    integer *ldb, floatreal *wr, floatreal *wi, floatreal *x, integer *ldx, floatreal *scale,
    floatreal *xnorm, integer *info);

/* Subroutine */ int F77NAME(slals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, floatreal *b, integer *ldb, floatreal *bx,
    integer *ldbx, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, floatreal *givnum, integer *ldgnum, floatreal *poles, floatreal *
    difl, floatreal *difr, floatreal *z__, integer *k, floatreal *c__, floatreal *s, floatreal *
    work, integer *info);

/* Subroutine */ int F77NAME(slalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, floatreal *b, integer *ldb, floatreal *bx, integer *ldbx, floatreal *
    u, integer *ldu, floatreal *vt, integer *k, floatreal *difl, floatreal *difr, floatreal *
    z__, floatreal *poles, integer *givptr, integer *givcol, integer *ldgcol,
    integer *perm, floatreal *givnum, floatreal *c__, floatreal *s, floatreal *work, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(slalsd)(char *uplo, integer *smlsiz, integer *n, integer
    *nrhs, floatreal *d__, floatreal *e, floatreal *b, integer *ldb, floatreal *rcond,
    integer *rank, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slamc1)(integer *beta, integer *t, logical *rnd, logical
    *ieee1);

/* Subroutine */ int F77NAME(slamc2)(integer *beta, integer *t, logical *rnd, floatreal *
    eps, integer *emin, floatreal *rmin, integer *emax, floatreal *rmax);

/* Subroutine */ int F77NAME(slamc4)(integer *emin, floatreal *start, integer *base);

/* Subroutine */ int F77NAME(slamc5)(integer *beta, integer *p, integer *emin,
    logical *ieee, integer *emax, floatreal *rmax);

/* Subroutine */ int F77NAME(slamrg)(integer *n1, integer *n2, floatreal *a, integer *
    strd1, integer *strd2, integer *index);

/* Subroutine */ int F77NAME(slanv2)(floatreal *a, floatreal *b, floatreal *c__, floatreal *d__, floatreal *
    rt1r, floatreal *rt1i, floatreal *rt2r, floatreal *rt2i, floatreal *cs, floatreal *sn);

/* Subroutine */ int F77NAME(slapll)(integer *n, floatreal *x, integer *incx, floatreal *y,
    integer *incy, floatreal *ssmin);

/* Subroutine */ int F77NAME(slapmt)(logical *forwrd, integer *m, integer *n, floatreal *x,
     integer *ldx, integer *k);

/* Subroutine */ int F77NAME(slaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     floatreal *ab, integer *ldab, floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal *
    colcnd, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqge)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *r__, floatreal *c__, floatreal *rowcnd, floatreal *colcnd, floatreal *amax, char *
    equed);

/* Subroutine */ int F77NAME(slaqp2)(integer *m, integer *n, integer *offset, floatreal *a,
     integer *lda, integer *jpvt, floatreal *tau, floatreal *vn1, floatreal *vn2, floatreal *
    work);

/* Subroutine */ int F77NAME(slaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, floatreal *a, integer *lda, integer *jpvt, floatreal *tau,
    floatreal *vn1, floatreal *vn2, floatreal *auxv, floatreal *f, integer *ldf);

/* Subroutine */ int F77NAME(slaqsb)(char *uplo, integer *n, integer *kd, floatreal *ab,
    integer *ldab, floatreal *s, floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqsp)(char *uplo, integer *n, floatreal *ap, floatreal *s, floatreal *
    scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqsy)(char *uplo, integer *n, floatreal *a, integer *lda,
    floatreal *s, floatreal *scond, floatreal *amax, char *equed);

/* Subroutine */ int F77NAME(slaqtr)(logical *ltran, logical *lfloatreal, integer *n, floatreal
    *t, integer *ldt, floatreal *b, floatreal *w, floatreal *scale, floatreal *x, floatreal *work,
    integer *info);

/* Subroutine */ int F77NAME(slar1v)(integer *n, integer *b1, integer *bn, floatreal *
    sigma, floatreal *d__, floatreal *l, floatreal *ld, floatreal *lld, floatreal *gersch, floatreal *
    z__, floatreal *ztz, floatreal *mingma, integer *r__, integer *isuppz, floatreal *
    work);

/* Subroutine */ int F77NAME(slar2v)(integer *n, floatreal *x, floatreal *y, floatreal *z__, integer
    *incx, floatreal *c__, floatreal *s, integer *incc);

/* Subroutine */ int F77NAME(slarf)(char *side, integer *m, integer *n, floatreal *v,
    integer *incv, floatreal *tau, floatreal *c__, integer *ldc, floatreal *work);

/* Subroutine */ int F77NAME(slarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, floatreal *v, integer *ldv,
    floatreal *t, integer *ldt, floatreal *c__, integer *ldc, floatreal *work, integer *
    ldwork);

/* Subroutine */ int F77NAME(slarfg)(integer *n, floatreal *alpha, floatreal *x, integer *incx,
    floatreal *tau);

/* Subroutine */ int F77NAME(slarft)(char *direct, char *storev, integer *n, integer *
    k, floatreal *v, integer *ldv, floatreal *tau, floatreal *t, integer *ldt);

/* Subroutine */ int F77NAME(slarfx)(char *side, integer *m, integer *n, floatreal *v,
    floatreal *tau, floatreal *c__, integer *ldc, floatreal *work);

/* Subroutine */ int F77NAME(slargv)(integer *n, floatreal *x, integer *incx, floatreal *y,
    integer *incy, floatreal *c__, integer *incc);

/* Subroutine */ int F77NAME(slarnv)(integer *idist, integer *iseed, integer *n, floatreal
    *x);

/* Subroutine */ int F77NAME(slarrb)(integer *n, floatreal *d__, floatreal *l, floatreal *ld, floatreal *
    lld, integer *ifirst, integer *ilast, floatreal *sigma, floatreal *reltol, floatreal
    *w, floatreal *wgap, floatreal *werr, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slarre)(integer *n, floatreal *d__, floatreal *e, floatreal *tol,
    integer *nsplit, integer *isplit, integer *m, floatreal *w, floatreal *woff,
    floatreal *gersch, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(slarrf)(integer *n, floatreal *d__, floatreal *l, floatreal *ld, floatreal *
    lld, integer *ifirst, integer *ilast, floatreal *w, floatreal *dplus, floatreal *
    lplus, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slarrv)(integer *n, floatreal *d__, floatreal *l, integer *isplit,
    integer *m, floatreal *w, integer *iblock, floatreal *gersch, floatreal *tol, floatreal *
    z__, integer *ldz, integer *isuppz, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(slartg)(floatreal *f, floatreal *g, floatreal *cs, floatreal *sn, floatreal *r__);

/* Subroutine */ int F77NAME(slartv)(integer *n, floatreal *x, integer *incx, floatreal *y,
    integer *incy, floatreal *c__, floatreal *s, integer *incc);

/* Subroutine */ int F77NAME(slaruv)(integer *iseed, integer *n, floatreal *x);

/* Subroutine */ int F77NAME(slarz)(char *side, integer *m, integer *n, integer *l,
    floatreal *v, integer *incv, floatreal *tau, floatreal *c__, integer *ldc, floatreal *
    work);

/* Subroutine */ int F77NAME(slarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, floatreal *v,
    integer *ldv, floatreal *t, integer *ldt, floatreal *c__, integer *ldc, floatreal *
    work, integer *ldwork);

/* Subroutine */ int F77NAME(slarzt)(char *direct, char *storev, integer *n, integer *
    k, floatreal *v, integer *ldv, floatreal *tau, floatreal *t, integer *ldt);

/* Subroutine */ int F77NAME(slas2)(floatreal *f, floatreal *g, floatreal *h__, floatreal *ssmin, floatreal *
    ssmax);

/* Subroutine */ int F77NAME(slascl)(char *type__, integer *kl, integer *ku, floatreal *
    cfrom, floatreal *cto, integer *m, integer *n, floatreal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(slasd0)(integer *n, integer *sqre, floatreal *d__, floatreal *e,
    floatreal *u, integer *ldu, floatreal *vt, integer *ldvt, integer *smlsiz,
    integer *iwork, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(slasd1)(integer *nl, integer *nr, integer *sqre, floatreal *
    d__, floatreal *alpha, floatreal *beta, floatreal *u, integer *ldu, floatreal *vt,
    integer *ldvt, integer *idxq, integer *iwork, floatreal *work, integer *
    info);

/* Subroutine */ int F77NAME(slasd2)(integer *nl, integer *nr, integer *sqre, integer
    *k, floatreal *d__, floatreal *z__, floatreal *alpha, floatreal *beta, floatreal *u, integer *
    ldu, floatreal *vt, integer *ldvt, floatreal *dsigma, floatreal *u2, integer *ldu2,
    floatreal *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc,
     integer *idxq, integer *coltyp, integer *info);

/* Subroutine */ int F77NAME(slasd3)(integer *nl, integer *nr, integer *sqre, integer
    *k, floatreal *d__, floatreal *q, integer *ldq, floatreal *dsigma, floatreal *u, integer *
    ldu, floatreal *u2, integer *ldu2, floatreal *vt, integer *ldvt, floatreal *vt2,
    integer *ldvt2, integer *idxc, integer *ctot, floatreal *z__, integer *
    info);

/* Subroutine */ int F77NAME(slasd4)(integer *n, integer *i__, floatreal *d__, floatreal *z__,
    floatreal *delta, floatreal *rho, floatreal *sigma, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(slasd5)(integer *i__, floatreal *d__, floatreal *z__, floatreal *delta,
    floatreal *rho, floatreal *dsigma, floatreal *work);

/* Subroutine */ int F77NAME(slasd6)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, floatreal *d__, floatreal *vf, floatreal *vl, floatreal *alpha, floatreal *beta,
     integer *idxq, integer *perm, integer *givptr, integer *givcol,
    integer *ldgcol, floatreal *givnum, integer *ldgnum, floatreal *poles, floatreal *
    difl, floatreal *difr, floatreal *z__, integer *k, floatreal *c__, floatreal *s, floatreal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slasd7)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *k, floatreal *d__, floatreal *z__, floatreal *zw, floatreal *vf,
    floatreal *vfw, floatreal *vl, floatreal *vlw, floatreal *alpha, floatreal *beta, floatreal *dsigma,
     integer *idx, integer *idxp, integer *idxq, integer *perm, integer *
    givptr, integer *givcol, integer *ldgcol, floatreal *givnum, integer *
    ldgnum, floatreal *c__, floatreal *s, integer *info);

/* Subroutine */ int F77NAME(slasd8)(integer *icompq, integer *k, floatreal *d__, floatreal *
    z__, floatreal *vf, floatreal *vl, floatreal *difl, floatreal *difr, integer *lddifr,
    floatreal *dsigma, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(slasd9)(integer *icompq, integer *ldu, integer *k, floatreal *
    d__, floatreal *z__, floatreal *vf, floatreal *vl, floatreal *difl, floatreal *difr, floatreal *
    dsigma, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(slasda)(integer *icompq, integer *smlsiz, integer *n,
    integer *sqre, floatreal *d__, floatreal *e, floatreal *u, integer *ldu, floatreal *vt,
    integer *k, floatreal *difl, floatreal *difr, floatreal *z__, floatreal *poles, integer *
    givptr, integer *givcol, integer *ldgcol, integer *perm, floatreal *givnum,
     floatreal *c__, floatreal *s, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(slasdq)(char *uplo, integer *sqre, integer *n, integer *
    ncvt, integer *nru, integer *ncc, floatreal *d__, floatreal *e, floatreal *vt,
    integer *ldvt, floatreal *u, integer *ldu, floatreal *c__, integer *ldc, floatreal *
    work, integer *info);

/* Subroutine */ int F77NAME(slasdt)(integer *n, integer *lvl, integer *nd, integer *
    inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ int F77NAME(slaset)(char *uplo, integer *m, integer *n, floatreal *alpha,
    floatreal *beta, floatreal *a, integer *lda);

/* Subroutine */ int F77NAME(slasq1)(integer *n, floatreal *d__, floatreal *e, floatreal *work,
    integer *info);

/* Subroutine */ int F77NAME(slasq2)(integer *n, floatreal *z__, integer *info);

/* Subroutine */ int F77NAME(slasq3)(integer *i0, integer *n0, floatreal *z__, integer *pp,
     floatreal *dmin__, floatreal *sigma, floatreal *desig, floatreal *qmax, integer *nfail,
    integer *iter, integer *ndiv, logical *ieee);

/* Subroutine */ int F77NAME(slasq4)(integer *i0, integer *n0, floatreal *z__, integer *pp,
     integer *n0in, floatreal *dmin__, floatreal *dmin1, floatreal *dmin2, floatreal *dn,
    floatreal *dn1, floatreal *dn2, floatreal *tau, integer *ttype);

/* Subroutine */ int F77NAME(slasq5)(integer *i0, integer *n0, floatreal *z__, integer *pp,
     floatreal *tau, floatreal *dmin__, floatreal *dmin1, floatreal *dmin2, floatreal *dn, floatreal *
    dnm1, floatreal *dnm2, logical *ieee);

/* Subroutine */ int F77NAME(slasq6)(integer *i0, integer *n0, floatreal *z__, integer *pp,
     floatreal *dmin__, floatreal *dmin1, floatreal *dmin2, floatreal *dn, floatreal *dnm1, floatreal *
    dnm2);

/* Subroutine */ int F77NAME(slasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, floatreal *c__, floatreal *s, floatreal *a, integer *lda);

/* Subroutine */ int F77NAME(slasrt)(char *id, integer *n, floatreal *d__, integer *info);

/* Subroutine */ int F77NAME(slassq)(integer *n, floatreal *x, integer *incx, floatreal *scale,
    floatreal *sumsq);

/* Subroutine */ int F77NAME(slasv2)(floatreal *f, floatreal *g, floatreal *h__, floatreal *ssmin, floatreal *
    ssmax, floatreal *snr, floatreal *csr, floatreal *snl, floatreal *csl);

/* Subroutine */ int F77NAME(slaswp)(integer *n, floatreal *a, integer *lda, integer *k1,
    integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(slasy2)(logical *ltranl, logical *ltranr, integer *isgn,
    integer *n1, integer *n2, floatreal *tl, integer *ldtl, floatreal *tr, integer *
    ldtr, floatreal *b, integer *ldb, floatreal *scale, floatreal *x, integer *ldx, floatreal
    *xnorm, integer *info);

/* Subroutine */ int F77NAME(slasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     floatreal *a, integer *lda, integer *ipiv, floatreal *w, integer *ldw, integer
    *info);

/* Subroutine */ int F77NAME(slatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, floatreal *ab, integer *ldab, floatreal *x,
    floatreal *scale, floatreal *cnorm, integer *info);

/* Subroutine */ int F77NAME(slatdf)(integer *ijob, integer *n, floatreal *z__, integer *
    ldz, floatreal *rhs, floatreal *rdsum, floatreal *rdscal, integer *ipiv, integer *
    jpiv);

/* Subroutine */ int F77NAME(slatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, floatreal *ap, floatreal *x, floatreal *scale, floatreal *cnorm,
    integer *info);

/* Subroutine */ int F77NAME(slatrd)(char *uplo, integer *n, integer *nb, floatreal *a,
    integer *lda, floatreal *e, floatreal *tau, floatreal *w, integer *ldw);

/* Subroutine */ int F77NAME(slatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, floatreal *a, integer *lda, floatreal *x, floatreal *scale, floatreal
    *cnorm, integer *info);

/* Subroutine */ int F77NAME(slatrz)(integer *m, integer *n, integer *l, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work);

/* Subroutine */ int F77NAME(slatzm)(char *side, integer *m, integer *n, floatreal *v,
    integer *incv, floatreal *tau, floatreal *c1, floatreal *c2, integer *ldc, floatreal *
    work);

/* Subroutine */ int F77NAME(slauu2)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(slauum)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(sopgtr)(char *uplo, integer *n, floatreal *ap, floatreal *tau,
    floatreal *q, integer *ldq, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sopmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, floatreal *ap, floatreal *tau, floatreal *c__, integer *ldc, floatreal *work,
    integer *info);

/* Subroutine */ int F77NAME(sorg2l)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sorg2r)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sorgbr)(char *vect, integer *m, integer *n, integer *k,
    floatreal *a, integer *lda, floatreal *tau, floatreal *work, integer *lwork, integer
    *info);

/* Subroutine */ int F77NAME(sorghr)(integer *n, integer *ilo, integer *ihi, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgl2)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sorglq)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgql)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgqr)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgr2)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sorgrq)(integer *m, integer *n, integer *k, floatreal *a,
    integer *lda, floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorgtr)(char *uplo, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sorm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sormbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__,
    integer *ldc, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, floatreal *a, integer *lda, floatreal *tau, floatreal *
    c__, integer *ldc, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sorml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sormlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormql)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sormr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, floatreal *a, integer *lda, floatreal *tau, floatreal *c__,
    integer *ldc, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sormrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, floatreal *a, integer *lda, floatreal *tau, floatreal *c__,
    integer *ldc, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(sormtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, floatreal *a, integer *lda, floatreal *tau, floatreal *c__, integer *ldc,
     floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(spbcon)(char *uplo, integer *n, integer *kd, floatreal *ab,
    integer *ldab, floatreal *anorm, floatreal *rcond, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(spbequ)(char *uplo, integer *n, integer *kd, floatreal *ab,
    integer *ldab, floatreal *s, floatreal *scond, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(spbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, floatreal *ab, integer *ldab, floatreal *afb, integer *ldafb, floatreal *b,
    integer *ldb, floatreal *x, integer *ldx, floatreal *ferr, floatreal *berr, floatreal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spbstf)(char *uplo, integer *n, integer *kd, floatreal *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, floatreal *ab, integer *ldab, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(spbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, floatreal *ab, integer *ldab, floatreal *afb, integer *ldafb,
    char *equed, floatreal *s, floatreal *b, integer *ldb, floatreal *x, integer *ldx,
    floatreal *rcond, floatreal *ferr, floatreal *berr, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(spbtf2)(char *uplo, integer *n, integer *kd, floatreal *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbtrf)(char *uplo, integer *n, integer *kd, floatreal *ab,
    integer *ldab, integer *info);

/* Subroutine */ int F77NAME(spbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, floatreal *ab, integer *ldab, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(spocon)(char *uplo, integer *n, floatreal *a, integer *lda,
    floatreal *anorm, floatreal *rcond, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spoequ)(integer *n, floatreal *a, integer *lda, floatreal *s, floatreal
    *scond, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(sporfs)(char *uplo, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *af, integer *ldaf, floatreal *b, integer *ldb, floatreal *x,
     integer *ldx, floatreal *ferr, floatreal *berr, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(sposv)(char *uplo, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatreal *a, integer *lda, floatreal *af, integer *ldaf, char *equed,
    floatreal *s, floatreal *b, integer *ldb, floatreal *x, integer *ldx, floatreal *rcond,
    floatreal *ferr, floatreal *berr, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spotf2)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotrf)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotri)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *info);

/* Subroutine */ int F77NAME(spotrs)(char *uplo, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sppcon)(char *uplo, integer *n, floatreal *ap, floatreal *anorm,
    floatreal *rcond, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sppequ)(char *uplo, integer *n, floatreal *ap, floatreal *s, floatreal *
    scond, floatreal *amax, integer *info);

/* Subroutine */ int F77NAME(spprfs)(char *uplo, integer *n, integer *nrhs, floatreal *ap,
    floatreal *afp, floatreal *b, integer *ldb, floatreal *x, integer *ldx, floatreal *ferr,
    floatreal *berr, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sppsv)(char *uplo, integer *n, integer *nrhs, floatreal *ap,
    floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatreal *ap, floatreal *afp, char *equed, floatreal *s, floatreal *b, integer *
    ldb, floatreal *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr, floatreal
    *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(spptrf)(char *uplo, integer *n, floatreal *ap, integer *info);

/* Subroutine */ int F77NAME(spptri)(char *uplo, integer *n, floatreal *ap, integer *info);

/* Subroutine */ int F77NAME(spptrs)(char *uplo, integer *n, integer *nrhs, floatreal *ap,
    floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptcon)(integer *n, floatreal *d__, floatreal *e, floatreal *anorm,
    floatreal *rcond, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(spteqr)(char *compz, integer *n, floatreal *d__, floatreal *e,
    floatreal *z__, integer *ldz, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sptrfs)(integer *n, integer *nrhs, floatreal *d__, floatreal *e,
    floatreal *df, floatreal *ef, floatreal *b, integer *ldb, floatreal *x, integer *ldx,
    floatreal *ferr, floatreal *berr, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sptsv)(integer *n, integer *nrhs, floatreal *d__, floatreal *e,
    floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptsvx)(char *fact, integer *n, integer *nrhs, floatreal *d__,
     floatreal *e, floatreal *df, floatreal *ef, floatreal *b, integer *ldb, floatreal *x, integer
    *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(spttrf)(integer *n, floatreal *d__, floatreal *e, integer *info);

/* Subroutine */ int F77NAME(spttrs)(integer *n, integer *nrhs, floatreal *d__, floatreal *e,
    floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sptts2)(integer *n, integer *nrhs, floatreal *d__, floatreal *e,
    floatreal *b, integer *ldb);

/* Subroutine */ int F77NAME(srscl)(integer *n, floatreal *sa, floatreal *sx, integer *incx);

/* Subroutine */ int F77NAME(ssbev)(char *jobz, char *uplo, integer *n, integer *kd,
    floatreal *ab, integer *ldab, floatreal *w, floatreal *z__, integer *ldz, floatreal *work,
     integer *info);

/* Subroutine */ int F77NAME(ssbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    floatreal *ab, integer *ldab, floatreal *w, floatreal *z__, integer *ldz, floatreal *work,
     integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, floatreal *ab, integer *ldab, floatreal *q, integer *ldq, floatreal *vl,
     floatreal *vu, integer *il, integer *iu, floatreal *abstol, integer *m, floatreal *
    w, floatreal *z__, integer *ldz, floatreal *work, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(ssbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, floatreal *ab, integer *ldab, floatreal *bb, integer *ldbb, floatreal *
    x, integer *ldx, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(ssbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, floatreal *ab, integer *ldab, floatreal *bb, integer *ldbb, floatreal *
    w, floatreal *z__, integer *ldz, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(ssbgvd)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, floatreal *ab, integer *ldab, floatreal *bb, integer *ldbb, floatreal *
    w, floatreal *z__, integer *ldz, floatreal *work, integer *lwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, floatreal *ab, integer *ldab, floatreal *bb, integer *
    ldbb, floatreal *q, integer *ldq, floatreal *vl, floatreal *vu, integer *il, integer
    *iu, floatreal *abstol, integer *m, floatreal *w, floatreal *z__, integer *ldz, floatreal
    *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    floatreal *ab, integer *ldab, floatreal *d__, floatreal *e, floatreal *q, integer *ldq,
    floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sspcon)(char *uplo, integer *n, floatreal *ap, integer *ipiv,
    floatreal *anorm, floatreal *rcond, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sspev)(char *jobz, char *uplo, integer *n, floatreal *ap,
    floatreal *w, floatreal *z__, integer *ldz, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sspevd)(char *jobz, char *uplo, integer *n, floatreal *ap,
    floatreal *w, floatreal *z__, integer *ldz, floatreal *work, integer *lwork, integer
    *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sspevx)(char *jobz, char *range, char *uplo, integer *n,
    floatreal *ap, floatreal *vl, floatreal *vu, integer *il, integer *iu, floatreal *abstol,
    integer *m, floatreal *w, floatreal *z__, integer *ldz, floatreal *work, integer *
    iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(sspgst)(integer *itype, char *uplo, integer *n, floatreal *ap,
     floatreal *bp, integer *info);

/* Subroutine */ int F77NAME(sspgv)(integer *itype, char *jobz, char *uplo, integer *
    n, floatreal *ap, floatreal *bp, floatreal *w, floatreal *z__, integer *ldz, floatreal *work,
    integer *info);

/* Subroutine */ int F77NAME(sspgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, floatreal *ap, floatreal *bp, floatreal *w, floatreal *z__, integer *ldz, floatreal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sspgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, floatreal *ap, floatreal *bp, floatreal *vl, floatreal *vu, integer *il,
     integer *iu, floatreal *abstol, integer *m, floatreal *w, floatreal *z__, integer *
    ldz, floatreal *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssprfs)(char *uplo, integer *n, integer *nrhs, floatreal *ap,
    floatreal *afp, integer *ipiv, floatreal *b, integer *ldb, floatreal *x, integer *
    ldx, floatreal *ferr, floatreal *berr, floatreal *work, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(sspsv)(char *uplo, integer *n, integer *nrhs, floatreal *ap,
    integer *ipiv, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatreal *ap, floatreal *afp, integer *ipiv, floatreal *b, integer *ldb, floatreal
    *x, integer *ldx, floatreal *rcond, floatreal *ferr, floatreal *berr, floatreal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ssptrd)(char *uplo, integer *n, floatreal *ap, floatreal *d__,
    floatreal *e, floatreal *tau, integer *info);

/* Subroutine */ int F77NAME(ssptrf)(char *uplo, integer *n, floatreal *ap, integer *ipiv,
    integer *info);

/* Subroutine */ int F77NAME(ssptri)(char *uplo, integer *n, floatreal *ap, integer *ipiv,
    floatreal *work, integer *info);

/* Subroutine */ int F77NAME(ssptrs)(char *uplo, integer *n, integer *nrhs, floatreal *ap,
    integer *ipiv, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(sstebz)(char *range, char *order, integer *n, floatreal *vl,
    floatreal *vu, integer *il, integer *iu, floatreal *abstol, floatreal *d__, floatreal *e,
    integer *m, integer *nsplit, floatreal *w, integer *iblock, integer *
    isplit, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(sstedc)(char *compz, integer *n, floatreal *d__, floatreal *e,
    floatreal *z__, integer *ldz, floatreal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstegr)(char *jobz, char *range, integer *n, floatreal *d__,
    floatreal *e, floatreal *vl, floatreal *vu, integer *il, integer *iu, floatreal *abstol,
    integer *m, floatreal *w, floatreal *z__, integer *ldz, integer *isuppz, floatreal *
    work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstein)(integer *n, floatreal *d__, floatreal *e, integer *m, floatreal
    *w, integer *iblock, integer *isplit, floatreal *z__, integer *ldz, floatreal *
    work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssteqr)(char *compz, integer *n, floatreal *d__, floatreal *e,
    floatreal *z__, integer *ldz, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(ssterf)(integer *n, floatreal *d__, floatreal *e, integer *info);

/* Subroutine */ int F77NAME(sstev)(char *jobz, integer *n, floatreal *d__, floatreal *e, floatreal *
    z__, integer *ldz, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(sstevd)(char *jobz, integer *n, floatreal *d__, floatreal *e, floatreal
    *z__, integer *ldz, floatreal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstevr)(char *jobz, char *range, integer *n, floatreal *d__,
    floatreal *e, floatreal *vl, floatreal *vu, integer *il, integer *iu, floatreal *abstol,
    integer *m, floatreal *w, floatreal *z__, integer *ldz, integer *isuppz, floatreal *
    work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(sstevx)(char *jobz, char *range, integer *n, floatreal *d__,
    floatreal *e, floatreal *vl, floatreal *vu, integer *il, integer *iu, floatreal *abstol,
    integer *m, floatreal *w, floatreal *z__, integer *ldz, floatreal *work, integer *
    iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssycon)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *ipiv, floatreal *anorm, floatreal *rcond, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(ssyev)(char *jobz, char *uplo, integer *n, floatreal *a,
    integer *lda, floatreal *w, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssyevd)(char *jobz, char *uplo, integer *n, floatreal *a,
    integer *lda, floatreal *w, floatreal *work, integer *lwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssyevr)(char *jobz, char *range, char *uplo, integer *n,
    floatreal *a, integer *lda, floatreal *vl, floatreal *vu, integer *il, integer *iu,
    floatreal *abstol, integer *m, floatreal *w, floatreal *z__, integer *ldz, integer *
    isuppz, floatreal *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ssyevx)(char *jobz, char *range, char *uplo, integer *n,
    floatreal *a, integer *lda, floatreal *vl, floatreal *vu, integer *il, integer *iu,
    floatreal *abstol, integer *m, floatreal *w, floatreal *z__, integer *ldz, floatreal *
    work, integer *lwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssygs2)(integer *itype, char *uplo, integer *n, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ssygst)(integer *itype, char *uplo, integer *n, floatreal *a,
    integer *lda, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ssygv)(integer *itype, char *jobz, char *uplo, integer *
    n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *w, floatreal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssygvd)(integer *itype, char *jobz, char *uplo, integer *
    n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *w, floatreal *work,
    integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(ssygvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *
    vl, floatreal *vu, integer *il, integer *iu, floatreal *abstol, integer *m,
    floatreal *w, floatreal *z__, integer *ldz, floatreal *work, integer *lwork, integer
    *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(ssyrfs)(char *uplo, integer *n, integer *nrhs, floatreal *a,
    integer *lda, floatreal *af, integer *ldaf, integer *ipiv, floatreal *b,
    integer *ldb, floatreal *x, integer *ldx, floatreal *ferr, floatreal *berr, floatreal *
    work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ssysv)(char *uplo, integer *n, integer *nrhs, floatreal *a,
    integer *lda, integer *ipiv, floatreal *b, integer *ldb, floatreal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, floatreal *a, integer *lda, floatreal *af, integer *ldaf, integer *ipiv,
    floatreal *b, integer *ldb, floatreal *x, integer *ldx, floatreal *rcond, floatreal *ferr,
     floatreal *berr, floatreal *work, integer *lwork, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(ssytd2)(char *uplo, integer *n, floatreal *a, integer *lda,
    floatreal *d__, floatreal *e, floatreal *tau, integer *info);

/* Subroutine */ int F77NAME(ssytf2)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(ssytrd)(char *uplo, integer *n, floatreal *a, integer *lda,
    floatreal *d__, floatreal *e, floatreal *tau, floatreal *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(ssytrf)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *ipiv, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ssytri)(char *uplo, integer *n, floatreal *a, integer *lda,
    integer *ipiv, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(ssytrs)(char *uplo, integer *n, integer *nrhs, floatreal *a,
    integer *lda, integer *ipiv, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(stbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, floatreal *ab, integer *ldab, floatreal *rcond, floatreal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, floatreal *ab, integer *ldab, floatreal *b, integer
    *ldb, floatreal *x, integer *ldx, floatreal *ferr, floatreal *berr, floatreal *work,
    integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, floatreal *ab, integer *ldab, floatreal *b, integer
    *ldb, integer *info);

/* Subroutine */ int F77NAME(stgevc)(char *side, char *howmny, logical *select,
    integer *n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *vl,
    integer *ldvl, floatreal *vr, integer *ldvr, integer *mm, integer *m, floatreal
    *work, integer *info);

/* Subroutine */ int F77NAME(stgex2)(logical *wantq, logical *wantz, integer *n, floatreal
    *a, integer *lda, floatreal *b, integer *ldb, floatreal *q, integer *ldq, floatreal *
    z__, integer *ldz, integer *j1, integer *n1, integer *n2, floatreal *work,
    integer *lwork, integer *info);

/* Subroutine */ int F77NAME(stgexc)(logical *wantq, logical *wantz, integer *n, floatreal
    *a, integer *lda, floatreal *b, integer *ldb, floatreal *q, integer *ldq, floatreal *
    z__, integer *ldz, integer *ifst, integer *ilst, floatreal *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(stgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, floatreal *a, integer *lda, floatreal *b, integer *
    ldb, floatreal *alphar, floatreal *alphai, floatreal *beta, floatreal *q, integer *ldq,
    floatreal *z__, integer *ldz, integer *m, floatreal *pl, floatreal *pr, floatreal *dif,
    floatreal *work, integer *lwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(stgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, floatreal *a, integer *lda,
     floatreal *b, integer *ldb, floatreal *tola, floatreal *tolb, floatreal *alpha, floatreal *
    beta, floatreal *u, integer *ldu, floatreal *v, integer *ldv, floatreal *q, integer *
    ldq, floatreal *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(stgsna)(char *job, char *howmny, logical *select,
    integer *n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *vl,
    integer *ldvl, floatreal *vr, integer *ldvr, floatreal *s, floatreal *dif, integer *
    mm, integer *m, floatreal *work, integer *lwork, integer *iwork, integer *
    info);

/* Subroutine */ int F77NAME(stgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *c__, integer *
    ldc, floatreal *d__, integer *ldd, floatreal *e, integer *lde, floatreal *f, integer
    *ldf, floatreal *scale, floatreal *rdsum, floatreal *rdscal, integer *iwork, integer
    *pq, integer *info);

/* Subroutine */ int F77NAME(stgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *c__, integer *
    ldc, floatreal *d__, integer *ldd, floatreal *e, integer *lde, floatreal *f, integer
    *ldf, floatreal *scale, floatreal *dif, floatreal *work, integer *lwork, integer *
    iwork, integer *info);

/* Subroutine */ int F77NAME(stpcon)(char *norm, char *uplo, char *diag, integer *n,
    floatreal *ap, floatreal *rcond, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatreal *ap, floatreal *b, integer *ldb, floatreal *x, integer *ldx,
     floatreal *ferr, floatreal *berr, floatreal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(stptri)(char *uplo, char *diag, integer *n, floatreal *ap,
    integer *info);

/* Subroutine */ int F77NAME(stptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatreal *ap, floatreal *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(strcon)(char *norm, char *uplo, char *diag, integer *n,
    floatreal *a, integer *lda, floatreal *rcond, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(strevc)(char *side, char *howmny, logical *select,
    integer *n, floatreal *t, integer *ldt, floatreal *vl, integer *ldvl, floatreal *vr,
    integer *ldvr, integer *mm, integer *m, floatreal *work, integer *info);

/* Subroutine */ int F77NAME(strexc)(char *compq, integer *n, floatreal *t, integer *ldt,
    floatreal *q, integer *ldq, integer *ifst, integer *ilst, floatreal *work,
    integer *info);

/* Subroutine */ int F77NAME(strrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *x,
    integer *ldx, floatreal *ferr, floatreal *berr, floatreal *work, integer *iwork,
    integer *info);

/* Subroutine */ int F77NAME(strsen)(char *job, char *compq, logical *select, integer
    *n, floatreal *t, integer *ldt, floatreal *q, integer *ldq, floatreal *wr, floatreal *wi,
    integer *m, floatreal *s, floatreal *sep, floatreal *work, integer *lwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(strsna)(char *job, char *howmny, logical *select,
    integer *n, floatreal *t, integer *ldt, floatreal *vl, integer *ldvl, floatreal *vr,
    integer *ldvr, floatreal *s, floatreal *sep, integer *mm, integer *m, floatreal *
    work, integer *ldwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(strsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, floatreal *a, integer *lda, floatreal *b, integer *ldb, floatreal *
    c__, integer *ldc, floatreal *scale, integer *info);

/* Subroutine */ int F77NAME(strti2)(char *uplo, char *diag, integer *n, floatreal *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(strtri)(char *uplo, char *diag, integer *n, floatreal *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(strtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, floatreal *a, integer *lda, floatreal *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(stzrqf)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, integer *info);

/* Subroutine */ int F77NAME(stzrzf)(integer *m, integer *n, floatreal *a, integer *lda,
    floatreal *tau, floatreal *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(xerbla)(char *srname, integer *info);

/* Subroutine */ int F77NAME(zbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
    nru, integer *ncc, doublereal *d__, doublereal *e, doublecomplex *vt,
    integer *ldvt, doublecomplex *u, integer *ldu, doublecomplex *c__,
    integer *ldc, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zdrot)(integer *n, doublecomplex *cx, integer *incx,
    doublecomplex *cy, integer *incy, doublereal *c__, doublereal *s);

/* Subroutine */ int F77NAME(zdrscl)(integer *n, doublereal *sa, doublecomplex *sx,
    integer *incx);

/* Subroutine */ int F77NAME(zgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
     integer *kl, integer *ku, doublecomplex *ab, integer *ldab,
    doublereal *d__, doublereal *e, doublecomplex *q, integer *ldq,
    doublecomplex *pt, integer *ldpt, doublecomplex *c__, integer *ldc,
    doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgbcon)(char *norm, integer *n, integer *kl, integer *ku,
     doublecomplex *ab, integer *ldab, integer *ipiv, doublereal *anorm,
    doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgbequ)(integer *m, integer *n, integer *kl, integer *ku,
     doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__,
    doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
    info);

/* Subroutine */ int F77NAME(zgbrfs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *
    afb, integer *ldafb, integer *ipiv, doublecomplex *b, integer *ldb,
    doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgbsv)(integer *n, integer *kl, integer *ku, integer *
    nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, doublecomplex *
    b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zgbsvx)(char *fact, char *trans, integer *n, integer *kl,
     integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab,
    doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed,
    doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb,
    doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr,
    doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
     doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
     doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgbtrs)(char *trans, integer *n, integer *kl, integer *
    ku, integer *nrhs, doublecomplex *ab, integer *ldab, integer *ipiv,
    doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zgebak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doublereal *scale, integer *m, doublecomplex *v,
    integer *ldv, integer *info);

/* Subroutine */ int F77NAME(zgebal)(char *job, integer *n, doublecomplex *a, integer
    *lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);

/* Subroutine */ int F77NAME(zgebd2)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq,
    doublecomplex *taup, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zgebrd)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq,
    doublecomplex *taup, doublecomplex *work, integer *lwork, integer *
    info);

/* Subroutine */ int F77NAME(zgecon)(char *norm, integer *n, doublecomplex *a,
    integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeequ)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd,
    doublereal *colcnd, doublereal *amax, integer *info);

/* Subroutine */ int F77NAME(zgees)(char *jobvs, char *sort, L_fp select, integer *n,
    doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w,
    doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork,
     doublereal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zgeesx)(char *jobvs, char *sort, L_fp select, char *
    sense, integer *n, doublecomplex *a, integer *lda, integer *sdim,
    doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *
    rconde, doublereal *rcondv, doublecomplex *work, integer *lwork,
    doublereal *rwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zgeev)(char *jobvl, char *jobvr, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl,
    integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work,
    integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *w,
    doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr,
    integer *ilo, integer *ihi, doublereal *scale, doublereal *abnrm,
    doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *
    lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgegs)(char *jobvsl, char *jobvsr, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl,
    integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *
    work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgegv)(char *jobvl, char *jobvr, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer
    *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer
    *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgehd2)(integer *n, integer *ilo, integer *ihi,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zgehrd)(integer *n, integer *ilo, integer *ihi,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgelq2)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zgelqf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgels)(char *trans, integer *m, integer *n, integer *
    nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgelsx)(integer *m, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work,
    doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgelsy)(integer *m, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work,
    integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeql2)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zgeqlf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgeqp3)(integer *m, integer *n, doublecomplex *a,
    integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work,
    integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeqpf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work,
    doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgeqr2)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zgeqrf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgerfs)(char *trans, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf,
    integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
    integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
     doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgerq2)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zgerqf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zgesc2)(integer *n, doublecomplex *a, integer *lda,
    doublecomplex *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

/* Subroutine */ int F77NAME(zgesv)(integer *n, integer *nrhs, doublecomplex *a,
    integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *
    info);

/* Subroutine */ int F77NAME(zgesvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
    ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__,
    doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgetc2)(integer *n, doublecomplex *a, integer *lda,
    integer *ipiv, integer *jpiv, integer *info);

/* Subroutine */ int F77NAME(zgetf2)(integer *m, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgetrf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zgetri)(integer *n, doublecomplex *a, integer *lda,
    integer *ipiv, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgetrs)(char *trans, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zggbak)(char *job, char *side, integer *n, integer *ilo,
    integer *ihi, doublereal *lscale, doublereal *rscale, integer *m,
    doublecomplex *v, integer *ldv, integer *info);

/* Subroutine */ int F77NAME(zggbal)(char *job, integer *n, doublecomplex *a, integer
    *lda, doublecomplex *b, integer *ldb, integer *ilo, integer *ihi,
    doublereal *lscale, doublereal *rscale, doublereal *work, integer *
    info);

/* Subroutine */ int F77NAME(zgges)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, integer *n, doublecomplex *a, integer *lda, doublecomplex *b,
    integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *
    beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer
    *ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork,
    logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
    delctg, char *sense, integer *n, doublecomplex *a, integer *lda,
    doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha,
    doublecomplex *beta, doublecomplex *vsl, integer *ldvsl,
    doublecomplex *vsr, integer *ldvsr, doublereal *rconde, doublereal *
    rcondv, doublecomplex *work, integer *lwork, doublereal *rwork,
    integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggev)(char *jobvl, char *jobvr, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer
    *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer
    *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zggevx)(char *balanc, char *jobvl, char *jobvr, char *
    sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *b,
    integer *ldb, doublecomplex *alpha, doublecomplex *beta,
    doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr,
    integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale,
    doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
    rcondv, doublecomplex *work, integer *lwork, doublereal *rwork,
    integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int F77NAME(zggglm)(integer *n, integer *m, integer *p,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *d__, doublecomplex *x, doublecomplex *y, doublecomplex
    *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zgghrd)(char *compq, char *compz, integer *n, integer *
    ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *b,
    integer *ldb, doublecomplex *q, integer *ldq, doublecomplex *z__,
    integer *ldz, integer *info);

/* Subroutine */ int F77NAME(zgglse)(integer *m, integer *n, integer *p,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *c__, doublecomplex *d__, doublecomplex *x,
    doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zggqrf)(integer *n, integer *m, integer *p,
    doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b,
     integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zggrqf)(integer *m, integer *p, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b,
     integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zggsvd)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *n, integer *p, integer *k, integer *l, doublecomplex *a,
    integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha,
    doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v,
    integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work,
    doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zggsvp)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex
    *b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k,
    integer *l, doublecomplex *u, integer *ldu, doublecomplex *v, integer
    *ldv, doublecomplex *q, integer *ldq, integer *iwork, doublereal *
    rwork, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zgtcon)(char *norm, integer *n, doublecomplex *dl,
    doublecomplex *d__, doublecomplex *du, doublecomplex *du2, integer *
    ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work,
    integer *info);

/* Subroutine */ int F77NAME(zgtrfs)(char *trans, integer *n, integer *nrhs,
    doublecomplex *dl, doublecomplex *d__, doublecomplex *du,
    doublecomplex *dlf, doublecomplex *df, doublecomplex *duf,
    doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb,
    doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zgtsv)(integer *n, integer *nrhs, doublecomplex *dl,
    doublecomplex *d__, doublecomplex *du, doublecomplex *b, integer *ldb,
     integer *info);

/* Subroutine */ int F77NAME(zgtsvx)(char *fact, char *trans, integer *n, integer *
    nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du,
    doublecomplex *dlf, doublecomplex *df, doublecomplex *duf,
    doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb,
    doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr,
    doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zgttrf)(integer *n, doublecomplex *dl, doublecomplex *
    d__, doublecomplex *du, doublecomplex *du2, integer *ipiv, integer *
    info);

/* Subroutine */ int F77NAME(zgttrs)(char *trans, integer *n, integer *nrhs,
    doublecomplex *dl, doublecomplex *d__, doublecomplex *du,
    doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zgtts2)(integer *itrans, integer *n, integer *nrhs,
    doublecomplex *dl, doublecomplex *d__, doublecomplex *du,
    doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zhbev)(char *jobz, char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__,
    integer *ldz, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbevd)(char *jobz, char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__,
    integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork,
    integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhbevx)(char *jobz, char *range, char *uplo, integer *n,
    integer *kd, doublecomplex *ab, integer *ldab, doublecomplex *q,
    integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *
    iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
     integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork,
     integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zhbgst)(char *vect, char *uplo, integer *n, integer *ka,
    integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb,
    integer *ldbb, doublecomplex *x, integer *ldx, doublecomplex *work,
    doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbgv)(char *jobz, char *uplo, integer *n, integer *ka,
    integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb,
    integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz,
    doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhbgvx)(char *jobz, char *range, char *uplo, integer *n,
    integer *ka, integer *kb, doublecomplex *ab, integer *ldab,
    doublecomplex *bb, integer *ldbb, doublecomplex *q, integer *ldq,
    doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *
    abstol, integer *m, doublereal *w, doublecomplex *z__, integer *ldz,
    doublecomplex *work, doublereal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(zhbtrd)(char *vect, char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, doublereal *d__, doublereal *e,
    doublecomplex *q, integer *ldq, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zhecon)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond,
    doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zheev)(char *jobz, char *uplo, integer *n, doublecomplex
    *a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork,
    doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zheevd)(char *jobz, char *uplo, integer *n,
    doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work,
    integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zheevr)(char *jobz, char *range, char *uplo, integer *n,
    doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu,
    integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
    w, doublecomplex *z__, integer *ldz, integer *isuppz, doublecomplex *
    work, integer *lwork, doublereal *rwork, integer *lrwork, integer *
    iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zheevx)(char *jobz, char *range, char *uplo, integer *n,
    doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu,
    integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
    w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *
    lwork, doublereal *rwork, integer *iwork, integer *ifail, integer *
    info);

/* Subroutine */ int F77NAME(zhegs2)(integer *itype, char *uplo, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhegst)(integer *itype, char *uplo, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhegv)(integer *itype, char *jobz, char *uplo, integer *
    n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
     integer *info);

/* Subroutine */ int F77NAME(zhegvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
     integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhegvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b,
    integer *ldb, doublereal *vl, doublereal *vu, integer *il, integer *
    iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
     integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork,
     integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zherfs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf,
    integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
    integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
     doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhesv)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b,
    integer *ldb, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zhesvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
    ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
     integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr,
    doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhetf2)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zhetrd)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau,
    doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zhetrf)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, doublecomplex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(zhetri)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zhetrs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zhgeqz)(char *job, char *compq, char *compz, integer *n,
    integer *ilo, integer *ihi, doublecomplex *a, integer *lda,
    doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
    beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
    ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpcon)(char *uplo, integer *n, doublecomplex *ap,
    integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zhpev)(char *jobz, char *uplo, integer *n, doublecomplex
    *ap, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhpevd)(char *jobz, char *uplo, integer *n,
    doublecomplex *ap, doublereal *w, doublecomplex *z__, integer *ldz,
    doublecomplex *work, integer *lwork, doublereal *rwork, integer *
    lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zhpevx)(char *jobz, char *range, char *uplo, integer *n,
    doublecomplex *ap, doublereal *vl, doublereal *vu, integer *il,
    integer *iu, doublereal *abstol, integer *m, doublereal *w,
    doublecomplex *z__, integer *ldz, doublecomplex *work, doublereal *
    rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zhpgst)(integer *itype, char *uplo, integer *n,
    doublecomplex *ap, doublecomplex *bp, integer *info);

/* Subroutine */ int F77NAME(zhpgv)(integer *itype, char *jobz, char *uplo, integer *
    n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex
    *z__, integer *ldz, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpgvd)(integer *itype, char *jobz, char *uplo, integer *
    n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex
    *z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *
    rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpgvx)(integer *itype, char *jobz, char *range, char *
    uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *
    vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol,
    integer *m, doublereal *w, doublecomplex *z__, integer *ldz,
    doublecomplex *work, doublereal *rwork, integer *iwork, integer *
    ifail, integer *info);

/* Subroutine */ int F77NAME(zhprfs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
    b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr,
    doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zhpsv)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhpsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv,
    doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zhptrd)(char *uplo, integer *n, doublecomplex *ap,
    doublereal *d__, doublereal *e, doublecomplex *tau, integer *info);

/* Subroutine */ int F77NAME(zhptrf)(char *uplo, integer *n, doublecomplex *ap,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zhptri)(char *uplo, integer *n, doublecomplex *ap,
    integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zhptrs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zhsein)(char *side, char *eigsrc, char *initv, logical *
    select, integer *n, doublecomplex *h__, integer *ldh, doublecomplex *
    w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr,
     integer *mm, integer *m, doublecomplex *work, doublereal *rwork,
    integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int F77NAME(zhseqr)(char *job, char *compz, integer *n, integer *ilo,
     integer *ihi, doublecomplex *h__, integer *ldh, doublecomplex *w,
    doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zlabrd)(integer *m, integer *n, integer *nb,
    doublecomplex *a, integer *lda, doublereal *d__, doublereal *e,
    doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, integer *
    ldx, doublecomplex *y, integer *ldy);

/* Subroutine */ int F77NAME(zlacgv)(integer *n, doublecomplex *x, integer *incx);

/* Subroutine */ int F77NAME(zlacon)(integer *n, doublecomplex *v, doublecomplex *x,
    doublereal *est, integer *kase);

/* Subroutine */ int F77NAME(zlacp2)(char *uplo, integer *m, integer *n, doublereal *
    a, integer *lda, doublecomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlacpy)(char *uplo, integer *m, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlacrm)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublereal *b, integer *ldb, doublecomplex *c__,
    integer *ldc, doublereal *rwork);

/* Subroutine */ int F77NAME(zlacrt)(integer *n, doublecomplex *cx, integer *incx,
    doublecomplex *cy, integer *incy, doublecomplex *c__, doublecomplex *
    s);

/* Subroutine */ int F77NAME(zlaed0)(integer *qsiz, integer *n, doublereal *d__,
    doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore,
    integer *ldqs, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlaed7)(integer *n, integer *cutpnt, integer *qsiz,
    integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__,
    doublecomplex *q, integer *ldq, doublereal *rho, integer *indxq,
    doublereal *qstore, integer *qptr, integer *prmptr, integer *perm,
    integer *givptr, integer *givcol, doublereal *givnum, doublecomplex *
    work, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlaed8)(integer *k, integer *n, integer *qsiz,
    doublecomplex *q, integer *ldq, doublereal *d__, doublereal *rho,
    integer *cutpnt, doublereal *z__, doublereal *dlamda, doublecomplex *
    q2, integer *ldq2, doublereal *w, integer *indxp, integer *indx,
    integer *indxq, integer *perm, integer *givptr, integer *givcol,
    doublereal *givnum, integer *info);

/* Subroutine */ int F77NAME(zlaein)(logical *rightv, logical *noinit, integer *n,
    doublecomplex *h__, integer *ldh, doublecomplex *w, doublecomplex *v,
    doublecomplex *b, integer *ldb, doublereal *rwork, doublereal *eps3,
    doublereal *smlnum, integer *info);

/* Subroutine */ int F77NAME(zlaesy)(doublecomplex *a, doublecomplex *b,
    doublecomplex *c__, doublecomplex *rt1, doublecomplex *rt2,
    doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);

/* Subroutine */ int F77NAME(zlaev2)(doublecomplex *a, doublecomplex *b,
    doublecomplex *c__, doublereal *rt1, doublereal *rt2, doublereal *cs1,
     doublecomplex *sn1);

/* Subroutine */ int F77NAME(zlags2)(logical *upper, doublereal *a1, doublecomplex *
    a2, doublereal *a3, doublereal *b1, doublecomplex *b2, doublereal *b3,
     doublereal *csu, doublecomplex *snu, doublereal *csv, doublecomplex *
    snv, doublereal *csq, doublecomplex *snq);

/* Subroutine */ int F77NAME(zlagtm)(char *trans, integer *n, integer *nrhs,
    doublereal *alpha, doublecomplex *dl, doublecomplex *d__,
    doublecomplex *du, doublecomplex *x, integer *ldx, doublereal *beta,
    doublecomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zlahef)(char *uplo, integer *n, integer *nb, integer *kb,
     doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w,
    integer *ldw, integer *info);

/* Subroutine */ int F77NAME(zlahqr)(logical *wantt, logical *wantz, integer *n,
    integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh,
    doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__,
    integer *ldz, integer *info);

/* Subroutine */ int F77NAME(zlahrd)(integer *n, integer *k, integer *nb,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *t,
    integer *ldt, doublecomplex *y, integer *ldy);

/* Subroutine */ int F77NAME(zlaic1)(integer *job, integer *j, doublecomplex *x,
    doublereal *sest, doublecomplex *w, doublecomplex *gamma, doublereal *
    sestpr, doublecomplex *s, doublecomplex *c__);

/* Subroutine */ int F77NAME(zlals0)(integer *icompq, integer *nl, integer *nr,
    integer *sqre, integer *nrhs, doublecomplex *b, integer *ldb,
    doublecomplex *bx, integer *ldbx, integer *perm, integer *givptr,
    integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum,
     doublereal *poles, doublereal *difl, doublereal *difr, doublereal *
    z__, integer *k, doublereal *c__, doublereal *s, doublereal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(zlalsa)(integer *icompq, integer *smlsiz, integer *n,
    integer *nrhs, doublecomplex *b, integer *ldb, doublecomplex *bx,
    integer *ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *
    k, doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
    poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
    perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
    rwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlapll)(integer *n, doublecomplex *x, integer *incx,
    doublecomplex *y, integer *incy, doublereal *ssmin);

/* Subroutine */ int F77NAME(zlapmt)(logical *forwrd, integer *m, integer *n,
    doublecomplex *x, integer *ldx, integer *k);

/* Subroutine */ int F77NAME(zlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
     doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__,
    doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqge)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd,
    doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqhb)(char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond,
    doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqhe)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, doublereal *s, doublereal *scond, doublereal *amax,
    char *equed);

/* Subroutine */ int F77NAME(zlaqhp)(char *uplo, integer *n, doublecomplex *ap,
    doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqp2)(integer *m, integer *n, integer *offset,
    doublecomplex *a, integer *lda, integer *jpvt, doublecomplex *tau,
    doublereal *vn1, doublereal *vn2, doublecomplex *work);

/* Subroutine */ int F77NAME(zlaqps)(integer *m, integer *n, integer *offset, integer
    *nb, integer *kb, doublecomplex *a, integer *lda, integer *jpvt,
    doublecomplex *tau, doublereal *vn1, doublereal *vn2, doublecomplex *
    auxv, doublecomplex *f, integer *ldf);

/* Subroutine */ int F77NAME(zlaqsb)(char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond,
    doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqsp)(char *uplo, integer *n, doublecomplex *ap,
    doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int F77NAME(zlaqsy)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, doublereal *s, doublereal *scond, doublereal *amax,
    char *equed);

/* Subroutine */ int F77NAME(zlar1v)(integer *n, integer *b1, integer *bn, doublereal
    *sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
    lld, doublereal *gersch, doublecomplex *z__, doublereal *ztz,
    doublereal *mingma, integer *r__, integer *isuppz, doublereal *work);

/* Subroutine */ int F77NAME(zlar2v)(integer *n, doublecomplex *x, doublecomplex *y,
    doublecomplex *z__, integer *incx, doublereal *c__, doublecomplex *s,
    integer *incc);

/* Subroutine */ int F77NAME(zlarcm)(integer *m, integer *n, doublereal *a, integer *
    lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc,
     doublereal *rwork);

/* Subroutine */ int F77NAME(zlarf)(char *side, integer *m, integer *n, doublecomplex
    *v, integer *incv, doublecomplex *tau, doublecomplex *c__, integer *
    ldc, doublecomplex *work);

/* Subroutine */ int F77NAME(zlarfb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, doublecomplex *v, integer
    *ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *
    ldc, doublecomplex *work, integer *ldwork);

/* Subroutine */ int F77NAME(zlarfg)(integer *n, doublecomplex *alpha, doublecomplex *
    x, integer *incx, doublecomplex *tau);

/* Subroutine */ int F77NAME(zlarft)(char *direct, char *storev, integer *n, integer *
    k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
    t, integer *ldt);

/* Subroutine */ int F77NAME(zlarfx)(char *side, integer *m, integer *n,
    doublecomplex *v, doublecomplex *tau, doublecomplex *c__, integer *
    ldc, doublecomplex *work);

/* Subroutine */ int F77NAME(zlargv)(integer *n, doublecomplex *x, integer *incx,
    doublecomplex *y, integer *incy, doublereal *c__, integer *incc);

/* Subroutine */ int F77NAME(zlarnv)(integer *idist, integer *iseed, integer *n,
    doublecomplex *x);

/* Subroutine */ int F77NAME(zlarrv)(integer *n, doublereal *d__, doublereal *l,
    integer *isplit, integer *m, doublereal *w, integer *iblock,
    doublereal *gersch, doublereal *tol, doublecomplex *z__, integer *ldz,
     integer *isuppz, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(zlartg)(doublecomplex *f, doublecomplex *g, doublereal *
    cs, doublecomplex *sn, doublecomplex *r__);

/* Subroutine */ int F77NAME(zlartv)(integer *n, doublecomplex *x, integer *incx,
    doublecomplex *y, integer *incy, doublereal *c__, doublecomplex *s,
    integer *incc);

/* Subroutine */ int F77NAME(zlarz)(char *side, integer *m, integer *n, integer *l,
    doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
    c__, integer *ldc, doublecomplex *work);

/* Subroutine */ int F77NAME(zlarzb)(char *side, char *trans, char *direct, char *
    storev, integer *m, integer *n, integer *k, integer *l, doublecomplex
    *v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *c__,
    integer *ldc, doublecomplex *work, integer *ldwork);

/* Subroutine */ int F77NAME(zlarzt)(char *direct, char *storev, integer *n, integer *
    k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
    t, integer *ldt);

/* Subroutine */ int F77NAME(zlascl)(char *type__, integer *kl, integer *ku,
    doublereal *cfrom, doublereal *cto, integer *m, integer *n,
    doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(zlaset)(char *uplo, integer *m, integer *n,
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, integer *
    lda);

/* Subroutine */ int F77NAME(zlasr)(char *side, char *pivot, char *direct, integer *m,
     integer *n, doublereal *c__, doublereal *s, doublecomplex *a,
    integer *lda);

/* Subroutine */ int F77NAME(zlassq)(integer *n, doublecomplex *x, integer *incx,
    doublereal *scale, doublereal *sumsq);

/* Subroutine */ int F77NAME(zlaswp)(integer *n, doublecomplex *a, integer *lda,
    integer *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int F77NAME(zlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
     doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w,
    integer *ldw, integer *info);

/* Subroutine */ int F77NAME(zlatbs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, integer *kd, doublecomplex *ab, integer *ldab,
    doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatdf)(integer *ijob, integer *n, doublecomplex *z__,
    integer *ldz, doublecomplex *rhs, doublereal *rdsum, doublereal *
    rdscal, integer *ipiv, integer *jpiv);

/* Subroutine */ int F77NAME(zlatps)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doublecomplex *ap, doublecomplex *x, doublereal *
    scale, doublereal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatrd)(char *uplo, integer *n, integer *nb,
    doublecomplex *a, integer *lda, doublereal *e, doublecomplex *tau,
    doublecomplex *w, integer *ldw);

/* Subroutine */ int F77NAME(zlatrs)(char *uplo, char *trans, char *diag, char *
    normin, integer *n, doublecomplex *a, integer *lda, doublecomplex *x,
    doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int F77NAME(zlatrz)(integer *m, integer *n, integer *l,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work);

/* Subroutine */ int F77NAME(zlatzm)(char *side, integer *m, integer *n,
    doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
    c1, doublecomplex *c2, integer *ldc, doublecomplex *work);

/* Subroutine */ int F77NAME(zlauu2)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zlauum)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpbcon)(char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, doublereal *anorm, doublereal *
    rcond, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpbequ)(char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond,
    doublereal *amax, integer *info);

/* Subroutine */ int F77NAME(zpbrfs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *
    ldafb, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
     doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zpbstf)(char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbsv)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(zpbsvx)(char *fact, char *uplo, integer *n, integer *kd,
    integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb,
    integer *ldafb, char *equed, doublereal *s, doublecomplex *b, integer
    *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *
    ferr, doublereal *berr, doublecomplex *work, doublereal *rwork,
    integer *info);

/* Subroutine */ int F77NAME(zpbtf2)(char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbtrf)(char *uplo, integer *n, integer *kd,
    doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int F77NAME(zpbtrs)(char *uplo, integer *n, integer *kd, integer *
    nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
    ldb, integer *info);

/* Subroutine */ int F77NAME(zpocon)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpoequ)(integer *n, doublecomplex *a, integer *lda,
    doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int F77NAME(zporfs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf,
    doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zposv)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zposvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
    ldaf, char *equed, doublereal *s, doublecomplex *b, integer *ldb,
    doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr,
    doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zpotf2)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpotrf)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpotri)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *info);

/* Subroutine */ int F77NAME(zpotrs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zppcon)(char *uplo, integer *n, doublecomplex *ap,
    doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(zppequ)(char *uplo, integer *n, doublecomplex *ap,
    doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int F77NAME(zpprfs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, doublecomplex *afp, doublecomplex *b, integer *ldb,
     doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zppsv)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zppsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublecomplex *ap, doublecomplex *afp, char *equed, doublereal *
    s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpptrf)(char *uplo, integer *n, doublecomplex *ap,
    integer *info);

/* Subroutine */ int F77NAME(zpptri)(char *uplo, integer *n, doublecomplex *ap,
    integer *info);

/* Subroutine */ int F77NAME(zpptrs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zptcon)(integer *n, doublereal *d__, doublecomplex *e,
    doublereal *anorm, doublereal *rcond, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zptrfs)(char *uplo, integer *n, integer *nrhs,
    doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef,
    doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(zptsv)(integer *n, integer *nrhs, doublereal *d__,
    doublecomplex *e, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(zptsvx)(char *fact, integer *n, integer *nrhs,
    doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef,
    doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zpttrf)(integer *n, doublereal *d__, doublecomplex *e,
    integer *info);

/* Subroutine */ int F77NAME(zpttrs)(char *uplo, integer *n, integer *nrhs,
    doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zptts2)(integer *iuplo, integer *n, integer *nrhs,
    doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb);

/* Subroutine */ int F77NAME(zrot)(integer *n, doublecomplex *cx, integer *incx,
    doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s);

/* Subroutine */ int F77NAME(zspcon)(char *uplo, integer *n, doublecomplex *ap,
    integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zspmv)(char *uplo, integer *n, doublecomplex *alpha,
    doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *
    beta, doublecomplex *y, integer *incy);

/* Subroutine */ int F77NAME(zspr)(char *uplo, integer *n, doublecomplex *alpha,
    doublecomplex *x, integer *incx, doublecomplex *ap);

/* Subroutine */ int F77NAME(zsprfs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
    b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr,
    doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(zspsv)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zspsvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv,
    doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsptrf)(char *uplo, integer *n, doublecomplex *ap,
    integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zsptri)(char *uplo, integer *n, doublecomplex *ap,
    integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zsptrs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(zstedc)(char *compz, integer *n, doublereal *d__,
    doublereal *e, doublecomplex *z__, integer *ldz, doublecomplex *work,
    integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork,
    integer *liwork, integer *info);

/* Subroutine */ int F77NAME(zstein)(integer *n, doublereal *d__, doublereal *e,
    integer *m, doublereal *w, integer *iblock, integer *isplit,
    doublecomplex *z__, integer *ldz, doublereal *work, integer *iwork,
    integer *ifail, integer *info);

/* Subroutine */ int F77NAME(zsteqr)(char *compz, integer *n, doublereal *d__,
    doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work,
    integer *info);

/* Subroutine */ int F77NAME(zsycon)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond,
    doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zsymv)(char *uplo, integer *n, doublecomplex *alpha,
    doublecomplex *a, integer *lda, doublecomplex *x, integer *incx,
    doublecomplex *beta, doublecomplex *y, integer *incy);

/* Subroutine */ int F77NAME(zsyr)(char *uplo, integer *n, doublecomplex *alpha,
    doublecomplex *x, integer *incx, doublecomplex *a, integer *lda);

/* Subroutine */ int F77NAME(zsyrfs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf,
    integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
    integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
     doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsysv)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b,
    integer *ldb, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zsysvx)(char *fact, char *uplo, integer *n, integer *
    nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
    ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
     integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr,
    doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(zsytf2)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int F77NAME(zsytrf)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, doublecomplex *work, integer *lwork,
    integer *info);

/* Subroutine */ int F77NAME(zsytri)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zsytrs)(char *uplo, integer *n, integer *nrhs,
    doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ztbcon)(char *norm, char *uplo, char *diag, integer *n,
    integer *kd, doublecomplex *ab, integer *ldab, doublereal *rcond,
    doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztbrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab,
    doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
    doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
    rwork, integer *info);

/* Subroutine */ int F77NAME(ztbtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab,
    doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ztgevc)(char *side, char *howmny, logical *select,
    integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer
    *ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
    ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork,
     integer *info);

/* Subroutine */ int F77NAME(ztgex2)(logical *wantq, logical *wantz, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz,
    integer *j1, integer *info);

/* Subroutine */ int F77NAME(ztgexc)(logical *wantq, logical *wantz, integer *n,
    doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz,
    integer *ifst, integer *ilst, integer *info);

/* Subroutine */ int F77NAME(ztgsen)(integer *ijob, logical *wantq, logical *wantz,
    logical *select, integer *n, doublecomplex *a, integer *lda,
    doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
    beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
    ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif,
    doublecomplex *work, integer *lwork, integer *iwork, integer *liwork,
    integer *info);

/* Subroutine */ int F77NAME(ztgsja)(char *jobu, char *jobv, char *jobq, integer *m,
    integer *p, integer *n, integer *k, integer *l, doublecomplex *a,
    integer *lda, doublecomplex *b, integer *ldb, doublereal *tola,
    doublereal *tolb, doublereal *alpha, doublereal *beta, doublecomplex *
    u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q,
    integer *ldq, doublecomplex *work, integer *ncycle, integer *info);

/* Subroutine */ int F77NAME(ztgsna)(char *job, char *howmny, logical *select,
    integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer
    *ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
    ldvr, doublereal *s, doublereal *dif, integer *mm, integer *m,
    doublecomplex *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ztgsy2)(char *trans, integer *ijob, integer *m, integer *
    n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd,
    doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf,
    doublereal *scale, doublereal *rdsum, doublereal *rdscal, integer *
    info);

/* Subroutine */ int F77NAME(ztgsyl)(char *trans, integer *ijob, integer *m, integer *
    n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
    doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd,
    doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf,
    doublereal *scale, doublereal *dif, doublecomplex *work, integer *
    lwork, integer *iwork, integer *info);

/* Subroutine */ int F77NAME(ztpcon)(char *norm, char *uplo, char *diag, integer *n,
    doublecomplex *ap, doublereal *rcond, doublecomplex *work, doublereal
    *rwork, integer *info);

/* Subroutine */ int F77NAME(ztprfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb,
    doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr,
    doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztptri)(char *uplo, char *diag, integer *n,
    doublecomplex *ap, integer *info);

/* Subroutine */ int F77NAME(ztptrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb,
    integer *info);

/* Subroutine */ int F77NAME(ztrcon)(char *norm, char *uplo, char *diag, integer *n,
    doublecomplex *a, integer *lda, doublereal *rcond, doublecomplex *
    work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrevc)(char *side, char *howmny, logical *select,
    integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl,
    integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer
    *m, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrexc)(char *compq, integer *n, doublecomplex *t,
    integer *ldt, doublecomplex *q, integer *ldq, integer *ifst, integer *
    ilst, integer *info);

/* Subroutine */ int F77NAME(ztrrfs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b,
    integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr,
    doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
    info);

/* Subroutine */ int F77NAME(ztrsen)(char *job, char *compq, logical *select, integer
    *n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq,
    doublecomplex *w, integer *m, doublereal *s, doublereal *sep,
    doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(ztrsna)(char *job, char *howmny, logical *select,
    integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl,
    integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s,
    doublereal *sep, integer *mm, integer *m, doublecomplex *work,
    integer *ldwork, doublereal *rwork, integer *info);

/* Subroutine */ int F77NAME(ztrsyl)(char *trana, char *tranb, integer *isgn, integer
    *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b,
    integer *ldb, doublecomplex *c__, integer *ldc, doublereal *scale,
    integer *info);

/* Subroutine */ int F77NAME(ztrti2)(char *uplo, char *diag, integer *n,
    doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(ztrtri)(char *uplo, char *diag, integer *n,
    doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int F77NAME(ztrtrs)(char *uplo, char *trans, char *diag, integer *n,
    integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b,
    integer *ldb, integer *info);

/* Subroutine */ int F77NAME(ztzrqf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, integer *info);

/* Subroutine */ int F77NAME(ztzrzf)(integer *m, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zung2l)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zung2r)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zungbr)(char *vect, integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zunghr)(integer *n, integer *ilo, integer *ihi,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungl2)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zunglq)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungql)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungqr)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungr2)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zungrq)(integer *m, integer *n, integer *k,
    doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zungtr)(char *uplo, integer *n, doublecomplex *a,
    integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunm2l)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zunm2r)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zunmbr)(char *vect, char *side, char *trans, integer *m,
    integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex
    *tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zunmhr)(char *side, char *trans, integer *m, integer *n,
    integer *ilo, integer *ihi, doublecomplex *a, integer *lda,
    doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *
    work, integer *lwork, integer *info);

/* Subroutine */ int F77NAME(zunml2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zunmlq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmql)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmqr)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmr2)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int F77NAME(zunmr3)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex
    *tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
    info);

/* Subroutine */ int F77NAME(zunmrq)(char *side, char *trans, integer *m, integer *n,
    integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zunmrz)(char *side, char *trans, integer *m, integer *n,
    integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex
    *tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
    lwork, integer *info);

/* Subroutine */ int F77NAME(zunmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doublecomplex *a, integer *lda, doublecomplex *tau,
    doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
     integer *info);

/* Subroutine */ int F77NAME(zupgtr)(char *uplo, integer *n, doublecomplex *ap,
    doublecomplex *tau, doublecomplex *q, integer *ldq, doublecomplex *
    work, integer *info);

/* Subroutine */ int F77NAME(zupmtr)(char *side, char *uplo, char *trans, integer *m,
    integer *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c__,
     integer *ldc, doublecomplex *work, integer *info);

} // extern "C"

#endif /* __CLAPACK_H */
