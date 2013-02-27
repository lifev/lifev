/*! Bessel function by CRBond, http://www.crbond.com/
  Algorithms and coefficient values from
  "Computation of Special Functions", Zhang and Jin, John
  Wiley and Sons, 1996.
*/
#include <complex>
#include <cmath>
#include <cstdlib>

#ifndef bessH
#define bessH


#define eps 1e-15
#define el 0.5772156649015329
namespace bessel
{
int msta1 (double x, int mp);
int msta2 (double x, int n, int mp);
int bessjy01a (double x, double& j0, double& j1, double& y0, double& y1,
               double& j0p, double& j1p, double& y0p, double& y1p);
int bessjy01b (double x, double& j0, double& j1, double& y0, double& y1,
               double& j0p, double& j1p, double& y0p, double& y1p);
int bessjyna (int n, double x, int& nm, double* jn, double* yn,
              double* jnp, double* ynp);
int bessjynb (int n, double x, int& nm, double* jn, double* yn,
              double* jnp, double* ynp);
int bessjyv (double v, double x, double& vm, double* jv, double* yv,
             double* jvp, double* yvp);
int bessik01a (double x, double& i0, double& i1, double& k0, double& k1,
               double& i0p, double& i1p, double& k0p, double& k1p);
int bessik01b (double x, double& i0, double& i1, double& k0, double& k1,
               double& i0p, double& i1p, double& k0p, double& k1p);
int bessikna (int n, double x, int& nm, double* in, double* kn,
              double* inp, double* knp);
int bessiknb (int n, double x, int& nm, double* in, double* kn,
              double* inp, double* knp);
int bessikv (double v, double x, double& vm, double* iv, double* kv,
             double* ivp, double* kvp);
int cbessjy01 (std::complex<double> z, std::complex<double>& cj0, std::complex<double>& cj1,
               std::complex<double>& cy0, std::complex<double>& cy1, std::complex<double>& cj0p,
               std::complex<double>& cj1p, std::complex<double>& cy0p, std::complex<double>& cy1p);
int cbessjyna (int n, std::complex<double> z, int& nm, std::complex<double>* cj,
               std::complex<double>* cy, std::complex<double>* cjp, std::complex<double>* cyp);
int cbessjynb (int n, std::complex<double> z, int& nm, std::complex<double>* cj,
               std::complex<double>* cy, std::complex<double>* cjp, std::complex<double>* cyp);
int cbessik01 (std::complex<double>z, std::complex<double>& ci0, std::complex<double>& ci1,
               std::complex<double>& ck0, std::complex<double>& ck1, std::complex<double>& ci0p,
               std::complex<double>& ci1p, std::complex<double>& ck0p, std::complex<double>& ck1p);
int cbessikna (int n, std::complex<double> z, int& nm, std::complex<double>* ci,
               std::complex<double>* ck, std::complex<double>* cip, std::complex<double>* ckp);
int cbessiknb (int n, std::complex<double> z, int& nm, std::complex<double>* ci,
               std::complex<double>* ck, std::complex<double>* cip, std::complex<double>* ckp);
int cbessjyva (double v, std::complex<double> z, double& vm, std::complex<double>* cjv,
               std::complex<double>* cyv, std::complex<double>* cjvp, std::complex<double>* cyvp);
int cbessikv (double v, std::complex<double>z, double& vm, std::complex<double>* civ,
              std::complex<double>* ckv, std::complex<double>* civp, std::complex<double>* ckvp);
}
#endif
