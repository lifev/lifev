#include "tab.hpp"
#include <cmath>

double gamma(double x)
{
  double angle_deg = 90;
  return ( angle_deg/180*M_PI);
}
double der_gamma(double x)
{
  return 0.;
}


KNM<double> fibrous_permea(double visc,KNM<double> perm_gen, double x)
{
  double inv_visc = 1./visc;
  KNM<double> Mik(3,3);
  double s=sin(gamma(x)),c=cos(gamma(x)),c2 = c*c, s2 = s*s,
    u11 = perm_gen(0,0),
    u21 = perm_gen(0,2),
    v   = perm_gen(1,1),
    u12 = perm_gen(2,0),
    u22 = perm_gen(2,2);
  //double res = c2*v + s2*u22;
  //return res;
  //cout << "Le resultat vaut" << res << endl;
  Mik(0,0) =   inv_visc     * u11;
  Mik(0,1) = - inv_visc * s * u21;
  Mik(0,2) =   inv_visc * c * u21;
  Mik(1,0) = - inv_visc * s * u12;
  Mik(1,1) =   inv_visc * (c2*v + s2*u22);
  Mik(1,2) =   inv_visc * c*s * (v - u22);
  Mik(2,0) =   inv_visc * c * u12;
  Mik(2,1) =   Mik(1,2);
  Mik(2,2) =   inv_visc * (c2*u22 + s2*v);
  return Mik;
}


double permeability_sd009(const double& x, const double& y, const double& z) {
  const double x1=790.;
  const double y1=1074.;
  const double x3=832.5;
  const double y3=1099; 
  double slope = (y3 - y1) / (x3 - x1);
  double position = y - y1 - slope * (x - x1);
  double perm;
  if ( position < 0 ) {
    perm = 1.e-12;
  }
  else {
    perm = 3e-4;
  }
  return perm ;
}

double permeability_sd010(const double& x, const double& y, const double& z) {
  const double x1=790.;
  const double y1=1074.;
  const double x3=832.5;
  const double y3=1099; 
  double slope = (y3 - y1) / (x3 - x1);
  double position = y - y1 - slope * (x - x1);
  double perm;
  if ( position < 0 ) {
    perm = 1.e-11;
  }
  else {
    perm = 3e-4;
  }
  return perm ;
}

