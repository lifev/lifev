#include "bdf.hpp"


Bdf::Bdf(const UInt n):
  _n(n),_s(0),_alpha(n+1),_beta(n)
{
  ASSERT(n >0 && n <= MAX_ORDER, "Error in bdf order specification");

  switch (n){ 
    case 1:
      _alpha[0] = 1.; // Backward Euler
      _alpha[1] = 1.;
      _beta[0] = 1.; // u^{n+1} \approx u^n
      break;
    case 2:
      _alpha[0] = 3./2.;
      _alpha[1] = 2.;
      _alpha[2] = -1./2.;
      _beta[0] = 2.;
      _beta[1] = -1.;
      break;
    case 3:
      _alpha[0] = 11./6.;
      _alpha[1] = 3.;
      _alpha[2] = -3./2.;
      _alpha[3] = 1./3.;
      _beta[0] = 3.;
      _beta[1] =-3.;
      _beta[2] = 1.;
      break;
 }   
  _unk.resize(n);
}



Bdf::~Bdf() {}


void Bdf::initialize_unk(Vector u0)
{
  vector< Vector >::iterator iter=_unk.begin();
  vector< Vector >::iterator iter_end=_unk.end();

  _s = u0.size();
  
  for (  iter=_unk.begin() ; iter != iter_end; iter++){
    *iter=u0;
  }

  return;
}


void Bdf::initialize_unk(vector<Vector> uv0)
{
  vector< Vector >::iterator iter=_unk.begin();
  vector< Vector >::iterator iter_end=_unk.end();
  
  vector< Vector >::iterator iter0 = uv0.begin();

  _s = iter0->size();

  UInt n0 = uv0.size();

  // Check if uv0 has the right dimensions
  ASSERT(n0<_n,"Initial data are not enough for the selected BDF")

  // if n0>n, only the first n unital data will be considered
  if (_n>n0) cout << "The initial data set is larger than needed by the BDF. Only the first n data will be considered. " << endl; 
  
  for (iter=_unk.begin() ; iter != iter_end; iter++){
    *iter=*iter0;
    iter0++;
  }

  return;
}

//
//
//





vector<Vector> Bdf::unk()
{
  return _unk;
}


double Bdf::coeff_der(UInt i)
{
  // Pay attention: i is c-based indexed 
  ASSERT(i>=0&i<_n+1,"Error in specification of the time derivative coefficient for the BDF formula (out of range error)");
  return _alpha[i];
}

double Bdf::coeff_ext(UInt i)
{
  // Pay attention: i is c-based indexed 
  ASSERT(i>=0&i<_n,"Error in specification of the time derivative coefficient for the BDF formula (out of range error)");
  return _beta[i];
}

void Bdf::showMe()
{
  cout << "*** BDF Time discretization of order " << _n << "***" << endl;
  cout << endl;
  cout << "*** Coefficients: " << endl;
  cout << endl;
  for (UInt i=0;i<_n+1;++i) 
    cout << "alpha("<<i<<") = "<<_alpha[i]<<endl;
  for (UInt i=0;i<_n;++i) 
    cout << "beta("<<i<<") = "<<_beta[i]<<endl;

  cout << "Length unknown vectors:" << _unk.size() << "," << _s << endl;

/*   vector< Vector >::iterator itg=_u_prec.begin(); */
/*   vector< Vector >::iterator itge=_u_prec.end(); */
/*   for ( ; itg!=itge;++itg){ */
/*     cout << "u_prec: " << itg-_u_prec.begin() << endl; */
/*     Vector::iterator itl=itg->begin(); */
/*     Vector::iterator itle=itg->end(); */
/*     cout << "[ "; */
/*     for ( ;itl!=itle; ++itl) */
/*       cout << *itl << ", " ;   */

/*     cout << "]" << endl; */
/*   }  */


  return;
}

void Bdf::shift_right(Vector unk_curr)
{
  vector< Vector >::iterator it=_unk.end()-1;  
  vector< Vector >::iterator itm1=_unk.end()-1;  
  vector< Vector >::iterator itb=_unk.begin();  

  for (; it!=itb; --it){
    itm1--;
    *it=*itm1; 
  }
  *itb=unk_curr;  
 
  return;
}


Vector Bdf::time_der(Real dt)
{
  // ! Compute the right hand side of the time derivative formula, i.e.
  // ! \alpha[0] u^{n+1} = right hand side
  Vector ut(_s);  

  for (UInt j=0;j<_s;++j) ut[j]=0.;

  for (UInt j=0;j<_s;++j){
   for (UInt i=0;i<_n;++i) 
     ut[j]=ut[j]+_alpha[i+1]/dt*_unk[i][j];
     }  

  return ut;
}

Vector Bdf::time_der()
{
  // ! Compute the right hand side of the time derivative formula, i.e.
  // ! \alpha[0] u^{n+1} = right hand side
  // ! In this case, the time step is considered elsewhere (e.g. included in the mass matrix)
  Vector ut(_s);  

  for (UInt j=0;j<_s;++j) ut[j]=0.;

  for (UInt j=0;j<_s;++j){
   for (UInt i=0;i<_n;++i) 
     ut[j]=ut[j]+_alpha[i+1]*_unk[i][j];
     }  

  return ut;
}

Vector Bdf::extrap()
{
  // ! Compute the approximation of u^{n+1}
  // ! extrapolated
  Vector ue(_s);  
  for (UInt j=0;j<_s;++j) ue[j]=0.;

  for (UInt j=0;j<_s;++j){
   for (UInt i=0;i<_n;++i) 
     ue[j]=ue[j]+_beta[i]*_unk[i][j];
     }  

  return ue;
}

