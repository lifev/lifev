#ifndef  RNM_tpl_
#define  RNM_tpl_

#include "RNM.hpp"

//   version du 22 nov 1999
//   Voila une debut de class vecteur + matrice 
//   les class avec termine avec un _ ne travail que sur 
//   un pointeur existant => pas de new et delete de pointeur
//   la class correspondant sans de _ genere les pointeurs
//  
//   avec ses classes on peut prendre une ligne 
//   ou une colonne d'une matrice
// -----------------------


template<class R>
void MatMul(KNM_<R> & ab, KNM_<R> &  a, KNM_<R> & b){
  // attention ne marche que si les adresses ne sont pas les memes
  int N= a.shapei.n;
  int M= a.shapej.n;
  K_assert(a.shapej.n == a.shapei.n);
  K_assert(a.shapei.n == ab.shapei.n);
  K_assert(b.shapej.n == ab.shapej.n);
  K_assert(b.v != ab.v);
  K_assert(a.v != ab.v);
  KN_<R> ai =a(0);
  for(int i=1;i<N;i++,++ai){
    KN_<R> bj=b[0];
    for(int j=1;j<M;j++,++bj)
      ab(i,j) = (ai , bj)  ;}
}



inline ostream & operator<<(ostream & f,const ShapeOfArray & s)
  { f << s.n ;
    if(s.step != 1)
      f << ":" << s.step ;
    if (s.step*s.n  !=  s.next ) 
      f << " n: " << setw(3) << s.next ;
    f << ",";
     return f;
   };

template<class R> ostream & operator<<(ostream & f,const KN_<const_R> & v)
  { //f <<  " KN_ : " << (ShapeOfArray) v << " "   <<  (const_R *) v << " :\n\t"  ;
    f << v.N() << "\t:\t" ;
    for (int i=0;i<v.N();i++)
      cout   << setw(3) << v[i] << ((i % 10) == 9 ? "\n\t" : "\t");
     return f;
   };

template<class R> ostream & operator<<(ostream & f,const KNM_<const_R> & v)
  {  //f << " KNM_ "<<v.N()<<"x"<<v.M()<< ": " << (ShapeOfArray) v  
     //<< " i "  << v.shapei
     // << " j "  << v.shapej
     // << " " << &v(0,0) << " :\n\t";
     f << v.N()<<'x'<<v.M() /*<< "  n" << v.next<<" :"<< v.shapei.next << "," << v.shapej.next */<< "\t:\n\t" ;
    for (int i=0;i<v.N();i++) {
      for (int j=0;j<v.M();j++) 
        cout << " " << setw(3) << v(i,j);
       cout << "\n\t";}
  return f;
    
   };

template<class R> ostream & operator<<(ostream & f,const KNMK_<const_R> & v)
  { //f << " KNM_" <<v.N()<<"x"<<v.M()<<"x"<<v.K()<< " : " << (ShapeOfArray) v  
    // << " i "  << v.shapei
    // << " j "  << v.shapej
    // << " k "  << v.shapek << endl;
    // << " " << (void *) & v(0,0,0) << "\n\t" ;
   f << v.N()<< 'x' <<v.M()<< 'x'<<v.K() << "\t:\n\t" ;    
  for (int i=0;i<v.shapei.n;i++){
    for (int j=0;j<v.shapej.n;j++){
      for (int k=0;k<v.shapek.n;k++)
	cout << " " << setw(3) << v(i,j,k);
      cout << "\n\t";}
    cout << "\n\t";}
  return f;
    
   };

template<class R>
 R  KN_<R>::operator,(const KN_<const_R> & u) const {
    K_assert(u.n == n);
    double  s=0; 
    R * l(v);
    R  *r(u.v);    
    for (int i=0;i<n;i++,l += step, r += u.step) s += *l * *r;
    return s;
  }


template<class R>
R  KN_<R>::KNMmin() const {
    R minv = v[index(0)];
    for (int i=1;i<n;i++)
      minv = minv < v[index(i)] ? minv : v[index(i)];
    return minv;
  }
template<class R>
R  KN_<R>::KNMmax() const {
    R maxv = v[index(0)];
    for (int i=1;i<n;i++)
      maxv = maxv > v[index(i)] ? maxv : v[index(i)];
    return maxv;
  }
  
template<class R>
R  KN_<R>::sum() const {
    R s = v[index(0)];
    for (int i=1;i<n;i++)
      s +=  v[index(i)];
    return s;
  }
template<class R>
 KN_<R>&  KN_<R>::map(R (*f)(R )) {
    for (int i=0;i<n;i++)
      {  R & x(v[index(i)]);
          x =  f(x);}
   
  }


///////////////// definition des operateurs d'affectation /////////////////////////
#define oper =
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper +=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper -=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper *=
#include "RNM_opc.hpp"
#define oper /=
#include "RNM_opc.hpp"

#endif
