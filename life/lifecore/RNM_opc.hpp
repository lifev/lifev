// Modif Miguel 02/12/2002
// "operator=" was changed to "operator oper" when 
// tab are consecutive in memory   

template<class R>
const KN_<R>&  KN_<R>::operator oper (const_R a)  {
    R * l(v);    
    for (int i=0;i<n;i++,l += step) 
      *l oper a;
    return *this;
  }
    
template<class R> 
inline   const KNM_<R> & KNM_<R>::operator oper (const_R a)    
{ 
  if(IsVector1() ) 
        KN_<R>::operator oper (a);
  else {  
          KN_<R>  lj(operator()('.',0)); //  (.,.,O)
          for (int  j=0;j<M();++j,++lj) 
             lj oper a;}       
  return *this;
}

template<class R> 
inline   const KNMK_<R> & KNMK_<R>::operator oper (const_R a)    
{ 
  if(IsVector1() ) 
        KN_<R>::operator oper (a);
  else {  
          KNM_<R>  lj(operator()('.','.',0)); //  (.,.,O)
          int j=K();
           while(j--)
            {lj oper a;++lj;}
       }
  return *this;
}

template<class R>
const KN_<R>&  KN_<R>::operator oper (const KN_<const_R> & u)   {
    K_assert(u.n == n);
    R * l(v);
    const R *r(u);    
    for (int i=0;i<n;i++,l += step, r += u.step) *l oper *r;
    return *this;
  }
  
template<class R> 
inline   const KNM_<R> & KNM_<R>::operator oper (const KNM_<const_R> & u)    
{ 
  if(IsVector1() && u.IsVector1() ) {
        KN_<R>::operator oper (u);
  }
  else {  
          KN_<R>  lj(operator()('.',0)); //  (.,O)
          KN_<const_R>  uj(u('.',0));
          int  j=M();
          while ( j--)
            { lj oper uj;++lj;++uj;} 
        }      
  return *this;
}


template<class R> 
inline   const KNMK_<R> & KNMK_<R>::operator oper (const KNMK_<const_R> & u)    
{ 
  if(IsVector1() && u.IsVector1() ) 
        KN_<R>::operator oper (u);
  else {  
          K_assert( K() == u.K());
          KNM_<R>  lj(operator()('.','.',0)); //  (.,O)
          KNM_<const_R>  uj(u('.','.',0));
          int j=K();
          while (j--)
           { lj oper uj;++lj;++uj;}
       }
  return *this;
}

#undef oper
