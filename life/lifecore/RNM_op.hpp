  
template<class R>
const KN_<R>& KN_<R>::operator oper (const Mul_KNM_KN_<R> & u)  {
    K_assert (SameShape(u.A.shapei) && !constant());
    R * l(v); KN_<const_R>  li(u.A(0,'.')); //  first line   
    for (int i=0;i<n;i++,l += step,++li)  
      *l oper (li,u.b); 
    return *this;}


template<class R>
const KN_<R>&  KN_<R>::operator oper (const Add_KN_<R> & u) {
    K_assert(u.a.N() == N()  );
    int stepa(u.a.step),stepb(u.b.step);
    R * l(v); const_R  *aa(u.a), *bb(u.b);    
    for (int i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper *aa+*bb;
    return *this;
  }

template<class R>
const KN_<R>&  KN_<R>::operator oper (const Sub_KN_<R> & u) {
    K_assert(u.a.N() == N()  );
    int stepa(u.a.step),stepb(u.b.step);
    R * l(v); const_R  *aa(u.a), *bb(u.b);    
    for (int i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper  *aa-*bb;
    return *this;
  }
  
template<class R>
const KN_<R>&  KN_<R>::operator oper (const Mulc_KN_<R> & u) {
    K_assert(u.a.N() == N()  );
    int stepa(u.a.step);
    R * l(v); const_R  *aa(u.a),bb(u.b)  ;
    for (int i=0;i<n;i++,l += step, aa +=stepa)
      *l oper *aa * bb;
    return *this;
  }
  
template<class R>
const KN_<R>&  KN_<R>::operator oper (const Add_Mulc_KN_<R> & u) {
    K_assert(u.a.N() == N()  );
    const int stepa(u.a.step),stepb(u.b.step);
    const R ca(u.ca),cb(u.cb);    
    R * l(v);
    const R *aa(u.a),*bb(u.b);    
    for (int i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper *aa*ca + *bb*cb;
    return *this;
  }
