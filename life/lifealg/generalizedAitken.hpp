
/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s):
Simone Deparis <simone.deparis@epfl.ch>
Gilles Fourestey <gilles.fourestey@epfl.ch>

Date: 2004-09-23

Copyright (C) 2004 EPFL

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _GENERALIZEDAITKEN_H
#define _GENERALIZEDAITKEN_H

namespace LifeV
{
    
    template<class Vector, class Real>
    class generalizedAitken
    {
        /*
          compute the acceleration with the vector variant of Aitken.
          lk(0,1)  are column vectors (lambda_k lambda_{k-1})
          lk(2) = lk(0)-lk(1)
          Rk(0,1)  are  column vectors of residuals(possibly more than one) .
          Rk(2) = Rk(0)-Rk(1)
          lk(M_k(0)): lambda at the iterate k
          lk(M_k(1)): lambda at the iterate k-1
          Rk: residual(s) associated to lk. in some cases = lk-f(lk)
          defaultOmega: defaultOmega can only be a scalar.
          If defaultOmega < 0, then do not use Aitken, but fixed relaxation
          parameter=abs(defaultOmega)
          If you want to assign a fixed parameter different from 
          the default one, you can use 
          Vector deltaLambda(lam, res, omega)
          be sure that omega has the size sizeOmega!
          If length(Rk1)==0, then Aitken acceleration is not applicable,
          hence we need a default value for omega.
          !: should be of the size 1,size(Rk,2)
          sizeOmega is the number of relaxation parameters to be computed.
          usually sizeOmega=1, but for example fo r a N-N preconditioner, =2.
        */
    public:

        // Constructors
        
        generalizedAitken(const int dof, const Real& defOmega=0.1, const int size=1);


        // Destructor

        ~generalizedAitken();

        // Member functions

        void   restart();

        inline void   setLambdaS  (const Vector&, const int);
        inline Real   getLambdaS  (const int, const int);
        inline void   setLambdaF  (const Vector&, const int);
        inline Real   getLambdaF  (const int, const int);

        inline void   setNuS      (const Vector&, const int);
        inline Real   getNuS      (const int, const int);
        inline void   setNuF      (const Vector&, const int);
        inline Real   getNuF      (const int, const int);

        
        Vector deltaLambda(const Vector& lam, const Vector res[]);
        Vector deltaLambda(const Vector& lam, const Vector res[], Real omega[]);

        // in this case, omega is taken as the default value

    private:

        //! array formed by the vectors [\lamda_s^k, \lamdba_s^{k-1}]
        Vector M_lambdaS;

        //! array formed by the vectors [\lamda_f^k, \lamdba_f^{k-1}]
        Vector M_lambdaF;

        //! array formed by the vectors [\nu_s^k, \nu_s^{k-1}]
        Vector M_nuS;

        //! array formed by the vectors [\nu_f^k, \nu_f^{k-1}]
        Vector M_nuF;

        //! new step coefficients
        Vector M_omega;
        
        Real   M_defaultOmega;

        int    M_sizeOmega;
        int    M_dof;

        int    M_k[2]; // to avoid useless copying, use this to say which is the 
        // k variable and which is the k-1

        Vector M_lk(const int k=0) const;
        Vector M_Rk(const int k=0, const int i=0) const;
    
        void   switchLamRes();
        void   computeDifferences();
    };


// Constructors

    template<class Vector, class Real>
    generalizedAitken< Vector,  Real>::
    generalizedAitken(const int dof, const Real& defOmega, const int size):
        M_dof         (dof),
        M_lambdaS     (dof*size),
        M_lambdaF     (dof*size),
        M_nuS         (dof*size),
        M_nuF         (dof*size),
        M_sizeOmega   (size),
        M_omega       (M_sizeOmega),
        M_defaultOmega(defOmega)
    {
        M_k[0]=-1;
        M_k[1]=-1;
    }
    
    //
    // Destructor
    //
    
    //
    // Member functions
    //


    // setter and getter for lambdaS and lamdbaF
    
    template<class Vector, class Real>
    inline void generalizedAitken<Vector, Real>::setLambdaS(const Vector& _lamdba, const int _pos)
    {
        int jt = _pos*M_dof;
        for (int it = 0; it < _lamdba.size(); ++it, ++jt)
            M_lambdaS[jt] = _lambda[it];
    }
    
    template<class Vector, class Real>
    inline Real generalizedAitken<Vector, Real>::getLambdaS(const int _i, const int _j)
    {
        return M_lambdaS[_j*M_dof + _i];
    }
        
    template<class Vector, class Real>
    inline void generalizedAitken<Vector, Real>::setLambdaF(const Vector& _lamdba, const int _pos)
    {
        int jt = _pos*M_dof;
        for (int it = 0; it < _lamdba.size(); ++it, ++jt)
            M_lambdaF[jt] = _lambda[it];
    }
    
    template<class Vector, class Real>
    inline Real generalizedAitken<Vector, Real>::getLambdaF(const int _i, const int _j)
    {
        return M_lambdaF[_j*M_dof + _i];
    }
    
    // setter and getter for nuS and nuF

    template<class Vector, class Real>
    inline void generalizedAitken<Vector, Real>::setNuS(const Vector& _nu, const int _pos)
    {
        int jt = _pos*M_dof;
        for (int it = 0; it < _lamdba.size(); ++it, ++jt)
            M_nuS[jt] = _nu[it];
    }
    
    template<class Vector, class Real>
    inline Real generalizedAitken<Vector, Real>::getNuS(const int _i, const int _j)
    {
        return M_nuS[_j*dof + _i];
    }
        
    template<class Vector, class Real>
    inline void generalizedAitken<Vector, Real>::setNuF(const Vector& _nu, const int _pos)
    {
        int jt = _pos*M_dof;
        for (int it = 0; it < _lamdba.size(); ++it, ++jt)
            M_nuF[jt] = _nu[it];
    }
    
    template<class Vector, class Real>
    inline Real generalizedAitken<Vector, Real>::getNuF(const int _i, const int _j)
    {
        return M_nuF[_j*dof + _i];
    }
    
    //

    
    
//     template<class Vector, class Real>
//     Vector generalizedAitken<Vector, Real>::lk(const int k) const
//     {
//         if (k<2) {
//             return(lamRes[M_k[k]]);
//         } else {
//             if (k==2) {
//                 return(lamRes[k]);
//             } else {	
//                 return(-1);
//             }
//         }
//     }

//     template<class Vector, class Real>
//     Vector generalizedAitken< Vector,  Real>::Rk(const int k, const int i) const
//     {
//         if (k<2 && i < sizeOmega) {
//             return(lamRes[3*(i+1)+M_k[k]]);
//         } else {
//             if (k==2) {
//                 return(lamRes[3*(i+1)+k]);
//             } else {	
//                 return(-1);
//             }
//         }
//     }

//     template<class Vector, class Real>
//     void generalizedAitken< Vector,  Real>::switchLamRes()
//     {
//         M_k[1] = M_k[0];
//         M_k[0] = 1-M_k[1];
//     }

//     template<class Vector, class Real>
//     void generalizedAitken< Vector,  Real>::computeDifferences()
//     {
//         for (int i(0); i<sizeOmega+1; i++)
//         {
//             lamRes[2+3*i] = lamRes[M_k[0]+3*i] -  lamRes[M_k[1]+3*i] ;
//         }
//     }

//     template<class Vector, class Real>
//     generalizedAitken< Vector, Real>::~generalizedAitken()
//     {
//     }

//     template<class Vector, class Real>
//     void generalizedAitken<Vector, Real>::restart()
//     {
//         M_k[0]-1;
//         M_k[1]-1;
//     }


//     template<class Vector, class Real>
//     Vector generalizedAitken<Vector,  Real>::deltaLambda(const Vector& lam, const Vector residue[], Real omega[])
//     {
//         Vector res(M_dof);
//         res = 0;
//         int i;

//         if (defaultOmega<=0)
//         {
//             for (i=0; i<sizeOmega; i++)
//                 res+= omega[i] * residue[i];
//             return(res);
//         }
  

//         if (Mk[0] == -1){ // back to default omega
//             Mk[0] = 0;
    
//             lamRes[K_k[0]] = lam;
//             for (i=0; i < sizeOmega; ++i)
//             {
//                 lamRes[3*(i+1)+M_k[0]] = residue[i];
//                 res+= omega[i] * Rk(0,i);
//             }
//             return(res);
//         }

//         switchLamRes(); // switch

//         lamRes[M_k[0]] = lam;
//         for (i=0; i<sizeOmega; i++){
//             lamRes[3*(i+1)+M_k[0]] = residue[i];
//         }
  
//         computeDifferences();

//         if (sizeOmega == 1) {
//             omega[0] = - (Rk(2,0)*lk(2)) / (Rk(2,0)*Rk(2,0)); // checks whether functions exist
//         } else {
//             cout << "not yet implemented" << endl;
//             //  *omega=-pinv(Rk-Rk1)*(lk-lk1);
//         }

//         for (i=0; i<sizeOmega; i++){
//             res+= omega[i] * Rk(0,i);
//         }
//         return(res);
//     }


//     template<class Vector, class Real>
//     Vector generalizedAitken< Vector,  Real>::deltaLambda(const Vector& lam, const Vector residue[])
//     {
//         Real* omega;
//         omega = new Real[sizeOmega];
//         int i;
//         for(i=0; i<sizeOmega; i++){
//             omega[i] = abs(defaultOmega);
//         }
//         return(deltaLambda(lam,residue,omega) );

//     }

}
#endif
