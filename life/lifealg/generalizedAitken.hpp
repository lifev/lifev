
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
        /*!
          Compute the acceleration with the vector variant of Aitken.
          See "A domain decompostion framework for fluid/structure interaction problems"
          by Simone Deparis, Marco Discacciati and Alfio Quarteroni for reference
        */
    public:

        // Constructors
        
        generalizedAitken(const int dof,
                          const Real defOmegaS = 0.1,
                          const Real defOmegaF = 0.1);


        // Destructor

        // ~generalizedAitken();

        // Member functions

        void   restart();

        inline void   setMuS      (const Vector&, const int);
        inline Real   getMuS      (const int    , const int);
        inline void   setMuF      (const Vector&, const int);
        inline Real   getMuF      (const int    , const int);

        inline Real   getNewOmegaS(){return M_omegaS;};
        inline Real   getNewOmegaF(){return M_omegaF;};
        
        Vector        computeLambda();
        
        //Vector deltaLambda(const Vector& lam, const Vector res[]);
        //Vector deltaLambda(const Vector& lam, const Vector res[], Real omega[]);

        // in this case, omega is taken as the default value

    private:

        //! M_lamdba0 = \lambda_{k} M_lambda1 = \lambda_{k - 1}
        Vector M_lambda0;
        Vector M_lambda1;

        //! array formed by the vectors [\mu_s^k, \mu_s^{k-1}]
        Vector M_muS;

        //! array formed by the vectors [\mu_f^k, \mu_f^{k-1}]
        Vector M_muF;

        //! new step coefficients
        Real   M_omegaS;
        Real   M_omegaF;

        //! fluid/structure interface dof count

        int    M_nDof;
    };


// Constructors

    template<class Vector, class Real>
    generalizedAitken<Vector,  Real>::
    generalizedAitken(const int dof, const Real defOmegaS, const Real defOmegaF):
        M_nDof         (  dof),
        M_lambda0      (  dof),
        M_lambda1      (  dof),
        M_muS          (2*dof),
        M_muF          (2*dof),
        M_omegaS       (defOmegaS),
        M_omegaF       (defOmegaS)
    {
    }
    
    //
    // Destructor
    //
    /*template<class Vector, class Real>
    generalizedAitken<Vector,  Real>::
    ~generalizedAitken()
    {
    }*/
    //
    // Member functions
    //


    // setter and getter for muS and muF

    template<class Vector, class Real>
    inline void generalizedAitken<Vector, Real>::setMuS(const Vector& _mu, const int _pos)
    {
        int jt = _pos*M_nDof;
        for (int it = 0; it < _mu.size(); ++it, ++jt)
            M_muS[jt] = _mu[it];
    }
    
    template<class Vector, class Real>
    inline Real generalizedAitken<Vector, Real>::getMuS(const int _i, const int _j)
    {
        return M_muS[_j*dof + _i];
    }
        
    template<class Vector, class Real>
    inline void generalizedAitken<Vector, Real>::setMuF(const Vector& _mu, const int _pos)
    {
        int jt = _pos*M_nDof;
        for (int it = 0; it < _mu.size(); ++it, ++jt)
            M_muF[jt] = _mu[it];
    }
    
    template<class Vector, class Real>
    inline Real generalizedAitken<Vector, Real>::getMuF(const int _i, const int _j)
    {
        return M_muF[_j*dof + _i];
    }
    
    //


    //! this functions computes omega and store the results in M_omegaS and M_omegaF
    template<class Vector, class Real>
    Vector generalizedAitken<Vector, Real>::computeLambda()
    {
        Vector lambda(M_nDof);
        
        lambda = M_lambda0 - M_lambda1;

        Real   a11 = 0.;
        Real   a12 = 0.;
        Real   a22 = 0.;
        Real   b1  = 0.;
        Real   b2  = 0.;

        Vector muS   (M_nDof);
        Vector muF   (M_nDof);

        for (int ii = 0; ii < M_nDof; ++ii)
        {
            muS[ii] = getMuS(ii, 0) - getMuS(ii, 1);
            muF[ii] = getMuF(ii, 0) - getMuF(ii, 1);

            a11    += muF[ii]*muF[ii];
            a21    += muF[ii]*muS[ii];
            a22    += muS[ii]*muS[ii];

            b1     += muF[ii]*lambda[ii];
            b2     += muS[ii]*lambda[ii];
        }
            
        Real omega1, omega2;

        Real det = a21*a21 - a22*a11;

        if (det != 0.)
        {
            omega1 = (a22*b1 - a21*b2)/det;
            omega2 = (a11*b2 - a21*b1)/det;
        }

        for (int ii = 0; ii < M_nDof; ++ii)
            lambda = lambda + omega1*getMuF[ii] + omega2*getMuS[ii];
                
        lambda1 = lambda0;
        lambda0 = lambda;

        return lambda;
        
            
    }


}

#endif
