
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

#ifndef _GENERALIZEDAITKEN_HPP
#define _GENERALIZEDAITKEN_HPP

#include "lifeV.hpp"

namespace LifeV
{
    template<class Vector, class Real>
    class generalizedAitken
    {
        /*!
          Compute the acceleration with the vector variant of Aitken.
          See "A domain decompostion framework for fluid/structure interaction problems"
          by Simone Deparis, Marco Discacciati and Alfio Quarteroni for references
        */
        
    public:
        
        // Constructors
        
        generalizedAitken(const int  _nDof,
                          const Real _defOmegaS = 0.1,
                          const Real _defOmegaF = 0.1);
        
        
        // Destructor
        
        ~generalizedAitken();
        
        // Member functions
        
        Vector        computeDeltaLambda(const Vector &,
                                         const Vector &,
                                         const Vector &);
        
    private:

        //! fluid/structure interface dof count
        UInt   M_nDof;

        //! \lambda^{k - 1}
        Vector M_lambda;

        //! \mu_s^{k - 1}
        Vector M_muS;

        //! \mu_f^{k - 1}
        Vector M_muF;

        //! defaults \omega_s and \omega_f
        Real   M_defOmegaS;
        Real   M_defOmegaF;

        //! first time call boolean

        bool   M_firstCall;
    };
    
    //
    // Constructors
    //
    
    template<class Vector, class Real>
    generalizedAitken<Vector, Real>::generalizedAitken(const int _nDof,
                                                       const Real _defOmegaS,
                                                       const Real _defOmegaF):
        M_nDof         (_nDof),
        M_lambda       (_nDof),
        M_muS          (_nDof),
        M_muF          (_nDof),
        M_defOmegaS    (_defOmegaS),
        M_defOmegaF    (_defOmegaF)
    {
        M_muS       = 0.;
        M_muF       = 0.;
        M_lambda    = 0.;
        M_firstCall = true;
    }
    
    //
    // Destructor
    //

    template<class Vector, class Real>
    generalizedAitken<Vector, Real>::~generalizedAitken()
    { //nothing needs to be done
    }
        
    //
    // Member functions
    //


    /*! this functions computes omega and return the new lambda
        _lambda is \lamdba^{k}
        _muS    is \mu_s^{k}
        _muF    is \mu_f^{k}
    */
    
    template<class Vector, class Real>
    Vector generalizedAitken<Vector, Real>::computeDeltaLambda(const Vector &_lambda, 
                                                 const Vector &_muS,
                                                 const Vector &_muF)
    {
        Vector    deltaLambda;

        if (!M_firstCall)
        {
            Real   a11 = 0.;
            Real   a21 = 0.;
            Real   a22 = 0.;
            Real   b1  = 0.;
            Real   b2  = 0.;
            
            Vector muS   (M_nDof);
            Vector muF   (M_nDof);
            
            /*! bulding the matrix and the right hand side
              see eq. (16) page 10
            */
            
            for (UInt ii = 0; ii < M_nDof; ++ii)
            {
                muS[ii] = _muS[ii] - M_muS[ii];
                muF[ii] = _muF[ii] - M_muF[ii];
                
                a11    += muF[ii]*muF[ii];
                a21    += muF[ii]*muS[ii];
                a22    += muS[ii]*muS[ii];
                
                b1     += muF[ii]*(_lambda[ii] - M_lambda[ii]);
                b2     += muS[ii]*(_lambda[ii] - M_lambda[ii]);
            }
            
            Real omegaS = M_defOmegaS;
            Real omegaF = M_defOmegaF;
            
            Real det    = a21*a21 - a22*a11;
            
            if (det != 0.)
            {
                omegaF = (a22*b1 - a21*b2)/det;
                omegaS = (a11*b2 - a21*b1)/det;
            }
            else if (a22 == 0)
            {
                //! eq. (12) page 8
                omegaS = 0.;
                omegaF = b1/a11*a11;
            }

            deltaLambda = omegaF*muF + omegaS*muS;
            
            M_lambda    = _lambda;
            M_muF       = _muF;
            M_muS       = _muS;
        }
        else
        {
            /*! first time aitken is called, the coefficients must be 
                set to their default values
            */

            M_firstCall = false;
                        
            deltaLambda = M_defOmegaF*_muF + M_defOmegaS*_muS;            

            M_lambda    = _lambda;
            M_muF       = _muF;
            M_muS       = _muS;
        }
        
        return deltaLambda;
    }

    
}

#endif
