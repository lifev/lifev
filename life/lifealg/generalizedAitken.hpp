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

#include <cstdlib>

#include <boost/shared_ptr.hpp>

namespace LifeV
{
template <typename VectorType, typename DataType>
class generalizedAitken
{
    /*!
      Compute the acceleration with the vector variant of Aitken.
      See "A domain decompostion framework for fluid/structure interaction problems"
      by Simone Deparis, Marco Discacciati and Alfio Quarteroni for references
    */
    typedef VectorType                     vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

public:

    // Constructors

    generalizedAitken():
        M_lambda    ( ),
        M_muS       ( ),
        M_muF       ( ),
        M_defOmegaS ( 0.1 ),
        M_defOmegaF ( 0.1 ),
        M_firstCall ( true ),
        M_issetup   ( false )
        {}

    generalizedAitken( const DataType _defOmegaS,
                       const DataType _defOmegaF = 0.1 );


    // Destructor

    ~generalizedAitken();

    // Member functions

    void setDefault( const DataType _defOmegaS = 0.1,
                     const DataType _defOmegaF = 0.1 );
    void restart();

    vector_type computeDeltaLambda( const vector_type &,
                                    const vector_type &,
                                   const vector_type & );

    vector_type computeDeltaLambda( const vector_type &,
                                   const vector_type & );


private:

    //! fluid/structure interface dof count

    //! \lambda^{k - 1}
    vector_ptrtype M_lambda;

    //! \mu_s^{k - 1}
    vector_ptrtype M_muS;

    //! \mu_f^{k - 1}
    vector_ptrtype M_muF;

    //! defaults \omega_s and \omega_f
    DataType M_defOmegaS;
    DataType M_defOmegaF;

    //! first time call boolean
    bool M_firstCall;

    bool M_issetup;

    //! If default omega is negative, then always use the
    // absolute value of the default omega. In this case
    //  M_usedefault=true
    bool M_useDefault;
};

//
// Constructors
//

template <class VectorType, class DataType>
generalizedAitken<VectorType, DataType>::generalizedAitken( const DataType _defOmegaF,
                                                            const DataType _defOmegaS ) :
    M_lambda    ( ),
    M_muS       ( ),
    M_muF       ( ),
    M_firstCall ( true )
{
    setDefault(_defOmegaF, _defOmegaS);
}

//
// Destructor
//

template <class VectorType, class DataType>
generalizedAitken<VectorType, DataType>::~generalizedAitken()
{ //nothing needs to be done
}

//
// Member functions
//
template <class VectorType, class DataType>
void generalizedAitken<VectorType, DataType>::
setDefault( const DataType _defOmegaS,
            const DataType _defOmegaF )
{
    if (( _defOmegaS < 0 ) || ( _defOmegaF< 0 ))
    {
        M_useDefault = true;
        M_defOmegaS = std::fabs(_defOmegaS);
        M_defOmegaF = std::fabs(_defOmegaF);
    } else {
        M_useDefault = false;
        M_defOmegaS = _defOmegaS;
        M_defOmegaF = _defOmegaF;
    }

    M_issetup = true;

}

template <class VectorType, class DataType>
void generalizedAitken<VectorType, DataType>::
restart()
{
    M_firstCall = true;
}

/*! this functions computes omega and return the new lambda
  _lambda is \lamdba^{k}
  _muS    is \mu_s^{k}
  _muF    is \mu_f^{k}
*/

template <class VectorType, class DataType>
typename generalizedAitken<VectorType, DataType>::vector_type generalizedAitken<VectorType, DataType>::
computeDeltaLambda( const vector_type &_lambda,
                    const vector_type &_muF,
                    const vector_type &_muS )
{
    VectorType deltaLambda(_lambda.Map());

    if (( !M_firstCall ) && ( !M_useDefault ))
    {
        DataType a11 = 0.;
        DataType a21 = 0.;
        DataType a22 = 0.;
        DataType b1  = 0.;
        DataType b2  = 0.;

        DataType muS( 0 );
        DataType muF( 0 );

        /*! bulding the matrix and the right hand side
          see eq. (16) page 10
        */

        muS = (_muS - M_muS);
        muF = (_muF - M_muF);

        a11 = muF*muF;
        a21 = muF*muS;
        a22 = muS*muS;

        b1 = muF * ( _lambda - *M_lambda);
        b2 = muS * ( _lambda - *M_lambda);

        /*
        for ( UInt ii = 0; ii < M_nDof; ++ii )
        {
            muS = _muS[ ii ] - *M_muS[ ii ];
            muF = _muF[ ii ] - *M_muF[ ii ];

            a11 += muF * muF;
            a21 += muF * muS;
            a22 += muS * muS;

            b1 += muF * ( _lambda[ ii ] - *M_lambda[ ii ] );
            b2 += muS * ( _lambda[ ii ] - *M_lambda[ ii ] );
        }
        */

        DataType omegaS ( M_defOmegaS );
        DataType omegaF ( M_defOmegaF );

        DataType det ( a22 * a11 - a21 * a21 );

        if ( std::fabs(det) != 0. )  //! eq. (12) page 8
        {
            omegaF = - ( a22 * b1 - a21 * b2 ) / det;
            omegaS = - ( a11 * b2 - a21 * b1 ) / det; // !

            if (omegaS == 0.) omegaS = M_defOmegaS;
            if (omegaF == 0.) omegaF = M_defOmegaF;
        }
        else if (  std::fabs(a22) == 0. )
        {
            std::cout << "generalizedAitken:  a22 = "
                      << std::fabs(a22) << std::endl;
            omegaS = 0.;
            omegaF = -b1 / a11;
        }
        else if (  std::fabs(a11) == 0. )
        {
            std::cout << "generalizedAitken:  a11 = "
                      << std::fabs(a11) << std::endl;
            omegaS = -b2 / a22;
            omegaF = 0.;
        }
        else
        {
            std::cout << "generalizedAitken: Failure: Det=0!!"
                      << fabs(det) << std::endl;
        }

        std::cout << " --------------- generalizedAitken: " << std::endl;
        std::cout << " omegaS = " << omegaS
                  << " omegaF = " << omegaF << std::endl;

        deltaLambda = omegaF * _muF + omegaS * _muS;

        *M_lambda = _lambda;
        *M_muF = _muF;
        *M_muS = _muS;
    }
    else
    {
        /*! first time aitken is called, the coefficients must be
          set to their default values
        */

        std::cout << " --------------- generalizedAitken: first call" << std::endl;
        M_firstCall = false;

        deltaLambda = M_defOmegaF * _muF + M_defOmegaS * _muS ;

        std::cout << "generalizedAitken: omegaS = " << M_defOmegaS
                  << " omegaF = " << M_defOmegaF << std::endl;

        M_lambda.reset( new vector_type(_lambda) );
        M_muF.reset( new vector_type(_muF) );
        M_muS.reset( new vector_type(_muS) );
    }


    return deltaLambda;
}


/*! one parameter version of the generalized aitken method. cf page 85 S. Deparis, PhD thesis */
template <class VectorType, class DataType>
typename generalizedAitken<VectorType, DataType>::vector_type generalizedAitken<VectorType, DataType>::
computeDeltaLambda( const vector_type &_lambda,
                    const vector_type &_mu )
{
    VectorType deltaLambda(_lambda.getMap());

    if (( !M_firstCall ) && ( !M_useDefault ))
    {
        VectorType deltaMu = _mu;
        deltaMu -= *M_muS;

        DataType omega = 0.;
        DataType norm  = 0.;

        deltaLambda  = _lambda;
        deltaLambda -= *M_lambda;

        omega = deltaLambda * deltaMu;
        norm  = deltaMu * deltaMu;

        /*
        for ( UInt ii = 0; ii < deltaLambda.size(); ++ii )
        {
            omega += deltaLambda[ ii ] * deltaMu[ ii ];
            norm  += deltaMu[ ii ] * deltaMu[ ii ];
        }
        */

        *M_lambda = _lambda;
        *M_muS    = _mu;

        omega    = - omega / norm ;

        if ( std::fabs(omega) < std::fabs(M_defOmegaS)/1024
             || std::fabs(omega) > std::fabs(M_defOmegaS)*1024 )
        {
            std::cout << "generalizedAitken: Failure: omega too small/big: "
                      << omega << std::endl;
            omega = M_defOmegaS;
        }

        deltaLambda = omega * _mu;

        std::cout << "generalizedAitken: omega = " << omega << std::endl;

    }
    else
    {
        M_firstCall = false;

        deltaLambda = M_defOmegaS * _mu;

        std::cout << "generalizedAitken: omega = " << M_defOmegaS << std::endl;

        M_lambda.reset( new vector_type( _lambda) );
        M_muS.reset( new vector_type(_mu) );
    }

    return deltaLambda;

}



}

#endif
