/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file generalizedAitken.hpp

  \version 1.0
  \date 23/09/2004
  \author Simone Deparis <simone.deparis@epfl.ch>
  \author Gilles Fourestey <gilles.fourestey@epfl.ch>

  \brief Compute the acceleration with the vector variant of Aitken.

  \version 1.48
  \date 15/10/2009
  \author Cristiano Malossi <cristiano.malossi@epfl.ch>

  - Added three new differen Aitken methods for Scalar (inverted omega),
    Vector, and Block relaxations;
  - Added some Get and Set Methods;
  - Added Doxygen to the class.

  TODO Modify computeDeltaLambdaFSI - we should have only one defaultOmega parameter
       to be more general, and the same for M_oldResidualS & M_oldResidualF.
*/

#ifndef _GENERALIZEDAITKEN_HPP
#define _GENERALIZEDAITKEN_HPP

#include <life/lifecore/life.hpp>

#include <cstdlib>

#include <boost/shared_ptr.hpp>

namespace LifeV {

//! generalizedAitken - LifeV class for Aitken method
/*
 *  @author Simone Deparis, Gilles Fourestey, Cristiano Malossi
 *  @see S. Deparis, M. Discacciati and A. Quarteroni, A domain decompostion framework for fluid/structure interaction problems.
 *
 *  Compute the acceleration with the vector variant of Aitken.
 */
template< typename VectorType, typename DataType = LifeV::Real >
class generalizedAitken
{
    typedef boost::shared_ptr<VectorType>      VectorType_ptr;

public:

    /** @name Constructors, destructor
     */
    //@{

    //! Constructor
    generalizedAitken();

    //! Constructor
	/*!
	 * \param defaultOmegaS - First  Omega (Solid for FSI problems)
	 * \param defaultOmegaF - Second Omega (Fluid for FSI problems)
	 */
    generalizedAitken( const DataType defaultOmegaS,
                       const DataType defaultOmegaF = 0.1 );

    ~generalizedAitken() {}

    //@}



    /** @name Set Methods
     */
    //@{

    //! Set starting values for Omega
    /*!
     * \param defaultOmegaS - First  Omega (Solid for FSI problems)
     * \param defaultOmegaF - Second Omega (Fluid for FSI problems)
     */
    void setDefault( const DataType& defaultOmegaS = 0.1,
                     const DataType& defaultOmegaF = 0.1 );

    //@}



    /** @name Get Methods
     */
    //@{

    //! Get the default value of OmegaS
    const DataType&         GetDefaultOmegaS()  const { return M_defaultOmegaS; }

    //! Get the default value of OmegaF
    const DataType&         GetDefaultOmegaF()  const { return M_defaultOmegaF; }

    //@}



    /** @name Methods
     */
    //@{

    //! Use default omega
    void UseDefaultOmega( const bool useDefaultOmega = true )   { M_useDefaultOmega = useDefaultOmega; }

    //! Restart using the default omega
    void restart( const bool restart = true )                   { M_restart = restart; }

    //! Compute OmegaS * muS + OmegaF * muF
	/*!
	 * \param _lambda - vector of unknown
	 * \param _muF - vector of residuals (Fluid for FSI problems)
	 * \param _muS - vector of residuals (Solid for FSI problems)
	 */
    VectorType computeDeltaLambdaFSI( const VectorType& _lambda,
                                      const VectorType& _muF,
                                      const VectorType& _muS);

    //! Compute Omega * residual - Paragraph 4.2.3 & 4.2.4 of S. Deparis PhD Thesis
	/*!
	 * \param solution - vector of unknown
	 * \param residual - vector of residuals
	 * \param invertedOmega - false (default): minimizing on omega; true: minimizing on omega^-1
	 */
    VectorType computeDeltaLambdaScalar( const VectorType& solution,
                                         const VectorType& residual,
                                         const bool        invertedOmega = false );

    //! Compute Omega * residual - Paragraph 4.2.6 of S. Deparis PhD Thesis
    /*!
     * \param solution - vector of unknown
     * \param residual - vector of residuals
     */
    VectorType computeDeltaLambdaVector( const VectorType& solution,
                                         const VectorType& residual,
                                         const bool        invertedOmega = false,
                                         const bool        independentOmega = false );

    //! Compute Omega * residual - Paragraph 4.2.6 of S. Deparis PhD Thesis
    /*!
     * \param solution - vector of unknown
     * \param residual - vector of residuals
     * \param blocksVector - vector identifying the different blocks ( ID start from 1 )
     * \param blocksNumber - number of different blocks == higher ID (if = 1, it equal to the scalar case)
     * \param invertedOmega - false (default): minimizing on omega; true: minimizing on omega^-1
     */
    VectorType computeDeltaLambdaVectorBlock( const VectorType& solution,
                                              const VectorType& residual,
                                              const VectorType& blocksVector,
                                              const UInt        blocksNumber  = 1,
                                              const bool        invertedOmega = false);

    //@}

private:

    // fluid/structure interface dof count
    VectorType_ptr M_oldSolution;   // \lambda^{k - 1}
    VectorType_ptr M_oldResidualS;  // \mu_s^{k - 1}
    VectorType_ptr M_oldResidualF;  // \mu_f^{k - 1}

    // defaults \omega_s and \omega_f
    DataType M_defaultOmegaS;
    DataType M_defaultOmegaF;

    // first time call boolean
    bool M_restart;

    // If default omega is negative, then always use the
    // absolute value of the default omega.
    // In this case M_usedefault=true
    bool M_useDefaultOmega;
};



// ===================================================
//! Constructors
// ===================================================
template <class VectorType, class DataType>
generalizedAitken<VectorType, DataType>::generalizedAitken() :
    M_oldSolution    ( ),
    M_oldResidualS   ( ),
    M_oldResidualF   ( ),
    M_defaultOmegaS  ( 0.1 ),
    M_defaultOmegaF  ( 0.1 ),
    M_restart        ( true ),
    M_useDefaultOmega( false )
{}

template <class VectorType, class DataType>
generalizedAitken<VectorType, DataType>::generalizedAitken( const DataType defaultOmegaS,
                                                            const DataType defaultOmegaF ) :
    M_oldSolution    ( ),
    M_oldResidualS   ( ),
    M_oldResidualF   ( ),
    M_defaultOmegaS  ( ),
    M_defaultOmegaF  ( ),
    M_restart        ( true ),
    M_useDefaultOmega( )
{
    setDefault( defaultOmegaF, defaultOmegaS );
}



// ===================================================
//! Methods
// ===================================================
template <class VectorType, class DataType>
void generalizedAitken<VectorType, DataType>::setDefault( const DataType& defaultOmegaS,
                                                          const DataType& defaultOmegaF )
{
    if (( defaultOmegaS < 0 ) || ( defaultOmegaF< 0 )) //If omega is < use always default value
    {
        M_useDefaultOmega = true;
        M_defaultOmegaS   = std::fabs( defaultOmegaS );
        M_defaultOmegaF   = std::fabs( defaultOmegaF );
    }
    else                                               //Else use Aitken method for compute omega
    {
        M_useDefaultOmega = false;
        M_defaultOmegaS   = defaultOmegaS;
        M_defaultOmegaF   = defaultOmegaF;
    }
}

template <class VectorType, class DataType>
VectorType
generalizedAitken<VectorType, DataType>::computeDeltaLambdaFSI( const VectorType &_lambda,
                                                                const VectorType &_muF,
                                                                const VectorType &_muS )
{
    VectorType deltaLambda(_lambda.Map());

    if (( !M_restart ) && ( !M_useDefaultOmega ))
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

        muS = (_muS - M_oldResidualS);
        muF = (_muF - M_oldResidualF);

        a11 = muF*muF;
        a21 = muF*muS;
        a22 = muS*muS;

        b1 = muF * ( _lambda - *M_oldSolution);
        b2 = muS * ( _lambda - *M_oldSolution);

        /*
        for ( UInt ii = 0; ii < M_nDof; ++ii )
        {
            muS = _muS[ ii ] - *M_oldResidualS[ ii ];
            muF = _muF[ ii ] - *M_oldResidualF[ ii ];

            a11 += muF * muF;
            a21 += muF * muS;
            a22 += muS * muS;

            b1 += muF * ( _lambda[ ii ] - *M_oldSolution[ ii ] );
            b2 += muS * ( _lambda[ ii ] - *M_oldSolution[ ii ] );
        }
        */

        DataType omegaS ( M_defaultOmegaS );
        DataType omegaF ( M_defaultOmegaF );

        DataType det ( a22 * a11 - a21 * a21 );

        if ( std::fabs(det) != 0. )  //! eq. (12) page 8
        {
            omegaF = - ( a22 * b1 - a21 * b2 ) / det;
            omegaS = - ( a11 * b2 - a21 * b1 ) / det; // !

            if (omegaS == 0.) omegaS = M_defaultOmegaS;
            if (omegaF == 0.) omegaF = M_defaultOmegaF;
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

        *M_oldSolution = _lambda;
        *M_oldResidualF = _muF;
        *M_oldResidualS = _muS;
    }
    else
    {
        /*! first time aitken is called, the coefficients must be
          set to their default values
        */

        std::cout << " --------------- generalizedAitken: first call" << std::endl;
        M_restart = false;

        deltaLambda = M_defaultOmegaF * _muF + M_defaultOmegaS * _muS ;

        std::cout << "generalizedAitken: omegaS = " << M_defaultOmegaS
                  << " omegaF = " << M_defaultOmegaF << std::endl;

        M_oldSolution.reset( new VectorType(_lambda) );
        M_oldResidualF.reset( new VectorType(_muF) );
        M_oldResidualS.reset( new VectorType(_muS) );
    }

    return deltaLambda;
}

/*! one parameter version of the generalized aitken method. cf page 85 S. Deparis, PhD thesis */
template <class VectorType, class DataType>
VectorType
generalizedAitken<VectorType, DataType>::computeDeltaLambdaScalar( const VectorType& solution,
                                                                   const VectorType& residual,
                                                                   const bool        invertedOmega )
{
	if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset ( new VectorType( solution ) );
        M_oldResidualS.reset( new VectorType( residual ) );

        std::cout << "generalizedAitken: omega = " << M_defaultOmegaS << std::endl;

        return M_defaultOmegaS * residual;
    }

	VectorType deltaX( solution );
	VectorType deltaR( residual );

	deltaX -= *M_oldSolution;
	deltaR -= *M_oldResidualS;

	*M_oldSolution  = solution;
	*M_oldResidualS = residual;

	DataType omega, norm;
	if ( invertedOmega )
	{
	    // Minimization of the inverse omega
	    omega = deltaX.Dot( deltaX );
	    norm  = deltaR.Dot( deltaX );
	}
	else
	{
        //Minimization of the original omega
        omega = deltaX.Dot( deltaR );
        norm  = deltaR.Dot( deltaR );
	}

	omega = - omega / norm ;

	if (    std::fabs( omega ) < std::fabs( M_defaultOmegaS )/1024
		 || std::fabs( omega ) > std::fabs( M_defaultOmegaS )*1024 )
	{
		std::cout << "generalizedAitken: Failure: omega too small/big: "
				  << omega << std::endl;
		omega = M_defaultOmegaS;
	}

	std::cout << "generalizedAitken: omega = " << omega << std::endl;

	return omega * residual;
}

template <class VectorType, class DataType>
VectorType
generalizedAitken<VectorType, DataType>::computeDeltaLambdaVector( const VectorType& solution,
                                                                   const VectorType& residual,
                                                                   const bool        invertedOmega,
                                                                   const bool        independentOmega )
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset ( new VectorType( solution ) );
        M_oldResidualS.reset( new VectorType( residual ) );

        std::cout << "generalizedAitken: omega = " << M_defaultOmegaS << std::endl;

        return M_defaultOmegaS * residual;
    }

    VectorType deltaX( solution );
    VectorType deltaR( residual );

    deltaX -= *M_oldSolution;
    deltaR -= *M_oldResidualS;

    *M_oldSolution  = solution;
    *M_oldResidualS = residual;

    VectorType omega( deltaX );
    DataType   norm = 1;
    if ( independentOmega )
    {
        deltaR += !deltaR; // add +1 where deltaR is equal to zero!
        omega /= deltaR;
    }
    else if ( invertedOmega )
    {
        omega *= deltaX;
        norm = deltaR.Dot( deltaX );
    }
    else
    {
        omega *= deltaR;
        norm = deltaR.Dot( deltaR );
    }

    omega /= -norm;

    //Check omega limits
    VectorType omegaAbs( omega ); omegaAbs.Abs();
    VectorType correction = ( omegaAbs < std::fabs( M_defaultOmegaS )/1024   ||
                              omegaAbs > std::fabs( M_defaultOmegaS )*1024 ) || ( !omegaAbs && residual );

    omega *= !correction;
    omega +=  M_defaultOmegaS * correction;

    std::cout << "generalizedAitken: omega = " << std::endl; omega.ShowMe();

    omega *= residual;

    return omega;
}

template <class VectorType, class DataType>
VectorType
generalizedAitken<VectorType, DataType>::computeDeltaLambdaVectorBlock( const VectorType& solution,
                                                                        const VectorType& residual,
                                                                        const VectorType& blocksVector,
                                                                        const UInt        blocksNumber,
                                                                        const bool        invertedOmega)
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset(  new VectorType( solution ) );
        M_oldResidualS.reset( new VectorType( residual ) );

        std::cout << "generalizedAitken: omega = " << M_defaultOmegaS << std::endl;

        return M_defaultOmegaS * residual;
    }

    VectorType deltaX( solution );
    VectorType deltaR( residual );

    deltaX -= *M_oldSolution;
    deltaR -= *M_oldResidualS;

    *M_oldSolution  = solution;
    *M_oldResidualS = residual;

    VectorType omega( deltaX ); omega = 0.0;
    DataType   tempOmega = 0;

    VectorType tempVector( blocksVector );
    VectorType tempDeltaX ( deltaX );
    VectorType tempDeltaR ( deltaR );

    for ( UInt i = 0 ; i < blocksNumber ; ++i )
    {
        tempVector = blocksVector - i; tempVector = !tempVector;

        tempDeltaX  = deltaX;
        tempDeltaX *= tempVector;

        tempDeltaR  = deltaR;
        tempDeltaR *= tempVector;

        if ( invertedOmega )
        {
            // Minimization of the inverse omega
            tempOmega = -( tempDeltaX.Dot( tempDeltaX ) ) / ( tempDeltaR.Dot( tempDeltaX ) );
        }
        else
        {
            //Minimization of the original omega
            tempOmega = -( tempDeltaX.Dot( tempDeltaR ) ) / ( tempDeltaR.Dot( tempDeltaR ) );
        }

        if (    std::fabs( tempOmega ) < std::fabs( M_defaultOmegaS )/1024
             || std::fabs( tempOmega ) > std::fabs( M_defaultOmegaS )*1024 )
        {
            std::cout << "generalizedAitken: Failure: omega too small/big: "
                      << tempOmega << std::endl;
            tempOmega = M_defaultOmegaS;
        }

        omega += tempOmega * tempVector;
    }

    std::cout << "generalizedAitken: omega = " << std::endl; omega.ShowMe();

    return omega * residual;
}

} // end namespace LifeV

#endif
