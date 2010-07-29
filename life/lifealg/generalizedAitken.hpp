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
  @file

  @version 1.0
  @date 23/09/2004
  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Gilles Fourestey <gilles.fourestey@epfl.ch>

  @brief Compute the acceleration with the vector variant of Aitken.

  @version 1.48
  @date 15/10/2009
  @author Cristiano Malossi <cristiano.malossi@epfl.ch>

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

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

namespace LifeV {

//! generalizedAitken - LifeV class for Aitken method
/*
 *  @author Simone Deparis, Gilles Fourestey, Cristiano Malossi
 *  @see S. Deparis, M. Discacciati and A. Quarteroni, A domain decompostion framework for fluid/structure interaction problems.
 *
 *  Compute the acceleration with the vector variant of Aitken.
 */
template< typename VectorType >
class generalizedAitken
{
    typedef boost::shared_ptr<VectorType>      Vector_PtrType;

public:

    /** @name Constructors, destructor
     */
    //@{

    //! Constructor
    generalizedAitken();

    ~generalizedAitken() {}

    //@}


    /** @name Methods
     */
    //@{

    //! Use default omega
    /*!
     * @param useDefaultOmega true: use always default Omega (relaxed fixed-point method)
     */
    void UseDefaultOmega( const bool useDefaultOmega = true );

    //! Reinitialize Omega to the default value
    void restart();

    //! Compute OmegaS * muS + OmegaF * muF
	/*!
	 * @param _lambda - vector of unknown
	 * @param _muF - vector of residuals (Fluid for FSI problems)
	 * @param _muS - vector of residuals (Solid for FSI problems)
	 */
    VectorType computeDeltaLambdaFSI( const VectorType& _lambda,
                                      const VectorType& _muF,
                                      const VectorType& _muS);

    //! Compute Omega * residual - Paragraph 4.2.3 & 4.2.4 of S. Deparis PhD Thesis
	/*!
	 * @param solution - vector of unknown
	 * @param residual - vector of residuals
	 * @param invertedOmega - false (default): minimizing on omega; true: minimizing on omega^-1
	 */
    VectorType computeDeltaLambdaScalar( const VectorType& solution,
                                         const VectorType& residual );

    //! Compute Omega * residual - Paragraph 4.2.6 of S. Deparis PhD Thesis
    /*!
     * @param solution - vector of unknown
     * @param residual - vector of residuals
     */
    VectorType computeDeltaLambdaVector( const VectorType& solution,
                                         const VectorType& residual,
                                         const bool        independentOmega = false );

    //! Compute Omega * residual - Paragraph 4.2.6 of S. Deparis PhD Thesis
    /*!
     * @param solution - vector of unknown
     * @param residual - vector of residuals
     * @param blocksVector - vector identifying the different blocks ( ID start from 1 )
     * @param blocksNumber - number of different blocks == higher ID (if = 1, it equal to the scalar case)
     */
    VectorType computeDeltaLambdaVectorBlock( const VectorType& solution,
                                              const VectorType& residual,
                                              const VectorType& blocksVector,
                                              const UInt        blocksNumber  = 1);

    //@}


    /** @name Set Methods
     */
    //@{

    //! Set starting values for Omega
    /*!
     * @param defaultOmegaS - First  Omega (Solid for FSI problems)
     * @param defaultOmegaF - Second Omega (Fluid for FSI problems)
     */
    void setDefaultOmega( const Real& defaultOmegaS = 0.1,
                          const Real& defaultOmegaF = 0.1 );

    //! Set the range of Omega.
    /*!
     * The range of Omega is defined as: OmegaMin < Omega < OmegaMax.
     *
     * @param OmegaRange array with the minimum and the maximum of Omega
     */
    void setOmegaRange( const boost::array< Real, 2 >& OmegaRange );

    //! Set the minimum of Omega.
    /*!
     * @param OmegaMin minimum of Omega
     */
    void setOmegaMin( const Real& OmegaMin );

    //! Set the maximum of Omega.
    /*!
     * @param OmegaMax maximum of Omega
     */
    void setOmegaMax( const Real& OmegaMax );

    //! Set minimization type.
    /*!
     * @param inverseOmega false: minimizing on omega; true: minimizing on omega^-1
     */
    void setMinimizationType( const bool& inverseOmega );

    //@}


    /** @name Get Methods
     */
    //@{

    //! Get the default value of OmegaS
    const Real& GetDefaultOmegaS()  const;

    //! Get the default value of OmegaF
    const Real& GetDefaultOmegaF()  const;

    //@}

private:

    /** @name Private Methods
     */
    //@{

    void checkRange( Real& Omega );

    //@}

    // fluid/structure interface dof count
    Vector_PtrType M_oldSolution;   // \lambda^{k - 1}
    Vector_PtrType M_oldResidualS;  // \mu_s^{k - 1}
    Vector_PtrType M_oldResidualF;  // \mu_f^{k - 1}

    // defaults \omega_s and \omega_f
    Real M_defaultOmegaS;
    Real M_defaultOmegaF;

    // first time call boolean
    bool M_restart;

    // If default omega is negative, then always use the
    // absolute value of the default omega.
    // In this case M_usedefault=true
    bool M_useDefaultOmega;

    // The max & min values for Omega
    boost::array< Real, 2 > M_rangeOmega;

    // Minimize on omega or omega^-1
    bool M_inverseOmega;
};



// ===================================================
// Constructors
// ===================================================
template <class VectorType>
generalizedAitken<VectorType>::generalizedAitken() :
    M_oldSolution      ( ),
    M_oldResidualS     ( ),
    M_oldResidualF     ( ),
    M_defaultOmegaS    ( 0.1 ),
    M_defaultOmegaF    ( 0.1 ),
    M_restart          ( true ),
    M_useDefaultOmega  ( false ),
    M_rangeOmega       ( ),
    M_inverseOmega     ( false )
{
    // Initializing the array
    M_rangeOmega[0] = 0.;
    M_rangeOmega[1] = 0.;
}

// ===================================================
// Methods
// ===================================================
template <class VectorType>
void
generalizedAitken<VectorType>::UseDefaultOmega( const bool useDefaultOmega )
{
    M_useDefaultOmega = useDefaultOmega;
}

template <class VectorType>
void
generalizedAitken<VectorType>::restart()
{
    M_restart = true;
}

template <class VectorType>
VectorType
generalizedAitken<VectorType>::computeDeltaLambdaFSI( const VectorType& _lambda,
                                                      const VectorType& _muF,
                                                      const VectorType& _muS )
{
    VectorType deltaLambda(_lambda.Map());

    if (( !M_restart ) && ( !M_useDefaultOmega ))
    {
        Real a11 = 0.;
        Real a21 = 0.;
        Real a22 = 0.;
        Real b1  = 0.;
        Real b2  = 0.;

        Real muS( 0 );
        Real muF( 0 );

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

        Real omegaS ( M_defaultOmegaS );
        Real omegaF ( M_defaultOmegaF );

        Real det ( a22 * a11 - a21 * a21 );

        if ( std::fabs(det) != 0. )  //! eq. (12) page 8
        {
            omegaF = - ( a22 * b1 - a21 * b2 ) / det;
            omegaS = - ( a11 * b2 - a21 * b1 ) / det; // !

            if (omegaS == 0.) omegaS = M_defaultOmegaS;
            if (omegaF == 0.) omegaF = M_defaultOmegaF;
        }
        else if (  std::fabs(a22) == 0. )
        {
            Debug(7020) << "generalizedAitken:  a22 = " << std::fabs(a22) << "\n";
            omegaS = 0.;
            omegaF = -b1 / a11;
        }
        else if (  std::fabs(a11) == 0. )
        {
            Debug(7020) << "generalizedAitken:  a11 = " << std::fabs(a11) << "\n";
            omegaS = -b2 / a22;
            omegaF = 0.;
        }
        else
        {
            Debug(7020) << "generalizedAitken: Failure: Det=0!!" << fabs(det) << "\n";
        }

        Debug(7020) << " --------------- generalizedAitken: \n";
        Debug(7020) << " omegaS = " << omegaS << " omegaF = " << omegaF << "\n";

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

        Debug(7020) << " --------------- generalizedAitken: first call\n";
        M_restart = false;

        deltaLambda = M_defaultOmegaF * _muF + M_defaultOmegaS * _muS ;

        Debug(7020) << "generalizedAitken: omegaS = " << M_defaultOmegaS << " omegaF = " << M_defaultOmegaF << "\n";

        M_oldSolution.reset( new VectorType(_lambda) );
        M_oldResidualF.reset( new VectorType(_muF) );
        M_oldResidualS.reset( new VectorType(_muS) );
    }

    return deltaLambda;
}

/*! one parameter version of the generalized aitken method. cf page 85 S. Deparis, PhD thesis */
template <class VectorType>
VectorType
generalizedAitken<VectorType>::computeDeltaLambdaScalar( const VectorType& solution,
                                                         const VectorType& residual )
{
	if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset ( new VectorType( solution ) );
        M_oldResidualS.reset( new VectorType( residual ) );

        Debug(7020) << "generalizedAitken: omega = " << M_defaultOmegaS << "\n";

        return M_defaultOmegaS * residual;
    }

	VectorType deltaX( solution );
	VectorType deltaR( residual );

	deltaX -= *M_oldSolution;
	deltaR -= *M_oldResidualS;

	*M_oldSolution  = solution;
	*M_oldResidualS = residual;

	Real omega, norm;
	if ( M_inverseOmega )
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

	//Check omega limits
	checkRange(omega);

	Debug(7020) << "generalizedAitken: omega = " << omega << "\n";
    //std::cout << "Omega" << std::endl;

	return omega * residual;
}

template <class VectorType>
VectorType
generalizedAitken<VectorType>::computeDeltaLambdaVector( const VectorType& solution,
                                                         const VectorType& residual,
                                                         const bool        independentOmega )
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset ( new VectorType( solution ) );
        M_oldResidualS.reset( new VectorType( residual ) );

        Debug(7020) << "generalizedAitken: omega = " << M_defaultOmegaS << "\n";

        return M_defaultOmegaS * residual;
    }

    VectorType deltaX( solution );
    VectorType deltaR( residual );

    deltaX -= *M_oldSolution;
    deltaR -= *M_oldResidualS;

    *M_oldSolution  = solution;
    *M_oldResidualS = residual;

    VectorType omega( deltaX );
    Real   norm = 1;
    if ( independentOmega )
    {
        deltaR += !deltaR; // add +1 where deltaR is equal to zero!
        omega /= deltaR;
    }
    else if ( M_inverseOmega )
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
    for ( UInt i(0) ; i < omega.size() ; ++i )
        checkRange(omega[i]);

    //Debug(7020) << "generalizedAitken: omega = " << "\n";
    //omega.ShowMe();

    omega *= residual;

    return omega;
}

template <class VectorType>
VectorType
generalizedAitken<VectorType>::computeDeltaLambdaVectorBlock( const VectorType& solution,
                                                              const VectorType& residual,
                                                              const VectorType& blocksVector,
                                                              const UInt        blocksNumber )
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset(  new VectorType( solution ) );
        M_oldResidualS.reset( new VectorType( residual ) );

        Debug(7020) << "generalizedAitken: omega = " << M_defaultOmegaS << "\n";

        return M_defaultOmegaS * residual;
    }

    VectorType deltaX( solution );
    VectorType deltaR( residual );

    deltaX -= *M_oldSolution;
    deltaR -= *M_oldResidualS;

    *M_oldSolution  = solution;
    *M_oldResidualS = residual;

    VectorType omega( deltaX ); omega = 0.0;
    Real   tempOmega = 0;

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

        if ( M_inverseOmega )
        {
            // Minimization of the inverse omega
            tempOmega = -( tempDeltaX.Dot( tempDeltaX ) ) / ( tempDeltaR.Dot( tempDeltaX ) );
        }
        else
        {
            //Minimization of the original omega
            tempOmega = -( tempDeltaX.Dot( tempDeltaR ) ) / ( tempDeltaR.Dot( tempDeltaR ) );
        }

        //Check omega limits
        checkRange(tempOmega);

        omega += tempOmega * tempVector;
    }

    //Debug(7020) << "generalizedAitken: omega = " << "\n";
    //omega.ShowMe();

    return omega * residual;
}

// ===================================================
// Set Methods
// ===================================================
template <class VectorType>
void generalizedAitken<VectorType>::setDefaultOmega( const Real& defaultOmegaS,
                                                     const Real& defaultOmegaF )
{
    M_defaultOmegaS = defaultOmegaS;
    M_defaultOmegaF = defaultOmegaF;
}

template <class VectorType>
void generalizedAitken<VectorType>::setOmegaRange( const boost::array< Real, 2 >& OmegaRange )
{
    M_rangeOmega = OmegaRange;
}

template <class VectorType>
void generalizedAitken<VectorType>::setOmegaMin( const Real& OmegaMin )
{
    M_rangeOmega[0] = std::abs( OmegaMin );
}

template <class VectorType>
void generalizedAitken<VectorType>::setOmegaMax( const Real& OmegaMax )
{
    M_rangeOmega[1] = std::abs( OmegaMax );
}

template <class VectorType>
void generalizedAitken<VectorType>::setMinimizationType( const bool& inverseOmega )
{
    M_inverseOmega = inverseOmega;
}

// ===================================================
// Get Methods
// ===================================================
template <class VectorType>
const Real&
generalizedAitken<VectorType>::GetDefaultOmegaS()  const
{
    return M_defaultOmegaS;
}

template <class VectorType>
const Real&
generalizedAitken<VectorType>::GetDefaultOmegaF()  const
{
    return M_defaultOmegaF;
}

// ===================================================
// Private Methods
// ===================================================
template <class VectorType>
void
generalizedAitken<VectorType>::checkRange( Real& Omega )
{
    if ( std::abs(Omega) < std::abs(M_rangeOmega[0]) )
    {
        if ( Omega < 0 )
            Omega = -M_rangeOmega[0];
        else
            Omega = M_rangeOmega[0];
    }
    else if ( std::abs(Omega) > std::abs(M_rangeOmega[1]) )
    {
        if ( Omega < 0 )
            Omega = -M_rangeOmega[1];
        else
            Omega = M_rangeOmega[1];
    }
}

} // end namespace LifeV

#endif
