//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the generalized Aitken algorithm
 *
 *  @date 23-09-2004
 *  @author Simone Deparis <simone.deparis@epfl.ch>
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef GENERALIZEDAITKEN_H
#define GENERALIZEDAITKEN_H

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <cstdlib>

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/life.hpp>

namespace LifeV
{

//! generalizedAitken - LifeV class for generalized Aitken algorithm
/*
 *  @author Simone Deparis, Gilles Fourestey, Cristiano Malossi
 *  @see S. Deparis, M. Discacciati and A. Quarteroni, A domain decompostion framework for fluid/structure interaction problems.
 *
 *  Compute the acceleration with the vector variant of Aitken.
 *
 *  TODO Modify computeDeltaLambdaFSI - we should have only one defaultOmega parameter
 *  to be more general, and the same for M_oldResidualSolid & M_oldResidualFluid.
 */
template< typename VectorType >
class generalizedAitken
{

public:

    //! @name Public Types
    //@{

    typedef VectorType                            vector_Type;
    typedef boost::shared_ptr< vector_Type >      vectorPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit generalizedAitken();

    //! Destructor
    virtual ~generalizedAitken() {}

    //@}


    //! @name Methods
    //@{

    //! Use default omega
    /*!
     * @param useDefaultOmega true: use always default Omega (relaxed fixed-point method)
     */
    void useDefaultOmega( const bool useDefaultOmega = true ) { M_useDefaultOmega = useDefaultOmega; }

    //! Reinitialize Omega to the default value
    void restart() { M_restart = true; }

    //! Compute OmegaS * deltaRSolid + OmegaF * deltaRFluid
    /*!
     * @param solution - vector of unknown
     * @param residualFluid - vector of residuals (Fluid for FSI problems)
     * @param residualSolid - vector of residuals (Solid for FSI problems)
     */
    vector_Type computeDeltaLambdaFSI( const vector_Type& solution,
                                       const vector_Type& residualFluid,
                                       const vector_Type& residualSolid );

    //! Compute Omega * residual - Paragraph 4.2.3 & 4.2.4 of S. Deparis PhD Thesis
    /*!
     * @param solution - vector of unknown
     * @param residual - vector of residuals
     * @param invertedOmega - false (default): minimizing on omega; true: minimizing on omega^-1
     */
    vector_Type computeDeltaLambdaScalar( const vector_Type& solution,
                                          const vector_Type& residual );

    //! Compute Omega * residual - Paragraph 4.2.6 of S. Deparis PhD Thesis
    /*!
     * @param solution - vector of unknown
     * @param residual - vector of residuals
     */
    vector_Type computeDeltaLambdaVector( const vector_Type& solution,
                                          const vector_Type& residual,
                                          const bool&        independentOmega = false );

    //! Compute Omega * residual - Paragraph 4.2.6 of S. Deparis PhD Thesis
    /*!
     * @param solution - vector of unknown
     * @param residual - vector of residuals
     * @param blocksVector - vector identifying the different blocks ( ID start from 1 )
     * @param blocksNumber - number of different blocks == higher ID (if = 1, it equal to the scalar case)
     */
    vector_Type computeDeltaLambdaVectorBlock( const vector_Type& solution,
                                               const vector_Type& residual,
                                               const vector_Type& blocksVector,
                                               const UInt&        blocksNumber = 1);

    //@}


    //! @name Set Methods
    //@{

    //! Set starting values for Omega
    /*!
     * @param defaultOmegaFluid default value for the omega fluid parameter
     * @param defaultOmegaSolid default value for the omega solid parameter
     */
    void setDefaultOmega( const Real& defaultOmegaFluid = 0.1, const Real& defaultOmegaSolid = 0.1 );

    //! Set the range of Omega.
    /*!
     * The range of Omega is defined as: OmegaMin < Omega < OmegaMax.
     *
     * @param omegaRange array with the minimum and the maximum of Omega
     */
    void setOmegaRange( const boost::array< Real, 2 >& omegaRange ) { M_rangeOmega = omegaRange; }

    //! Set the minimum of Omega.
    /*!
     * @param omegaMin minimum of Omega
     */
    void setOmegaMin( const Real& omegaMin ) { M_rangeOmega[0] = std::abs( omegaMin ); }

    //! Set the maximum of Omega.
    /*!
     * @param omegaMax maximum of Omega
     */
    void setOmegaMax( const Real& omegaMax ) { M_rangeOmega[1] = std::abs( omegaMax ); }

    //! Set minimization type.
    /*!
     * @param inverseOmega false: minimizing on omega; true: minimizing on omega^-1
     */
    void setMinimizationType( const bool& inverseOmega ) { M_inverseOmega = inverseOmega; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the default value of omega fluid.
    /*!
     * @return default value of omega fluid
     */
    const Real& defaultOmegaFluid()  const { return M_defaultOmegaFluid; }

    //! Get the default value of omega solid.
    /*!
     * @return default value of omega solid
     */
    const Real& defaultOmegaSolid()  const { return M_defaultOmegaSolid; }

    //@}

private:

    //! @name Private unimplemented Methods
    //@{

    generalizedAitken( const generalizedAitken& aitken );

    generalizedAitken& operator=( const generalizedAitken& aitken );

    //@}


    //! @name Private Methods
    //@{

    void checkRange( Real& omega );

    //@}

    // fluid/structure interface dof count
    vectorPtr_Type M_oldSolution;   // \lambda^{k - 1}
    vectorPtr_Type M_oldResidualFluid;  // \mu_f^{k - 1}
    vectorPtr_Type M_oldResidualSolid;  // \mu_s^{k - 1}

    // defaults \omega_s and \omega_f
    Real M_defaultOmegaFluid;
    Real M_defaultOmegaSolid;

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
template < class VectorType >
generalizedAitken< VectorType >::generalizedAitken() :
        M_oldSolution      ( ),
        M_oldResidualFluid ( ),
        M_oldResidualSolid ( ),
        M_defaultOmegaFluid( 0.1 ),
        M_defaultOmegaSolid( 0.1 ),
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
template < class VectorType >
typename generalizedAitken< VectorType >::vector_Type
generalizedAitken< VectorType >::computeDeltaLambdaFSI( const vector_Type& solution,
                                                        const vector_Type& residualFluid,
                                                        const vector_Type& residualSolid )
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset ( new vector_Type( solution ) );
        M_oldResidualSolid.reset( new vector_Type( residualSolid ) );
        M_oldResidualFluid.reset( new vector_Type( residualFluid) );

#ifdef HAVE_LIFEV_DEBUG
        Debug(7020) << "generalizedAitken: omegaFluid = " << M_defaultOmegaFluid << " omegaSolid = " << M_defaultOmegaSolid << "\n";
#endif
        return M_defaultOmegaFluid * residualFluid + M_defaultOmegaSolid * residualSolid;
    }
    else
    {
        Real deltaRSolid ( residualSolid - M_oldResidualSolid );
        Real deltaRFluid ( residualFluid - M_oldResidualFluid );

        Real a11 = deltaRFluid * deltaRFluid;
        Real a21 = deltaRFluid * deltaRSolid;
        Real a22 = deltaRSolid * deltaRSolid;

        Real b1  = deltaRFluid * ( solution - *M_oldSolution);
        Real b2  = deltaRSolid * ( solution - *M_oldSolution);

        Real det ( a22 * a11 - a21 * a21 );

        Real omegaFluid ( M_defaultOmegaFluid );
        Real omegaSolid ( M_defaultOmegaSolid );

        if ( std::fabs( det ) != 0. )  //! eq. (12) page 8
        {
            omegaFluid = - ( a22 * b1 - a21 * b2 ) / det;
            omegaSolid = - ( a11 * b2 - a21 * b1 ) / det;

            if ( omegaFluid == 0. )
                omegaFluid = M_defaultOmegaFluid;

            if ( omegaSolid == 0. )
                omegaSolid = M_defaultOmegaSolid;
        }
        else if ( std::fabs( a22 ) == 0. )
        {
#ifdef HAVE_LIFEV_DEBUG
            Debug(7020) << "generalizedAitken:  a22 = " << std::fabs(a22) << "\n";
#endif
            omegaFluid = -b1 / a11;
            omegaSolid = 0.;
        }
        else if ( std::fabs(a11) == 0. )
        {
#ifdef HAVE_LIFEV_DEBUG
            Debug(7020) << "generalizedAitken:  a11 = " << std::fabs(a11) << "\n";
#endif
            omegaFluid = 0.;
            omegaSolid = -b2 / a22;
        }
#ifdef HAVE_LIFEV_DEBUG
        else
            Debug(7020) << "generalizedAitken: Failure: Det=0!!" << fabs(det) << "\n";

        Debug(7020) << " --------------- generalizedAitken: \n";
        Debug(7020) << " omegaSolid = " << omegaSolid << " omegaFluid = " << omegaFluid << "\n";
#endif

        *M_oldSolution = solution;
        *M_oldResidualFluid = residualFluid;
        *M_oldResidualSolid = residualSolid;

        return omegaFluid * residualFluid + omegaSolid * residualSolid;
    }
}

/*! one parameter version of the generalized aitken method. cf page 85 S. Deparis, PhD thesis */
template < class VectorType >
typename generalizedAitken< VectorType >::vector_Type
generalizedAitken< VectorType >::computeDeltaLambdaScalar( const vector_Type& solution,
                                                           const vector_Type& residual )
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset ( new vector_Type( solution ) );
        M_oldResidualFluid.reset( new vector_Type( residual ) );

#ifdef HAVE_LIFEV_DEBUG
        Debug(7020) << "generalizedAitken: omega = " << M_defaultOmegaFluid << "\n";
#endif

        return M_defaultOmegaFluid * residual;
    }

    vector_Type deltaX( solution );
    vector_Type deltaR( residual );

    deltaX -= *M_oldSolution;
    deltaR -= *M_oldResidualFluid;

    *M_oldSolution  = solution;
    *M_oldResidualFluid = residual;

    Real omega, norm;
    if ( M_inverseOmega )
    {
        // Minimization of the inverse omega
        omega = deltaX.dot( deltaX );
        norm  = deltaR.dot( deltaX );
    }
    else
    {
        //Minimization of the original omega
        omega = deltaX.dot( deltaR );
        norm  = deltaR.dot( deltaR );
    }

    omega = - omega / norm ;

    //Check omega limits
    checkRange(omega);

#ifdef HAVE_LIFEV_DEBUG
    Debug(7020) << "generalizedAitken: omega = " << omega << "\n";
#endif

    return omega * residual;
}

template < class VectorType >
typename generalizedAitken< VectorType >::vector_Type
generalizedAitken< VectorType >::computeDeltaLambdaVector( const vector_Type& solution,
                                                           const vector_Type& residual,
                                                           const bool&        independentOmega )
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset ( new vector_Type( solution ) );
        M_oldResidualFluid.reset( new vector_Type( residual ) );

#ifdef HAVE_LIFEV_DEBUG
        Debug(7020) << "generalizedAitken: omega = " << M_defaultOmegaFluid << "\n";
#endif

        return M_defaultOmegaFluid * residual;
    }

    vector_Type deltaX( solution );
    vector_Type deltaR( residual );

    deltaX -= *M_oldSolution;
    deltaR -= *M_oldResidualFluid;

    *M_oldSolution  = solution;
    *M_oldResidualFluid = residual;

    vector_Type omega( deltaX );
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

#ifdef HAVE_LIFEV_DEBUG
    Debug(7020) << "generalizedAitken: omega = " << "\n";
    omega.ShowMe();
#endif

    omega *= residual;

    return omega;
}

template < class VectorType >
typename generalizedAitken< VectorType >::vector_Type
generalizedAitken< VectorType >::computeDeltaLambdaVectorBlock( const vector_Type& solution,
                                                                const vector_Type& residual,
                                                                const vector_Type& blocksVector,
                                                                const UInt&        blocksNumber )
{
    if ( M_restart || M_useDefaultOmega )
    {
        M_restart = false;

        M_oldSolution.reset(  new vector_Type( solution ) );
        M_oldResidualFluid.reset( new vector_Type( residual ) );

#ifdef HAVE_LIFEV_DEBUG
        Debug(7020) << "generalizedAitken: omega = " << M_defaultOmegaFluid << "\n";
#endif

        return M_defaultOmegaFluid * residual;
    }

    vector_Type deltaX( solution );
    vector_Type deltaR( residual );

    deltaX -= *M_oldSolution;
    deltaR -= *M_oldResidualFluid;

    *M_oldSolution  = solution;
    *M_oldResidualFluid = residual;

    vector_Type omega( deltaX );
    omega = 0.0;
    Real   tempOmega = 0;

    vector_Type tempVector( blocksVector );
    vector_Type tempDeltaX ( deltaX );
    vector_Type tempDeltaR ( deltaR );

    for ( UInt i = 0 ; i < blocksNumber ; ++i )
    {
        tempVector = blocksVector - i;
        tempVector = !tempVector;

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

#ifdef HAVE_LIFEV_DEBUG
    Debug(7020) << "generalizedAitken: omega = " << "\n";
    omega.ShowMe();
#endif

    return omega * residual;
}

// ===================================================
// Set Methods
// ===================================================
template < class VectorType >
inline void
generalizedAitken< VectorType >::setDefaultOmega( const Real& defaultOmegaFluid, const Real& defaultOmegaSolid )
{
    M_defaultOmegaFluid = defaultOmegaFluid;
    M_defaultOmegaSolid = defaultOmegaSolid;
}

// ===================================================
// Private Methods
// ===================================================
template < class VectorType >
inline void
generalizedAitken< VectorType >::checkRange( Real& omega )
{
    if ( std::abs(omega) < std::abs(M_rangeOmega[0]) )
    {
        if ( omega < 0 )
            omega = -M_rangeOmega[0];
        else
            omega = M_rangeOmega[0];
    }
    else if ( std::abs(omega) > std::abs(M_rangeOmega[1]) )
    {
        if ( omega < 0 )
            omega = -M_rangeOmega[1];
        else
            omega = M_rangeOmega[1];
    }
}

} // end namespace LifeV

#endif // GENERALIZEDAITKEN_H
