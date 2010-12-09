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
    @file
    @brief File containing a base class providing physical operations for the 1D model data.

    @version 1.0
    @date 01-07-2004
    @author Vincent Martin

    @version 2.0
    @date 13-04-2010
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>

    @contributor Simone Rossi <simone.rossi@epfl.ch>

    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/OneDimensionalModel_Physics.hpp>

namespace LifeV
{

std::map< std::string, OneDimensionalModel_PhysicsTypes > OneDimensionalModel_PhysicsMap;

// ===================================================
// Constructors
// ===================================================
OneDimensionalModel_Physics::OneDimensionalModel_Physics() :
    M_data      ()
{
}

OneDimensionalModel_Physics::OneDimensionalModel_Physics( const dataPtr_Type data ) :
    M_data      ( data )
{
}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_Physics::celerity0( const UInt& i ) const
{
    return std::sqrt( M_data -> beta0(i) * M_data -> beta1(i) / M_data -> densityRho() );
}

Real
OneDimensionalModel_Physics::elasticPressure( const Real& A, const UInt& i ) const
{
    return ( M_data -> beta0(i) * ( std::pow( A/M_data -> area0(i), M_data -> beta1(i) ) - 1 ) ) + M_data -> externalPressure();
}

Real
OneDimensionalModel_Physics::viscoelasticPressure( const Real& Anp1, const Real& An, const Real& Anm1, const Real& timeStep, const UInt& i ) const
{
    Real area(Anp1);
    if ( M_data -> linearizeStringModel() )
        area = M_data -> area0(i);

    return M_data -> viscoelasticModulus() / ( 2*sqrt( Pi*area ) ) * dAdt(Anp1, An, Anm1, timeStep);
}

Real
OneDimensionalModel_Physics::dAdt( const Real& Anp1, const Real& An, const Real& Anm1, const Real& timeStep ) const
{
    if ( M_data -> DPdtSteps() == 0 )
        return ( Anp1 - An ) / timeStep;
    else
        return ( 3 / 2 * Anp1 - 2 * An + 1 / 2 * Anm1 ) / timeStep;
}

Real
OneDimensionalModel_Physics::dPdA( const Real& A, const UInt& i ) const
{
    return M_data -> beta0(i) * M_data -> beta1(i) * std::pow( A / M_data -> area0(i), M_data -> beta1(i) ) / A;
}

Real
OneDimensionalModel_Physics::dAdP( const Real& P, const UInt& i ) const
{
    return M_data -> area0(i) / ( M_data -> beta0(i) * M_data -> beta1(i) ) * std::pow( 1 + ( P - M_data -> externalPressure() ) / M_data -> beta0(i), 1 / M_data -> beta1(i) - 1 );
}

Real
OneDimensionalModel_Physics::fromPToA( const Real& P, const UInt& i ) const
{
    return ( M_data -> area0(i) * std::pow( ( P - M_data -> externalPressure() ) / M_data -> beta0(i) + 1, 1/M_data -> beta1(i) )  );
}

Real
OneDimensionalModel_Physics::totalPressure( const Real& A, const Real& Q, const UInt& i ) const
{
    return elasticPressure( A, i ) + M_data -> densityRho() / 2 * Q * Q / ( A * A );
}

Real
OneDimensionalModel_Physics::totalPressureDiff( const Real& A, const Real& Q, const ID& id, const UInt& i) const
{
    if ( id == 1 ) // dPt/dA
        return dPdA( A, i ) - M_data -> densityRho() * Q * Q / ( A * A * A );

    if ( id == 2 ) // dPt/dQ
        return M_data -> densityRho() * Q / ( A * A);

    ERROR_MSG("Total pressure's differential function has only 2 components.");
    return -1.;
}

void
OneDimensionalModel_Physics::stiffenVesselLeft( const Real& xl,         const Real& xr,
                                                const Real& factor,     const Real& alpha,
                                                const Real& delta,      const Real& n,
                                                const Real& minDeltaX, const UInt& yesAdaptive )
{
    /* Stiffen Left boundary with a fifth order polynomial law
       if (alpha-delta/2) <= x < alpha
       coeff = ( (alpha + delta/2) - x )^5 * 2^4 / delta^5;

       if alpha <= x <= alpha + delta/2
       coeff = 1 - ( x - (alpha - delta/2))^5 * 2^4 / delta^5;
    */
    if (yesAdaptive)
    {
        Real ratio, n_elem_delta,n_elem_r,n_elem_l;

        UInt iz=0, alpha_iz;

        //      alpha_iz = static_cast<UInt>( alpha / (xr-xl) * static_cast<Real>( M_data -> numberOfElements()-1 ) );
        alpha_iz = static_cast<int>( std::floor( (alpha - delta / 2 ) / minDeltaX + 0.5 ) ) +
                   ( ( M_data -> numberOfElements() - 1 ) -
                     static_cast<int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) ) -
                     static_cast<int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) ) ) / 2;

        //      n_elem_r = static_cast<Real>( (M_data -> numberOfElements()-1) - alpha_iz );
        n_elem_r = ( ( M_data -> numberOfElements() - 1 ) - alpha_iz ) -
                   static_cast<int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) );

        //      n_elem_l = static_cast<Real>( alpha_iz );
        n_elem_l = alpha_iz -
                   static_cast<int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) );

        n_elem_delta = static_cast<Real>( M_data -> numberOfElements() - 1 ) / ( xr - xl ) * delta;

        //      n_elem_delta = n_elem_r + n_elem_l;
        Real x_current,deltax,deltax_adaptive,deltax_uniform;

        x_current = alpha;

        do
        {
            //! beta0
            // fifth order
            ratio=( ( ( alpha + delta / 2 ) - x_current ) / delta);

            M_data -> setdBeta0dz( M_data -> beta0(alpha_iz + iz) *
                                 ( factor * (- n / delta) * ( pow(2,(n-1)) * pow(ratio, (n-1)) ) ), alpha_iz + iz );
            M_data -> setdBeta0dz( M_data -> dBeta0dz(alpha_iz + iz), alpha_iz - iz );

            M_data -> setBeta0( M_data -> beta0(alpha_iz + iz) * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n) ) ), alpha_iz + iz );
            M_data -> setBeta0( M_data -> beta0(alpha_iz + iz) / ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n) ) )
                              * ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) ), alpha_iz - iz );

            // first order
            //        M_dPressbeta0dz[iz] = M_Pressbeta0[iz] * ( -factor * (n / delta) *
            //                       ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
            //M_Pressbeta0[iz] = M_Pressbeta0[iz] * ( 1 + factor * ratio );


            deltax_adaptive = ( -1 / n_elem_delta ) * ( 1 / ( ( -n / delta) * pow(2,(n-1) ) *
                                                            pow( ratio , (n-1) ) ) );

            deltax_uniform = ( ( alpha + delta / 2) - x_current ) / ( n_elem_r - iz );

            iz++;

            deltax = ( ( deltax_adaptive < deltax_uniform ) && ( iz < n_elem_r ) )
                     ? deltax_adaptive : deltax_uniform;

            //( xr - xl ) / M_edgeList.size();
            ASSERT_PRE( deltax > 0 , "The left point is on the right..." );

            x_current += deltax;

        }
        while ( ( x_current < ( alpha + delta/2 ) ) && ( ( alpha_iz - ( iz - 1 ) ) > 0) );

        if ( ( alpha_iz - ( iz - 1 ) ) > 0)
        {
            do
            {
                M_data -> setBeta0( M_data -> beta0(alpha_iz - iz) * ( 1 + factor ), alpha_iz - iz );
                iz++;
            }
            while ( (alpha_iz - ( iz - 1 ) ) > 0 );

            //      M_PressBeta0[0] = M_PressBeta0[0] *
            //  ( 1 + factor );
        }
        else
            std::cout << "[stiffenVesselRight] error! out of left boundary" << std::endl;
    }

    else
    {
        UInt iz=0;

        Real ratio, x_current=xl, deltax;

        deltax=( xr - xl ) / static_cast<Real>(M_data -> numberOfElements() - 1 );

        while ( ( x_current < ( alpha - delta / 2 ) ) && ( iz < M_data -> numberOfElements() ) )
        {
            M_data -> setBeta0( M_data->beta0(iz) * ( 1 + factor ), iz );
            iz++;
            x_current+=deltax;
        }

        while ( (x_current < alpha) && (iz < M_data -> numberOfElements()) )
        {
            ratio=(( x_current - (alpha-delta/2) ) / delta);

            M_data -> setdBeta0dz( M_data -> beta0(iz) * ( factor * (- n / delta ) * ( pow(2,(n-1)) * pow(ratio,(n-1) ) ) ), iz );

            M_data -> setBeta0( M_data -> beta0(iz) * ( 1 + factor * ( 1 - pow(2,(n-1)) * pow(ratio,n) ) ), iz );
            iz++;
            x_current+=deltax;
        }

        while ( ( x_current < ( alpha + delta / 2 ) ) && (iz < M_data -> numberOfElements()) )
        {
            ratio=( ( ( alpha + delta / 2 ) - x_current ) / delta );

            M_data -> setdBeta0dz( M_data -> beta0(iz) * ( factor * ( -n / delta) * ( pow(2,(n-1)) * pow(ratio,(n-1) ) ) ), iz );

            M_data -> setBeta0( M_data -> beta0(iz) * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n) ) ), iz );
            iz++;
            x_current += deltax;
        }
    }
}

void
OneDimensionalModel_Physics::stiffenVesselRight( const Real& xl,     const Real& xr,
                                                 const Real& factor, const Real& alpha,
                                                 const Real& delta,  const Real& n,
                                                 const Real& minDeltaX, const UInt& yesAdaptive )
{
    Debug( 6320 ) << "stiffenVesselright ...\n";
    /* Stiffen Left boundary with a fifth order polynomial law

       if (alpha-delta/2) <= x < alpha
       coeff = ( x - (alpha - delta/2) )^5 * 2^4 / delta^5;

       if alpha <= x <= alpha + delta/2
       coeff = 1 + ( x - (alpha + delta/2))^5 * 2^4 / delta^5;
    */
    if ( yesAdaptive )
    {
        Real ratio, n_elem_delta,n_elem_r,n_elem_l;

        UInt iz=0, alpha_iz;

        //      alpha_iz = static_cast<UInt>( alpha / (xr-xl) * ( static_cast<Real>( M_data -> numberOfElements()-1 ) ) );
        alpha_iz = static_cast<int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) ) +
                   ( (M_data -> numberOfElements() - 1 ) -
                     static_cast<int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) ) -
                     static_cast<int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) ) ) / 2;

        n_elem_delta = static_cast<Real>( M_data -> numberOfElements() - 1 ) / ( xr - xl ) * delta;

        //      n_elem_r = static_cast<Real>( (M_data -> numberOfElements()-1) - alpha_iz );
        n_elem_r = ( ( M_data -> numberOfElements() - 1 ) - alpha_iz ) -
                   static_cast<int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) );

        //      n_elem_l = static_cast<Real>( alpha_iz );
        n_elem_l = alpha_iz -
                   static_cast<int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) );

        Real x_current,deltax,deltax_adaptive,deltax_uniform;

        x_current = alpha;

        do
        {
            //! beta0
            // fifth order
            ratio=( ( ( alpha + delta / 2 ) - x_current ) / delta );

            M_data -> setdBeta0dz( M_data -> beta0( alpha_iz + iz ) * ( factor * ( n / delta) *
                                                                  ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), alpha_iz + iz );

            M_data -> setdBeta0dz( M_data -> dBeta0dz(alpha_iz + iz), alpha_iz - iz );

            M_data -> setBeta0( M_data -> beta0(alpha_iz + iz) *
                              ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) ), (alpha_iz + iz) );

            M_data -> setBeta0( M_data -> beta0(alpha_iz + iz) /
                              ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) )
                              * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n) ) ), alpha_iz - iz );

            // first order
            //        M_data -> dBeta0dz(iz) = M_data -> Beta0(iz) * ( -factor * (n / delta) *
            //                       ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
            //M_data -> Beta0(iz) = M_data -> Beta0(iz) * ( 1 + factor * ratio );

            deltax_adaptive = ( -1 / n_elem_delta ) *
                              ( 1 / ( (-n / delta ) * pow(2,(n-1)) *
                                      pow( ratio , (n-1) )
                                    )
                              );

            deltax_uniform = ( ( alpha + delta / 2) - x_current ) / ( n_elem_r - iz );

            iz++;

            deltax = ( ( deltax_adaptive < deltax_uniform ) && ( iz < n_elem_r) )
                     ? deltax_adaptive : deltax_uniform;

            //( xr - xl ) / M_edgeList.size();
            ASSERT_PRE( deltax > 0 ,
                        "The left point is on the right..." );

            x_current += deltax;

        }
        while ( x_current < ( alpha + delta / 2 ) && ( ( alpha_iz - ( iz - 1 ) ) > 0) );

        if ( ( alpha_iz + iz ) <= (M_data -> numberOfElements() -1 ) )
        {
            do
            {
                M_data -> setBeta0( M_data -> beta0(alpha_iz + iz ) * ( 1 + factor ), alpha_iz + iz );
                iz++;
            }
            while ( ( alpha_iz + iz - 1 ) < ( M_data -> numberOfElements() -1 ) );

            //      M_PressBeta0[0] = M_PressBeta0[0] *
            //  ( 1 + factor );
        }
        else
            std::cout << "\n[stiffenVesselRight] error! out of right boundary" << std::endl;
    }
    else
    {
        UInt iz = M_data -> numberOfElements()-1;

        Real ratio, x_current=xr, deltax;

        deltax = ( xr - xl ) / static_cast<Real>(M_data -> numberOfElements() - 1 );

        while ( ( x_current > ( alpha + delta / 2 ) ) && ( ( iz + 1 ) > 0 ) )
        {
            M_data -> setBeta0( M_data -> beta0(iz) * ( 1 + factor ), iz );
            iz--;
            x_current -= deltax;
        }

        while ( ( x_current > alpha ) && ( ( iz + 1 ) > 0 ) )
        {
            ratio=( ( ( alpha + delta / 2 ) - x_current ) / delta );

            M_data -> setdBeta0dz( M_data -> beta0(iz) * ( factor * ( n / delta) *  ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), iz );

            M_data -> setBeta0( M_data -> beta0(iz) * ( 1 + factor * ( 1 - pow(2,(n-1)) * pow(ratio,n) ) ), iz );
            iz--;
            x_current -= deltax;
        }

        while ( ( x_current > ( alpha - delta / 2 ) ) && ( ( iz + 1 ) > 0 ) )
        {
            ratio = ( ( x_current - ( alpha - delta / 2 ) ) / delta );

            M_data -> setdBeta0dz( M_data -> beta0(iz) * ( factor * ( n / delta) * ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), iz );

            M_data -> setBeta0( M_data -> beta0(iz) * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n) ) ), iz );
            iz--;
            x_current -= deltax;
        }
    }
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Physics::setData( const dataPtr_Type& Data )
{
    M_data = Data;
}

// ===================================================
// Get Methods
// ===================================================
OneDimensionalModel_Physics::dataPtr_Type
OneDimensionalModel_Physics::data() const
{
    return M_data;
}

}
