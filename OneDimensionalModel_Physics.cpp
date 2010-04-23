//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a base class providing physical operations for the 1D model data.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date 01-07-2004
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 13-04-2010
 */

#include <lifemc/lifesolver/OneDimensionalModel_Physics.hpp>

namespace LifeV {

std::map< std::string, OneDimensionalModel_PhysicsTypes > OneDimensionalModel_PhysicsMap;

// ===================================================
// Constructors
// ===================================================
OneDimensionalModel_Physics::OneDimensionalModel_Physics() :
    M_Data      ()
{
}

OneDimensionalModel_Physics::OneDimensionalModel_Physics( const Data_PtrType Data ) :
    M_Data      ( Data )
{
}

// ===================================================
// Methods
// ===================================================
Real
OneDimensionalModel_Physics::Celerity0( const UInt& i ) const
{
    return std::sqrt( M_Data->Beta0(i) * M_Data->Beta1(i) / M_Data->DensityRho() );
}

Real
OneDimensionalModel_Physics::Length() const
{
    return M_Data->xRight() - M_Data->xLeft();
}

ScalVec
OneDimensionalModel_Physics::pressure( const Real& _A,
                                             const Real& _A_n,
                                             const Real& _A_nm1,
                                             const Real& dt,
                                             const UInt& indz,
                                             const UInt& steps,
                                             const bool& visco,
                                             const bool& linearized ) const
{
    ScalVec _a(2), _b(2), _c(2), area(2), result(4);

    _a(0) =  1.;
    _b(0) = -1.;
    _c(0) =  0.;

    _a(1) =  3./2.;
    _b(1) = -2.;
    _c(1) =  1./2.;

    area(0) = _A;
    area(1) = M_Data->Area0(indz);

    Real _pi( 4*std::atan(1) );

    result(3) = ( _a(steps) * _A + _b(steps) * _A_n + _c(steps) * _A_nm1 ) / dt; //> dA/dt

    result(2) = M_Data->Gamma() * result(3) / ( 2*sqrt(_pi*area(linearized)) );               //> visc_component

    result(1) = pressure( _A, indz );                                                  //> elast_component

    result(0) = result(1) + visco * result(2);

    return result;
}

Real
OneDimensionalModel_Physics::pressure( const Real& _A, const UInt& indz ) const
{
    return ( M_Data->Beta0(indz) * ( std::pow( _A/M_Data->Area0(indz), M_Data->Beta1(indz) ) - 1 ) );
}

Real
OneDimensionalModel_Physics::pressureDiff( const Real& _A, const UInt& indz ) const
{
    //Real AoverA0POWbeta1( std::pow( _A / M_Area0[indz], M_PressBeta1[indz] ) );
    //std::cout << indz << " -> " <<  M_PressBeta0[indz] << " " <<  M_PressBeta1[indz] << " " <<  AoverA0POWbeta1 << " " << _A << std::endl;
    return M_Data->Beta0(indz) * M_Data->Beta1(indz) * std::pow( _A / M_Data->Area0(indz), M_Data->Beta1(indz) ) / _A;
}

Real
OneDimensionalModel_Physics::totalPressure( const Real& _A, const Real& _Q, const UInt& indz ) const
{
    return ( pressure( _A, indz ) + (M_Data->DensityRho()/2) * _Q * _Q / ( _A * _A ) );
}

Real
OneDimensionalModel_Physics::totalPressureDiff( const Real& _A, const Real& _Q, const ID& i, const UInt& indz) const
{
    Real vel( _Q / _A );

    if( i == 1 ) //! dPt/dA
    {
        Real dPtdA( pressureDiff( _A, indz ) - M_Data->DensityRho() * vel * vel / _A );

        return dPtdA;
    }

    if( i == 2 ) //! dPt/dQ
    {
        return ( M_Data->DensityRho() * vel / _A );
    }

    ERROR_MSG("Total pressure's differential function has only 2 components.");
    return -1.;
}

Real
OneDimensionalModel_Physics::A_from_P( const Real& _P, const UInt& indz ) const
{
    return ( M_Data->Area0(indz) * std::pow( _P / M_Data->Beta0(indz) + 1, 1/M_Data->Beta1(indz) )  );
}

void
OneDimensionalModel_Physics::stiffenVesselLeft( const Real& xl,         const Real& xr,
                                                const Real& factor,     const Real& alpha,
                                                const Real& delta,      const Real& n,
                                                const Real& min_deltax, const UInt& yesAdaptive )
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

        //      alpha_iz = static_cast<UInt>( alpha / (xr-xl) * static_cast<Real>( M_Data->nbElem()-1 ) );
        alpha_iz = static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) +
            ( (M_Data->nbElem() - 1) -
              static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) ) -
              static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) ) / 2;

        //      n_elem_r = static_cast<Real>( (M_Data->nbElem()-1) - alpha_iz );
        n_elem_r = ( (M_Data->nbElem()-1) - alpha_iz ) -
            static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) );

        //      n_elem_l = static_cast<Real>( alpha_iz );
        n_elem_l = alpha_iz -
            static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) );

        n_elem_delta = static_cast<Real>(M_Data->nbElem() - 1) / (xr - xl) * delta;

        //      n_elem_delta = n_elem_r + n_elem_l;


        Real x_current,deltax,deltax_adaptive,deltax_uniform;

        x_current = alpha;

        do
        {
            //! beta0
            // fifth order
            ratio=(( (alpha + delta/2) - x_current ) / delta);

            M_Data->setdBeta0dz( M_Data->Beta0(alpha_iz+iz) *
                ( factor * (- n / delta) * ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), alpha_iz+iz );
            M_Data->setdBeta0dz( M_Data->dBeta0dz(alpha_iz+iz), alpha_iz-iz );

            M_Data->setBeta0( M_Data->Beta0(alpha_iz+iz) * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n))), alpha_iz+iz );
            M_Data->setBeta0( M_Data->Beta0(alpha_iz+iz) / ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n)))
                              * ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) ), alpha_iz-iz );

            // first order
            //        M_dPressBeta0dz[iz] = M_PressBeta0[iz] * ( -factor * (n / delta) *
            //                       ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
            //M_PressBeta0[iz] = M_PressBeta0[iz] * ( 1 + factor * ratio );


            deltax_adaptive = ( -1/n_elem_delta ) * ( 1 / ( (-n/delta) * pow(2,(n-1)) *
                              pow( ratio , (n-1) ) ) );

            deltax_uniform = ( (alpha+delta/2) - x_current) / ( n_elem_r - iz );

            iz++;

            deltax = ( ( deltax_adaptive < deltax_uniform ) && (iz < n_elem_r) )
                ? deltax_adaptive : deltax_uniform;

            //( xr - xl ) / M_edgeList.size();
            ASSERT_PRE( deltax > 0 , "The left point is on the right..." );

            x_current += deltax;

        }
        while ( ( x_current < ( alpha + delta/2 ) ) && ( (alpha_iz - (iz - 1)) > 0) );

        if ( ( alpha_iz - (iz - 1)) > 0)
        {
            do
            {
                M_Data->setBeta0( M_Data->Beta0(alpha_iz-iz) * ( 1 + factor ), alpha_iz-iz );
                iz++;
            }
            while ( (alpha_iz - (iz - 1)) > 0 );

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

        deltax=(xr-xl)/static_cast<Real>(M_Data->nbElem()-1);

        while ( (x_current < (alpha - delta/2)) && (iz < M_Data->nbElem()) )
        {
            M_Data->setBeta0( M_Data->Beta0(iz) * ( 1 + factor ), iz );
            iz++;
            x_current+=deltax;
        }

        while ( (x_current < alpha) && (iz < M_Data->nbElem()) )
        {
            ratio=(( x_current - (alpha-delta/2) ) / delta);

            M_Data->setdBeta0dz( M_Data->Beta0(iz) * ( factor * (- n / delta) * ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), iz );

            M_Data->setBeta0( M_Data->Beta0(iz) * ( 1 + factor * ( 1 - pow(2,(n-1)) * pow(ratio,n)) ), iz );
            iz++;
            x_current+=deltax;
        }

        while ( (x_current < (alpha+delta/2)) && (iz < M_Data->nbElem()) )
        {
            ratio=(( (alpha+delta/2) - x_current ) / delta);

            M_Data->setdBeta0dz( M_Data->Beta0(iz) * ( factor * ( -n / delta) * ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), iz );

            M_Data->setBeta0( M_Data->Beta0(iz) * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n)) ), iz );
            iz++;
            x_current+=deltax;
        }
    }
}

void
OneDimensionalModel_Physics::stiffenVesselRight( const Real& xl,     const Real& xr,
                                                 const Real& factor, const Real& alpha,
                                                 const Real& delta,  const Real& n,
                                                 const Real& min_deltax, const UInt& yesAdaptive )
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

        //      alpha_iz = static_cast<UInt>( alpha / (xr-xl) * ( static_cast<Real>( M_Data->nbElem()-1 ) ) );
        alpha_iz = static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) +
            ( (M_Data->nbElem() - 1) -
              static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) ) -
              static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) ) / 2;

        n_elem_delta = static_cast<Real>(M_Data->nbElem() - 1) / (xr - xl) * delta;

        //      n_elem_r = static_cast<Real>( (M_Data->nbElem()-1) - alpha_iz );
        n_elem_r = ( (M_Data->nbElem()-1) - alpha_iz ) -
            static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) );

        //      n_elem_l = static_cast<Real>( alpha_iz );
        n_elem_l = alpha_iz -
            static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) );

        Real x_current,deltax,deltax_adaptive,deltax_uniform;

        x_current = alpha;

        do
        {
            //! beta0
            // fifth order
            ratio=(( (alpha + delta/2) - x_current) / delta);

            M_Data->setdBeta0dz( M_Data->Beta0( alpha_iz+iz ) * ( factor * ( n / delta) *
                                 ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), alpha_iz+iz );

            M_Data->setdBeta0dz( M_Data->dBeta0dz(alpha_iz+iz), alpha_iz-iz );

            M_Data->setBeta0( M_Data->Beta0(alpha_iz+iz) *
                              ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) ), (alpha_iz+iz) );

            M_Data->setBeta0( M_Data->Beta0(alpha_iz+iz) /
                                ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) )
                              * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n) ) ), alpha_iz-iz );

            // first order
            //        M_Data->dBeta0dz(iz) = M_Data->Beta0(iz) * ( -factor * (n / delta) *
            //                       ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
            //M_Data->Beta0(iz) = M_Data->Beta0(iz) * ( 1 + factor * ratio );

            deltax_adaptive = ( -1/n_elem_delta ) *
                ( 1 / ( (-n/delta) * pow(2,(n-1)) *
                        pow( ratio , (n-1) )
                        )
                  );

            deltax_uniform = ( (alpha+delta/2) - x_current) / ( n_elem_r - iz );

            iz++;

            deltax = ( ( deltax_adaptive < deltax_uniform ) && (iz < n_elem_r) )
                ? deltax_adaptive : deltax_uniform;

            //( xr - xl ) / M_edgeList.size();
            ASSERT_PRE( deltax > 0 ,
                        "The left point is on the right..." );

            x_current += deltax;

        }
        while ( x_current < ( alpha + delta/2 ) && ( (alpha_iz - (iz - 1)) > 0) );

        if ( ( alpha_iz + iz ) <= (M_Data->nbElem() -1) )
        {
            do
            {
                M_Data->setBeta0( M_Data->Beta0(alpha_iz+iz) * ( 1 + factor ), alpha_iz+iz );
                iz++;
            }
            while ( (alpha_iz + iz - 1) < (M_Data->nbElem() -1) );

            //      M_PressBeta0[0] = M_PressBeta0[0] *
            //  ( 1 + factor );
        }
        else
            std::cout << "\n[stiffenVesselRight] error! out of right boundary" << std::endl;
    }
    else
    {
        UInt iz=M_Data->nbElem()-1;

        Real ratio, x_current=xr, deltax;

        deltax=(xr-xl)/static_cast<Real>(M_Data->nbElem()-1);

        while ( (x_current > (alpha+delta/2)) && ((iz+1) > 0) )
        {
            M_Data->setBeta0( M_Data->Beta0(iz) * ( 1 + factor ), iz );
            iz--;
            x_current-=deltax;
        }

        while ( (x_current > alpha) && ((iz+1) > 0 ) )
        {
            ratio=(( (alpha+delta/2) - x_current ) / delta);

            M_Data->setdBeta0dz( M_Data->Beta0(iz) * ( factor * ( n / delta) *  ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), iz );

            M_Data->setBeta0( M_Data->Beta0(iz) * ( 1 + factor * ( 1 - pow(2,(n-1)) * pow(ratio,n)) ), iz );
            iz--;
            x_current-=deltax;
        }

        while ( (x_current > (alpha-delta/2)) && ((iz+1) > 0) )
        {
            ratio=(( x_current - (alpha-delta/2) ) / delta);

            M_Data->setdBeta0dz( M_Data->Beta0(iz) * ( factor * ( n / delta) * ( pow(2,(n-1)) * pow(ratio,(n-1)) ) ), iz );

            M_Data->setBeta0( M_Data->Beta0(iz) * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n)) ), iz );
            iz--;
            x_current-=deltax;
        }
    }
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Physics::SetData( const Data_PtrType Data )
{
    M_Data = Data;
}

// ===================================================
// Get Methods
// ===================================================
OneDimensionalModel_Physics::Data_PtrType
OneDimensionalModel_Physics::Data() const
{
    return M_Data;
}

}
