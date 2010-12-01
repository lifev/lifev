/*
 This file is part of the LifeV library

 Authors: A. Veneziani
          C. Prud'homme <christophe.prudhomme@epfl.ch>
          C. Winkelmann
          M. Perego
          U. Villa

 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
/*!
  \file bdf.hpp
  \author A. Veneziani
  \author C. Prud'homme
  \author C. Winkelmann
  \author M. Perego
  \author U. Villa

  File containing a class for an easy handling of different order time
  discretizations/extrapolations BDF based

*/
#ifndef _BDF_VARIABLE_TIMESTEP_H
#define _BDF_VARIABLE_TIMESTEP_H
#include <string>
#include <iostream>
#include <numeric>
#include <life/lifearray/EpetraVector.hpp>


namespace LifeV
{
const uint BDF_MAX_ORDER = 6;


/*!
  \class Bdf
  \brief Backward differencing formula time discretization

  A differential equation of the form

  \f$ M u' = A u + f \f$

  is discretized in time as

  \f$ M p'(t_{k+1}) = A u_{k+1} + f_{k+1} \f$

  where p denotes the polynomial of order n in t that interpolates
  (t_i,u_i) for i = k-n+1,...,k+1.

  The approximative time derivative \f$ p'(t_{k+1}) \f$ is a linear
  combination of state vectors u_i:

  \f$ p'(t_{k+1}) = \frac{1}{\Delta t} (\alpha_0 u_{k+1} - \sum_{i=0}^n \alpha_i u_{k+1-i} )\f$

  Thus we have

  \f$ \frac{\alpha_0}{\Delta t} M u_{k+1} = A u_{k+1} + f + M \bar{p} \f$

  with

  \f$ \bar{p} = \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i u_{k+1-i} \f$

  This class stores the n last state vectors in order to be able to
  calculate \f$ \bar{p} \f$. It also provides alpha_i
  and can extrapolate the the new state from the n last states with a
  polynomial of order n-1:

  \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$

  The class allows to change the time-steps, recomputing the coefficients.
  to change \Delta_t at each time-step,
  use shift_right(u_curr, deltat_new) instead of shift_right(u_curr).

  Other methods that could be helpful while adapting in time are:
  set_deltat, store_unk and restore_unk
*/
template<typename FeVector = EpetraVector>
class BdfVS
{
public:

    typedef boost::shared_ptr<FeVector> FeVector_ptr;
    typedef std::vector<FeVector_ptr>   BdfContainer;

    //! Empty constructor
    BdfVS();

    /*! Constructor
     *  @param n order of the BDF
     */
    BdfVS( const UInt n );

    ~BdfVS();

    void setup( const UInt n );

    //! Initialize
    //@{
    //! Initialize all the entries of the unknown vector to be derived with the
    //! vector u0 (duplicated if startup=0)
    /*! If startup = 0, it initializes all the entries of the unknown vectors with a given function
        The array of initial conditions needed by the selected BDF is
        initialized as follows: _unk=[ u0(t0), u0(t0-dt), u0(t0-2*dt), ...]

        When startup = true, only _unk[0] is initialized (with u_0).
    	At the first step, class Bfd computes the coefficients of BDF1,
    	at the second step, the coefficient of BDF2, and so on,
    	until it reaches the wanted method.
    	Second order of accuracy is obtained using this procedure to startup BDF2.
    	Using the startup procedure with BDF3, the accuracy order is limited to the second order;
    	anyway it is better to use the startup procedure instead of taking u_{-2} = u_{-1} = u_0.
    */
    void initialize_unk( FeVector u0, Real const dt, bool startup = 0  );

    //! Initialize all the entries of the unknown vector to be derived with a
    //! set of vectors uv0
    void initialize_unk( std::vector<FeVector> uv0, Real const dt);

    /* Initialize all the entries of the unknonwn vectors with a given function
        The array of initial conditions needed by the selected BDF is
        initialized as follows: _unk=[ u0(t0), u0(t0-dt), u0(t0-2*dt), ...]
        For the space dependence of the initial conditions we need informations
        on:
        -# the function we want to interpolate
        -# a vector type u_vect, u_vect is a corrected mapped vector in input, and contains u0(t0) in output
    	-# a finite element space
        -# which is the initial time (t0) and the time step (for solutions
           before the initial instant)
        Based on the NavierStokesHandler::initialize by M. Fernandez
    */
    template<typename FunctClass, typename FESpaceClass>
    void initialize_unk(  const FunctClass& u0, FeVector& u_vect,
                          FESpaceClass & feSpace, Real t0, Real dt  );

    //@}

    //! Core methods
    //@{

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula divided by dt
    FeVector time_der_dt( ) const;

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula divided by dt (backward compatibility version, will be discontinued)
    FeVector time_der( Real dt = 1  ) const;


    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    FeVector extrap() const;

    /*! Update the vectors of the previous time steps by shifting on the right
     *  the old values.
     *  @param u_curr current (new) value of the state vector
     */
    void shift_right( FeVector const& u_curr );


    /*! Update the vectors of the previous time steps by shifting on the right
     *  the old values. Set the the new time step and shift right the old time steps.
     *  @param u_curr current (new) value of the state vector
     *  @param deltat_new new value of the the time step
     */
    void shift_right( FeVector const& u_curr, Real deltat_new);

    //@}

    //! Getters
    //@{
    //! Return the i-th coefficient of the time derivative alpha_i
    double coeff_der( UInt i ) const;

    //! Return the i-th coefficient of the time derivative alpha_i divided by dt
    double coeff_der_dt(UInt i) const;

    //! Return the i-th coefficient of the time extrapolation beta_i
    double coeff_ext( UInt i ) const;

    //! Return the vector of the time steps, ordered starting from the most recent one.
    inline Vector vec_deltat() {return _M_deltat;};

    //! Return a vector with the last n state vectors
    inline const BdfContainer& unk() const {return _M_unknowns;}
    inline const BdfContainer& unk_bk() const {return _M_unknowns_bk;}
    //@}

    //! Setters
    //@{
    //! Replace the current time-step with deltat and computes the coefficients a_i and beta_i as functions of _M_delta_t.
    void set_deltat(Real deltat);
    //@}

    //! Other methods
    //@{
    //! save the current vector _M_unknowns and the current vector _M_deltat
    void store_unk();

    //! restore the vector _M_unknowns and the vector _M_deltat with the ones saved with store_unk()
    void restore_unk();

    //! restore_unk() + set_deltat()
    void restore_unk(Real delta_t);

    void showMe() const;
    //@}
private:

    //! Private methods
    //@{
    /*! Computes the coefficients a_i and beta_i as functions of _M_delta_t.*/
    void comput_coeff();
    //@}

    //! Order of the BDF derivative/extrapolation: the time-derivative
    //! coefficients vector has size n+1, the extrapolation vector has size n
    UInt _M_order;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    Vector _M_alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    Vector _M_beta;

    //! Vector \f$ \delta_t \f$ of time intervals
    Vector _M_deltat;
    Vector _M_deltat_bk;


    //! Last n state vectors
    BdfContainer _M_unknowns;
    BdfContainer _M_unknowns_bk;

};

//============================================================================
// Empty constructor
template<typename FeVector>
BdfVS<FeVector>::BdfVS()
        :
        _M_order( 0 ),
        _M_alpha( _M_order + 1 ),
        _M_beta( _M_order ),
        _M_deltat(ScalarVector(_M_order, 1.))
{}


template<typename FeVector>
BdfVS<FeVector>::BdfVS( const UInt n )
{
    setup( n );
}


template<typename FeVector>
BdfVS<FeVector>::~BdfVS()
{}

//============================================================================
// Initialize methods
template<typename FeVector>
void BdfVS<FeVector>::initialize_unk( FeVector u0, Real const dt, bool startup )
{
    _M_unknowns.resize(0);
    for (UInt i=0; i<_M_order; i++)
    {
        FeVector_ptr feVectorPtr(new FeVector(u0));
        _M_unknowns.push_back(feVectorPtr);
        _M_deltat[i] = dt;
    }

    if (startup)
    {
        for (UInt i=1; i<_M_order; i++)
            _M_deltat[i] = 1e20;
    }

    comput_coeff();

    return ;
}

template<typename FeVector>
void BdfVS<FeVector>::initialize_unk(std::vector<FeVector> uv0, Real const dt)
{
    UInt n0(uv0.size());
    // Check if uv0 has the right dimensions
    ASSERT( n0 >= _M_order, "Initial data are not enough for the selected BDF" );

    _M_unknowns.resize(0);
    for (UInt i=0; i<_M_order; i++)
    {
        FeVector_ptr tmp(new FeVector(uv0[i]));
        _M_unknowns.push_back(tmp);
        _M_deltat[i] = dt;
    }

    comput_coeff();

    // if n0>n, only the first n inital data will be considered
    if ( n0 > _M_order )
    {
        std::cout << "The initial data set is larger than needed by the BDF."
                  << std::endl;
        std::cout << "Only the first " << _M_order << " data will be considered. "
                  << std::endl;
    }

    return ;
}

template<typename FeVector>
template<typename FunctClass, typename FESpaceClass>
void BdfVS<FeVector>::initialize_unk(  const FunctClass& u0, FeVector& u_vect, FESpaceClass & feSpace, Real t0, Real dt  )
{

    _M_unknowns.resize(0);

    for (UInt i = 0 ; i < _M_order; ++i )
    {
        FeVector_ptr tmp(new FeVector(u_vect));
        _M_unknowns.push_back(tmp);
        feSpace.interpolate( u0, *_M_unknowns[i], t0-i*dt);
        _M_deltat[i] = dt;
    }
    u_vect = *_M_unknowns[0];
    return ;
}

//==============================================================================================
// Core Methods

template<typename FeVector>
void BdfVS<FeVector>::setup( const UInt n )
{
    _M_order = n;
    _M_alpha.resize( _M_order + 1 );
    _M_beta.resize( _M_order );
    _M_deltat = ScalarVector(_M_order, 1.);

    if ( n <= 0 || n > BDF_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
        << " you want to use BDF order " << n << "\n"
        << " we support BDF order from 1 to " << BDF_MAX_ORDER << "\n";
        throw std::invalid_argument( __ex.str() );
    }
    comput_coeff();
    _M_unknowns.reserve( n );
}


template<typename FeVector>
FeVector
BdfVS<FeVector>::time_der_dt() const
{
    FeVector ut(*_M_unknowns[0]);
    ut *= _M_alpha[ 1 ]/_M_deltat[0];

    for ( UInt i = 1; i < _M_order; ++i )
        ut += (_M_alpha[ i + 1 ]/_M_deltat[0]) * *_M_unknowns[ i ];

    return ut;
}

template<typename FeVector>
FeVector
BdfVS<FeVector>::time_der(Real dt) const
{
    FeVector ut(*_M_unknowns[0]);
    ut *= _M_alpha[ 1 ]/dt;

    for ( UInt i = 1; i < _M_order; ++i )
        ut += (_M_alpha[ i + 1 ]/dt) * *_M_unknowns[ i ];

    return ut;
}


template<typename FeVector>
FeVector
BdfVS<FeVector>::extrap() const
{
    FeVector ue(*_M_unknowns[0]);
    ue *= _M_beta[ 0 ];

    for ( UInt i = 1; i < _M_order; ++i )
        ue += _M_beta[ i ] * *_M_unknowns[ i ];

    return ue;
}

template<typename FeVector>
void
BdfVS<FeVector>::shift_right( FeVector const& unk_curr )
{
    typedef typename BdfContainer::iterator BdfContIt;
    BdfContIt it = _M_unknowns.end() - 1;
    BdfContIt itm1 = _M_unknowns.end() - 1;
    BdfContIt itb = _M_unknowns.begin();

    for ( ; it != itb; --it )
    {
        itm1--;
        *it = *itm1;
    }
    FeVector_ptr tmp(new FeVector(unk_curr));
    *itb = tmp;

    for (UInt i=_M_order-1; i>0; i--)
        _M_deltat[i] = _M_deltat[i-1];
    comput_coeff();
}


template<typename FeVector>
void
BdfVS<FeVector>::shift_right( FeVector const& unk_curr, Real deltat_new)
{

    typedef typename BdfContainer::iterator BdfContIt;
    BdfContIt it = _M_unknowns.end() - 1;
    BdfContIt itm1 = _M_unknowns.end() - 1;
    BdfContIt itb = _M_unknowns.begin();

    for ( ; it != itb; --it )
    {
        itm1--;
        *it = *itm1;
    }

    FeVector_ptr tmp(new FeVector(unk_curr));
    *itb = tmp;

    for (UInt i=_M_order-1; i>0; i--)
        _M_deltat[i] = _M_deltat[i-1];
    _M_deltat[0] = deltat_new;
    comput_coeff();
}


//===============================================================================================
// Getters

template<typename FeVector>
Real
BdfVS<FeVector>::coeff_der( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i < _M_order + 1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_alpha[ i ];
}

template<typename FeVector>
Real
BdfVS<FeVector>::coeff_der_dt( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i < _M_order + 1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_alpha[ i ]/_M_deltat[0];
}

template<typename FeVector>
Real
BdfVS<FeVector>::coeff_ext( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT(  i < _M_order,
             "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_beta[ i ];
}

//====================================================================================
//Setters
template<typename FeVector>
void BdfVS<FeVector>::set_deltat(Real delta_t)
{
    _M_deltat[0] = delta_t;
    comput_coeff();
}

//====================================================================================
//Other methods
template<typename FeVector>
void BdfVS<FeVector>::store_unk()
{
    typedef typename BdfContainer::iterator BdfContIt;

    _M_unknowns_bk.resize(0);
    for (BdfContIt it = _M_unknowns.begin(); it < _M_unknowns.end(); ++it)
        _M_unknowns_bk.push_back(new FeVector_ptr(**it));

    _M_deltat_bk = _M_deltat;
}

template<typename FeVector>
void BdfVS<FeVector>::restore_unk()
{
    typedef typename BdfContainer::iterator BdfContIt;

    _M_unknowns.resize(0);
    for (BdfContIt it = _M_unknowns_bk.begin(); it < _M_unknowns_bk.end(); ++it)
        _M_unknowns.push_back(new FeVector_ptr(**it));


    _M_deltat = _M_deltat_bk;
    comput_coeff();
}

template<typename FeVector>
void BdfVS<FeVector>::restore_unk(Real delta_t)
{
    typedef typename BdfContainer::iterator BdfContIt;

    _M_unknowns.resize(0);
    for (BdfContIt it = _M_unknowns_bk.begin(); it < _M_unknowns_bk.end(); ++it)
        _M_unknowns.push_back(new FeVector_ptr(**it));

    _M_deltat = _M_deltat_bk;
    _M_deltat[0] = delta_t;
    comput_coeff();
}


template<typename FeVector>
void
BdfVS<FeVector>::showMe() const
{
    std::cout << "*** BDF Time discretization of order " << _M_order << " ***"
              << std::endl;
    std::cout << "    Coefficients: " << std::endl;
    for ( UInt i = 0; i < _M_order + 1; ++i )
        std::cout << "       alpha(" << i << ") = " << _M_alpha[ i ]
                  << std::endl;
    for ( UInt i = 0; i < _M_order; ++i )
        std::cout << "       beta (" << i << ") = " << _M_beta[ i ]
                  << std::endl;

    return ;
}


//=============================================================================================
// Private Methods

//Arbitrary order, variable time step BDF coefficients.  (Umberto Villa <uvilla@emory.edu>)
//Reference: Comincioli pag. 618-619
//Matlab implementation:
/*
    function [alpha0, alpha, beta] = compute_coef(dts, order)
    assert( (length(dts) == order) );
    rho = dts(1)./cumsum(dts);

    alpha0 = sum(rho);
    alpha = zeros(order,1);
    beta = zeros(order,1);

    for j = 1:order
        ind = setdiff([1:order],j);
        beta(j) = 1/prod(1-rho(ind)/rho(j));
        alpha(j) = rho(j)*beta(j);
    end

end
 */

template<typename FeVector>
void BdfVS<FeVector>::comput_coeff()
{
    Vector rho(ScalarVector(_M_order, _M_deltat[0]));

    Vector cumsum_dt(_M_order);
    std::partial_sum(_M_deltat.begin(), _M_deltat.end(), cumsum_dt.begin());

    // rho = _M_deltat[0]./cumsum_M_deltat;
    std::transform(rho.begin(), rho.end(), cumsum_dt.begin(), rho.begin(), std::divides<double>());

    _M_alpha[0] = boost::numeric::ublas::sum(rho);

    for (UInt j=0; j<_M_order; ++j)
    {
        double tmp(1.0);

        if ( rho[j]!=0 ) //rho[j] may be 0 with the start-up procedure
        {
            for (UInt kk=0; kk<_M_order; ++kk)
            {
                if (kk!=j)
                {
                    tmp *= (1-rho[kk]/rho[j]);
                }
            }
            _M_alpha[j+1] = rho[j]/tmp;
            _M_beta[j]    = 1./tmp;
        }
        else //we don't have enough initial data to start a higher order BDF.
        {
            _M_alpha[j+1] = 0;
            _M_beta[j] = 0;
        }
    }
}


} /*namespace LifeV*/
#endif
