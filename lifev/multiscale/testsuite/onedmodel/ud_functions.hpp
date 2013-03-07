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
 *  @brief File containing function for imposing boundary conditions of the 1D model.
 *
 *  @version 1.0
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *  @date 01-04-2007
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef ONED_FUNCTIONS_1D_H
#define ONED_FUNCTIONS_1D_H

#include <lifev/one_d_fsi/function/OneDFSIFunction.hpp>

namespace LifeV
{

//! Const - Base class for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class Constant
{
public:

    explicit Constant ( const Real& value ) :
        M_value ( value )
    {}

    virtual ~Constant() {}

    Real operator() ( const Real& /*time*/ )
    {
        return M_value;
    }

private:

    Real M_value;
};



//! Sin - Sinusoidal wave.
/*!
 *  @author Lucia Mirabella
 */
class Sin
{
public:

    /*!
      \brief The constructor

      \param[in] mean the time average of the sinus \f$ B \f$
      \param[in] scale the sinus amplitude \f$ A \f$
      \param[in] period the sinus period \f$ T \f$
      \param[in] phase the sinus phase \f$ \phi \f$
    */
    explicit Sin ( const Real mean   = 0, const Real scale  = 10, const Real period = .01, const Real phase  = 0. ) :
        M_mean (mean),
        M_scale (scale),
        M_period (period),
        M_phase (phase)
    {}

    virtual ~Sin() {}

    /*!
      \brief Compute the sinus at the specified time

      \param[in] time the time
      \return \f$ B + A sin( \frac{2 \pi t}{T} + \phi ) \f$
    */
    Real operator() ( const Real& time )
    {
        if (time < M_period)
        {
            std::cout << time << " Flux BC = " << M_mean + M_scale* std::sin (M_phase + 2 * M_PI * time / M_period) << std::endl;
            return - (M_mean + M_scale * std::sin (M_phase + 2 * M_PI * time / M_period) );
        }
        else
        {
            return 0.;
        }
    }

private:

    Real M_mean;
    Real M_scale;
    Real M_period;
    Real M_phase;
};


//! Cos_min_Sin - A superimposition of a sinusoidal and a cosinusoidal waves, whose amplitude is damped by exponential terms
/*!
 *  @author Lucia Mirabella
 */
class Cos_min_Sin
{
public:

    /*!
      \brief The constructor

      \param[in] coeff_exp_t_cos the coefficient multiplying time in the
      exponential damping the cosinus amplitude \f$ \lambda_c \f$
      \param[in] mean_cos the time average of the cosinus \f$ B_c \f$
      \param[in] amplitude_cos the cosinus amplitude \f$ A_c \f$
      \param[in] frequency_cos the cosinus period \f$ T_c \f$
      \param[in] phase_cos the cosinus phase \f$ \phi_c \f$
      \param[in] coeff_exp_t_sin the coefficient multiplying time in the
      exponential damping the sinus amplitude \f$ \lambda_s \f$
      \param[in] mean_sin the time average of the sinus \f$ B_s \f$
      \param[in] amplitude_sin the sinus amplitude \f$ A_s \f$
      \param[in] frequency_sin the sinus period \f$ T_s \f$
      \param[in] phase_sin the sinus phase \f$ \phi_s \f$
    */
    explicit Cos_min_Sin (const Real coeff_exp_t_cos = 0, const Real mean_cos = 0,
                          const Real amplitude_cos = 10, const Real frequency_cos = 8.*atan (1.),
                          const Real phase_cos = 0.,
                          const Real coeff_exp_t_sin = 0, const Real mean_sin = 0,
                          const Real amplitude_sin = 10, const Real frequency_sin = 8.*atan (1.),
                          const Real phase_sin = 0. ) :
        M_coeff_exp_t_cos (coeff_exp_t_cos),
        M_mean_cos (mean_cos),
        M_amplitude_cos (amplitude_cos),
        M_frequency_cos (frequency_cos),
        M_phase_cos (phase_cos),
        M_coeff_exp_t_sin (coeff_exp_t_sin),
        M_mean_sin (mean_sin),
        M_amplitude_sin (amplitude_sin),
        M_frequency_sin (frequency_sin),
        M_phase_sin (phase_sin)
    {}

    virtual ~Cos_min_Sin() {}

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
      ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    Real operator() ( const Real& time )
    {
        Real result = M_mean_cos
                      + M_amplitude_cos * std::cos (M_phase_cos + time * M_frequency_cos)
                      - ( M_mean_sin + M_amplitude_sin * std::sin (M_phase_sin + time * M_frequency_sin) );

        return result;
    }

    /*! \name Getters
     */
    //@{

    Real& coeff_exp_t_cos()
    {
        return M_coeff_exp_t_cos;
    }
    Real& mean_cos()
    {
        return M_mean_cos;
    }
    Real& amplitude_cos()
    {
        return M_amplitude_cos;
    }
    Real& frequency_cos()
    {
        return M_frequency_cos;
    }
    Real& phase_cos()
    {
        return M_phase_cos;
    }
    Real& coeff_exp_t_sin()
    {
        return M_coeff_exp_t_sin;
    }
    Real& mean_sin()
    {
        return M_mean_sin;
    }
    Real& amplitude_sin()
    {
        return M_amplitude_sin;
    }
    Real& frequency_sin()
    {
        return M_frequency_sin;
    }
    Real& phase_sin()
    {
        return M_phase_sin;
    }

    //@}

private:

    Real M_coeff_exp_t_cos;
    Real M_mean_cos;
    Real M_amplitude_cos;
    Real M_frequency_cos;
    Real M_phase_cos;
    Real M_coeff_exp_t_sin;
    Real M_mean_sin;
    Real M_amplitude_sin;
    Real M_frequency_sin;
    Real M_phase_sin;
};


//! Analytical_Solution - A particular case of Cos_min_Sin
/*!
 *  @author Lucia Mirabella
 *
 *  This analytical solution is in the form
 *
 *  U = ( Re(U) + i Im(U) ) * exp( Im(k) x ) *
 *  exp[ i ( \omega t - Re(k) x ) ] .
 *
 *  Its real part is
 *  Re(U) = [ Re(U) * exp( Im(k) x ) ] * cos( \omega t - Re(k) x ) +
 *  - [ Im(U) * exp( Im(k) x ) ] * sin( \omega t - Re(k) x ) +
 */
class Analytical_Solution : public Cos_min_Sin
{
public:

    /*!
      \brief The constructor

      \param[in] sol_amplitude_Re the real part of the wave amplitude
      \param[in] sol_amplitude_Im the imaginary part of the wave amplitude
      \param[in] kappa_Re the real part of the wave number \f$ k \f$
      \param[in] kappa_Im the imaginary part of the wave number \f$ k \f$
      \param[in] omega the wave time frequency \f$ \omega \f$ (supposed to be real)
    */
    explicit Analytical_Solution ( Real const& sol_amplitude_Re, Real const& sol_amplitude_Im,
                                   Real const& kappa_Re, Real const& kappa_Im,
                                   Real const& omega ) :
        Cos_min_Sin (0., 0., sol_amplitude_Re, omega, 0.,
                     0., 0., sol_amplitude_Im, omega, 0.),
        M_sol_amplitude_Re (sol_amplitude_Re),
        M_sol_amplitude_Im (sol_amplitude_Im),
        M_kappa_Re (kappa_Re),
        M_kappa_Im (kappa_Im)
    {}

    virtual ~Analytical_Solution() {}

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
      ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    Real operator() ( const Real& time )
    {
        return Cos_min_Sin::operator() (time);
    }

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
      ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    void update_x ( const Real& _x )
    {
        this->amplitude_cos() = M_sol_amplitude_Re * std::exp ( M_kappa_Im * _x );
        this->phase_cos() = - M_kappa_Re * _x;
        this->amplitude_sin() = M_sol_amplitude_Im * std::exp ( M_kappa_Im * _x );
        this->phase_sin() = - M_kappa_Re * _x;
    }

private:

    Real M_sol_amplitude_Re;
    Real M_sol_amplitude_Im;
    Real M_kappa_Re;
    Real M_kappa_Im;
};


//! PhysiologicalFlux - Base class for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class PhysiologicalFlux
{
public:

    explicit PhysiologicalFlux() :
        M_rampT(),
        M_time_step(),
        M_scale()
    {}

    virtual ~PhysiologicalFlux() {}

    Real operator() ( const Real& time );

private:

    Real M_rampT;
    Real M_time_step;
    Real M_scale;
};

Real
PhysiologicalFlux::operator() ( const Real& t )
{
    double time = t;
    int    numData = 100;
    double flux[101] = {  0.55457447187998,
                          0.55312181720914,
                          0.55299302643153,
                          0.55302818124406,
                          0.55317321131557,
                          0.55353364652152,
                          0.55374878962962,
                          0.55406829977313,
                          0.55585881887584,
                          0.55879983633299,
                          0.56387572718194,
                          0.57079488161935,
                          0.58817359652018,
                          0.61511673048210,
                          0.65025652077432,
                          0.79143597227331,
                          1.06959564837144,
                          1.33301653745648,
                          1.55916094914244,
                          1.69807981838757,
                          1.73941326221337,
                          1.69691789994768,
                          1.61546505715344,
                          1.51298277554169,
                          1.35636910481872,
                          1.18958468647846,
                          1.02381943290960,
                          0.89472539111700,
                          0.79994401665900,
                          0.74338513206540,
                          0.72984443940672,
                          0.74311502021666,
                          0.77155202899612,
                          0.80920592343822,
                          0.84528951107944,
                          0.87014986946118,
                          0.89016029509068,
                          0.89999452080415,
                          0.89511807514404,
                          0.87823357219620,
                          0.85663326250089,
                          0.82858260857792,
                          0.79671836139948,
                          0.75291106131325,
                          0.70711033067577,
                          0.65740190152584,
                          0.61066125251620,
                          0.56488426267136,
                          0.51713402331108,
                          0.46453623816504,
                          0.41513731950517,
                          0.47836113912116,
                          0.56452765299777,
                          0.62096051336166,
                          0.66202024502726,
                          0.69173157064612,
                          0.71835003021294,
                          0.74183126309604,
                          0.75295645424862,
                          0.75292455314576,
                          0.74514787317446,
                          0.72467414023271,
                          0.70473486985061,
                          0.68057326827129,
                          0.66194232245132,
                          0.64681425465222,
                          0.63714376881254,
                          0.62991615896879,
                          0.62662699778909,
                          0.62724985397200,
                          0.62674770176751,
                          0.62666043242736,
                          0.62617509524360,
                          0.62556258658310,
                          0.62581913341632,
                          0.62604032520998,
                          0.62582937168093,
                          0.62404471163034,
                          0.61923663804136,
                          0.61378537728592,
                          0.60976137345625,
                          0.60596975158344,
                          0.60144172708524,
                          0.59702451106965,
                          0.59319136754468,
                          0.58982107329344,
                          0.58718879911670,
                          0.58474181066352,
                          0.58208280034445,
                          0.57913818429409,
                          0.57588144074776,
                          0.57289019558638,
                          0.57076133371909,
                          0.56912637026578,
                          0.56776096894206,
                          0.56622327633393,
                          0.56376396210446,
                          0.56132345661888,
                          0.55914786876437,
                          0.55714215894326,
                          0.55457447187998
                       };

    double timescale = M_scale / numData;

    for (;;)
    {
        if (time < M_scale)
        {
            break;
        }
        time = time - M_scale;
    }

    int     ipos     = time / timescale;
    double t2 = timescale * (ipos + 1);

    double a = (flux[ipos + 1] - flux[ipos]) / timescale;
    double b = flux[ipos + 1] - a * t2;

    std::cout << "BC: Flux = " << time* a + b
              << " period = " << M_scale << " pos = " << ipos << std::endl;
    return time * a + b;
}

//! PressureRamp
/*!
 *  @author Lucia Mirabella
 */
class PressureRamp
{
public:
    explicit PressureRamp ( const Real& startT   = .001,
                            const Real& duration = 0.7,
                            const Real& endvalue = 106400 );

    virtual ~PressureRamp() {}

    Real operator() ( const Real& time );

private:

    Real M_startT;
    Real M_duration;
    Real M_endvalue;
};

PressureRamp::PressureRamp ( const Real& startT,
                             const Real& duration,
                             const Real& endvalue ) :
    M_startT    ( startT ),
    M_duration  ( duration ),
    M_endvalue  ( endvalue )
{}

Real
PressureRamp::operator() ( const Real& time )
{
    Real t = time;

    Int  numData = 80;
    Real pressure[81] = { 110170,
                          109540,
                          108930,
                          108320,
                          107710,
                          107120,
                          106530,
                          111130,
                          115440,
                          118690,
                          121460,
                          123940,
                          126350,
                          128890,
                          131510,
                          133980,
                          136200,
                          138330,
                          140350,
                          142290,
                          144360,
                          146130,
                          147530,
                          148780,
                          149740,
                          150320,
                          150470,
                          150250,
                          149750,
                          148990,
                          148220,
                          147210,
                          145940,
                          144960,
                          143750,
                          141980,
                          139900,
                          137260,
                          133970,
                          131670,
                          131320,
                          133150,
                          132710,
                          131570,
                          130280,
                          129750,
                          129330,
                          128910,
                          128360,
                          127680,
                          127000,
                          126410,
                          125920,
                          125480,
                          125040,
                          124560,
                          124050,
                          123530,
                          123000,
                          122440,
                          121840,
                          121220,
                          120580,
                          119950,
                          119330,
                          118710,
                          118100,
                          117470,
                          116840,
                          116200,
                          115560,
                          114920,
                          114280,
                          113650,
                          113020,
                          112400,
                          111790,
                          111200,
                          110620,
                          110060,
                          110170
                        };

    Real P = 0;
    if (t < 0)
    {
        P = t / M_startT * pressure[0];
    }
    else
    {
        Real timescale = M_duration / numData;

        for (;;)
        {
            if (t < M_duration)
            {
                break;
            }
            t -= M_duration;
        }

        Int  ipos = t / timescale;
        Real t2   = timescale * (ipos + 1);

        Real a = ( pressure[ipos + 1] - pressure[ipos] ) / timescale;
        Real b =   pressure[ipos + 1] - a * t2;

        P = t * a + b;

        std::cout << "BC: Pressure = " << P
                  << " period = " << M_duration << " pos = " << ipos << std::endl;
    }

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 6030 ) << "[PressureRamp::evaluate] imposed pressure = " << P << "\n";
#endif

    return P;
}

}

#endif
