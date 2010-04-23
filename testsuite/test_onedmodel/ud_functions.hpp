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
  \file ud_functions.hpp
  \author Tiziano Passerini
  \author Lucia Mirabella *queen*

  \date 04/2007

  \brief This file contains functions to impose boundary conditions for 1D tubes.

*/

#ifndef _ONED_FUNCTIONS_1D_
#define _ONED_FUNCTIONS_1D_
#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>

namespace LifeV {

typedef OneDimensionalModel_Physics                          Physics_Type;
typedef boost::shared_ptr< Physics_Type >                    Physics_PtrType;

typedef OneDimensionalModel_BCFunction::Flux_Type            Flux_Type;
typedef OneDimensionalModel_BCFunction::Flux_PtrType         Flux_PtrType;

typedef OneDimensionalModel_BCFunction::Source_Type          Source_Type;
typedef OneDimensionalModel_BCFunction::Source_PtrType       Source_PtrType;

typedef OneDimensionalModel_BCFunction::Mesh_Type            Mesh_Type;
typedef OneDimensionalModel_BCFunction::FESpace_Type         FESpace_Type;
typedef OneDimensionalModel_BCFunction::Vector_Type          Vector_Type;

//! Const - Base class for One Dimensional BC Functions
/*!
 *  @author Lucia Mirabella
 */
class Const : public OneDimensionalModel_BCFunction
{
public:
    Const(const Real val) : M_val(val) {}
    Real evaluate( const Real& /*time*/ ) {return M_val;}
    ~Const() {}
private:
    Real M_val;
};

/*!
  \brief A sinusoidal wave
*/
class Sin : public OneDimensionalModel_BCFunction
{
public:
    /*!
      \brief The constructor

      \param[in] mean the time average of the sinus \f$ B \f$
      \param[in] scale the sinus amplitude \f$ A \f$
      \param[in] period the sinus period \f$ T \f$
      \param[in] phase the sinus phase \f$ \phi \f$
    */
    Sin(const Real mean   = 0,
        const Real scale  = 10,
        const Real period = .01,
        const Real phase  = 0. ) :
            M_mean(mean), M_scale(scale), M_period(period), M_phase(phase) {}

    /*!
      \brief Compute the sinus at the specified time

      \param[in] time the time
      \return \f$ B + A sin( \frac{2 \pi t}{T} + \phi ) \f$
    */
    Real evaluate( const Real& time )
    {
        //         std::cout << "mean  " << M_mean << std::endl;
        //         std::cout << "scale  " << M_scale << std::endl;
        //         std::cout << "period " << M_period << std::endl;
        //         std::cout << "phase  " << M_phase << std::endl;


        if (time < M_period)
            {
                std::cout << time << " Flux BC = " << M_mean + M_scale*std::sin(M_phase+2*M_PI*time/M_period) << std::endl;
                return M_mean + M_scale*std::sin(M_phase+2*M_PI*time/M_period);
            }
        else
            {
                return 0.;
            }
//         std::cout << "BC imposed at time  " << time << " = " << M_mean + M_scale*std::sin(M_phase + 2*M_PI*time/M_period) << std::endl;
//         return M_mean + M_scale*std::sin(M_phase + 2*M_PI*time/M_period);
    };

    ~Sin() {}
private:
    Real M_mean;
    Real M_scale;
    Real M_period;
    Real M_phase;
};


/*!
  \brief A superimposition of a sinusoidal and a cosinusoidal waves,
  whose amplitude is damped by exponential terms
*/
class Cos_min_Sin : public OneDimensionalModel_BCFunction
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
    Cos_min_Sin(const Real coeff_exp_t_cos = 0, const Real mean_cos = 0,
                const Real amplitude_cos = 10, const Real frequency_cos = 8.*atan(1.),
                const Real phase_cos = 0.,
                const Real coeff_exp_t_sin = 0, const Real mean_sin = 0,
                const Real amplitude_sin = 10, const Real frequency_sin = 8.*atan(1.),
                const Real phase_sin = 0. ) :
            M_coeff_exp_t_cos(coeff_exp_t_cos), M_mean_cos(mean_cos),
            M_amplitude_cos(amplitude_cos), M_frequency_cos(frequency_cos), M_phase_cos(phase_cos),
            M_coeff_exp_t_sin(coeff_exp_t_sin), M_mean_sin(mean_sin),
            M_amplitude_sin(amplitude_sin), M_frequency_sin(frequency_sin), M_phase_sin(phase_sin) {}

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
      ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    Real evaluate( const Real& time )
    {
        Real result = M_mean_cos
            + M_amplitude_cos * std::cos(M_phase_cos+time*M_frequency_cos)
            - ( M_mean_sin
                + M_amplitude_sin * std::sin(M_phase_sin+time*M_frequency_sin) );

        return result;
    };

    /*! \name Getters
     */
    //@{
    Real& coeff_exp_t_cos(){ return M_coeff_exp_t_cos;}
    Real& mean_cos(){ return M_mean_cos;}
    Real& amplitude_cos(){ return M_amplitude_cos;}
    Real& frequency_cos(){ return M_frequency_cos;}
    Real& phase_cos(){ return M_phase_cos;}
    Real& coeff_exp_t_sin(){ return M_coeff_exp_t_sin;}
    Real& mean_sin(){ return M_mean_sin;}
    Real& amplitude_sin(){ return M_amplitude_sin;}
    Real& frequency_sin(){ return M_frequency_sin;}
    Real& phase_sin(){ return M_phase_sin;}
    //@}

    ~Cos_min_Sin() {}
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


/*!
  \brief A particular case of Cos_min_Sin

  This analytical solution is in the form

  U = ( Re(U) + i Im(U) ) * exp( Im(k) x ) *
  exp[ i ( \omega t - Re(k) x ) ] .

  Its real part is
  Re(U) = [ Re(U) * exp( Im(k) x ) ] * cos( \omega t - Re(k) x ) +
  - [ Im(U) * exp( Im(k) x ) ] * sin( \omega t - Re(k) x ) +
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
    Analytical_Solution( Real const& sol_amplitude_Re, Real const& sol_amplitude_Im,
                         Real const& kappa_Re, Real const& kappa_Im,
                         Real const& omega ) :
            Cos_min_Sin(0., 0., sol_amplitude_Re, omega, 0.,
                        0., 0., sol_amplitude_Im, omega, 0.),
            M_sol_amplitude_Re(sol_amplitude_Re),
            M_sol_amplitude_Im(sol_amplitude_Im),
            M_kappa_Re(kappa_Re),
            M_kappa_Im(kappa_Im)
    {}

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
      ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    Real evaluate( const Real& time )
    {
        return Cos_min_Sin::evaluate(time);
    };

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
      ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    void update_x( const Real& _x )
    {
        this->amplitude_cos() = M_sol_amplitude_Re * std::exp( M_kappa_Im * _x );
        this->phase_cos() = - M_kappa_Re * _x;
        this->amplitude_sin() = M_sol_amplitude_Im * std::exp( M_kappa_Im * _x );
        this->phase_sin() = - M_kappa_Re * _x;
    };

    ~Analytical_Solution() {}

private:
    Real M_sol_amplitude_Re;
    Real M_sol_amplitude_Im;
    Real M_kappa_Re;
    Real M_kappa_Im;
};


class PhysiologicalFlux : public OneDimensionalModel_BCFunction
{
public:
    PhysiologicalFlux( GetPot const& data_file );

    Real evaluate( const Real& time );

    ~PhysiologicalFlux() {}
private:
	Real M_rampT;
	Real M_time_step;
    Real M_scale;
};






Real
PhysiologicalFlux::evaluate( const Real& t )
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
                            0.55457447187998 };


      double timescale = M_scale/numData;

      for (;;)
      {
          if (time < M_scale) break;
          time = time - M_scale;
      }

      int     ipos     = time/timescale;
      double t2 =timescale*(ipos + 1);

      double a = (flux[ipos + 1] - flux[ipos])/timescale;
      double b = flux[ipos + 1] - a*t2;

      double slope = ipos -  time;
      std::cout << "BC: Flux = " << time*a + b
                << " period = " << M_scale << " pos = " << ipos << std::endl;
      return time*a + b;
}





class PressureRamp : public Compatibility
{
public:
    PressureRamp( const FESpace_Type&                            fespace,
                  const Flux_PtrType                             fluxFun,
                  const Source_PtrType                           sourceFun,
                  const std::vector<Vector_Type>& U_thistime,
                  const Real&                                    dt,
                  const std::string&                             border,
                  const std::string&                             var,
                  const Physics_PtrType                          onedparam,
                  const Real&                                    startT     = .001,
                  const Real&                                    duration   = 0.7,
                  const Real&                                    endvalue   = 106400 );
    // const RegionMesh& mesh,
    //                  const FLUX& fluxFun,
    //                  const SOURCE& sourceFun,
    //                  const std::vector<DataOneDModel::vector_type>& U_thistime,
    //                  const Real& dt,
    //                  const std::string& border,
    //                  const std::string & var,
    //                  const PARAM& onedparam,
    //                  const Real& startT = 0.,
    //                  const Real& duration = 0.05,
    //                  const Real& endvalue = 106400 );

    Real evaluate( const Real& time );

    ~PressureRamp() {}
private:
    //! Reference to the solver parameters
    Physics_PtrType M_onedparam;

    Real M_startT;
    Real M_duration;
    Real M_endvalue;
    Real M_dt;
};


PressureRamp::PressureRamp (const FESpace_Type&          fespace,
                                                 const Flux_PtrType                             fluxFun,
                                                 const Source_PtrType                           sourceFun,
                                                 const std::vector<Vector_Type>& U_thistime,
                                                 const Real&                                    dt,
                                                 const std::string&                             border,
                                                 const std::string&                             var,
                                                 const Physics_PtrType                          onedparam,
                                                 const Real&                                    startT,
                                                 const Real&                                    duration,
                                                 const Real&                                    endvalue):
        Compatibility             (fespace, fluxFun, sourceFun, U_thistime,/* W_thistime,*/ dt, border, var),
        M_onedparam               (onedparam),
        M_startT                  (startT),
        M_duration                (duration),
        M_endvalue                (endvalue),
        M_dt                      (dt)
{}



Real PressureRamp::evaluate( const Real& time )
{
    Real W_out, result(0.);

//     Real _P = ( time < (M_startT + M_duration) ) ?
//         ( ( time - M_startT ) / M_duration ) : 1;
    //Real _P = M_endvalue;
    Real t = time;

    int    numData = 80;
    double pressure[81] =
        {
            110170,
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

    Real _P = 0;
//     Real _P = ( time < (M_startT + M_duration) ) ?
//         ( ( time - M_startT ) / M_duration ) : 1
    if (t < 0)
        _P = t/M_startT*pressure[0];
    else
    {

        double timescale = M_duration/numData;

        for (;;)
        {
            if (t < M_duration) break;
            t = t - M_duration;
        }

        int     ipos     = t/timescale;
        double t2        = timescale*(ipos + 1);

        double a = (pressure[ipos + 1] - pressure[ipos])/timescale;
        double b = pressure[ipos + 1] - a*t2;

        _P = t*a + b;

        double slope = ipos -  t;
        std::cout << "BC: Pressure = " << _P
                  << " period = " << M_duration << " pos = " << ipos << std::endl;

        //Real _P = 10000;
    }

    Debug( 6030 ) << "[PressureRamp::evaluate] imposed pressure = " << _P << "\n";
    //    Debug( 6030 ) << "[PressureRamp::evaluate] target pressure = " << M_endvalue << "\n";

    switch( this->M_oneDBCFunctionsMapStringValues[this->M_var] )
        {
        case OneDBCW1:
            W_out = this->extrapolate_W( OneDBCW2 );
            result = M_onedparam->W_from_P( _P, W_out, 2, this->M_boundaryDof);
            break;
        case OneDBCW2:
            W_out = this->extrapolate_W( OneDBCW1 );
            result = M_onedparam->W_from_P( _P, W_out, 1, this->M_boundaryDof);
            break;
        default:
            std::cout << "\n[PressureRamp::evaluate] incorrect variable identifier: " << this->M_var << std::endl;
        }

    Debug( 6030 ) << "[PressureRamp::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

    return result;

}




class Resi : public Compatibility
{
public:
    Resi(  const Real & resistance, // const GetPot& data_file,
           const Physics_PtrType onedparam,
           const FESpace_Type& FESpace,
           const Flux_PtrType fluxFun,
           const Source_PtrType sourceFun,
           const std::vector<Vector_Type>& U_thistime,
           const Real& dt,
           const std::string& border,
           const std::string & var,
           const bool& absorbing = false);

    Real evaluate( const Real& time );

    ~Resi() {}
private:
    //! Resi value
    Real M_resistance;
    //! Reference to the solver parameters
    Physics_PtrType M_onedparam;
    //! to impose absorbing conditions
    bool M_absorbing;
};


Resi::Resi(const Real &                                   resistance,
                                const Physics_PtrType                          onedparam, // const GetPot& data_file,
                                const FESpace_Type&          fespace,
                                const Flux_PtrType                                    fluxFun,
                                const Source_PtrType                                  sourceFun,
                                const std::vector<Vector_Type>& U_thistime,
                                const Real&                                    dt,
                                const std::string&                             border,
                                const std::string&                             var,
                                const bool&                                    absorbing ):
        Compatibility( fespace, fluxFun, sourceFun, U_thistime, /*W_thistime,*/ dt, border, var),
        M_resistance(resistance),
        M_onedparam( onedparam ),
        M_absorbing( absorbing )
{
    Debug( 6030 ) << "[Resi::Resi] resistance = " << M_resistance << "\n";
}



Real Resi::evaluate( const Real& /*time*/ )
{
    //! Coefficients
    Real W_out(0.), result;
    Real a1, a2, a11, a22, b1, b2, c1, c2;

    this->update_U_boundary();
    this->update_U_internalBd();

    Debug( 6030 ) << "[Resi::Resi] at node " << this->M_boundaryDof
                  << ", A  = " << this->M_U_boundary[0] << "( " << this->M_U_thistime[0][1] << " ) "
                  << ", Q  = " << this->M_U_boundary[1]
                  << ", W1 = " << this->M_W_boundary[0]
                  << ", W2 = " << this->M_W_boundary[1]
                  << "\n";

    this->computeEigenValuesVectors();

    a1 = M_onedparam->pressure(this->M_U_boundary[0], this->M_boundaryDof - 1); // pressure at previous time step

    a2 = this->M_U_boundary[1]; // flux at previous time step

    b1 = M_onedparam->pressure_WDiff( this->M_W_boundary[0], this->M_W_boundary[1], 1, this->M_boundaryDof - 1);  // dP / dW1

    b2 = this->M_U_boundary[0] / 2; // dQ / dW1

    c1 = M_onedparam->pressure_WDiff( this->M_W_boundary[0], this->M_W_boundary[1], 2, this->M_boundaryDof - 1);  // dP / dW2

    c2 = b2; // dQ / dW2

    Debug( 6030 ) << "[Resi::evaluate] P(A) = " << a1 << "\n";
    Debug( 6030 ) << "[Resi::evaluate] P(W1,W2) = "
                  << M_onedparam->pressure_W(this->M_W_boundary[0], this->M_W_boundary[1], this->M_boundaryDof - 1) << "\n";

    a11 = a1 - b1*this->M_W_boundary[0] - c1*this->M_W_boundary[1];
    a22 = a2 - b2*this->M_W_boundary[0] - c2*this->M_W_boundary[1];

    switch( this->M_oneDBCFunctionsMapStringValues[this->M_var] )
        {
        case OneDBCW1:
            W_out = this->extrapolate_L_dot_U(this->M_eigval2, this->M_left_eigvec2)
                - dot( this->M_left_eigvec2, this->M_U_boundary ) + this->M_W_boundary[1];

            break;
        case OneDBCW2:
            W_out = this->extrapolate_L_dot_U(this->M_eigval1, this->M_left_eigvec1)
                - dot( this->M_left_eigvec1, this->M_U_boundary ) + this->M_W_boundary[0];

            break;
        default:
            std::cout << "\n[Resi::evaluate] incorrect variable identifier: " << this->M_border << std::endl;
        }

    Debug( 6030 ) << "[Resi::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

    if( M_absorbing ) {
        M_resistance = b1 / b2;
        Debug( 6030 ) << "[Resi::evaluate] imposing absorbing condition, R = " << M_resistance << "\n";
    }

    result = W_out * ((b2*M_resistance-b1)/(c1-c2*M_resistance))
        + ((a22*M_resistance-a11)/(c1-c2*M_resistance));

    Debug( 6030 ) << "[Resi::evaluate] a1 = " << a1 << "\n";
    Debug( 6030 ) << "[Resi::evaluate] b1 = " << b1 << "\n";
    Debug( 6030 ) << "[Resi::evaluate] c1 = " << c1 << "\n";
    Debug( 6030 ) << "[Resi::evaluate] a2 = " << a2 << "\n";
    Debug( 6030 ) << "[Resi::evaluate] b2 = " << b2 << "\n";
    Debug( 6030 ) << "[Resi::evaluate] c2 = " << c2 << "\n";

    return result;
}
}



#endif
