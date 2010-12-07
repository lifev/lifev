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
    @brief  Heart Functors for the Luo-Rudy Kinetics

    @date 04−2010
    @author

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */


#ifndef _HEARTFUNCTORS_H_
#define _HEARTFUNCTORS_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/EpetraMatrix.hpp>

namespace LifeV
{

class HeartFunctors
{
public:

    //! @name Public Types
    //@{

    typedef boost::function<Real ( Real const& x, Real const& y, Real const& z, Real const&, ID const& id , Real const&)> region_Type;
    typedef boost::function<Real ( Real const& x, Real const& y, Real const& z, Real const&, ID const& id)> region1_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    HeartFunctors();

    HeartFunctors( GetPot& dataFile );

    //@}


    //! @name Set Methods
    //@{

    /**
     * current volume source
     *
     * Define the stimulation current
     */
    Real Iapp ( const Real& x, const Real& y, const Real& z, const Real& t, const EntityFlag& id ) const;

    Real IappZygote(const double& t, const double& x, const double& y, const double& z, const ID& i, const EntityFlag& ref);

    Real stim ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   id) const;





    /**
     *
     * Reduces the conductivity in a sphere
     *
     */
    Real reduced_sigma_sphere( const Real& x, const Real& y, const Real& z, const Real& t, const ID&   id, const Real& sigma) const;

    /**
     *
     * Reduces the conductivity in a cylinder
     *
     */
    Real reduced_sigma_cylinder( const Real& x, const Real& y, const Real& z, const Real& t, const ID&   id, const Real& sigma ) const;

    Real reduced_sigma_box( const Real& x, const Real& y, const Real& z, const Real& t, const ID& id, const Real& sigma ) const;

    //@}


    //! @name Get Methods
    //@{

    Real initial_scalar( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );

    Real zero_scalar( const Real& t, const Real& x, const Real& y, const Real& z , const ID& id );


    region1_Type get_Iapp();

    region1_Type get_stim();

    const region_Type get_reduced_sigma_sphere();

    const region_Type get_reduced_sigma_cylinder();

    const region_Type get_reduced_sigma_box();

    const region1_Type get_initial_scalar();

    const region1_Type get_zero_scalar();

    //@}

    GetPot M_dataFile;

    boost::shared_ptr<Epetra_Comm>   M_comm;

    Int  M_stimulusSource;
    Real M_stimulusPeriod1;
    Real M_stimulusPeriod2;
    Real M_stimulusPeriod3;
    Real M_stimulusPeriod4;
    Real M_stimulusPeriod5;
    Real M_stimulusPeriod6;

    Real M_stimulusStart1;
    Real M_stimulusStop1;
    Real M_stimulusValue1;
    Real M_stimulusRadius1;
    KN<Real> M_stimulusCenter1;

    Real M_stimulusStart2;
    Real M_stimulusStop2;
    Real M_stimulusValue2;
    Real M_stimulusRadius2;
    KN<Real> M_stimulusCenter2;

    Real M_stimulusStart3;
    Real M_stimulusStop3;
    Real M_stimulusValue3;
    Real M_stimulusRadius3;
    KN<Real> M_stimulusCenter3;

    Real M_stimulusStart4;
    Real M_stimulusStop4;
    Real M_stimulusValue4;
    Real M_stimulusRadius4;
    KN<Real> M_stimulusCenter4;

    Real M_stimulusStart5;
    Real M_stimulusStop5;
    Real M_stimulusValue5;
    Real M_stimulusRadius5;
    KN<Real> M_stimulusCenter5;

    Real M_stimulusStart6;
    Real M_stimulusStop6;
    Real M_stimulusValue6;
    Real M_stimulusRadius6;
    KN<Real> M_stimulusCenter6;

    Real M_sphereX;
    Real M_sphereY;
    Real M_sphereZ;
    Real M_sphereR;
    KN<Real> M_sigmaReduction;

    Real M_cylinderX;
    Real M_cylinderY;
    Real M_cylinderZ;
    Real M_cylinderA;
    Real M_cylinderB;
    Real M_cylinderC;
    Real M_cylinderR;
    Real M_minimumCylinderX;
    Real M_maximumCylinderX;

    Real M_minimumBoxX;
    Real M_minimumBoxY;
    Real M_minimumBoxZ;

    Real M_maximumBoxX;
    Real M_maximumBoxY;
    Real M_maximumBoxZ;

    //parametres Iapp IappZygote REO
    Real M_timePeriod;
    Real M_appliedCurrentRightVentriculeAngle;
    Real M_appliedCurrentLeftVentriculeAngle;
    Real M_appliedCurrentStimulusTimeRightVentricule;
    Real M_appliedCurrentStimulusTimeLeftVentrivcule;
    Real M_ventricularFibrillation;
    Real M_nb_fibrillationSources;
    Real M_fibrillationSources;
    Real M_restPotential;

private:

    //! @name Unimplemented Methods
    //@{



    HeartFunctors( const HeartFunctors& heartFunctors );

    HeartFunctors& operator=( const HeartFunctors& heartFunctors );

    //@}

};
}

#endif
