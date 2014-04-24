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
 @brief Class for applying cardiac stimulus at Purkinje-muscle junctions

 @date 11-2013
 @author Toni Lassila <toni.lassila@epfl.ch>

 @last update 11-2013
 */

#ifndef STIMULUSPMJ_HPP_
#define STIMULUSPMJ_HPP_

#include <lifev/electrophysiology/stimulus/ElectroStimulus.hpp>

namespace LifeV
{

class StimulusPMJ : public ElectroStimulus
{

public:

    //! @name Type definitions
    //@{

    /*! @struct StimulusPMJ_Activation
     */
    struct StimulusPMJ_Activation
    {
        Real x;
        Real y;
        Real z;
        Real time;
        Real duration;
    };

    typedef std::vector<StimulusPMJ_Activation >  activationData_type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    StimulusPMJ();

    //! Destructor
    virtual ~StimulusPMJ() {}

    //@}

    //! @name Get Methods
    //@{
    inline activationData_type activationData( )
    {
        return M_activationData;
    }

    //@}

    //! @name Set Methods
    //@{
    void setPMJFromFile ( std::string fileName );
    void setPMJAddJunction ( Real x, Real y, Real z, Real time, Real duration );
    inline void setProblemFolder ( std::string problemFolder )
    {
        M_problemFolder=problemFolder;
    }

    inline void setRadius ( Real r )
    {
        ASSERT (r > 0, "Invalid radius value.");
        M_radius = r;
    }

    inline void setTotalCurrent ( Real I )
    {
        ASSERT (I >= 0, "Invalid current value.");
        M_totalCurrent = I;
    }

	void setParameters (list_Type&  list)
	{
		this->setRadius( list.get ("applied_current_radius", 0.2) );
		this->setTotalCurrent( list.get ("applied_total_current", 1.0) );
		this->setProblemFolder ( list.get ("PMJ_activation_folder", "./") );
		this->setPMJFromFile( M_problemFolder + list.get ("PMJ_activation_file", "PMJ_activation.txt") );
	}

    //@}

    //! @name Copy Methods
    //@{

    //@}

    //! @name Methods
    //@{
    Real appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i );
    void showMe ();
    //@}

private:

    activationData_type  M_activationData;
    Real                 M_radius;
    Real                 M_totalCurrent;
    std::string          M_problemFolder;

};

} // namespace LifeV

#endif /* STIMULUSPMJ_HPP_ */
