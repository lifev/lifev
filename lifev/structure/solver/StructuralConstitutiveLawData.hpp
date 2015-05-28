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
 *  @brief DataElasticStructure - File containing a data container for solid problems with elastic structure
 *
 *  @version 1.0
 *  @date 01-10-2003
 *  @author M.A. Fernandez
 *
 *  @date 10-06-2010
 *  @author Cristiano Malossi
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef StructuralConstitutiveLawData_H
#define StructuralConstitutiveLawData_H

#include <string>
#include <iostream>
#include <map>

#include <boost/shared_ptr.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/StringUtility.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/TimeAdvanceData.hpp>

namespace LifeV
{

//! DataElasticStructure - Data container for solid problems with elastic structure
class StructuralConstitutiveLawData
{
public:

    //! @name Type definitions
    //@{

    typedef TimeData                               time_Type;
    typedef boost::shared_ptr< time_Type >         timePtr_Type;

    typedef TimeAdvanceData                        timeAdvance_Type;
    typedef boost::shared_ptr<timeAdvance_Type>    timeAdvancePtr_Type;

    typedef std::map<UInt, Real>                   materialContainer_Type;
    typedef materialContainer_Type::const_iterator materialContainerIterator_Type;
    typedef std::vector<UInt>                      vectorFlags_Type;
    typedef std::vector<Real>                      vectorParameters_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    StructuralConstitutiveLawData();

    //! Copy constructor
    /*!
     * @param StructuralConstitutiveLawData - StructuralConstitutiveLawData
     */
    StructuralConstitutiveLawData ( const StructuralConstitutiveLawData& structuralConstitutiveLawData );

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param StructuralConstitutiveLawData - StructuralConstitutiveLawData
     */
    StructuralConstitutiveLawData& operator= ( const StructuralConstitutiveLawData& structuralConstitutiveLawData );

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup ( const GetPot& dataFile, const std::string& section = "solid" );

    //! Display the values
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set methods
    //@{

    //! Set data time container
    /*!
     * @param TimeData shared_ptr to TimeData container
     */
    void setTimeData ( const timePtr_Type timeData )
    {
        M_time = timeData;
    }

    //! Set data time advance container
    /*!
     * @param timeAdvanceData shared_ptr to TimeAdvanceData container
     */
    void setTimeAdvanceData ( const timeAdvancePtr_Type timeAdvanceData )
    {
        M_timeAdvance = timeAdvanceData;
    }

    //! Set data external pressure for the external surface of the solid
    /*!
     * @param externalPressure external pressure value
     */
    void setExternalPressure ( const Real& externalPressure )
    {
        M_externalPressure = externalPressure;
    }

    //! Set density
    /*!
     * @param density solid density value
     */
    void setDensity ( const Real& density )
    {
        M_density = density;
    }

    //! Set thickness
    /*!
     * @param thickness solid thickness value
     */
    void setThickness ( const Real& thickness )
    {
        M_thickness = thickness;
    }

    //! Set poisson
    /*!
     * @param poisson solid poisson value
     * @param material material ID (1 by default)
     */
    void setPoisson ( const Real& poisson, const UInt& material )
    {
        M_materialsFlagSet = true;
        M_poisson[material] = poisson;
    }

    //! Set Young modulus
    /*!
     * @param Solid Young modulus value
     * @param material material ID (1 by default)
     */
    void setYoung ( const Real& young, const UInt& material )
    {
        M_materialsFlagSet = true;
        M_young[material] = young;
    }

    //! Set bulk modulus (nearly incompressible materials)
    /*!
     * @param bulk modulus value
     * @param material material ID (1 by default)
     */
    void setBulk ( const Real& bulk, const UInt& material )
    {
        M_materialsFlagSet = true;
        M_bulk[material] = bulk;
    }

    //! Set Alfa modulus (nearly incompressible materials)
    /*!
     * @param Alfa modulus value
     * @param material material ID (1 by default)
     */
    void setAlpha ( const Real& alpha, const UInt& material )
    {
        M_materialsFlagSet = true;
        M_alpha[material] = alpha;
    }

    //! Set Gamma (nearly imcompressible materials)
    /*!
     * @param Gamma modulus value
     * @param material material ID (1 by default)
     */
    void setGamma ( const Real& gamma, const UInt& material )
    {
        M_materialsFlagSet = true;
        M_gamma[material] = gamma;
    }

    //@}


    //! @name Get methods
    //@{

    //! Get data time container
    /*!
     * @return shared_ptr to TimeData container
     */
    timePtr_Type dataTime() const
    {
        return M_time;
    }

    //! Get data time advance container
    /*!
     * @return shared_ptr to TimeAdvanceData container
     */
    timeAdvancePtr_Type dataTimeAdvance() const
    {
        return M_timeAdvance;
    }

    //! Get the external pressure to be applied to the external surface of the solid
    /*!
     * @return the value of the external pressure
     */
    const Real& externalPressure() const
    {
        return M_externalPressure;
    }

    //! Get solid density
    /*!
     * @return Solid density
     */
    const Real& rho() const
    {
        return M_density;
    }

    //! Get solid thickness
    /*!
     * @return Solid thickness
     */
    const Real& thickness() const
    {
        return M_thickness;
    }

    //! Get solid thickness
    /*!
     * @return Solid thickness
     */
    const UInt maxSubIterationNumber() const
    {
        return M_maxSubIterationNumber;
    }
    //! Get solid thickness
    /*!
     * @return Solid thickness
     */
    const Real absoluteTolerance() const
    {
        return M_absoluteTolerance;
    }
    //! Get solid reltol newton
    /*!
     * @return Solid reltol newton
     */
    const Real relativeTolerance() const
    {
        return M_relativeTolerance;
    }
    //! Get solid etamax
    /*!
     * @return Solid etamax
     */
    const Real errorTolerance() const
    {
        return M_errorTolerance;
    }
    //! Get solid nonLinear Line Search ofr Newton
    /*!
     * @return nonLinearLineSearch
     */
    const UInt NonLinearLineSearch() const
    {
        return M_NonLinearLineSearch;
    }


    //! Get solid poisson coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid poisson coefficient
     */
    Real poisson ( const UInt& material ) const;

    //! Get solid young modulus
    /*!
     * @param material material ID (1 by default)
     * @return Solid young modulus
     */
    Real young ( const UInt& material ) const;

    //! Get solid first lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid first Lame coefficient
     */
    Real lambda ( const UInt& material ) const;

    //! Get solid second Lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid second Lame coefficient
     */
    Real mu ( const UInt& material ) const;

    //! Get bulk modulus (nearly incompressible materials)
    /*!
     * @param material material ID (1 by default)
     * @return bulk modulus
     */
    Real bulk ( const UInt& material = 1 ) const;

    //! Get alpha parameter (nearly incompressible materials)
    /*!
     * @param material material ID (1 by default)
     * @return alpha parameter (neraly incompressible materials)
     */
    Real alpha ( const UInt& material = 1 ) const;

    //! Get gamma parameter (nearly incompressible materials)
    /*!
     * @param material material ID (1 by default)
     * @return gamma parameter (nearly incompressible materials)
     */
    Real gamma ( const UInt& material = 1 ) const;

    //! Get FE order
    /*!
     * @return FE order
     */
    const std::string& order() const
    {
        return M_order;
    }

    //! Get verbose level
    /*!
     * @return verbose level
     */
    const UInt& verbose() const
    {
        return M_verbose;
    }

    //! Get isotropic law
    /*!
     * @return solid type isotropic
     */
    const std::string& solidTypeIsotropic()
    {
        return M_solidTypeIsotropic;
    }

    //! Get anisotropic law
    /*!
     * @return solid type anisotropic
     */

    const std::string& constitutiveLaw()
    {
        return M_constitutiveLaw;
    }

    const std::string& solidTypeAnisotropic()
    {
        return M_solidTypeAnisotropic;
    }

    const UInt numberFibersFamilies()
    {
        return M_numberFibers;
    }

    const Real ithStiffnessFibers( const UInt i )
    {
        return M_stiffnessParametersFibers[ i ];
    }

    const Real ithNonlinearityFibers( const UInt i )
    {
        return M_nonlinearityParametersFibers[ i ];
    }

    const Real ithCharacteristicStretch( const UInt i )
    {
        return M_characteristicStretch[ i ];
    }

    const Real ithDistributionFibers( const UInt i )
    {
        return M_distributionParametersFibers[ i ];
    }

    const Real smoothness( void )
    {
        return M_epsilon;
    }

    const std::string fiberActivation( void )
    {
        return M_fiberActivation;
    }

    const Real toleranceActivation( void )
    {
        return M_toleranceActivation;
    }



    //! Get law type
    /*!
     * @return law type
     */
    const std::string& lawType()
    {
        return M_lawType;
    }

    //! Get whether to use or not exact Jacobian
    /*!
     * @return true: if using exact Jacobian, false: otherwise
     */
    const bool& getUseExactJacobian() const
    {
        return M_useExactJacobian;
    }


    //! Get the vector of the set material_flags
    /*!
     * @return the vector of the material_flags set in the data file
     */
    const vectorFlags_Type& vectorFlags() const
    {
        return M_vectorMaterialFlags;
    }

    //@}

private:

    //! Data containers for time and mesh
    timePtr_Type           M_time;
    timeAdvancePtr_Type    M_timeAdvance;

    //! Physics
    Real                   M_density;
    Real                   M_thickness;
    Real                   M_externalPressure;

    bool                   M_materialsFlagSet;

    //! Young Modulus and Poisson ratio
    materialContainer_Type M_poisson;
    materialContainer_Type M_young;

    //! Bulk modulus k, alpha, gamma
    materialContainer_Type M_bulk;
    materialContainer_Type M_alpha;
    materialContainer_Type M_gamma;

    //! Space discretization
    std::string            M_order;

    //! Miscellaneous
    UInt                   M_verbose; // temporal output verbose

    std::string            M_solidTypeIsotropic;
    std::string            M_constitutiveLaw;
    std::string            M_solidTypeAnisotropic;
    UInt                   M_numberFibers;
    vectorParameters_Type  M_stiffnessParametersFibers;
    vectorParameters_Type  M_nonlinearityParametersFibers;
    vectorParameters_Type  M_characteristicStretch;
    vectorParameters_Type  M_distributionParametersFibers;
    Real                   M_epsilon;
    std::string            M_fiberActivation;
    Real                   M_toleranceActivation;
    std::string            M_lawType;
    bool                   M_useExactJacobian;

    vectorFlags_Type       M_vectorMaterialFlags;

    UInt                   M_maxSubIterationNumber;
    Real                   M_absoluteTolerance;
    Real                   M_relativeTolerance;
    Real                   M_errorTolerance;
    UInt                   M_NonLinearLineSearch;
};

} // end namespace LifeV

#endif // StructuralConstitutiveLawData_H
