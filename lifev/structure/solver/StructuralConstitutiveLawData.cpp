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
 *  @version 1.18
 *  @date 10-06-2010
 *  @author Cristiano Malossi
 *
 *  @version 1.19
 *  @author Gilles Fourestey
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

namespace LifeV
{

//=====================================================
// Constructors
//=====================================================
StructuralConstitutiveLawData::StructuralConstitutiveLawData() :
    M_time                             ( ),
    M_timeAdvance                      ( ),
    M_density                          ( ),
    M_thickness                        ( ),
    M_externalPressure                 ( ),
    M_materialsFlagSet                 ( false ),
    M_poisson                          ( ),
    M_young                            ( ),
    M_bulk                             ( ),
    M_alpha                            ( ),
    M_gamma                            ( ),
    M_order                            ( ),
    M_verbose                          ( ),
    M_solidTypeIsotropic               ( ),
    M_constitutiveLaw                  ( ),
    M_solidTypeAnisotropic             ( ),
    M_numberFibers                     ( 0 ),
    M_stiffnessParametersFibers        ( ),
    M_nonlinearityParametersFibers     ( ),
    M_characteristicStretch            ( ),
    M_distributionParametersFibers     ( ),
    M_epsilon                          ( 0 ),
    M_fiberActivation                  ( ),
    M_toleranceActivation              ( ),
    M_lawType                          ( ),
    M_useExactJacobian                 ( false ),
    M_vectorMaterialFlags              ( ),
    M_maxSubIterationNumber            ( ),
    M_absoluteTolerance                ( ),
    M_relativeTolerance                ( ),
    M_errorTolerance                   ( ),
    M_NonLinearLineSearch              ( ),
    M_thinLayer						   ( false ),
    M_thinLayerThickness			   ( ),
    M_thinLayerDensity				   ( ),
    M_thinLayerLameI				   ( ),
    M_thinLayerLameII				   ( ),
    M_interfaceFlag					   ( ),
    M_LameThickByFunctors			   ( false )
{
}

StructuralConstitutiveLawData::StructuralConstitutiveLawData ( const StructuralConstitutiveLawData& structuralConstitutiveLawData ) :
    M_time                             ( structuralConstitutiveLawData.M_time ),
    M_timeAdvance                      ( structuralConstitutiveLawData.M_timeAdvance ),
    M_density                          ( structuralConstitutiveLawData.M_density ),
    M_thickness                        ( structuralConstitutiveLawData.M_thickness ),
    M_externalPressure                 ( structuralConstitutiveLawData.M_externalPressure ),
    M_materialsFlagSet                 ( structuralConstitutiveLawData.M_materialsFlagSet ),
    M_poisson                          ( structuralConstitutiveLawData.M_poisson ),
    M_young                            ( structuralConstitutiveLawData.M_young ),
    M_bulk                             ( structuralConstitutiveLawData.M_bulk ),
    M_alpha                            ( structuralConstitutiveLawData.M_alpha ),
    M_gamma                            ( structuralConstitutiveLawData.M_gamma ),
    M_order                            ( structuralConstitutiveLawData.M_order ),
    M_verbose                          ( structuralConstitutiveLawData.M_verbose ),
    M_solidTypeIsotropic               ( structuralConstitutiveLawData.M_solidTypeIsotropic ),
    M_constitutiveLaw                  ( structuralConstitutiveLawData.M_constitutiveLaw ),
    M_solidTypeAnisotropic             ( structuralConstitutiveLawData.M_solidTypeAnisotropic ),
    M_numberFibers                     ( structuralConstitutiveLawData.M_numberFibers ),
    M_stiffnessParametersFibers        ( structuralConstitutiveLawData.M_stiffnessParametersFibers ),
    M_nonlinearityParametersFibers     ( structuralConstitutiveLawData.M_nonlinearityParametersFibers ),
    M_characteristicStretch            ( structuralConstitutiveLawData.M_characteristicStretch ),
    M_distributionParametersFibers     ( structuralConstitutiveLawData.M_distributionParametersFibers ),
    M_epsilon                          ( structuralConstitutiveLawData.M_epsilon ),
    M_fiberActivation                  ( structuralConstitutiveLawData.M_fiberActivation ),
    M_toleranceActivation              ( structuralConstitutiveLawData.M_toleranceActivation ),
    M_lawType                          ( structuralConstitutiveLawData.M_lawType ),
    M_useExactJacobian                 ( structuralConstitutiveLawData.M_useExactJacobian ),
    M_vectorMaterialFlags              ( structuralConstitutiveLawData.M_vectorMaterialFlags ),
    M_maxSubIterationNumber            ( structuralConstitutiveLawData.M_maxSubIterationNumber ),
    M_absoluteTolerance                ( structuralConstitutiveLawData.M_absoluteTolerance ),
    M_relativeTolerance                ( structuralConstitutiveLawData.M_relativeTolerance ),
    M_errorTolerance                   ( structuralConstitutiveLawData.M_errorTolerance ),
    M_NonLinearLineSearch              ( structuralConstitutiveLawData.M_NonLinearLineSearch ),
    M_thinLayer						   ( structuralConstitutiveLawData.M_thinLayer ),
    M_thinLayerThickness			   ( structuralConstitutiveLawData.M_thinLayerThickness ),
    M_thinLayerDensity				   ( structuralConstitutiveLawData.M_thinLayerDensity ),
    M_thinLayerLameI				   ( structuralConstitutiveLawData.M_thinLayerLameI ),
    M_thinLayerLameII				   ( structuralConstitutiveLawData.M_thinLayerLameII ),
    M_interfaceFlag                    ( structuralConstitutiveLawData.M_interfaceFlag),
    M_LameThickByFunctors              ( structuralConstitutiveLawData.M_LameThickByFunctors)
{
}

// ===================================================
// Operators
// ===================================================
StructuralConstitutiveLawData&
StructuralConstitutiveLawData::operator= ( const StructuralConstitutiveLawData& structuralConstitutiveLawData )
{
    if ( this != &structuralConstitutiveLawData )
    {
        M_time                             = structuralConstitutiveLawData.M_time;
        M_timeAdvance                      = structuralConstitutiveLawData.M_timeAdvance;
        M_density                          = structuralConstitutiveLawData.M_density;
        M_thickness                        = structuralConstitutiveLawData.M_thickness;
        M_externalPressure                 = structuralConstitutiveLawData.M_externalPressure;
        M_materialsFlagSet                 = structuralConstitutiveLawData.M_materialsFlagSet;
        M_poisson                          = structuralConstitutiveLawData.M_poisson;
        M_young                            = structuralConstitutiveLawData.M_young;
        M_bulk                             = structuralConstitutiveLawData.M_bulk;
        M_alpha                            = structuralConstitutiveLawData.M_alpha;
        M_gamma                            = structuralConstitutiveLawData.M_gamma;
        M_order                            = structuralConstitutiveLawData.M_order;
        M_verbose                          = structuralConstitutiveLawData.M_verbose;
        M_solidTypeIsotropic               = structuralConstitutiveLawData.M_solidTypeIsotropic;
        M_constitutiveLaw                  = structuralConstitutiveLawData.M_constitutiveLaw;
        M_solidTypeAnisotropic             = structuralConstitutiveLawData.M_solidTypeAnisotropic;
        M_numberFibers                     = structuralConstitutiveLawData.M_numberFibers;
        M_stiffnessParametersFibers        = structuralConstitutiveLawData.M_stiffnessParametersFibers;
        M_nonlinearityParametersFibers     = structuralConstitutiveLawData.M_nonlinearityParametersFibers;
        M_characteristicStretch            = structuralConstitutiveLawData.M_characteristicStretch;
        M_distributionParametersFibers     = structuralConstitutiveLawData.M_distributionParametersFibers;
        M_epsilon                          = structuralConstitutiveLawData.M_epsilon;
        M_fiberActivation                  = structuralConstitutiveLawData.M_fiberActivation;
        M_toleranceActivation              = structuralConstitutiveLawData.M_toleranceActivation;
        M_lawType                          = structuralConstitutiveLawData.M_lawType;
        M_useExactJacobian                 = structuralConstitutiveLawData.M_useExactJacobian;
        M_vectorMaterialFlags              = structuralConstitutiveLawData.M_vectorMaterialFlags;
        M_maxSubIterationNumber            = structuralConstitutiveLawData.M_maxSubIterationNumber;
        M_absoluteTolerance                = structuralConstitutiveLawData.M_absoluteTolerance;
        M_relativeTolerance                = structuralConstitutiveLawData.M_relativeTolerance;
        M_errorTolerance                   = structuralConstitutiveLawData.M_errorTolerance;
        M_NonLinearLineSearch              = structuralConstitutiveLawData.M_NonLinearLineSearch;
        M_thinLayer						   = structuralConstitutiveLawData.M_thinLayer;
        M_thinLayerThickness			   = structuralConstitutiveLawData.M_thinLayerThickness;
        M_thinLayerDensity				   = structuralConstitutiveLawData.M_thinLayerDensity;
        M_thinLayerLameI				   = structuralConstitutiveLawData.M_thinLayerLameI;
        M_thinLayerLameII				   = structuralConstitutiveLawData.M_thinLayerLameII;
        M_interfaceFlag					   = structuralConstitutiveLawData.M_interfaceFlag;
        M_LameThickByFunctors			   = structuralConstitutiveLawData.M_LameThickByFunctors;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
StructuralConstitutiveLawData::setup ( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
    {
        M_time.reset ( new time_Type ( dataFile, section + "/time_discretization" ) );
    }

    if ( !M_timeAdvance.get() )
    {
        M_timeAdvance.reset ( new timeAdvance_Type ( dataFile, section + "/time_discretization" ) );
    }

    // physics
    M_solidTypeIsotropic = dataFile ( ( section + "/model/solidTypeIsotropic" ).data(), "NO_DEFAULT_SOLID_TYPE_ISOTROPIC" );

    // Reading the type of anisotropic part and the number of fibers
    M_constitutiveLaw = dataFile ( ( section + "/model/constitutiveLaw" ).data(), "isotropic" );

    if(  !M_constitutiveLaw.compare("anisotropic") ) // anisotropic laws
    {
        M_solidTypeAnisotropic = dataFile ( ( section + "/model/solidTypeAnisotropic" ).data(), "NO_DEFAULT_SOLID_TYPE_ANISOTROPIC" );
        M_numberFibers = dataFile ( ( section + "/model/fibers/numberFamilies" ).data(), 0 );

        ASSERT( M_numberFibers, " The number of fibers of the anisotropic law has to be different from 0!" );
    }

    // The check can be done on the isotropic part since the anisotropic is for sure nonlinear
    if( !M_solidTypeIsotropic.compare("linearVenantKirchhoff") )
      {
           M_lawType = "linear";
      }
    else
      {
          M_lawType = "nonlinear";
      }

    if(  !M_constitutiveLaw.compare("anisotropic") ) // anisotropic laws => no LE
      {
          ASSERT ( M_lawType.compare ("linear"), "The Linear Elastic law cannot be used with anisotropic laws");
      }

    M_externalPressure = dataFile ( ( section + "/physics/externalPressure" ).data(), 0. );
    M_density   = dataFile ( ( section + "/physics/density"   ).data(), 1. );
    M_thickness = dataFile ( ( section + "/physics/thickness" ).data(), 0.1 );


    UInt materialsNumber = dataFile.vector_variable_size ( ( section + "/physics/material_flag" ).data() );

    if( !M_constitutiveLaw.compare("anisotropic") )
      {
          UInt numberOfStiffnesses = dataFile.vector_variable_size ( ( section + "/model/fibers/stiffness" ).data() );
          UInt numberOfNonlinearities = dataFile.vector_variable_size ( ( section + "/model/fibers/nonlinearity" ).data() );


          ASSERT( M_numberFibers  , " The number of fiber families is equal to zero, change the variable constitutiveLaw from anisotropic to isotropic " );
          ASSERT( ( M_numberFibers == numberOfStiffnesses ) && ( M_numberFibers == numberOfNonlinearities ), " Inconsistency in the set up of the fiber parameters" );

          if( !M_solidTypeAnisotropic.compare("multimechanism") ||
              !M_solidTypeAnisotropic.compare("holzapfelGeneralized") )
          {
              UInt numberStretches = dataFile.vector_variable_size ( ( section + "/model/fibers/stretch" ).data() );
              ASSERT( ( M_numberFibers == numberOfStiffnesses ) &&
                      ( M_numberFibers == numberOfNonlinearities ) &&
                      ( M_numberFibers == numberStretches ) , " Inconsistency in the set up of the fiber parameters" );
          }

          if( !M_solidTypeAnisotropic.compare("distributedHolzapfel") )
          {
              UInt numberOfDistribution   = dataFile.vector_variable_size ( ( section + "/model/fibers/distribution" ).data() );
              ASSERT( ( M_numberFibers == numberOfStiffnesses ) && ( M_numberFibers == numberOfNonlinearities )
                      && ( M_numberFibers == numberOfDistribution ), " Inconsistency in the set up of the fiber parameters" );
          }
      }

    // Reading the material for isotropic laws
    ASSERT ( materialsNumber, "Set the materrial_flag variable in [solid]/physics");

    if ( materialsNumber == 0 )
    {
        // If no material is specified in the data file the code assume that there is just one material
        // and by default it is memorized with ID 1. Getters and Setters have been designed to deal with thic choice.
        M_materialsFlagSet = false;
        M_vectorMaterialFlags.resize (1);

        M_vectorMaterialFlags[0] = 1;
        M_young[1]   = dataFile ( ( section + "/model/young"   ).data(), 0. );
        M_poisson[1] = dataFile ( ( section + "/model/poisson" ).data(), 0. );

        M_bulk[1] = dataFile ( ( section + "/model/bulk"   ).data(), 0. );
        M_alpha[1] = dataFile ( ( section + "/model/alpha" ).data(), 0. );
        M_gamma[1] = dataFile ( ( section + "/model/gamma" ).data(), 0. );
    }
    else
    {
        M_materialsFlagSet = true;

        // These asserts are commented because in some cases we need to initialize the materials with default Young and Poisson, setting the correct values a posteriori.
        //ASSERT( M_materialsFlagSet == dataFile.vector_variable_size( ( section + "/physics/young"   ).data()),  "!!! ERROR: Inconsistent size for Young Modulus !!!");
        //ASSERT( M_materialsFlagSet == dataFile.vector_variable_size( ( section + "/physics/poisson" ).data() ), "!!! ERROR: Inconsistent size for Poisson Coef. !!!");

        UInt material (0);
        for ( UInt i (0) ; i < materialsNumber ; ++i )
        {

            M_vectorMaterialFlags.resize ( materialsNumber );
            material            = dataFile ( ( section + "/physics/material_flag" ).data(), 0., i );

            M_vectorMaterialFlags[i] = material;
            M_young[material]   = dataFile ( ( section + "/model/young"         ).data(), 0., i );
            M_poisson[material] = dataFile ( ( section + "/model/poisson"       ).data(), 0., i );

            M_bulk[material] = dataFile ( ( section + "/model/bulk"         ).data(), 0.0, i );
            M_alpha[material] = dataFile ( ( section + "/model/alpha"       ).data(), 0.0, i );
            M_gamma[material] = dataFile ( ( section + "/model/gamma"       ).data(), 0.0, i );
        }
    }

    if( !M_constitutiveLaw.compare("anisotropic") )
      {
          M_stiffnessParametersFibers .resize ( M_numberFibers  );
          M_nonlinearityParametersFibers .resize ( M_numberFibers  );
          M_distributionParametersFibers .resize ( M_numberFibers  );

          if( !M_solidTypeAnisotropic.compare("multimechanism") ||
              !M_solidTypeAnisotropic.compare("holzapfelGeneralized") )
          {
              M_characteristicStretch .resize ( M_numberFibers  );
          }

          for ( UInt i (0) ; i < M_numberFibers ; ++i )
          {
              M_stiffnessParametersFibers[ i ]      = dataFile ( ( section + "/model/fibers/stiffness"    ).data(), 0., i );
              M_nonlinearityParametersFibers[ i ]   = dataFile ( ( section + "/model/fibers/nonlinearity" ).data(), 0., i );

              if( !M_solidTypeAnisotropic.compare("multimechanism") ||
                  !M_solidTypeAnisotropic.compare("holzapfelGeneralized") )
              {
                  M_characteristicStretch[ i ]   = dataFile ( ( section + "/model/fibers/stretch" ).data(), 0., i );
              }

              if( !M_solidTypeAnisotropic.compare("distributedHolzapfel") )
              {
                  M_distributionParametersFibers[ i ]   = dataFile ( ( section + "/model/fibers/distribution" ).data(), 0., i );
              }
          }
          M_epsilon = dataFile ( ( section + "/model/fibers/smoothness"   ).data(), 0. );
          M_fiberActivation = dataFile ( ( section + "/model/fiberActivation" ).data(), "implicit" );
          M_toleranceActivation = dataFile ( ( section + "/model/fibers/tolActivation"   ).data(), 0.001 );
      }

    // space_discretization
    M_order            = dataFile ( ( section + "/space_discretization/order" ).data(), "P1" );

    // miscellaneous
    M_verbose          = dataFile ( ( section + "/miscellaneous/verbose" ).data(), 0 );
    M_useExactJacobian = dataFile ( ( section + "/useExactJacobian"      ).data(), false );

    // Problem - Non Linear Richardson Parameters
    M_maxSubIterationNumber = dataFile ( ( section + "/newton/maxSubIter" ).data(), 300 );
    M_absoluteTolerance = dataFile ( ( section + "/newton/abstol" ).data(), 1.e-07 );
    M_relativeTolerance = dataFile ( ( section + "/newton/reltol" ).data(), 1.e-07 );
    M_errorTolerance = dataFile ( ( section + "/newton/etamax" ).data(), 1.e-03 );
    M_NonLinearLineSearch = static_cast<Int> ( dataFile ( ( section + "/newton/NonLinearLineSearch" ).data(), 0 ) );

    M_LameThickByFunctors = dataFile ( ( section + "/physics/LameThickByFunctor" ).data(), false );

    // thin layer parameters
    M_thinLayer = dataFile ( ( section + "/physics/use_thin" ).data(), false );

    if ( M_thinLayer )
    {
		M_thinLayerThickness = dataFile ( ( section + "/physics/h_thin" ).data(), 0.001 );
		M_thinLayerDensity = dataFile ( ( section + "/physics/rho_thin" ).data(), 1.2 );
		M_interfaceFlag = dataFile ( ( section + "/physics/interface" ).data(), 1 );

		Real poisson_thin = dataFile ( ( section + "/physics/poisson_thin" ).data(), 0.3 );
		Real young_thin = dataFile ( ( section + "/physics/young_thin" ).data(), 3.0e+6 );

		M_thinLayerLameI = young_thin / (2*(1.0+poisson_thin));
		M_thinLayerLameII = (young_thin*poisson_thin) / ( (1.0+poisson_thin)*(1.0-poisson_thin) );

//		std::cout << "\n\nM_thinLayer = "          << M_thinLayer << "\n"
//				  << "\n\nM_thinLayerThickness = " << M_thinLayerThickness << "\n"
//				  << "\n\nM_thinLayerDensity = "   << M_thinLayerDensity << "\n"
//				  << "\n\nM_thinLayerLameI = "     << M_thinLayerLameI << "\n"
//				  << "\n\nM_thinLayerLameII = "    << M_thinLayerLameII << "\n"
//				  << "\n\nM_interfaceFlag = "      << M_interfaceFlag << "\n";
    }
}

void
StructuralConstitutiveLawData::showMe ( std::ostream& output ) const
{
    // physics
    output << "\n*** Values for data [solid/physics]\n\n";
    output << "external pressure                = " << M_externalPressure << std::endl;
    output << "density                          = " << M_density << std::endl;
    output << "thickness                        = " << M_thickness << std::endl;
    for ( materialContainerIterator_Type i = M_young.begin() ; i != M_young.end() ; ++i )
    {
        output << "young[" << i->first << "]                         = " << i->second << std::endl;
    }
    for ( materialContainerIterator_Type i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
    {
        output << "poisson[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_bulk.begin() ; i != M_bulk.end() ; ++i )
    {
        output << "bulk[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_alpha.begin() ; i != M_alpha.end() ; ++i )
    {
        output << "alpha[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_gamma.begin() ; i != M_gamma.end() ; ++i )
    {
        output << "gamma[" << i->first << "]                       = " << i->second << std::endl;
    }

    for ( materialContainerIterator_Type i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
    {
        output << "Lame - lambda[" << i->first << "]                 = " << lambda ( i->first ) << std::endl;
        output << "Lame - mu[" << i->first << "]                     = " << mu ( i->first ) << std::endl;
    }

    for ( UInt i (0); i < M_vectorMaterialFlags.size(); i++ )
    {
        output << "Position:" << i << " -> Material Flag:            = " << M_vectorMaterialFlags[i] << std::endl;
    }


    output << "\n*** Values for data [solid/miscellaneous]\n\n";
    output << "verbose                          = " << M_verbose << std::endl;

    output << "\n*** Values for data [solid/space_discretization]\n\n";
    output << "FE order                         = " << M_order << std::endl;

    output << "\n*** Values for data [solid/time_discretization]\n\n";
    M_time->showMe ( output );
    M_timeAdvance->showMe ( output );

    output << " Informations on the constitutive law " << std::endl;
    output << " Type of constitutive law " << M_constitutiveLaw << std::endl;
    output << " Isotropic Part:  " << M_solidTypeIsotropic << std::endl;

    if( !M_constitutiveLaw.compare("anisotropic") )
      {
          output << " Anisotropic Part:  " << M_solidTypeAnisotropic << std::endl;

          for ( UInt i (0) ; i < M_numberFibers ; ++i )
          {
              std::cout << i + 1 << "-th coupled of parameters ( stiffness, nonlinearity ) : ( " << M_stiffnessParametersFibers[ i ]
                        << ", " << M_nonlinearityParametersFibers[ i ] << " ); " << std::endl;
          }

          if( !M_solidTypeAnisotropic.compare("multimechanism") ||
              !M_solidTypeAnisotropic.compare("holzapfelGeneralized"))
          {
              for ( UInt i (0) ; i < M_numberFibers ; ++i )
              {
                  std::cout << i + 1
                            << "-th characteristic stretch : " << M_characteristicStretch[ i ]
                            << std::endl;
              }
          }

          if( !M_solidTypeAnisotropic.compare("distributedHolzapfel") )
          {
              for ( UInt i (0) ; i < M_numberFibers ; ++i )
              {
                  std::cout << i + 1
                            << "-th distribution: " << M_distributionParametersFibers[ i ]
                            << std::endl;
              }
          }
      }
}

// ===================================================
// Get Method
// ===================================================
Real
StructuralConstitutiveLawData::poisson ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;

    if ( M_materialsFlagSet )
    {
        IT = M_poisson.find ( material );
    }
    else
    {
        IT = M_poisson.find ( 1 );
    }

    if ( IT != M_poisson.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the Poisson modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::young ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_young.find ( material );
    }
    else
    {
        IT = M_young.find ( 1 );
    }

    if ( IT != M_young.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the Young modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::bulk ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_bulk.find ( material );
    }
    else
    {
        IT = M_bulk.find ( 1 );
    }

    if ( IT != M_bulk.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the bulk modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::alpha ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_alpha.find ( material );
    }
    else
    {
        IT = M_alpha.find ( 1 );
    }

    if ( IT != M_alpha.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the alpha modulus has not been set !!!" << std::endl;
        return 0;
    }
}

Real
StructuralConstitutiveLawData::gamma ( const UInt& material ) const
{
    materialContainer_Type::const_iterator IT;
    if ( M_materialsFlagSet )
    {
        IT = M_gamma.find ( material );
    }
    else
    {
        IT = M_gamma.find ( 1 );
    }

    if ( IT != M_gamma.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " !!! Warning: the gamma modulus has not been set !!!" << std::endl;
        return 0;
    }
}


Real
StructuralConstitutiveLawData::lambda ( const UInt& material ) const
{
    Real youngC, poissonC;

    youngC   = this->young ( material );
    poissonC = this->poisson ( material );

    return youngC * poissonC / ( ( 1.0 + poissonC ) * ( 1.0 - 2.0 * poissonC ) );
}

Real
StructuralConstitutiveLawData::mu ( const UInt& material ) const
{
    Real youngC, poissonC;

    youngC   = this->young ( material );
    poissonC = this->poisson ( material );

    return youngC / ( 2.0 * ( 1.0 + poissonC ) );
}

}
