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
    @file
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 25 Nov 2010

    A more detailed description of the file (if necessary)
 */

#include <lifev/core/fem/QuadratureRuleProvider.hpp>

namespace LifeV
{

QuadratureRuleProvider::NoPreciseExactness QuadratureRuleProvider::S_BehaviorNoPreciseExactness = QuadratureRuleProvider::ReturnSup;
QuadratureRuleProvider::TooHighExactness QuadratureRuleProvider::S_BehaviorTooHighExactness = QuadratureRuleProvider::ErrorTooHigh;
QuadratureRuleProvider::NegativeWeight QuadratureRuleProvider::S_BehaviorNegativeWeight = QuadratureRuleProvider::Accept;

// ===================================================
// Methods
// ===================================================

QuadratureRule
QuadratureRuleProvider::
provideExactness (const ReferenceShapes& shape, const UInt& exactness)
{
    switch (shape)
    {
        case TETRA:
        {
            if (S_BehaviorNegativeWeight == Reject)
            {
                return provideExactnessTetraNoNeg (exactness);
            }
            return provideExactnessTetra (exactness);
            break;
        }

        case PRISM:
            return provideExactnessPrism (exactness);
            break;

        case HEXA:
            return provideExactnessHexa (exactness);
            break;

        case QUAD:
            return provideExactnessQuad (exactness);
            break;

        case TRIANGLE:
            if (S_BehaviorNegativeWeight == Reject)
            {
                return provideExactnessTriangleNoNeg (exactness);
            }
            return provideExactnessTriangle (exactness);
            break;

        case LINE:
            return provideExactnessLine (exactness);
            break;

        case POINT:
            return provideExactnessPoint (exactness);
            break;

        case NONE:
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature for this shape! " << std::endl;
            std::abort();
    };

    // Impossible case, but avoids a warning
    return QuadratureRule();
}


QuadratureRule
QuadratureRuleProvider::
provideMaximal (const ReferenceShapes& shape)
{
    switch (shape)
    {
        case TETRA:
        {
            if (S_BehaviorNegativeWeight == Reject)
            {
                return quadRuleTetra64pt;
            }
            manageWarningNegativeWeight();
            QuadratureRule qr;
            qr.import ( QRKeast<7>() );
            return qr;
            break;
        }

        case HEXA:
            return quadRuleHexa8pt;
            break;

        case QUAD:
            return quadRuleQuad9pt;
            break;

        case TRIANGLE:
            return quadRuleTria7pt;
            break;

        case LINE:
            return quadRuleTria3pt;
            break;

        case POINT:
            return quadRuleNode1pt;
            break;

        case PRISM:
        case NONE:
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature for this shape! " << std::endl;
            std::abort();
    };

    // Impossible case, but avoids a warning
    return QuadratureRule();
}

// ===================================================
// Private Methods
// ===================================================

QuadratureRule
QuadratureRuleProvider::
provideExactnessTetra (const UInt& exactness)
{
    switch (exactness)
    {
        case 0:
        case 1:
            return quadRuleTetra1pt;
            break;
        case 2:
            return quadRuleTetra4pt;
            break;
        case 3:
            manageWarningNegativeWeight();
            return quadRuleTetra5pt;
            break;
        case 4:
        {
            manageWarningNegativeWeight();
            QuadratureRule qr;
            qr.import ( QRKeast<4>() );
            return qr;
            break;
        }
        case 5:
            return quadRuleTetra15pt;
            break;
        case 6:
        {
            QuadratureRule qr;
            qr.import ( QRKeast<6>() );
            return qr;
            break;
        }

        case 7:
        {
            manageWarningNegativeWeight();
            QuadratureRule qr;
            qr.import ( QRKeast<7>() );
            return qr;
            break;
        }
        default:

            manageTooHighExactnessCase();
            manageWarningNegativeWeight();
            QuadratureRule qr;
            qr.import ( QRKeast<7>() );
            return qr;
    };

    return QuadratureRule();
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessPrism (const UInt& exactness)
{
    switch (exactness)
    {
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature for this exactness (prism) ";
            std::cerr << std::endl;
            std::abort();
    };

    /*
     * Fix to remove warning
     */
    return QuadratureRule();
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessHexa (const UInt& exactness)
{
    switch (exactness)
    {
        case 0:

        case 1:
            return quadRuleHexa1pt;
            break;

        case 2:
            manageNoPreciseExactnessCase();
            // No break here!

        case 3:
            return quadRuleHexa8pt;
            break;

        default:

            manageTooHighExactnessCase();
            return quadRuleHexa8pt;

    };

    return QuadratureRule();
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessQuad (const UInt& exactness)
{
    switch (exactness)
    {
        case 0:

        case 1:
            return quadRuleQuad1pt;
            break;

        case 2:
            manageNoPreciseExactnessCase();
            // No break!

        case 3:
            return quadRuleQuad4pt;
            break;

        case 4:
            manageNoPreciseExactnessCase();
            // No break!

        case 5:
            return quadRuleQuad9pt;
            break;

        default:

            manageTooHighExactnessCase();
            return quadRuleQuad9pt;
    };

    return QuadratureRule();
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessTriangle (const UInt& exactness)
{
    switch (exactness)
    {
        case 0:

        case 1:
            return quadRuleTria1pt;
            break;

        case 2:
            return quadRuleTria3pt;
            break;

        case 3:
            manageWarningNegativeWeight();
            return quadRuleTria4pt;
            break;

        case 4:
            return quadRuleTria6pt;
            break;

        case 5:
            return quadRuleTria7pt;
            break;

        default:
            manageTooHighExactnessCase();
            return quadRuleTria7pt;

    };

    return QuadratureRule();
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessLine (const UInt& exactness)
{
    switch (exactness)
    {
        case 0:

        case 1:
            return quadRuleSeg1pt;
            break;

        case 2:
            return quadRuleSeg2pt;
            break;

        case 3:
            return quadRuleTria3pt;
            break;

        default:

            manageTooHighExactnessCase();
            return quadRuleTria3pt;
    };

    return QuadratureRule();
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessPoint (const UInt& exactness)
{
    switch (exactness)
    {
        default:
            return quadRuleNode1pt;
    };
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessTetraNoNeg (const UInt& exactness)
{
    switch (exactness)
    {
        case 0:
        case 1:
            return quadRuleTetra1pt;
            break;
        case 2:
            return quadRuleTetra4pt;
            break;
        case 3:
            manageNoPreciseExactnessCase();
            // No break!
        case 4:
            manageNoPreciseExactnessCase();
            // No break!
        case 5:
            return quadRuleTetra15pt;
            break;
        case 6:
        {
            QuadratureRule qr;
            qr.import ( QRKeast<6>() );
            return qr;
            break;
        }

        case 7:
            return quadRuleTetra64pt;
            break;

        default:

            manageTooHighExactnessCase();
            return quadRuleTetra64pt;
    };

    return QuadratureRule();
}

QuadratureRule
QuadratureRuleProvider::
provideExactnessTriangleNoNeg (const UInt& exactness)
{
    switch (exactness)
    {
        case 0:

        case 1:
            return quadRuleTria1pt;
            break;

        case 2:
            return quadRuleTria3pt;
            break;

        case 3:
            manageNoPreciseExactnessCase();
            // No break!

        case 4:
            return quadRuleTria6pt;
            break;

        case 5:
            return quadRuleTria7pt;
            break;

        default:
            manageTooHighExactnessCase();
            return quadRuleTria7pt;

    };

    return QuadratureRule();
}



void
QuadratureRuleProvider::
manageNoPreciseExactnessCase()
{
    if ( S_BehaviorNoPreciseExactness == ErrorNoPrecise )
    {
        std::cerr << "QuadratureRuleProvider: Error: required degree does not exist" << std::endl;
        std::abort();
    }
    else if (S_BehaviorNoPreciseExactness == WarningAndReturnSup)
    {
        std::cerr << "QuadratureRuleProvider: Warning: required degree does not exist" << std::endl;
        std::cerr << "                        => attempting to find with degree+1" << std::endl;

    }
}

void
QuadratureRuleProvider::
manageTooHighExactnessCase()
{
    if (S_BehaviorTooHighExactness == ErrorTooHigh)
    {
        std::cerr << "QuadratureRuleProvider: Error: required degree too high. " << std::endl;
        std::abort();
    }
    else if (S_BehaviorTooHighExactness == WarningAndReturnMax)
    {
        std::cerr << "QuadratureRuleProvider: Warning: required degree too high. " << std::endl;
        std::cerr << "                                 => returning highest degree" << std::endl;
    }
}

void
QuadratureRuleProvider::
manageWarningNegativeWeight()
{
    if (S_BehaviorNegativeWeight == WarningAndAccept)
    {
        std::cerr << "QuadratureRuleProvider: Warning: negative weights. " << std::endl;
    }
}



} // Namespace LifeV
