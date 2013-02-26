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

const UInt QuadratureRuleProvider::S_maxExactnessTetra = 7;
const UInt QuadratureRuleProvider::S_maxExactnessPrism = 0;
const UInt QuadratureRuleProvider::S_maxExactnessHexa = 3;
const UInt QuadratureRuleProvider::S_maxExactnessQuad = 5;
const UInt QuadratureRuleProvider::S_maxExactnessTriangle = 5;
const UInt QuadratureRuleProvider::S_maxExactnessLine = 3;

// ===================================================
// Constructors & Destructor
// ===================================================

// ===================================================
// Operators
// ===================================================

// ===================================================
// Methods
// ===================================================

const QuadratureRule&
QuadratureRuleProvider::
provideExactness (const ReferenceShapes& shape, const UInt& exactness)
{
    switch (shape)
    {
        case TETRA:
            return provideExactnessTetra (exactness);
            break;

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
            std::cerr << " QuadratureRuleProvider: No quadrature can be furnished for this shape! " << std::endl;
            abort();
    };

    // In case you have found nothing, return the maximal
    // quadrature rule.
    return provideMaximal (shape);
}

const QuadratureRule&
QuadratureRuleProvider::
provideExactnessMax (const ReferenceShapes& shape, const UInt& exactness)
{
    switch (shape)
    {
        case TETRA:
            if (exactness <= S_maxExactnessTetra)
            {
                return provideExactnessTetra (exactness);
            }

        case PRISM:
            if (exactness <= S_maxExactnessPrism)
            {
                return provideExactnessPrism (exactness);
            }
            break;

        case HEXA:
            if (exactness <= S_maxExactnessHexa)
            {
                return provideExactnessHexa (exactness);
            }
            break;

        case QUAD:
            if (exactness <= S_maxExactnessQuad)
            {
                return provideExactnessQuad (exactness);
            }
            break;

        case TRIANGLE:
            if (exactness <= S_maxExactnessTriangle)
            {
                return provideExactnessTriangle (exactness);
            }
            break;

        case LINE:
            if (exactness <= S_maxExactnessLine)
            {
                return provideExactnessLine (exactness);
            }
            break;

        case POINT:
            // No matter the exactness, this is always exact!
            return provideExactnessPoint (exactness);
            break;

        case NONE:
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature can be furnished for this shape! " << std::endl;
            abort();
    };

    // In case you have found nothing, return the maximal
    // quadrature rule.
    return provideMaximal (shape);
}

const QuadratureRule&
QuadratureRuleProvider::
provideMaximal (const ReferenceShapes& shape)
{
    switch (shape)
    {
        case TETRA:
            return quadRuleTetra64pt;
            break;

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
            std::cerr << " QuadratureRuleProvider: No quadrature can be furnished for this shape! " << std::endl;
            abort();
    };

    return quadRuleTetra64pt;
}

// ===================================================
// Set Methods
// ===================================================

// ===================================================
// Get Methods
// ===================================================

// ===================================================
// Private Methods
// ===================================================

const QuadratureRule&
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
            return quadRuleTetra5pt;
            break;
        case 4:
        case 5:
            return quadRuleTetra15pt;
            break;
        case 6:
        case 7:
            return quadRuleTetra64pt;
            break;
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature rule can be furnished with such an exactness (tetra) ";
            std::cerr << std::endl;
            abort();
    };

    return quadRuleTetra64pt;
}

const QuadratureRule&
QuadratureRuleProvider::
provideExactnessPrism (const UInt& exactness)
{
    switch (exactness)
    {
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature rule can be furnished with such an exactness (prism) ";
            std::cerr << std::endl;
            abort();
    };

    /*
     * Fix to remove warning
     * This line should be changed when a QuadratureRule object for
     * the prism will be available!
     */
    return quadRuleTetra64pt;
}

const QuadratureRule&
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
        case 3:
            return quadRuleHexa8pt;
            break;
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature rule can be furnished with such an exactness (hexa) ";
            std::cerr << std::endl;
            abort();
    };

    return quadRuleHexa8pt;
}

const QuadratureRule&
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
        case 3:
            return quadRuleQuad4pt;
            break;
        case 4:
        case 5:
            return quadRuleQuad9pt;
            break;
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature rule can be furnished with such an exactness (quad) ";
            std::cerr << std::endl;
            abort();
    };

    return quadRuleQuad9pt;
}

const QuadratureRule&
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
            return quadRuleTria4pt;
            break;
        case 4:
            return quadRuleTria6pt;
            break;
        case 5:
            return quadRuleTria7pt;
            break;
        default:
            std::cerr << " QuadratureRuleProvider: No quadrature rule can be furnished with such an exactness (triangle) ";
            std::cerr << std::endl;
            abort();
    };

    return quadRuleTria7pt;
}

const QuadratureRule&
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
            std::cerr << " QuadratureRuleProvider: No quadrature rule can be furnished with such an exactness (line) ";
            std::cerr << std::endl;
            abort();
    };

    return quadRuleTria3pt;
}

const QuadratureRule&
QuadratureRuleProvider::
provideExactnessPoint (const UInt& exactness)
{
    switch (exactness)
    {
        default:
            return quadRuleNode1pt;
    };
}




} // Namespace LifeV
