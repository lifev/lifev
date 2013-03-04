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

#ifndef QUADRATURE_RULE_PROVIDER_H
#define QUADRATURE_RULE_PROVIDER_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/core/fem/QRKeast.hpp>

namespace LifeV
{

//! QuadratureRuleProvider - This class is used to generate quadrature rules
/*!
    @author Samuel Quinodoz

    Thanks to this class, it is possible to adopt the best suitable quadrature
    rule for each case, depending on the degree of the expression to be
    integrated.

    Two important methods are provided:
    <ul>
    <li> provideExactness that attempts to give, for a given shape and exactness
    the best suitable quadrature rule (i.e. with the lowest possible number of points)
    <li> provideMaximal that returns the highest order quadrature available for the
    given shape.
    </ul>

    These two methods are influenced by the 3 behaviors available, which
    give the way to handle special situations:
    <ul>
        <li> NoPreciseExactness: if the exactness required is not available,
        but more precise quadrature rule exist. There are 3 possible behaviors:
        <ol>
             <li> ErrorNoPrecise: display an error and abort.
             <li> WarningAndReturnSup: display a warning and return a more precise quadrature.
             <li> ReturnSup: (default) return a more precise quadrature.
        </ol>
        <li> TooHighExactness: if there is no quadrature precise enough available.
        3 behaviors are possible:
        <ol>
             <li> ErrorTooHigh: (default) display an error and abort.
             <li> WarningAndReturnMax: display a warning and return the best quadrature available.
             <li> ReturnMax: return the best quadrature available.
        </ol>
        <li> NegativeWeight: sometimes, the quadrature with the least number of
        points contains negative weights. 3 possible choices are available in that
        case:
        <ol>
             <li> Accept: (default) return the quadrature with negative weight.
             <li> WarningAndAccept: issue a warning and return the quadrature with negative weight.
             <li> Reject: search the quadrature only among the ones with positive weights only.
        </ol>
    </ul>

    These 3 behaviors can be changed via the ad-hoc setter.

    <b> Remark </b> When running in parallel, it is a good practice
    to issue warning only with the leader process. Do not "Reject"
    negative quadratures with only one process (different processes
    would have different quadratures!).

 */
class QuadratureRuleProvider
{
public:

    //! @name Public Types
    //@{

    enum NoPreciseExactness { ErrorNoPrecise, WarningAndReturnSup, ReturnSup };
    enum TooHighExactness { ErrorTooHigh, WarningAndReturnMax, ReturnMax };
    enum NegativeWeight { Accept, WarningAndAccept, Reject};

    //@}


    //! @name Methods
    //@{

    //! Provide a quadrature rule with the given exactness
    /*!
      Given a shape, this method will try to return a quadrature rule that has the
      given exactness. If such a quadrature rule is not defined, the program will
      abort.
     */
    static QuadratureRule provideExactness (const ReferenceShapes& shape, const UInt& exactness);

    //! Provide the quadrature rule with the highest exactness available.
    static QuadratureRule provideMaximal (const ReferenceShapes& shape);

    //@}



    //! @name Operators
    //@{

    //@}


    //! @name Set Methods
    //@{


    /*!
      Setter for the behavior in case a quadrature rule with the precise exactness
      required could not be found.
     */
    static void setBehaviorNoPreciseExactness ( const NoPreciseExactness& behavior)
    {
        S_BehaviorNoPreciseExactness = behavior;
    }

    /*!
      Setter for the behavior in case the exactness required cannot be achieved
      by a known quadrature rule.
     */
    static void setBehaviorTooHighExactness ( const TooHighExactness& behavior)
    {
        S_BehaviorTooHighExactness = behavior;
    }

    /*!
      Setter for the behavior in case the quadrature rule asked for has
      negative weights for some quadrature nodes.
     */
    static void setBehaviorNegativeWeight ( const NegativeWeight& behavior)
    {
        S_BehaviorNegativeWeight = behavior;
    }


    //@}


private:

    //! @name Private Methods
    //@{

    //! Empty Constructor
    QuadratureRuleProvider();

    //! Copy Constructor
    QuadratureRuleProvider ( const QuadratureRuleProvider&);

    //! Destructor
    virtual ~QuadratureRuleProvider();

    //! Method for the differentShapes

    static QuadratureRule provideExactnessTetra (const UInt& exactness);
    static QuadratureRule provideExactnessPrism (const UInt& exactness);
    static QuadratureRule provideExactnessHexa (const UInt& exactness);
    static QuadratureRule provideExactnessQuad (const UInt& exactness);
    static QuadratureRule provideExactnessTriangle (const UInt& exactness);
    static QuadratureRule provideExactnessLine (const UInt& exactness);
    static QuadratureRule provideExactnessPoint (const UInt& exactness);

    static QuadratureRule provideExactnessTetraNoNeg (const UInt& exactness);
    static QuadratureRule provideExactnessTriangleNoNeg (const UInt& exactness);

    static void manageNoPreciseExactnessCase();
    static void manageTooHighExactnessCase();
    static void manageWarningNegativeWeight();

    //@}

    static NoPreciseExactness S_BehaviorNoPreciseExactness;
    static TooHighExactness S_BehaviorTooHighExactness;
    static NegativeWeight S_BehaviorNegativeWeight;
};

} // Namespace LifeV

#endif /* QUADRATURERULEPROVIDER_H */
