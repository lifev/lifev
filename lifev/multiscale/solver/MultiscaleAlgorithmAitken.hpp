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
 *  @brief File containing the Multiscale Aitken Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleAlgorithmAitken_H
#define MultiscaleAlgorithmAitken_H 1

#include <lifev/core/algorithm/NonLinearAitken.hpp>

#include <lifev/multiscale/solver/MultiscaleAlgorithm.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleAlgorithmAitken - The Multiscale Algorithm implementation of Aitken
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleAlgorithmAitken is an implementation of multiscaleAlgorithm_Type
 *  which implements the Aitken method.
 */
class MultiscaleAlgorithmAitken : public virtual multiscaleAlgorithm_Type
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleAlgorithmAitken();

    //! Destructor
    virtual ~MultiscaleAlgorithmAitken() {}

    //@}


    //! @name Multiscale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param FileName Name of the data file.
     */
    void setupData ( const std::string& fileName );

    //! Perform sub-iteration on the coupling variables
    void subIterate();

    //! Display some information about the algorithm
    void showMe();

    //@}

    //! @name Set Methods
    //@{

    //! Set the the main parameters of the algorithm (tolerance, maximum number of subiterations, etc.)
    /*!
     * @param parameterList teuchos list of parameters
     */
    void setAlgorithmParameters ( const multiscaleParameterList_Type& parameterList );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleAlgorithmAitken ( const MultiscaleAlgorithmAitken& algorithm );

    MultiscaleAlgorithmAitken& operator= ( const MultiscaleAlgorithmAitken& algorithm );

    //@}

    enum methodType
    {
        Scalar, Vectorial, VectorialBlock
    };

    std::map< std::string, methodType >            M_methodMap;
    methodType                                     M_method;
    NonLinearAitken< multiscaleVector_Type >       M_generalizedAitken;

};

//! Factory create function
inline multiscaleAlgorithm_Type* createMultiscaleAlgorithmAitken()
{
    return new MultiscaleAlgorithmAitken();
}

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleAlgorithmAitken_H */
