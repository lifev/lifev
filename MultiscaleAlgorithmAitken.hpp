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
 *  @brief File containing the MultiScale Aitken Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MS_Algorithm_Aitken_H
#define MS_Algorithm_Aitken_H 1

#include <life/lifealg/generalizedAitken.hpp>

#include <lifemc/lifesolver/MultiscaleAlgorithm.hpp>

namespace LifeV
{

//! MS_Algorithm_Aitken - The MultiScale Algorithm implementation of Aitken
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Algorithm_Aitken is an implementation of MS_Algorithm
 *  which implements the Aitken method.
 */
class MS_Algorithm_Aitken : public virtual MS_Algorithm
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm_Aitken();

    //! Destructor
    virtual ~MS_Algorithm_Aitken() {}

    //@}


    //! @name MultiScale Algorithm Virtual Methods
    //@{

    //! Setup the data of the algorithm using a data file
    /*!
     * @param FileName Name of the data file.
     */
    void setupData( const std::string& fileName );

    //! Perform sub-iteration on the coupling variables
    void subIterate();

    //! Display some information about the algorithm
    void showMe();

    //@}

protected:

    //! @name Protected Methods
    //@{

    //@}

    enum methodType
    {
        Scalar, Vectorial, VectorialBlock
    };

    std::map< std::string, methodType >            M_methodMap;
    methodType                                     M_method;
    generalizedAitken< MS_Vector_Type >            M_generalizedAitken;

private:

    //! @name Unimplemented Methods
    //@{

    MS_Algorithm_Aitken( const MS_Algorithm_Aitken& algorithm );

    MS_Algorithm_Aitken& operator=( const MS_Algorithm_Aitken& algorithm );

    //@}
};

//! Factory create function
inline MS_Algorithm* createMultiscaleAlgorithmAitken()
{
    return new MS_Algorithm_Aitken();
}

} // Namespace LifeV

#endif /* MS_Algorithm_Aitken_H */
