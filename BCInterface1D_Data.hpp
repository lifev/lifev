//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @file
 *  @brief BCInterface_Data
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 10-05-2010
 */

#ifndef BCInterface_Data1D_H
#define BCInterface_Data1D_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>

namespace LifeV
{

//! BCInterface1D_Data - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface1D_Data class provides a general container to pass information
 *  to all the BCInterface functions.
 */
class BCInterface1D_Data
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_BaseList                                                     BCBaseList_Type;
    typedef std::vector< Real >                                                        ResistanceContainer_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface1D_Data();

    //! Copy constructor
    /*!
     * @param data BCInterface1D_Data
     */
    BCInterface1D_Data( const BCInterface1D_Data& data );

    //! Destructor
    ~BCInterface1D_Data() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterface1D_Data
     * @return reference to a copy of the class
     */
    BCInterface1D_Data& operator=( const BCInterface1D_Data& data );

    //@}


    //! @name Methods
    //@{
    /*!
     * @param FileName Name of the data file.
     * @param dataSection BC section
     * @param name name of the boundary condition
     */
    void ReadBC( const std::string& FileName, const std::string& dataSection, const BCName& name );

    //@}


    //! @name Set Methods
    //@{

    //! Set the side of the boundary condition
    /*!
     * @param flag Boundary condition side
     */
    void SetSide( const OneD_BCSide& side );

    //! Set the line of the boundary condition
    /*!
     * @param line Boundary condition line
     */
    void SetLine( const OneD_BCLine& line );

    //! Set the quantity of the boundary condition
    /*!
     * @param quantity Boundary condition quantity
     */
    void SetQuantity( const OneD_BC& quantity );

    //! Set the base string of the boundary condition
    /*!
     * @param baseString Boundary condition base string
     */
    void SetBaseString( const std::string& baseString );

    //! Set the base type of the boundary condition
    /*!
     * @param base Boundary condition base type
     */
    void SetBase( const std::pair< std::string, BCBaseList_Type >& base );

    //@}


    //! @name Get Methods
    //@{

    //! Get the flag of the boundary condition
    const OneD_BCSide& GetSide() const;

    //! Get the mode of the boundary condition
    const OneD_BCLine& GetLine() const;

    //! Get the quantity of the boundary condition
    const OneD_BC& GetQuantity() const;

    //! Get the base string of the boundary condition
    const std::string& GetBaseString() const;

    //! Get the base type of the boundary condition
    const std::pair< std::string, BCBaseList_Type >& GetBase() const;

    //! Get the resistance vector {R1, R2, R3 ...}
    const ResistanceContainer_Type& GetResistance() const;

    //@}

private:

    //! @name Private Methods
    //@{

    inline void ReadSide( const std::string& FileName, const char* side );

    inline void ReadLine( const std::string& FileName, const char* line );

    inline void ReadQuantity( const std::string& FileName, const char* quantity );

    inline void ReadBase( const std::string& FileName, const std::string& base );

    inline bool IsBase( const std::string& FileName, const char* base );

    inline void ReadResistance( const std::string& FileName, const char* resistance );

    //@}

    OneD_BCSide                               M_side;
    OneD_BCLine                               M_line;
    OneD_BC                                   M_quantity;
    std::string                               M_baseString;
    std::pair< std::string, BCBaseList_Type > M_base;

    ResistanceContainer_Type                  M_resistance;

    // Maps
    std::map< std::string, OneD_BCSide >      M_mapSide;
    std::map< std::string, OneD_BC >          M_mapQuantity;
    std::map< std::string, OneD_BCLine >      M_mapLine;
    std::map< std::string, BCBaseList_Type >  M_mapBase;
};

} // Namespace LifeV

#endif /* BCInterface_Data1D_H */
